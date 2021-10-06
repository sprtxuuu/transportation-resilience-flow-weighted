# version 1.5 updated 30 Sep 2021
# travel distance calculation added

# version 1.4 updated
# add unweighted relocation analysis (node level)

# version 1.3 updated
# add parameter trip_edge= {trip: [edge,...]}
# (need manual read load_pet)
# can load edge list from trip_edge dict

# version 1.2 updated
# add capacity-weighted model based on trips/routes
# add capacity-related parameters

# version 1.1 updated
# can combine two ClassResilience providing cross-layer edges

import csv
import numpy as np
import networkx as nx
import random
from collections import defaultdict
from collections import OrderedDict
import copy
from multiprocessing import Pool
import matplotlib.pyplot as plt

import pandas as pd
import math
from math import radians, cos, sin, asin, sqrt
from networkx.algorithms import approximation as approx
from tqdm import tqdm
from itertools import permutations
import json
import itertools
import multiprocessing
from itertools import combinations
from itertools import product
import matplotlib


class Resilience:
    def __init__(self, graph_name, indexing=False):
        self.name = graph_name
        self.G = nx.DiGraph()
        # geospatial data
        self.node_coordinates = {}  # (lat, lon)
        # flow-weighted model
        # self.flow_matrix = None
        self.od_flow = {}
        self.node_flow = {}
        self.node_flow_centrality = {}
        self.indexing = indexing
        if self.indexing:
            self.node2index = None
            self.index2node = None
        self._matrix_header = None  # load from adjacency matrix
        self._edge_dict = None
        self._relocation_edge_dict = None
        self._relocation_edge_weight = None
        self._restoration_edge_dict = None
        self._restoration_node_weight = None
        # capacity-weighted model based on trips/routes (GTFS...)
        self.edge_trip = {}
        self.trip_edge = {}  # {trip: list of edges}
        self.edge_param = defaultdict(dict)
        self.node_param = defaultdict(dict)
        self.edge_capacity = {}
        self.node_capacity = {}
        # multi-processing
        self.core_num = 8

    def save_weights_to_json(self):
        res = {'gps': self.node_coordinates,
               'edge_trip': self.edge_trip,
               'trip_edge': self.trip_edge,
               'edge_param': self.edge_param,
               'node_param': self.node_param,
               'edge_capacity': self.edge_capacity,
               'node_capacity': self.node_capacity}
        save_pet(res, f'xc_{self.name}_weights')

    def load_weights_from_save(self):
        res = load_pet(f'xc_{self.name}_weights')
        self.node_coordinates = res['gps']
        self.edge_trip = res['edge_trip']
        self.trip_edge = res['trip_edge']
        self.edge_param = res['edge_param']
        self.node_param = res['node_param']
        self.edge_capacity = res['edge_capacity']
        self.node_capacity = res['node_capacity']

    def get_node_list(self):
        return list(self.G.nodes)

    def get_edge_list(self):
        return list(self.G.edges)

    def get_edge_dict(self, update=True):
        if update:
            edge_dict = defaultdict(list)
            for edge in self.G.edges():
                x, y = edge[0], edge[1]
                edge_dict[x].append(y)
            self._edge_dict = edge_dict
        return self._edge_dict

    def load_adjacency_matrix(self, file_path, contain_header=True):
        node_list = []
        with open(file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            matrix = []
            for row in csv_reader:
                matrix.append(row)
        if contain_header:
            for i in range(len(matrix[0])):
                if matrix[0][i] != matrix[i][0]:
                    print('error: adjacency matrix has asymmetric headers')
            # delete headers
            for row in matrix:
                del row[0]
            header = matrix.pop(0)
            # use (network_name, node_name) represents node
            self._matrix_header = header
            node_list = [(self.name, node) for node in header]
        self.G.add_nodes_from(node_list)
        if self.indexing:
            self.node2index = {node: index for index, node in enumerate(node_list)}
            self.index2node = {index: node for index, node in enumerate(node_list)}
        for idx, x in enumerate(node_list):
            for idy, y in enumerate(node_list):
                if int(matrix[idx][idy]) > 0:
                    self.G.add_edge(x, y)
        print('\nnetwork created:',
              f'name = {self.name}, '
              f'number of nodes = {self.G.number_of_nodes()}, '
              f'number of edges = {self.G.number_of_edges()}')

    def load_edge_list(self, file_path, contain_header=True, u_col=1, v_col=2):
        with open(file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            edge_list = []
            for row in csv_reader:
                edge_list.append(row)
        if contain_header:
            del edge_list[0]  # delete headers
        node_list = []
        for row in edge_list:
            node_list.extend([row[u_col], row[v_col]])
        node_list = [(self.name, node) for node in list(set(node_list))]
        self.G.add_nodes_from(node_list)
        for row in edge_list:
            u, v = (self.name, row[u_col]), (self.name, row[v_col])
            self.G.add_edge(u, v)
        print(f'{self.name}, '
              f'number of nodes = {self.G.number_of_nodes()}, '
              f'number of edges = {self.G.number_of_edges()}')

    def load_edges_from_trip_edge(self):
        """
        use if self.trip_edge already defined
        self.trip_edge = {trip: list of edges}
        :return:
        """
        edge_pool = [edge for edge_seq in list(self.trip_edge.values())
                     for edge in edge_seq]
        self.G.add_edges_from(remove_duplicate(edge_pool))
        print(f'{self.name}, '
              f'number of nodes = {self.G.number_of_nodes()}, '
              f'number of edges = {self.G.number_of_edges()}')

    def load_edge_parameter(self, file_path, u_col=None, v_col=None, edge_col=None):
        """
        dict of dict format {(u,v):{param:value}}
        dataframe read csv; need header
        :param file_path: csv file
        :param u_col: column index for the first nodes of edges
        :param v_col: column index for the second nodes of edges
        :param edge_col: use it if edges are writen as (u, v) in single column
        :return:
        """
        if edge_col is not None:
            param = pd.read_csv(file_path, index_col=edge_col)
        else:
            param = pd.read_csv(file_path, index_col=[u_col, v_col])
        # eval("[]")
        for column in param.columns:
            for row in param.iterrows():
                param.loc[row[0], column] = eval(param.loc[row[0], column])
        param = param.to_dict(orient='index')
        self.edge_param = {edge: param[int(edge[0][1]), int(edge[1][1])]
                           for edge in self.get_edge_list()}

    def load_node_parameter(self, file_path, node_col):
        """
        dict of dict format {node:{param:value}}
        :param file_path: csv file
        :param node_col: column index for nodes' ids or names
        :return:
        """
        # dataframe read csv; need header.
        param = pd.read_csv(file_path, index_col=node_col).to_dict(orient='index')
        self.node_param = {node: param[node[1]]
                           for node in self.get_node_list()}

    def load_flow_matrix(self, file_path, contain_header=True):
        # same format as adjacency matrix
        flow_matrix = []
        with open(file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                flow_matrix.append(row)
        if contain_header:
            for i in range(1, len(flow_matrix[0])):
                if flow_matrix[0][i] != flow_matrix[i][0]:
                    print('Error: flow matrix has asymmetric headers')
                if flow_matrix[0][i].strip() != self.get_node_list()[i - 1][1].strip():
                    print(f'Warning: flow matrix has headers with unidentical station name: '
                          f'\"{flow_matrix[0][i]}\", \"{self.get_node_list()[i - 1][1]}\"')
                # if flow_matrix[0][i].strip() != self.matrix_header[i-1].strip():
                #     print(f'Warning: flow matrix has headers with unidentical station name: '
                #           f'\"{flow_matrix[0][i]}\", \"{self.matrix_header[i-1]}\"')
            for row in flow_matrix:
                del row[0]
            del flow_matrix[0]
        for x in range(len(flow_matrix)):
            for y in range(len(flow_matrix[x])):
                if len(flow_matrix[x][y].strip()) < 1 or flow_matrix[x][y] == '/':
                    # print(f"data missing at position({x, y}), header excluded..")
                    flow_matrix[x][y] = 0
                else:
                    flow_matrix[x][y] = float(flow_matrix[x][y])
        for idx in range(len(flow_matrix)):
            for idy in range(len(flow_matrix[idx])):
                x, y = self.get_node_list()[idx], self.get_node_list()[idy]
                self.od_flow[x, y] = flow_matrix[idx][idy]

    def load_gps_coordinates(self, file_path, contain_header=True,
                             node_col=0, lat_col=1, lon_col=2):
        # # table: Node, Lat, Lon
        coordinates = []
        with open(file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                coordinates.append(row)
        # print(coordinates)
        if contain_header:
            del coordinates[0]
        for row in coordinates:
            node, lat, lon = row[node_col], eval(row[lat_col]), eval(row[lon_col])
            self.node_coordinates[(self.name, node)] = (lat, lon)

    def reachable_nodes(self, origin, edge_dict=None):
        current_node = origin
        visited = []
        path = [current_node]
        if edge_dict is None:
            edge_dict = self.get_edge_dict()
        while path:
            current_node = path.pop(0)
            visited.append(current_node)
            destinations = edge_dict[current_node]
            for next_node in destinations:
                if next_node not in visited and next_node not in path:
                    path.append(next_node)
        return visited

    def get_node_degree(self):
        return dict(nx.degree(self.G))

    def get_node_betweenness_centrality(self):
        return dict(nx.betweenness_centrality(self.G))

    def get_local_node_connectivity(self):
        od_pairs = permutations(self.get_node_list(), 2)
        return {od: approx.local_node_connectivity(G=self.G, source=od[0], target=od[1]) for od in od_pairs}

    def get_node_flow(self):
        if not self.od_flow:
            return None
        self.node_flow = {}
        for u in self.get_node_list():
            self.node_flow[u] = 0
        for u in self.get_node_list():
            # for v in self.reachable_nodes(u):
            for v in self.get_node_list():
                if v != u:
                    self.node_flow[u] += self.od_flow[u, v]
        return self.node_flow

    def get_node_flow_centrality(self):
        node_list, queue = self.get_node_list(), self.get_node_list()
        if not self.od_flow:
            return None
        self.node_flow_centrality = {}
        total_flow = np.sum([value for value in self.od_flow.values()])
        sp = {}
        while queue:
            v = queue.pop(0)
            od_flow_v = 0
            for e1 in node_list:
                for e2 in node_list:
                    od_flow = self.od_flow[e1, e2]
                    if e1 != e2 and v not in [e1, e2] and od_flow > 0:
                        if (e1, e2) not in sp.keys():
                            try:
                                paths = list(nx.all_shortest_paths(self.G, e1, e2))
                                sp[e1, e2] = paths
                            except:
                                paths = None
                                sp[e1, e2] = paths
                        else:
                            paths = sp[e1, e2]
                        if paths:
                            nosp = len(paths)
                            for path in paths:
                                if v in path:
                                    od_flow_v += od_flow / nosp
            if total_flow == 0:
                self.node_flow_centrality[v] = 0
            else:
                self.node_flow_centrality[v] = od_flow_v / total_flow
        return self.node_flow_centrality

    def get_travel_distance(self, mean_value=False):
        spl = dict(nx.all_pairs_shortest_path_length(self.G))
        travel_distance_distribution = {}
        for od_pair, flow in self.od_flow.items():
            trip_len = spl[od_pair[0]][od_pair[1]]
            if trip_len in travel_distance_distribution.keys():
                travel_distance_distribution[trip_len] += flow
            else:
                travel_distance_distribution[trip_len] = flow
        # print(travel_distance_distribution)
        if mean_value:
            total_trip_len, total_flow = 0.0, 0.0
            for trip_len, flow in travel_distance_distribution.items():
                total_trip_len += trip_len * flow
                total_flow += flow
            mean_trip_len = total_trip_len / total_flow
            return round(mean_trip_len, 3)
        else:
            return travel_distance_distribution

    def preparedness_node_degree_gini(self):
        node_degree = self.get_node_degree()
        return gini([value for value in node_degree.values()])

    def preparedness_node_betweenness_centrality_gini(self):
        node_betweenness_centrality = self.get_node_betweenness_centrality()
        return gini([value for value in node_betweenness_centrality.values()])

    def preparedness_node_flow_gini(self):
        if self.od_flow:
            node_flow = self.get_node_flow()
            return gini([value for value in node_flow.values()])
        else:
            return None

    def preparedness_node_flow_centrality_gini(self):
        if self.od_flow:
            node_flow_centrality = self.get_node_flow_centrality()
            return gini([value for value in node_flow_centrality.values()])
        else:
            return None

    def preparedness_gini(self):
        return {'Gini_ND': self.preparedness_node_degree_gini(),
                'Gini_BC': self.preparedness_node_betweenness_centrality_gini(),
                'Gini_NF': self.preparedness_node_flow_gini(),
                'Gini_FC': self.preparedness_node_flow_centrality_gini()}

    def _unweighted_sequential_attack(self, strategy=None, multiple_removal=1):
        """
        single scenario test
        :param strategy: 'node_degree' or 'node_betweenness_centrality'
        :param multiple_removal: step size, increase it to reduce computational time
        :return: list of y value in a curve
        """
        temp_g = copy.deepcopy(self.G)
        node_list, n0 = list(temp_g.nodes()), temp_g.number_of_nodes()
        removed_list, degradation_curve = [], [1.0]
        step = 0
        while node_list:
            if strategy == 'node_degree':
                nd_dic = dict(nx.degree(temp_g))
                i = search_for_max(nd_dic, multiple_search=multiple_removal)
            elif strategy == 'node_betweenness_centrality':
                bc_dic = nx.betweenness_centrality(temp_g)
                i = search_for_max(bc_dic, multiple_search=multiple_removal)
            else:
                if multiple_removal > len(node_list):
                    i = copy.deepcopy(node_list)
                else:
                    i = random.sample(node_list, k=multiple_removal)
            if type(i) is list:
                removed_list.extend(i)
                temp_g.remove_nodes_from(i)
            else:
                removed_list.append(i)
                temp_g.remove_node(i)
            node_list = list(temp_g.nodes())
            step += 1
            if node_list:  # identify giant components
                gcc = sorted(nx.strongly_connected_components(temp_g), key=len, reverse=True)
                n_g0 = temp_g.subgraph(gcc[0]).number_of_nodes()
                degradation_curve.append(n_g0 / n0)
            else:
                degradation_curve.append(0)
        return degradation_curve

    def _flow_weighted_sequential_attack(self, strategy=None, multiple_removal=1):
        """
        single scenario test
        :param strategy: 'node_degree' or 'node_betweenness_centrality'
        :param multiple_removal: step size, increase it to reduce computational time
        :return: list of y value in a curve
        """
        temp = copy.deepcopy(self)
        node_list, n0 = list(temp.G.nodes()), temp.G.number_of_nodes()
        removed_list, degradation_curve = [], [1.0]
        step = 0
        total_flow = np.sum([value for value in self.od_flow.values()])
        while node_list:
            if strategy == 'node_degree':
                nd_dic = dict(nx.degree(temp.G))
                i = search_for_max(nd_dic, multiple_search=multiple_removal)
            elif strategy == 'node_betweenness_centrality':
                bc_dic = nx.betweenness_centrality(temp.G)
                i = search_for_max(bc_dic, multiple_search=multiple_removal)
            elif strategy == 'node_flow':
                nf_dic = temp.get_node_flow()
                i = search_for_max(nf_dic, multiple_search=multiple_removal)
            elif strategy == 'node_flow_centrality':
                fc_dic = temp.get_node_flow_centrality()
                i = search_for_max(fc_dic, multiple_search=multiple_removal)
            else:
                if multiple_removal > len(node_list):
                    i = copy.deepcopy(node_list)
                else:
                    i = random.sample(node_list, k=multiple_removal)
            if type(i) is list:
                removed_list.extend(i)
                temp.G.remove_nodes_from(i)
            else:
                removed_list.append(i)
                temp.G.remove_node(i)
            node_list = list(temp.G.nodes())
            step += 1
            if node_list:  # identify remaining flow
                remaining_flow = 0
                for origin in node_list:
                    reached = temp.reachable_nodes(origin)
                    for destination in reached:
                        remaining_flow += temp.od_flow[origin, destination]
                degradation_curve.append(remaining_flow / total_flow)
            else:
                degradation_curve.append(0)
        return degradation_curve

    def _robustness_repeated_tests(self, strategy, weight,
                                   number_of_tests=100,
                                   multiple_removal=1,
                                   multi_processing=False):
        if multi_processing:
            curves = []
            pool = Pool(processes=self.core_num)
            if weight:
                pool_result = [pool.apply_async(self._flow_weighted_sequential_attack,
                                                args=(strategy, multiple_removal)) for test in
                               tqdm(range(number_of_tests))]
            else:
                pool_result = [pool.apply_async(self._unweighted_sequential_attack,
                                                args=(strategy, multiple_removal)) for test in range(number_of_tests)]
            pool.close()
            pool.join()
            for nt in pool_result:
                nt_curve = nt.get()
                curves.append(nt_curve)
            average_curve = np.mean(curves, axis=0)
            non = len(self.get_node_list())
            xs = [x / non for x in range(0, non, multiple_removal)]
            xs.append(1.0)
            print(f'rb_{strategy} =', numerical_integral_nml(average_curve, xs=xs))
        else:
            curves = []
            for test in range(number_of_tests):
                if weight:
                    curves.append(self._flow_weighted_sequential_attack(
                        strategy=strategy, multiple_removal=multiple_removal))
                else:
                    curves.append(self._unweighted_sequential_attack(
                        strategy=strategy, multiple_removal=multiple_removal))
            average_curve = np.mean(curves, axis=0)
            non = len(self.get_node_list())
            xs = [x / non for x in range(0, non, multiple_removal)]
            xs.append(1.0)
            print(f'rb_{strategy} =', numerical_integral_nml(average_curve, xs=xs))
        return [list(average_curve), xs]

    def robustness_unweighted_degree_based_attack(self, number_of_tests=100,
                                                  multiple_removal=1,
                                                  multi_processing=False):
        strategy, weight = 'node_degree', False
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_flow_weighted_degree_based_attack(self, number_of_tests=100,
                                                     multiple_removal=1,
                                                     multi_processing=False):
        strategy, weight = 'node_degree', True
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_unweighted_betweenness_based_attack(self, number_of_tests=100,
                                                       multiple_removal=1,
                                                       multi_processing=False):
        strategy, weight = 'node_betweenness_centrality', False
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_flow_weighted_betweenness_based_attack(self, number_of_tests=100,
                                                          multiple_removal=1,
                                                          multi_processing=False):
        strategy, weight = 'node_betweenness_centrality', True
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_unweighted_random_attack(self, number_of_tests=1000,
                                            multiple_removal=1,
                                            multi_processing=False):
        strategy, weight = None, False
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_flow_weighted_random_attack(self, number_of_tests=1000,
                                               multiple_removal=1,
                                               multi_processing=False):
        strategy, weight = None, True
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_flow_weighted_node_flow_based_attack(self, number_of_tests=100,
                                                        multiple_removal=1,
                                                        multi_processing=False):
        strategy, weight = 'node_flow', True
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def robustness_flow_weighted_flow_centrality_based_attack(self, number_of_tests=100,
                                                              multiple_removal=1,
                                                              multi_processing=False):
        strategy, weight = 'node_flow_centrality', True
        return self._robustness_repeated_tests(strategy=strategy, weight=weight,
                                               number_of_tests=number_of_tests,
                                               multiple_removal=multiple_removal,
                                               multi_processing=multi_processing)

    def _attack_sequence_generation(self, strategy, multiple_removal=1):
        temp = copy.deepcopy(self)
        node_list = copy.deepcopy(self.get_node_list())
        attack_sequence = []
        while node_list:
            if strategy == 'node_degree':
                nd_dic = dict(nx.degree(temp.G))
                i = search_for_max(nd_dic, multiple_search=multiple_removal)
            elif strategy == 'node_betweenness_centrality':
                bc_dic = nx.betweenness_centrality(temp.G)
                i = search_for_max(bc_dic, multiple_search=multiple_removal)
            elif strategy == 'node_flow':
                nf_dic = temp.get_node_flow()
                i = search_for_max(nf_dic, multiple_search=multiple_removal)
            elif strategy == 'node_flow_centrality':
                fc_dic = temp.get_node_flow_centrality()
                i = search_for_max(fc_dic, multiple_search=multiple_removal)
            else:
                if multiple_removal > len(node_list):
                    i = copy.deepcopy(node_list)
                else:
                    i = random.sample(node_list, k=multiple_removal)
            if type(i) is list:
                attack_sequence.extend(i)
                temp.G.remove_nodes_from(i)
            else:
                attack_sequence.append(i)
                temp.G.remove_node(i)
            node_list = temp.get_node_list()
        return attack_sequence

    def _relocation_dijsktra_weighted(self, initial, end):
        # shortest paths is a dict of nodes
        # whose value is a tuple of (previous node, weight)
        edge_dict = self._relocation_edge_dict
        edge_weight = self._relocation_edge_weight
        shortest_paths = {initial: (None, 0)}
        current_node = initial
        visited = set()
        while current_node != end:
            visited.add(current_node)
            destinations = edge_dict[current_node]
            weight_to_current_node = shortest_paths[current_node][1]
            for next_node in destinations:
                weight = edge_weight[(current_node, next_node)] + weight_to_current_node
                # weight = 1 + weight_to_current_node
                if next_node not in shortest_paths:
                    shortest_paths[next_node] = (current_node, weight)
                else:
                    current_shortest_weight = shortest_paths[next_node][1]
                    if current_shortest_weight > weight:
                        shortest_paths[next_node] = (current_node, weight)
            next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
            if not next_destinations:
                return math.inf
            # next node is the destination with the lowest weight
            current_node = min(next_destinations, key=lambda k: next_destinations[k][1])
        return shortest_paths[end][1]

    def _haversine_distance_between_nodes(self, node_x, node_y):
        lat1, lon1 = self.node_coordinates[node_x]
        lat2, lon2 = self.node_coordinates[node_y]
        return haversine(lat1, lon1, lat2, lon2)

    def unweighted_relocation(self):
        """
        :return: dict of node relocation
        """

        def df(d):
            will = 1 - d / 1600
            if 0 <= will <= 1:
                return will
            else:
                return 0

        od_pairs = permutations(self.get_node_list(), 2)
        distance_matrix = {od: self._haversine_distance_between_nodes(od[0], od[1]) for od in od_pairs}
        relocation_matrix = {od: df(dst) for od, dst in distance_matrix.items()}
        relocation_potential = {}
        for node in self.get_node_list():
            relocation_potential[node] = np.sum(
                [relocation_matrix[node, other] for other in self.get_node_list() if node != other])
        return relocation_potential

    def _flow_weighted_relocation(self, v_a, v_b, cum_dst_prev=None):
        """
        used in recovery()
        functions: distance(), dijsktra_weighted(), df().
        :return:
        """

        def df(d):
            will = 1 - d / 1600
            if 0 <= will <= 1:
                return will
            else:
                return 0

        cum_dst = {}  # cumulative walking distance
        if cum_dst_prev is not None:
            for od_pair, cd_value in cum_dst_prev.items():
                # in restoration scenario, we keep value for intact/recovered path
                if cd_value == 0:
                    cum_dst[od_pair] = 0
        self._relocation_edge_dict = defaultdict(list)  # edge_dict[current_node]=next_nodes`
        self._relocation_edge_weight = {}  # edge_weight[node1,node2]=walking_distance
        max_flows = np.sum([value for value in self.od_flow.values()])
        for x in v_a:
            for y in v_a:
                if y in self._edge_dict[x]:  # connected by operational metro line
                    self._relocation_edge_dict[x].append(y)
                    self._relocation_edge_weight[x, y] = 0  # no walking distance
            for z in v_b:
                dst = self._haversine_distance_between_nodes(x, z)  # walking distance
                if dst <= 1600:  # introduce a relocation edge if no farther than 1600 metres
                    self._relocation_edge_dict[x].append(z)
                    self._relocation_edge_dict[z].append(x)
                    self._relocation_edge_weight[x, z] = dst
                    self._relocation_edge_weight[z, x] = dst
        total_relocation = 0
        for od_pair, flow in self.od_flow.items():
            # print('Debug: cd =', cd)
            if od_pair not in cum_dst.keys():  # cd=cumulative walking distance for the shortest path
                cum_dst[od_pair] = self._relocation_dijsktra_weighted(od_pair[0], od_pair[1])
            elif cum_dst[od_pair] != 0:  # keep the 0 value for intact/recovered path
                # if not 0 value, recalculate: (inf) previously disconnected or (>0) need walking
                cum_dst[od_pair] = self._relocation_dijsktra_weighted(od_pair[0], od_pair[1])
            if cum_dst[od_pair] > 0:
                will = df(cum_dst[od_pair])  # relocation rate based on cd
                relocation = will * flow  # relocated flow for the od_pair
                total_relocation += relocation
        return total_relocation / max_flows, cum_dst

    def _restoration_dijsktra_weighted(self, initial, end):
        # shortest paths is a dict of nodes
        # whose value is a tuple of (previous node, weight)
        node_weight = self._restoration_node_weight
        edge_dict = self._restoration_edge_dict
        shortest_paths = {initial: (None, node_weight[initial])}
        current_node = initial
        visited = set()
        while current_node != end:
            visited.add(current_node)
            # destinations = graph.edges[current_node]
            destinations = edge_dict[current_node]
            weight_to_current_node = shortest_paths[current_node][1]
            for next_node in destinations:
                weight = node_weight[next_node] + weight_to_current_node
                # weight = 1 + weight_to_current_node
                if next_node not in shortest_paths:
                    shortest_paths[next_node] = (current_node, weight)
                else:
                    current_shortest_weight = shortest_paths[next_node][1]
                    if current_shortest_weight > weight:
                        shortest_paths[next_node] = (current_node, weight)
            next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
            if not next_destinations:
                print(f'error:{initial}---{end}')
                return False
            # next node is the destination with the lowest weight
            current_node = min(next_destinations, key=lambda k: next_destinations[k][1])
        # Work back through destinations in shortest path
        path = []
        while current_node is not None:
            path.append(current_node)
            next_node = shortest_paths[current_node][0]
            current_node = next_node
        # Reverse path
        path = path[::-1]
        # return path
        return shortest_paths[end][1], path

    def flow_weighted_restoration(self, disruption_k, disruption_strategy, simulation_num):
        def rafr_total_flow_calculator(flow_dict, od_cost, node_list):
            max_flows = 0
            for key, item in flow_dict.items():
                max_flows += item
            total_flows = 0
            for x in node_list:
                for y in node_list:
                    if od_cost[x, y] < 1:  # edges are bi-directional
                        total_flows += flow_dict[x, y]
            return total_flows / max_flows  # normalized flow

        def linear_interpolation(curve):
            pn = 0
            while pn < len(curve) - 1:
                x1, x2 = curve[pn][0], curve[pn + 1][0]
                dx = x2 - x1
                if dx > 1:
                    y1, y2, z1, z2 = curve[pn][1], curve[pn + 1][1], curve[pn][2], curve[pn + 1][2]
                    dy = (y2 - y1) / dx
                    dz = (z2 - z1) / dx
                    for m in range(1, dx):
                        curve.insert(pn + m, (x1 + m, y1 + m * dy, z1 + m * dz))
                    pn += dx
                    # print(f'Warning: linear interpolation applied between point {x1} and {x1 + dx}')
                else:
                    pn += 1
            return curve

        curves = []
        k_nodes = int(self.G.number_of_nodes() * disruption_k)
        for simulation_no in range(simulation_num):
            print(f'Simulation No. {simulation_no + 1} / {simulation_num}')
            v_b = self._attack_sequence_generation(strategy=disruption_strategy, multiple_removal=1)[0:k_nodes]
            v_a = []
            self._restoration_node_weight = {}
            # used as node_weight for path finding; 0 for operational, 1 for disrupted.
            for node in self.get_node_list():
                if node not in v_b:
                    v_a.append(node)
                    self._restoration_node_weight[node] = 0
                else:
                    self._restoration_node_weight[node] = 1
            print(f'{len(v_a)} operational nodes, {len(v_b)} disrupted nodes')
            pbar = tqdm(total=len(v_b), desc=f"Simulation No.{simulation_no + 1}/{simulation_num}")
            self._restoration_edge_dict = self.get_edge_dict()
            od_cost, paths = {}, {}  # initialize OD path recovery cost; represent the number of disrupted nodes on path
            curve = []
            steps = 0  # x value for curve; cumulative od_cost in steps
            cum_dst = None  # cumulative walking distance for relocation
            while v_b:  # recovery cycle
                for x in self.get_node_list():  # od_cost computation
                    for y in self.get_node_list():
                        # initialize or update od_cost and paths
                        if (x, y) not in od_cost.keys():
                            # od pairs are kept which have been connected in previous cycles
                            od_cost[x, y], paths[x, y] = self._restoration_dijsktra_weighted(x, y)
                            if od_cost[x, y] == math.inf:
                                print(f'error: path ({x},{y}) not found; ad_matrix may not be complete')
                            else:
                                od_cost[y, x] = od_cost[x, y]  # if all edges are bi-directional
                                paths[y, x] = paths[x, y]
                # evaluate and save initial/current recovery step (determined by previous loop)
                rl, cum_dst = self._flow_weighted_relocation(v_a, v_b, cum_dst)
                rt = rafr_total_flow_calculator(self.od_flow, od_cost, self.get_node_list())
                curve.append([steps, rt, rt + rl])  # update curve
                # determine the next step based on cost_benefit
                cost_benefit = {}
                for x in self.get_node_list():
                    for y in self.get_node_list():
                        if od_cost[x, y] < 1:  # i.e. path connected, no need for recovery
                            cost_benefit[x, y] = -1
                        else:
                            cost_benefit[x, y] = (self.od_flow[x, y] + self.od_flow[y, x]) / od_cost[x, y]
                # identify prior od pair to connect
                best_od = max(cost_benefit, key=lambda p: cost_benefit[p])
                # print(prior_od, cost_benefit[prior_od], od_cost[prior_od])
                step_length = od_cost[best_od]
                steps += step_length  # update steps; increased value = od_cost
                pbar.update(step_length)
                path = paths[best_od]  # the least cost path
                for node in path:  # update node sets, node state and edge_dict
                    if node in v_b:
                        v_a.append(node)
                        v_b.remove(node)
                        self._restoration_node_weight[node] = 0
                od_cost_prev = od_cost  # keep copy of od_cost
                od_cost = {}  # reset od_cost
                for x in self.get_node_list():  # keep od pairs that has already been connected
                    for y in self.get_node_list():
                        if od_cost_prev[x, y] == 0:
                            od_cost[x, y] = 0
            # evaluate and save the last step
            for x in self.get_node_list():  # od_cost computation for od pairs previously unconnected
                for y in self.get_node_list():
                    # od pairs are kept which have been connected in previous cycles
                    if (x, y) not in od_cost.keys():
                        od_cost[x, y], paths[x, y] = self._restoration_dijsktra_weighted(x, y)
                        if od_cost[x, y] is False:
                            print(f'error: path ({x},{y}) not found; ad_matrix may not be complete')
                        else:
                            od_cost[y, x] = od_cost[x, y]  # if all edges are bi-directional
                            paths[y, x] = paths[x, y]
            rl, cum_dst = self._flow_weighted_relocation(v_a, v_b, cum_dst)
            rt = rafr_total_flow_calculator(self.od_flow, od_cost, self.get_node_list())
            curve.append((steps, rt, rt + rl))  # update curve
            pbar.close()
            curves.append(curve)
        for cn in range(len(curves)):
            curves[cn] = linear_interpolation(curves[cn])
        curves = np.mean(curves, axis=0)
        return curves

    def plot(self, show=False, save=False, save_path=None,
             with_labels=False, caption=False, resize=1.0, legend=False,
             alpha=1,
             node_size=3,
             node_color='#1f78b4',
             node_shape='o',
             linewidths=0.2,
             width=0.5,
             edge_color='dimgray',
             style='solid',
             arrowsize=0.5,
             font_size=6,
             edge_cmap=None,
             edge_vmin=None,
             edge_vmax=None):
        plt.figure(1, figsize=(12, 12))
        pos = {key: (value[1], value[0]) for key, value in self.node_coordinates.items()}
        nx.draw_networkx(self.G, pos=pos,
                         with_labels=with_labels,
                         alpha=alpha,
                         node_size=node_size,
                         node_color=node_color,
                         node_shape=node_shape,
                         linewidths=linewidths,
                         width=width,
                         edge_color=edge_color,
                         style=style,
                         arrowsize=arrowsize,
                         font_size=font_size,
                         font_family='Times',
                         label=self.name,
                         edge_cmap=edge_cmap,
                         edge_vmin=edge_vmin,
                         edge_vmax=edge_vmax)
        if legend:
            plt.legend()
        if caption:
            plt.title(self.name)
        if save_path is not None:
            plt.savefig(save_path, transparent=True, dpi=300)
        elif save:
            plt.savefig(f'fig_{self.name}.png', transparent=True, dpi=300 * resize)
        if show:
            plt.show()


def search_for_max(dic, multiple_search=1):
    # list(dict(sorted(dic.items(), key=lambda item: item[1])).keys())[-1]
    result = []
    temp_dic = copy.deepcopy(dic)
    while temp_dic and len(result) < multiple_search:
        max_value = max(temp_dic.values())
        max_value_keys = []
        for key, value in temp_dic.items():
            if value == max_value:
                max_value_keys.append(key)
        the_key = random.choice(max_value_keys)
        result.append(the_key)
        temp_dic.pop(the_key)
    return result


def remove_duplicate(a):
    return list(OrderedDict.fromkeys(a))


def get_edge_dict_from(G):
    edge_dict = defaultdict(list)
    for edge in G.edges():
        x, y = edge[0], edge[1]
        edge_dict[x].append(y)
    return edge_dict


def print_dict(dic, num):
    inum = 0
    for key, values in dic.items():
        if inum > num:
            break
        print(key, values)
        inum += 1


def numerical_integral_nml(curve, xs=None, dec=None):
    """
    Riemann sum Midpoint rule
    :param curve:
    :param xs:
    :param dec:
    :return:
    """
    y = np.asarray(curve)
    if xs is not None:
        if len(curve) != len(xs):
            print('error: curve and xs have different length')
            return
        else:
            x = np.asarray(xs)
            nml = 0.5 * ((x[1] - x[0]) * y[0] + np.sum((x[2:] - x[:-2]) * y[1:-1]) + (x[-1] - x[-2]) * y[-1])
    else:
        n = len(y) - 1
        nml = (0.5 * (y[0] + y[-1]) + np.sum(y[1:-1])) / n
    if dec:
        return round(nml, dec)
    else:
        return nml


def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    array = np.asarray(array).astype(float).flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1, array.shape[0] + 1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))


def dijsktra_weighted(edge_dict, edge_weight, initial, end):
    # shortest paths is a dict of nodes
    # whose value is a tuple of (previous node, weight)
    shortest_paths = {initial: (None, 0)}
    current_node = initial
    visited = set()
    while current_node != end:
        visited.add(current_node)
        # destinations = graph.edges[current_node]
        destinations = edge_dict[current_node]
        weight_to_current_node = shortest_paths[current_node][1]
        for next_node in destinations:
            weight = edge_weight[(current_node, next_node)] + weight_to_current_node
            # weight = 1 + weight_to_current_node
            if next_node not in shortest_paths:
                shortest_paths[next_node] = (current_node, weight)
            else:
                current_shortest_weight = shortest_paths[next_node][1]
                if current_shortest_weight > weight:
                    shortest_paths[next_node] = (current_node, weight)
        next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
        if not next_destinations:
            return math.inf
        # next node is the destination with the lowest weight
        current_node = min(next_destinations, key=lambda k: next_destinations[k][1])
    # Work back through destinations in shortest path
    # path = []
    # while current_node is not None:
    #     path.append(current_node)
    #     next_node = shortest_paths[current_node][0]
    #     current_node = next_node
    # Reverse path
    # path = path[::-1]
    # return path
    return shortest_paths[end][1]


def haversine(lat1, lon1, lat2, lon2):  # (decimal）
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6371
    dst = c * r * 1000
    return dst


def combine_network_weighted(XC1, XC2, cross_layer_edges, new_name='Combined_network'):
    # not support flow
    CN = Resilience(new_name)
    CN.G.add_edges_from(list(XC1.G.edges()))
    CN.G.add_edges_from(list(XC2.G.edges()))
    CN.G.add_edges_from(cross_layer_edges)
    CN.node_coordinates.update(XC1.node_coordinates)
    CN.node_coordinates.update(XC2.node_coordinates)
    print(f'{CN.name}, '
          f'number of nodes = {CN.G.number_of_nodes()}, '
          f'number of edges = {CN.G.number_of_edges()}')
    return CN


def export_list(xp_list, filename, mode='w'):
    with open(filename, mode, encoding='utf-8-sig', newline='') as output_file:
        print(xp_list)
        csv_writer = csv.writer(output_file, delimiter=',')
        for x in xp_list:
            csv_writer.writerow(x)
    print('file:', filename, 'created')


def save_pet(pet, filename='temporary file'):
    with open(filename, 'w') as f:
        f.write(json.dumps(str(pet)))


def load_pet(filename):
    with open(filename) as f:
        pet = json.loads(f.read())
    return eval(pet)


def revert_dict_of_list(original_dict):
    """
    :param original_dict: a dict of lists
    """
    new_dict = defaultdict(list)
    for key, value in original_dict.items():
        for v in value:
            new_dict[v].append(key)
    return dict(new_dict)
