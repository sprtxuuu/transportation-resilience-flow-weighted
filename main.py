import xc_resilience as xc

if __name__ == '__main__':
    # # paths to files
    ad_path, flow_path, gps_path = {}, {}, {}
    ad_path['LDN'] = 'Dataset/London Adjacency Matrix.csv'
    ad_path['BA'] = 'Dataset/BART Adjacency Matrix 2017.csv'
    ad_path['DC'] = 'Dataset/Washington DC Adjacency Matrix (1).csv'
    ad_path['QLD'] = 'Dataset/Queensland Adjacency Matrix.csv'
    ad_path['SG'] = 'Dataset/Singapore_Adjacency_Merged revised.csv'
    flow_path['LDN'] = 'Dataset/London_OD_Total_Matrix.csv'
    flow_path['BA'] = 'Dataset/BART 2017 Combined.csv'
    flow_path['DC'] = 'Dataset/WashingtonDC_OD_Matrix_Final (1).csv'
    flow_path['QLD'] = 'Dataset/Queensland_Total_Matrix.csv'
    flow_path['SG'] = 'Dataset/Singapore_OD_Matrix_Merged.csv'
    gps_path['LDN'] = "Dataset/London GPS.csv"
    gps_path['BA'] = "Dataset/BART GPS.csv"
    gps_path['DC'] = "Dataset/WashingtonDC GPS.csv"
    gps_path['QLD'] = "Dataset/Queensland GPS.csv"
    gps_path['SG'] = "Dataset/Singapore GPS March 2020 Final.csv"
    metros = ['LDN', 'BA', 'DC', 'QLD', 'SG']

    # # create networks and load dataset
    LDN = xc.Resilience('LDN')
    LDN.load_adjacency_matrix(file_path=ad_path['LDN'])
    LDN.load_flow_matrix(file_path=flow_path['LDN'])
    LDN.load_gps_coordinates(gps_path['LDN'])

    BA = xc.Resilience('BA')
    BA.load_adjacency_matrix(file_path=ad_path['BA'])
    BA.load_flow_matrix(file_path=flow_path['BA'])
    BA.load_gps_coordinates(gps_path['BA'])

    DC = xc.Resilience('DC')
    DC.load_adjacency_matrix(file_path=ad_path['DC'])
    DC.load_flow_matrix(file_path=flow_path['DC'])
    DC.load_gps_coordinates(gps_path['DC'])

    QLD = xc.Resilience('QLD')
    QLD.load_adjacency_matrix(file_path=ad_path['QLD'])
    QLD.load_flow_matrix(file_path=flow_path['QLD'])
    QLD.load_gps_coordinates(gps_path['QLD'])

    SG = xc.Resilience('SG')
    SG.load_adjacency_matrix(file_path=ad_path['SG'])
    SG.load_flow_matrix(file_path=flow_path['SG'])
    SG.load_gps_coordinates(gps_path['SG'])

    networks = [LDN, BA, DC, QLD, SG]

    # # preparedness analysis
    # # four Gini coefficient: node degree, node flow, node betweenness, node flow betweenness centrality.
    for network in networks:
        print(network.name, network.preparedness_gini())

    # # robustness analysis
    # # use multi-processing to speed up
    for network in networks:
        network.core_num = 4  # number of CPU cores for multi-processing, default value = 8
        network.robustness_flow_weighted_degree_based_attack(multi_processing=True, number_of_tests=1000)
        network.robustness_flow_weighted_node_flow_based_attack(multi_processing=True, number_of_tests=100)
        network.robustness_flow_weighted_betweenness_based_attack(multi_processing=False, number_of_tests=10)
        network.robustness_flow_weighted_flow_centrality_based_attack(multi_processing=False, number_of_tests=10)
        network.robustness_flow_weighted_random_attack(multi_processing=False)

    # # recovery analysis
    # # disruption_k: damage level (fraction of nodes removed)
    # # disruption_strategy: include 'node_degree', 'node_flow', 'node_betweenness_centrality', 'node_flow_centrality'
    # # simulation_num: number of repeated simulations
    for network in networks:
        network.flow_weighted_restoration(disruption_k=0.5, disruption_strategy='node_degree', simulation_num=10)
        network.flow_weighted_restoration(disruption_k=0.5, disruption_strategy='node_flow', simulation_num=2)
        network.flow_weighted_restoration(disruption_k=0.5, disruption_strategy='node_betweenness_centrality',
                                          simulation_num=10)
        network.flow_weighted_restoration(disruption_k=0.5, disruption_strategy='node_flow_centrality',
                                          simulation_num=2)

    # BA.robustness_flow_weighted_betweenness_based_attack()

    # adaptation analysis
    # BA.robustness_flow_weighted_degree_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_node_flow_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_betweenness_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_flow_centrality_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_random_attack(multi_processing=False)
    #
