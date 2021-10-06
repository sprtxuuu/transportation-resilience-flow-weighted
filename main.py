import xc_resilience as xc

if __name__ == '__main__':
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

    # create networks and load dataset
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

    # preparedness analysis
    for network in networks:
        print(network.name, network.preparedness_gini())
    # print('LDN', BA.preparedness_gini())
    # print('LDN', DC.preparedness_gini())
    # print('LDN', QLD.preparedness_gini())
    # print('LDN', SG.preparedness_gini())

    # robustness analysis
    # BA.robustness_flow_weighted_degree_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_node_flow_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_betweenness_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_flow_centrality_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_random_attack(multi_processing=False)


    # print(BA.get_node_list())
    # print(BA.node_coordinates)
    # print(BA.od_flow)
    # print(BA.flow_weighted_restoration(0.5, 'node_flow', 2))
    # BA.robustness_flow_weighted_betweenness_based_attack()

    # BA.robustness_flow_weighted_degree_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_node_flow_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_betweenness_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_flow_centrality_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_random_attack(multi_processing=False)
    #

