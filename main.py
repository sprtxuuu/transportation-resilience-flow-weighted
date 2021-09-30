import xc_resilience as xc

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    ad_path, f_path, gps_path = {}, {}, {}
    ad_path['LDN'] = 'Metro networks project/London Adjacency Matrix.csv'
    ad_path['BA'] = 'Metro networks project/BART Adjacency Matrix 2017.csv'
    ad_path['DC'] = 'Metro networks project/Washington DC Adjacency Matrix (1).csv'
    ad_path['QLD'] = 'Metro networks project/Queensland Adjacency Matrix.csv'
    ad_path['SG'] = 'Metro networks project/Singapore_Adjacency_Merged revised.csv'
    f_path['LDN'] = 'Metro networks project/London_OD_Total_Matrix.csv'
    f_path['BA'] = 'Metro networks project/BART 2017 Combined.csv'
    f_path['DC'] = 'Metro networks project/WashingtonDC_OD_Matrix_Final (1).csv'
    f_path['QLD'] = 'Metro networks project/Queensland_Total_Matrix.csv'
    f_path['SG'] = 'Metro networks project/Singapore_OD_Matrix_Merged.csv'
    gps_path['LDN'] = "Metro networks project/London GPS.csv"
    gps_path['BA'] = "Metro networks project/BART GPS.csv"
    gps_path['DC'] = "Metro networks project/WashingtonDC GPS.csv"
    gps_path['QLD'] = "Metro networks project/Queensland GPS.csv"
    gps_path['SG'] = "Metro networks project/Singapore GPS March 2020 Final.csv"
    metros = ['LDN', 'BA', 'DC', 'QLD', 'SG']

    # LDN = Resilience('LDN')
    # LDN.load_adjacency_matrix(file_path=ad_path['LDN'])
    # LDN.load_flow_matrix(file_path=f_path['LDN'])

    BA = xc.Resilience('BA')
    BA.load_adjacency_matrix(file_path=ad_path['BA'])
    BA.load_flow_matrix(file_path=f_path['BA'])
    # print(BA.preparedness_gini())
    BA.load_gps_coordinates(gps_path['BA'], contain_header=True)
    print(BA.get_node_list())
    print(BA.node_coordinates)
    print(BA.od_flow)
    print(BA.flow_weighted_restoration(0.5, 'node_flow', 2))
    # BA.robustness_flow_weighted_betweenness_based_attack()

    # BA.robustness_flow_weighted_degree_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_node_flow_based_attack(multi_processing=True)
    # BA.robustness_flow_weighted_betweenness_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_flow_centrality_based_attack(multi_processing=False, number_of_tests=1)
    # BA.robustness_flow_weighted_random_attack(multi_processing=False)
    #
    # DC = Resilience('DC')
    # DC.load_adjacency_matrix(file_path=ad_path['DC'])
    # DC.load_flow_matrix(file_path=f_path['DC'])
    #
    # QLD = Resilience('QLD')
    # QLD.load_adjacency_matrix(file_path=ad_path['QLD'])
    # QLD.load_flow_matrix(file_path=f_path['QLD'])
    #
    # SG = Resilience('SG')
    # SG.load_adjacency_matrix(file_path=ad_path['SG'])
    # SG.load_flow_matrix(file_path=f_path['SG'])

