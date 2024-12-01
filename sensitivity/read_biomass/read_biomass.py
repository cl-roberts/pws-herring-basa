import pandas as pd
def read_biomass_base():
    biomass_df = pd.read_csv('../data_outputs/outputs-for-management_base.csv')
    return biomass_df

def read_biomass_sensitivity(dir):
    filepath = '../data_outputs/' + dir + '/outputs-for-management.csv'
    sensitivity_df = pd.read_csv(filepath)
    return sensitivity_df