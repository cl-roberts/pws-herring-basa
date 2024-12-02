"""
Read biomass time series
"""

import pandas as pd

def read_biomass(dir_outputs):
    """
    This function reads the estimated biomass time series from a BASA model run

    :param str dir_outputs: Path to BASA outputs directory
    """

    biomass_df = pd.read_csv(dir_outputs + '/outputs-for-management.csv')
    return biomass_df
