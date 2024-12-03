"""This module defines a function for calculating error statistics"""
import pandas as pd

def sens_comparison(sens_run):
    """This function reads csv, calculates error metrics, returns stats"""
    # Read base model outputs from CSV
    base = pd.read_csv("../../data_outputs/outputs-for-management_base.csv")
    # Read the sensitivity output
    sens = pd.read_csv(sens_run)
    # Pull the median biomass from both datasets
    biomass_base = base['Median Pre-fishery biomass (in 1000s metric tons)']
    biomass_sens = sens['Median Pre-fishery biomass (in 1000s metric tons)']
    # Calculate percentage error per year (unit is percent)
    biomass_perc_error = 100*(biomass_base - biomass_sens).abs()/biomass_base
    # Calculate mean percent error
    mean_perc_error = biomass_perc_error.mean()
    # Calculate minimum percent error
    min_perc_error = biomass_perc_error.min()
    # Calculate maximum percent error
    max_perc_error = biomass_perc_error.max()
    # Write the percentage errors to a csv file
    biomass_perc_error.to_csv('../../data_outputs/sensitivity_comparison.csv')
    # Prepare the special values to be returned to user
    sens_metrics = {'Mean Perc. Error': mean_perc_error, \
    'Min Perc. Error': min_perc_error, 'Max Perc. Error': max_perc_error}
    sens_df = pd.DataFrame(sens_metrics, index = [0])
    sens_df.index.name = "Values"
    # Return the error metrics
    return(sens_df)

# To test, use "../../data_outputs/outputs-for-management.csv" for sens_run
