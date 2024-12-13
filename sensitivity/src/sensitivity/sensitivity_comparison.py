"""This module defines a function for calculating error statistics"""
import pandas as pd

def sens_comparison(dir_outputs):
    """This function reads csv, calculates error metrics, returns statistics
    
    :param str dir_outputs: path to model outputs

    """
    if isinstance(dir_outputs, str) == False:
        raise ValueError("The path should be a string.")
    
    # Read base model outputs from CSV
    base = pd.read_csv(dir_outputs + "/outputs-for-management_base.csv")
    # Read the sensitivity output
    sens = pd.read_csv(dir_outputs + "/outputs-for-management.csv")
    # Pull the median biomass from both datasets
    biomass_base = base['Median Pre-fishery biomass (in 1000s metric tons)']
    biomass_sens = sens['Median Pre-fishery biomass (in 1000s metric tons)']
    biomass_years = base['Years']
    # Calculate percentage error per year (unit is percent)
    biomass_perc_error = 100*(biomass_base - biomass_sens).abs()/biomass_base
    # Combine percentage error and years
    biomass_perc_df = pd.concat([biomass_perc_error, biomass_years], axis = 1)
    # Update column name for percentage error
    biomass_perc_df.columns = ['Percentage Error', 'Year']
    # Save as csv
    biomass_perc_df.to_csv(dir_outputs + '/sensitivity_comparison.csv')
    # Calculate mean percent error
    mean_perc_error = biomass_perc_error.mean()
    # Calculate minimum percent error
    min_perc_error = biomass_perc_error.min()
    # Calculate maximum percent error
    max_perc_error = biomass_perc_error.max()
    # Prepare the special values to be returned to user
    sens_metrics = {'Mean Perc. Error': mean_perc_error, \
    'Min Perc. Error': min_perc_error, 'Max Perc. Error': max_perc_error}
    sens_df = pd.DataFrame(sens_metrics, index = [0])
    sens_df.index.name = "Values"
    # Return the error metrics
    return sens_df

# To test, use "sensitivity/data-outputs" for dir_outputs
