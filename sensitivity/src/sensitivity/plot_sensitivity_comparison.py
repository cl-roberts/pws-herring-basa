"""This module defines a function for creating a plot of error metrics"""
import pandas as pd
import matplotlib.pyplot as plt

def plot_sens_comparison(directory, output_path):
    """This function reads a csv file and return a plot"""
    # Read in the csv file
    perc_df = pd.read_csv(directory)
    # Drop the first column (unnecessary index) if present
    if perc_df.shape[1] == 3:
        perc_df = perc_df.drop(perc_df.columns[0], axis = 1)
    # Create plot
    plt.figure(figsize = (10,6))
    plt.plot(perc_df["Year"], perc_df["Percentage Error"])
    plt.xlabel("Year")
    plt.ylabel("Percentage Error")
    plt.title("Percentage Error of Biomass: Sensitivity vs Base Model")
    plt.legend
    plt.savefig(output_path, dpi = 300)
    plt.close


# To test, use directory "../../data_outputs/sensitivity_comparison.csv"
# And for output_path, use "../../data_outputs/sensitivity_plot.png"