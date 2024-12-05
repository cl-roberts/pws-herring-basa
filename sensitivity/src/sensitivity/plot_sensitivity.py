"""
Creates a time series plot for sensitivity data and base data.
"""

import matplotlib.pyplot as plt
import pandas as pd

def plot_sensitivity(input_directory, output_path):
    """
    Plots the sensitivity model biomass time series alongsize the base model.
    
    Parameters:
    input_directory: string, input_directory containing the sensitivity and outputs
    output_path: string, the file path to save the plot

    Return:
    none
    """

    # Define the paths for the files; Q: Not sure how I should read in the data? Does this right?
    sensitivity_file = f"{input_directory}/outputs-for-management.csv"
    base_file = f"{input_directory}/outputs-for-management-base.csv"

    # Read in the csv data
    try:
        sensitivity_data = pd.read_csv(sensitivity_file)
        base_data = pd.read_csv(base_file)
    except FileNotFoundError as err:
        raise FileNotFoundError(f"Could not find file: {err.filename}") from err

    # Verifying necessary columns are present
    if not {"Years", "Median Pre-fishery biomass (in 1000s metric tons)"
            }.issubset(sensitivity_data.columns):
        raise ValueError("Sensitivity data must have 'Years' and 'Median Pre-fishery biomass (in 1000s metric tons)' columns.") # Qs: IDK how to shorten this line to make pylint happy (so I'm making it longer with this comment lol)
    if not {"Years", "Median Pre-fishery biomass (in 1000s metric tons)"
            }.issubset(base_data.columns):
        raise ValueError("Base data must have 'Years' and 'Median Pre-fishery biomass (in 1000s metric tons)ss' coulmns.")

    # Plots
    plt.figure(figsize=(10,6))
    plt.plot(sensitivity_data["Years"],
             sensitivity_data["Median Pre-fishery biomass (in 1000s metric tons)"],
             label = "Sensitivity Model",
             linestyle = "-",
             color = "blue")

    plt.plot(sensitivity_data["Years"],
             sensitivity_data["Median Pre-fishery biomass (in 1000s metric tons)"],
             label = "Base Model",
             linestyle = "--",
             color = "orange")

    plt.xlabel("Years")
    plt.ylabel("Median Biomass")
    plt.title("Biomass Time Series: Sensitivity vs Base Models")
    plt.legend()
    plt.grid(True)

    plt.savefig(output_path, dpi=300) # Q: should I have some sort of return failsafe?
    plt.close()
