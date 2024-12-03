import matplotlib.pyplot as plt
import pandas as pd

def plot_sensitivity(directory, output_path):
    """
    Plots the sensitvity model biomass time series alongsize the base model.
    
    Parameters:
    directory: string, directory containing the sensitivity and outputs
    output_path: string, the file path to save the plot

    Return:
    none
    """

    # Define the paths for the files
    sensitivity_file = f"{directory}/outputs-for-management.csv"
    base_file = f"{directory}/outputs-for-management-base.csv"

    # Read in the csv data
    try:
        sensitivity_data = pd.read_csv(sensitivity_file)
        base_data = pd.read_csv(base_file)
    except FileNotFoundError as err:
        raise FileNotFoundError(f"Could not find file: {err.filename}") from err

    # Verifying necessary columns are present (I'm assuming year ?)
    if not {"Years", "Median Pre-fishery biomass (in 1000s metric tons)"}.issubset(sensitivity_data.columns):
        raise ValueError("Sensitivity data must have 'Years' and 'Median Pre-fishery biomass (in 1000s metric tons)' columns.")
    if not {"Years", "Median Pre-fishery biomass (in 1000s metric tons)"}.issubset(base_data.columns):
        raise ValueError("Base data must have 'Years' and 'biomaMedian Pre-fishery biomass (in 1000s metric tons)ss' coulmns.")

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

    plt.savefig(output_path, dpi=300)
    plt.close()
