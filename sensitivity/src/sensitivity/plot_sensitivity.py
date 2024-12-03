import matplotlib.pyplot as plt
import pandas as pd

def plot_sensitivity(directory, output_path):
    """
    Plots the sensitvity model biomass time series alongsize the base model.
    
    Parameters:
    directory: string, directory containing the sensitivty and outputs
    output_path: string, the file path to save the plot

    Return:
    save plt.Figure: matplotlib figure object
    """

    # Define the paths for the files
    sensitivity_file = f"{directory}/outputs-for-management.csv"
    base_file = f"{directory}/outputs-for-management-base.csv"

    # Read in the csv data
    try:
        sensitivity_data = pd.read_csv(sensitivity_file)
        base_data = pd.read_csv(base_file)
    except FileNotFoundError as err:
        raise FileNotFoundError(f"Could not find file: {err.filename}")
    
    # Verifying necessary columns are present (I'm assuming year ?)
    if not {"year", "biomass"}.issubset(sensitivity_data.columns):
        raise ValueError("Sensitivity data must have 'year' and 'biomass' columns.")
    if not {"year", "biomass"}.issubset(base_data.columns):
        raise ValueError("Base data must have 'year' and 'biomass' coulmns.")

    # Plots
    plt.figure(figsize=(10,6))
    plt.plot(sensitivity_data["year"], 
             sensitivity_data["biomass"], 
             label = "Sensitivity Model",
             linestyle = "-",
             color = "blue")
    
    plt.plot(sensitivity_data["year"], 
             sensitivity_data["biomass"], 
             label = "Base Model",
             linestyle = "--",
             color = "orange")

    plt.xlabel("Year")
    plt.ylabel("Biomass")
    plt.title("Biomass Time Series: Sensitivity vs Base Models")
    plt.legend
    plt.grid(True)

    if not output_path:
        raise ValueError("Output path needs to be specified when save = True")
    else:
        plt.savefig(output_path, dpi=300)
        plt.close()
