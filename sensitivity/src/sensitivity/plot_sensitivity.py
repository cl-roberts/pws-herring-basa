import matplotlib.pyplot as plt

def plot_sensitivity(directory, output_path):
    """
    Plots the sensitvity model biomass time series alongsize the base model.
    
    Parameters:
    directory: string, directory containing the sensitivty and outputs
    save: boolean, if True the file will be saved instead of returned (save will default to False)
    output_path: string, the file path to save the plot when save = True

    Return:
    plt.Figure: matplotlib figure object (only if save = False)
    """

    sensitivity_data = read_biomass_sensitivity(directory)
    base_data = read_biomass_base(directory)

    # Verifying necessary columns are present (I'm assuming year ?)
    if not {"year", "biomass"}.issubset(sensitivity_data.columns):
        raise ValueError("Sensitivity data must have 'year' and 'biomass' columns.")
    if not {"year", "biomass"}.issubset(base_data.columns):
        raise ValueError("Base data must have 'year' and 'biomass' coulmns.")

    # Plots
    plt.figure(figsize=(10,6))
    plt.plot(sensitivity_data["year"], sensitivity_data["biomass"], label = "Sensitivity Model")
    plt.plot(sensitivity_data["year"], sensitivity_data["biomass"], label = "Base Model")

    plt.xlabel("Year")
    plt.ylabel("Biomass")
    plt.title("Biomass Time Series: Sensitivity vs Base Models")
    plt.legend
    plt.grid(True)

    return plt.gcf()