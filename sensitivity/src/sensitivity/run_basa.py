"""
Execute BASA
"""

import subprocess

def run_basa(dir_sensitivity):
    """
    This is a wrapper function for running the Bayesian age-structured assessment
    (BASA) model. It simply is a python wrapper for calling an R script which 
    runs the model.

    :param str dir_sensitivity: Path to run_basa.r script
    """

    cmd_run_basa = "Rscript " + dir_sensitivity + "/run_basa.r"
    cmd_save_outputs = "Rscript " + dir_sensitivity + "/plotting/plot_management_outputs.R"

    try:
        subprocess.check_output(cmd_run_basa)
        subprocess.call(cmd_save_outputs)
    except subprocess.CalledProcessError:
        print("Model did not converge, try another value for M")
        return 1

    return 0
