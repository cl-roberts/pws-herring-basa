"""
Execute BASA
"""

import subprocess
import os

def run_basa(dir_sensitivity): 
    """
    This is a wrapper function for running the Bayesian age-structured assessment
    (BASA) model. It simply is a python wrapper for calling an R script which 
    runs the model.

    :param str dir_sensitivity: Path to run_basa.r script
    """

    cmd = "Rscript " + dir_sensitivity + "/run_basa.r"


    try: 
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError:
        print("Model did not converge, try another value for M")
        return(1)

    return(0)
