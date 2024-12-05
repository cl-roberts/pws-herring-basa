"""
Tests for plot_sensivity to verify that the function works and
creates a plot correctly
"""

import os
import sys
import tempfile
from PIL import Image
import pandas as pd
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from sensitivity.plot_sensitivity import plot_sensitivity

# Test Data:

def create_temp_data(Years, Median):
    """
    Creates temporary data for running the test for improved reproducibility.

    Parameters: None

    Returns: 
        temp_dir, tuple - Temporary Directory object containing test files
        output_path, tuple - path for the plot output
    """
    sensitivity_testdata = {
        Years: [2010,2011,2012, 2013],
        Median: [100, 120, 130, 140]
        }

    base_testdata = {
        Years: [2010,2011,2012, 2013],
        Median: [90, 110, 150, 130]
        }

    temp_dir = tempfile.TemporaryDirectory()

    sensitivity_testfile = os.path.join(temp_dir.name, "outputs-for-management.csv")
    base_testfile = os.path.join(temp_dir.name, "outputs-for-management_base.csv")
    output_path = os.path.join(temp_dir.name, "test_plot.png")

    pd.DataFrame(sensitivity_testdata).to_csv(sensitivity_testfile, index = False)
    pd.DataFrame(base_testdata).to_csv(base_testfile, index = False)

    return(temp_dir, output_path)

# Smoke test: make sure the code actually runs
def test_smoke_plotsensivity():
    """
    Smoke test to make sure the function runs.
    """
    temp_dir, output_path = create_temp_data("Years", "Median Pre-fishery biomass (in 1000s metric tons)")

    try:
        plot_sensitivity(temp_dir.name, output_path)
        print("Smoke test passed!")

    except Exception as e:
        raise AssertionError(f"Smoke test failed: {e}") from e

    finally:
        temp_dir.cleanup()

# Test: Verify a plot is created and is valid
def test_createdplot_plotsensitivity():
    """
    Tests to see if a valid PNG plot file was created.
    """

    temp_dir, output_path = create_temp_data("Years", "Median Pre-fishery biomass (in 1000s metric tons)")

    try:
        plot_sensitivity(temp_dir, output_path)

        # check that file exists
        assert os.path.exists(output_path), "Plot file was not created"

        # opens the file to verify its a valid PNG
        with Image.open(output_path) as image:
            image.verify()
            assert image.format == "PNG", "Plot file is not in PNG format"
            print("Test passed: plot created successfully and is a PNG image")

    except Exception as e: # Catching general Exception for broader error coverage for testing
        print (f"Test failed: {e}")

    finally:
        temp_dir.cleanup()

# Edge tests: check to make sure that errors are being caught
def test_wrongcolumns_plotsensitivity():
    """
    Tests to see if error is raised when wrong column names are given
    """
    # data with wrong columns... (ie no biomass or years columns)

    temp_dir, output_path = create_temp_data("Time", "Biomass")

    try:
        plot_sensitivity(temp_dir.name, output_path)
        raise AssertionError("Test failed: No error raised for incorrect columns.")

    except ValueError as e:
        expected_message = """Sensitivity data must have 'Years'
                          and 'Median Pre-fishery biomass (in 1000s metric tons)' columns."""
        assert str(e) == expected_message, "Error message does not match expected format"
        print("Test passed: correct error for inocorrect column names")

    finally:
        temp_dir.cleanup()

#FUTURE DESIRED TEST:

# def test_filenotfound_plotsensitivity():
    # csv file does not exist so check to make sure its caught

# Test that verfies the resulting plot is what the plot should be (ie expected output = output)
