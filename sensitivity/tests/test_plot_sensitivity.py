import os
import tempfile
from PIL import Image
import pandas as pd
from plot_sensitivity import plot_sensitivity

# Test Data:
sensitivity_testdata = {
    "Years": [2010,2011,2012, 2013],
    "Median Pre-fishery biomass (in 1000s metric tons)": [100, 120, 130, 140]
    }

base_testdata = {
    "Years": [2010,2011,2012, 2013],
    "Median Pre-fishery biomass (in 1000s metric tons)": [90, 110, 150, 130]
    }

# Smoke test: make sure the code actually runs
def test_smoke_plotsensivity():
    with tempfile.TemporaryDirectory() as tempdir:
        sensitivity_testfile = os.path.join(tempdir, "outputs-for-management.csv")
        base_testfile = os.path.join(temp_dir, "outputs-for-management-base.csv")

        pd.DataFrame(sensitivity_testdata, index = False)
        pd.DataFrame(base_testdata, index = False)

        output_path = os.path.join(tempdir, "test_plot.png")

        try:
            plot_sensitivity(tempdir, output_path)
            print("Smoke test passed!")
        except Exception as e:
            raise AssertionError(f"Smoke test failed: {e}")
        
# Test: Verify a plot is created and is valid
def test_plotcreated_plotsensitivity():
    # verify a plot is created
    
    # some code here (test inputs or something)

    try:
        plot_sensitivity(input_directory, output_path)

        # check that file exists
        assert os.path.exists(output_path), "Plot file was not created"

        # opens the file to verify its a valid PNG
        with Image.open(output_path) as image:
            image.verify()
            assert image.format == "PNG", "Plot file is not in PNG format"
            print("Test passed: plot created successfully and is a PNG image")

    except Exception as e:
        print (f"Test failed: {e}")

    finally: # cleanup and remove test output
        if os.path.exists(output_path):
            os.remove(output_path)

#FUTURE DESIRED TEST: 

# Edge tests: check to make sure that errors are being caught
# def test_filenotfound_plotsensitivity():
    # csv file does not exist so check to make sure its caught

# def test_wrongcolumns_plotsensitivity():
    # data with wrong columns... (ie no biomass or years columns)

# Test that verfies the resulting plot is what the plot should be (ie expected output = output)
