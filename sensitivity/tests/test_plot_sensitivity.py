import os
from PIL import Image

# Smoke test: make sure the code actually runs
def test_smoke_plotsensivity():
    plot_sensitivity()

def test_plotcreated_plotsensitivity():
    # verify a plot is created
    
    # some code here (test inputs or something)

    try:
        plot_sensitivity(directory, output_path)

        # check that file exists
        assert os.path.exists(output_path), "Plot file was no created"

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

# Edge tests: check to make sure that errors are being caught
def test_filenotfound_plotsensitivity():
    # csv file does note exist so check to make sure its caught

def test_wrongcolumns_plotsensitivity():
    # data with wrong columns... (ie no biomass or years)