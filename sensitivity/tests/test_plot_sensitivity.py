from PIL import Image

# Smoke test: make sure the code actually runs
def test_smoke_plotsensivity():
    plot_sensitivity()

def test_plotcreated_plotsensitivity():
    # verify a plot is created
    
# Edge tests: check to make sure that errors are being caught
def test_filenotfound_plotsensitivity():
    # csv file does note exist so check to make sure its caught

def test_wrongcolumns_plotsensitivity():
    # data with wrong columns... (ie no biomass or years)