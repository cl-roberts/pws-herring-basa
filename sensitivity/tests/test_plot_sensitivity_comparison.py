"""This is a test module for plot_sensitivity_comparison.py"""
import sys
sys.path.append('../src/sensitivity')
import plot_sensitivity_comparison as psc

# Created a smoke test to see if the function runs

def test_smoke_psc():
    """This is an smoke test for the sens_comparison function"""
    psc.plot_sens_comparison("../data_outputs/", "../data_outputs/")

# In the future, should add edge tests for invalid data paths or issues with image generation

# For this test, run the pytest commmand line prompt in sensitivity\tests
