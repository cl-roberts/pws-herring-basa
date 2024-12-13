"""This is a test module for plot_sensitivity_comparison.py"""
import sensitivity

# Created a smoke test to see if the function runs

def test_smoke_psc():
    """This is an smoke test for the sens_comparison function"""
    # need to run sensitivity comparison to generate sensitivity_comparison.csv
    sensitivity.sens_comparison("../data_outputs/")
    sensitivity.plot_sens_comparison("../data_outputs", "../data_outputs/sensitivity_comparison_plot")

def test_string_1():
    dir_outputs = 1
    dir_plot = "../data_outputs/sensitivity_comparison_plot"
    try:
        sensitivity.plot_sens_comparison(dir_outputs, dir_plot)
    except (ValueError) as err:
        print("The path should be a string.")


def test_string_2():
    dir_outputs = "../data_outputs/"
    dir_plot = 1
    try:
        sensitivity.plot_sens_comparison(dir_outputs, dir_plot)
    except (ValueError) as err:
        print("The path should be a string.")
# In the future, should add edge tests for invalid data paths or issues with image generation

# For this test, run the pytest commmand line prompt in sensitivity\tests
