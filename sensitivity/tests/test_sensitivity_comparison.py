"""This is a test module for sensitivity_comparison.py"""
# import module
import sensitivity

# Created a smoke test to see if the function runs
# Due to randomness, a one-shot test won't work

def test_smoke_sc():
    """This is an smoke test for the sens_comparison function"""
    sensitivity.sens_comparison("../data_outputs/")

# In the future, should add edge tests for invalid data paths or issues with the provided csv files

# For this test, run the pytest commmand line prompt in sensitivity\tests
