"""This is a test module for modify_mortality.py"""
import pandas as pd
import pytest


import sensitivity

# edge test checking that having a mortality rate lower than 0 will raise an error
def test_lower_zero():
    new_value = -0.5
    dir_model = "../model"
    try:
        sensitivity.modify_mortality(new_value, dir_model)
    except (ValueError) as err:
        print("Instantaneous natural mortality cannot be negative. Change new_value")

# edge test checking that "dir_model" is a string
def test_string():
    new_value = 2
    dir_model = 3
    try:
        sensitivity.modify_mortality(new_value, dir_model)
    except (ValueError) as err:
        print("The path to BASA model directory should be a string.")


# smoke test
def test_smoke():
    new_value = 2
    dir_model = "../model"
    sensitivity.modify_mortality(new_value, dir_model)

