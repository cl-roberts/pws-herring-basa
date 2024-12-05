"""This is a test module for modify_mortality.py"""
import pandas as pd
import pytest

# add file path to modify_mortality module
import sys
import os
sys.path.insert(0, os.path.abspath('../src/sensitivity'))
import modify_mortality as mm

# edge test checking that having a mortality rate lower than 0 will raise an error
def test_lower_zero():
    new_value = -0.5
    dir_model = "../model"
    try:
        mm.modify_mortality(new_value, dir_model)
    except (ValueError) as err:
        print("Instantaneous natural mortality cannot be negative. Change new_value")

# smoke test
def test_smoke():
    new_value = 2
    dir_model = "../model"
    mm.modify_mortality(new_value, dir_model)

