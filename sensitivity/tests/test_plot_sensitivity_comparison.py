"""This is a test module for sensitivity_comparison.py"""
import pandas as pd
import pytest
from sensitivity_comparison import sc

# Change this to a smoke test -- just to see that the function runs
# Randomness means that a one-shot test won't work
def test_oneshot():
    round(sc.sens_comparison("../../data_outputs/outputs-for-management.csv").iloc[0,0], 10) == 13.0621471667