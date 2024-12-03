"""This is a test module for sensitivity_comparison.py"""
import pandas as pd
import pytest
from sensitivity_comparison import sc

def test_oneshot():
    round(sc.sens_comparison("../../data_outputs/outputs-for-management.csv").iloc[0,0], 10) == 13.0621471667