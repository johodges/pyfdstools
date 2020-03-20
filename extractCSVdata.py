# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:24:57 2020

@author: jhodges
"""

import pandas as pd

def load_devc(file):
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data

def load_hrr(file):
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data
