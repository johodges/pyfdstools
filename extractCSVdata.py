#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This software is part of a python library to assist in developing and
# analyzing simulation results from Fire Dynamics Simulator (FDS).
# FDS is an open source software package developed by NIST. The source
# code is available at: https://github.com/firemodels/fds
#
# EXAMPLES:
# See the examples subroutine for example operation.
#
#=======================================================================
# # IMPORTS
#=======================================================================
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
