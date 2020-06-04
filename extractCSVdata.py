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
import os
from .utilities import zopen, getFileList

def load_csv(modeldir, chid, suffix='_devc'):
    if 'zip' in modeldir:
        csv_files = getFileList(modeldir, chid, 'csv')
        suff_files = [x for x in csv_files if suffix in x]
        print(modeldir, chid, suff_files)
        f = zopen(suff_files[0])
        line = f.readline()
        line = f.readline()
        lines = f.readlines()
        header = (line.decode('utf-8')).replace('\n','').split(',')
        lines = [[float(y) for y in (x.decode('utf-8')).split(',')] for x in lines]
        data = pd.DataFrame(lines, columns=header,)
        f.close()
    else:
        file = "%s%s%s%s.csv"%(modeldir, os.sep, chid, suffix)
        with open(file, 'r') as f:
            line = f.readline()
            line = f.readline()
        header = line.replace('\n','').split(',')
        data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data

def load_hrr(file):
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.replace('\n','').split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data
