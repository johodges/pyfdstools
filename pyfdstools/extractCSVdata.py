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

def findHeaderLength(lines):
    replacements =[
        [b'kW/m\xb2',b'kW/m2'],
        [b'\xb0C',b'C']
        ]
    counter = 0
    headerCheck = True
    while headerCheck and counter < 100:
        linetmp = lines[counter]
        for replacement in replacements:
            linetmp = linetmp.replace(replacement[0],replacement[1])
        lines[counter] = linetmp
        line = (linetmp.decode('utf-8',errors='ignore')).replace('\r\n','')
        while line[-1] == ',': line = line[:-1]
        try:
            [float(y) for y in line.split(',')]
            counter = counter - 1
            headerCheck = False
        except:
            counter = counter + 1
    if counter < 100:
        return counter
    else:
        print("Unable to find header length, returning 0")
        return 0

def cleanDataLines(lines2, headerLines):
    lines = lines2[headerLines+1:]
    for i in range(0, len(lines)):
        line = (lines[i].decode('utf-8')).replace('\r\n','')
        while line[-1] == ',': line = line[:-1]
        lines[i] = [float(y) for y in line.split(',')]
    return lines

def load_csv(modeldir, chid, suffix='_devc', labelRow=-1):
    if 'zip' in modeldir:
        csv_files = getFileList(modeldir, chid, 'csv')
        suff_files = [x for x in csv_files if suffix in x]
        print(modeldir, chid, suff_files)
        f = zopen(suff_files[0])
    else:
        file = "%s%s%s%s.csv"%(modeldir, os.sep, chid, suffix)
        f = open(file, 'rb')
    lines = f.readlines()
    f.close()
    headerLines = findHeaderLength(lines)
    if labelRow == -1:
        header = (lines[headerLines].decode('utf-8')).replace('\r\n','').replace('\n','').split(',')
    else:
        header = (lines[labelRow].decode('utf-8')).replace('\r\n','').replace('\n','').split(',')
    dataLines = cleanDataLines(lines, headerLines)
    data = pd.DataFrame(dataLines, columns=header,)
    return data

def load_hrr(file):
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.replace('\n','').split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data
