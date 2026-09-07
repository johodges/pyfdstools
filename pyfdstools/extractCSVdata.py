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
    """Finds the index of the last header line in an FDS csv file

    FDS csv output starts with one or more header lines (units, then
    labels) before the numeric data. The header ends at the last line
    which cannot be parsed as a row of floats.

    Parameters
    ----------
    lines : list
        List of the raw bytes lines of the file. Non-ASCII unit symbols
        FDS writes (degree signs, superscript twos) are replaced in
        place so that the caller decodes cleanly

    Returns
    -------
    int
        Index of the label line, or 0 if no data line was found within
        the first 100 lines
    """

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

def cleanDataLines(lines2, headerLines, skipcols=None):
    """Parses the numeric rows of an FDS csv file

    Parameters
    ----------
    lines2 : list
        List of the raw bytes lines of the file
    headerLines : int
        Index of the label line, from findHeaderLength
    skipcols : list, optional
        Indices of columns to drop

    Returns
    -------
    list
        List of rows, each a list of floats
    """

    lines = lines2[headerLines+1:]
    for i in range(0, len(lines)):
        line = (lines[i].decode('utf-8')).replace('\r\n','')
        while line[-1] == ',': line = line[:-1]
        if skipcols is None:
            lines[i] = [float(y) for y in line.split(',')]
        else:
            split = line.split(',')
            lines[i] = [float(split[i]) for i in range(0, len(split)) if i not in skipcols]
    return lines

def load_csv(modeldir, chid, suffix='_devc', labelRow=-1, skipcols=None):
    """Reads one of the csv files FDS writes for a case

    Parameters
    ----------
    modeldir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case
    suffix : str, optional
        Output suffix identifying which csv to read, for example
        '_devc', '_hrr' or '_ctrl' (default '_devc')
    labelRow : int, optional
        Index of the row holding the column labels. The header length is
        detected automatically when this is -1 (default -1)
    skipcols : list, optional
        Indices of columns to drop

    Returns
    -------
    pandas.DataFrame
        Table of the csv contents, with a column per device

    Raises
    ------
    FileNotFoundError
        If no csv with that suffix exists for the case
    """

    if 'zip' in modeldir:
        csv_files = getFileList(modeldir, chid, 'csv')
        suff_files = [x for x in csv_files if suffix in x]
        if len(suff_files) == 0:
            raise FileNotFoundError(
                "No csv file matching '%s' for chid %s was found in %s."
                % (suffix, chid, modeldir))
        f = zopen(suff_files[0])
    else:
        file = "%s%s%s%s.csv"%(modeldir, os.sep, chid, suffix)
        f = open(file, 'rb')
    lines = f.readlines()
    f.close()
    if labelRow == -1:
        headerLines = findHeaderLength(lines)
        header = (lines[headerLines].decode('utf-8')).replace('\r\n','').replace('\n','').split(',')
    else:
        headerLines = labelRow
        header = (lines[labelRow].decode('utf-8')).replace('\r\n','').replace('\n','').split(',')
    if skipcols is not None:
        header1 = [header[i] for i in range(0, len(header)) if i not in skipcols]
        header = header1
    dataLines = cleanDataLines(lines, headerLines, skipcols)
    data = pd.DataFrame(dataLines, columns=header)
    return data

def load_hrr(file):
    """Reads an FDS _hrr.csv file directly by path

    Parameters
    ----------
    file : str
        Path to the csv file

    Returns
    -------
    pandas.DataFrame
        Table of the csv contents

    See Also
    --------
    load_csv : reads the same file given a directory and CHID, and works
        for archives and for every other csv suffix
    """

    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.replace('\n','').split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data
