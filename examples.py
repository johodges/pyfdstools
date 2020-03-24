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
#=======================================================================
# # IMPORTS
#=======================================================================
import sys

sys.path.append('E:\\projects\\customPythonModules\\')

import fdsTools as fds
import os

def exampleImportFile():
    fdsPath = os.path.abspath("examples%s%s.fds"%(os.sep, "case001"))
    fdsFile = fds.fdsFileOperations()
    print("\tStarting to import file from %s"%(fdsPath))
    fdsFile.importFile(fdsPath)
    print("\tFinished importing file.")
    return fdsFile
    
def exampleSaveFile(file=None):
    if file is None:
        file = exampleImportFile()
    location = os.path.abspath("generated%s%s.fds"%(os.sep, file.head['ID']['CHID']))
    print("\tStarting to save model.")
    file.saveModel(1, location, allowMeshSplitting=False)
    print("\tFinished saving model.")
    return location

def exampleErrorCalculation():
    data, keys = fds.readErrorTable()
    
    errorvalues = fds.calculatePercentile([100, 200, 300, 400, 500, 600], 'Plume Temperature', 0.95, fdsVersion='6.7.1')
    fig, ax1 = fds.plotPercentile(500, 'Plume Temperature', fdsVersion='6.7.1')
    return errorvalues

if __name__ == '__main__':
    
    outDirName = "E:\\1MJP00014.002\\"
    inDirName = "E:\\1MJP00014.002\\"
    
    chid = "case001"
    orientations = [0]
    vName = 'WALL TEMPERATURE'
    defaultThreshold = 330
    meshes = [-1]
    tStart= 0
    finalTime = 3600
    intervalTime = 10
    bandTime = 10
    
    print("Importing model example")
    file = exampleImportFile()
    print("Saving model example")
    exampleSaveFile(file)
    print("Plotting error example")
    exampleErrorCalculation()
    

