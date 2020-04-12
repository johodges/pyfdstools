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

sys.path.append('C:\\Users\\Jonathan\\Desktop\\')

import pyfdstools as fds
import os
from collections import defaultdict
import numpy as np

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

def exampleReadSlcf2dResults(resultDir=None, chid=None,
                             fdsQuantities = ['TEMPERATURE'],
                             tStart=0, tEnd=120):
    if chid is None: chid = "case001"
    if resultDir is None: resultDir = os.path.abspath("examples%s%s.zip"%(os.sep, chid))
    
    datas = defaultdict(bool)
    
    for qty in fdsQuantities:
        grid, data, times = fds.readSLCF2Ddata(chid, resultDir, qty)
        tStartInd = np.argwhere(times >= tStart)[0][0]
        tEndInd = np.argwhere(times <= tEnd)[-1][0]
        data_tavg = np.nanmean(data[:, :, :, tStartInd:tEndInd], axis=3)
        datas[qty] = data_tavg.copy()
    datas['GRID'] = grid
    datas['TIMES'] = times
    
    return datas

def exampleExtract2dFromSlcf2d(datas,
                               fdsQuantities = ['TEMPERATURE'],
                               fdsUnits = ['C'],
                               axis=1, value=2.5,
                               qnty_mn=20, qnty_mx=150):
    datas2D = defaultdict(bool)    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        x, z, data_slc = fds.findSliceLocation(datas['GRID'], datas[qty], axis, value)
        datas2D[qty] = data_slc
        fig = fds.plotSlice(x, z, data_slc, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    datas2D['X'] = x
    datas2D['Z'] = z
    
    return datas2D, fig

def exampleReadSlcf3dResults(resultDir=None, chid=None,
                             fdsQuantities = ['TEMPERATURE'],
                             tStart=0, tEnd=120):
    if chid is None: chid = "case001"
    if resultDir is None: resultDir = os.path.abspath("examples%s%s.zip"%(os.sep, chid))
    
    datas = defaultdict(bool)
    
    for qty in fdsQuantities:
        grid, data, times = fds.readSLCF3Ddata(chid, resultDir, qty)
        tStartInd = np.argwhere(times >= tStart)[0][0]
        tEndInd = np.argwhere(times <= tEnd)[-1][0]
        data_tavg = np.nanmean(data[:, :, :, tStartInd:tEndInd], axis=3)
        datas[qty] = data_tavg.copy()
    datas['GRID'] = grid
    datas['TIMES'] = times
    
    return datas

def exampleExtract2dFromSlcf3d(datas,
                               fdsQuantities = ['TEMPERATURE'],
                               fdsUnits = ['C'],
                               axis=2, value=4.4,
                               qnty_mn=20, qnty_mx=150):
    datas2D = defaultdict(bool)    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        x, z, data_slc = fds.findSliceLocation(datas['GRID'], datas[qty], axis, value)
        datas2D[qty] = data_slc
        fig = fds.plotSlice(x, z, data_slc, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    datas2D['X'] = x
    datas2D['Z'] = z
    
    return datas2D, fig

def exampleImportBndf(resultDir=None, chid=None, fdsFile=None,
                       fdsQuantities = ['WALL TEMPERATURE'],
                       fdsUnits = ['C'],
                       tStart=0, tEnd=120,
                       axis=-2, value=4.4,
                       qnty_mn=20, qnty_mx=100):
    if chid is None: chid = "case001"
    if resultDir is None: resultDir = os.path.abspath("examples%s%s.zip"%(os.sep, chid))
    if fdsFile is None: fdsFilePath = fds.getFileList(resultDir, chid, 'fds')[0]
    
    datas, times = fds.queryBndf(resultDir, chid, fdsFilePath, fdsQuantities, fdsUnits, axis, value)
    tStartInd = np.argwhere(times >= tStart)[0][0]
    tEndInd = np.argwhere(times <= tEnd)[-1][0]
    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        for mesh in list(datas[qty].keys()):
            data = datas[qty][mesh]['DATA']
            x = datas[qty][mesh]['X']
            z = datas[qty][mesh]['Z']
        meanData = np.mean(data[:, :, tStartInd:tEndInd], axis=2)
        fig = fds.plotSlice(x, z, meanData, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    return datas, fig

def exampleExtractBndfMax(resultDir=None, chid=None, fdsFile=None, smvFile=None,
                          outputName=None,
                          fdsQuantities = ['WALL TEMPERATURE'], fdsUnits = ['C'],
                          tStart=0, tEnd=120, tInt=1, tBand=3, orientations=[0],
                          axis=-2, value=4.4,
                          qnty_mn=20, qnty_mx=100):
    """This class defines operations that can be used to read, modify,
    and generate FDS input files.
    
    This is the further elaboration of the docstring. Within this section,
    you can elaborate further on details as appropriate for the situation.
    Notice that the summary and the elaboration is separated by a blank new
    line.
    """
    if chid is None: chid = "case001"
    if resultDir is None: resultDir = os.path.abspath("examples%s%s.zip"%(os.sep, chid))
    if fdsFile is None: fdsFilePath = fds.getFileList(resultDir, chid, 'fds')[0]
    if smvFile is None: smvFilePath = fds.getFileList(resultDir, chid, 'smv')[0]
    if outputName is None: outputName = os.path.abspath("generated%s%s_max_"%(os.sep, chid))
    
    datas = fds.extractMaxBndfValues(fdsFilePath, smvFilePath, resultDir, chid, fdsQuantities,
                    tStart=tStart, tEnd=tEnd, tInt=tInt, tBand=tBand, orientations=orientations)
    
    figs = []
    for qty in fdsQuantities:
        times = datas[qty]['TIMES']
        mPts = datas[qty]['DATA']
        names = datas[qty]['NAMES']
        outName = "%s%s"%(outputName, qty)
        fds.maxValueCSV(times, mPts, names, outName)
        fig = fds.maxValuePlot(times, mPts, names, outName, vName=qty, yticks=[20, 50, 100, 150, 200, 250, 300, 350])
        figs.append(fig)
    return datas, figs

if __name__ == '__main__':
    
    print("Importing model example")
    file = exampleImportFile()
    print("Saving model example")
    exampleSaveFile(file)
    print("Plotting error example")
    exampleErrorCalculation()
    print("Importing SLCF2D results example")
    datas = exampleReadSlcf2dResults()
    print("Extracting slice from SLCF2D results example")
    datas2D, fig = exampleExtract2dFromSlcf2d(datas)
    print("Importing SLCF3D results example")
    datas = exampleReadSlcf3dResults()
    print("Extracting 2-D slice from SLCF3D results example")
    datas2D, fig = exampleExtract2dFromSlcf3d(datas)
    print("Importing BNDF results example")
    datas, fig = exampleImportBndf()
    print("Extracting max value from BNDF results example")
    datas, figs = exampleExtractBndfMax()
