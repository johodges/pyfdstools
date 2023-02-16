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
import pyfdstools as fds
import os
from collections import defaultdict
import numpy as np


def exampleImportFile(fdsPath=None):
    if fdsPath is None:
        systemPath = os.path.dirname(os.path.abspath(__file__))
        fdsPath = os.path.join(systemPath, "examples", "case001.fds")
    fdsFile = fds.fdsFileOperations()
    print("\tStarting to import file from %s"%(fdsPath))
    fdsFile.importFile(fdsPath)
    print("\tFinished importing file.")
    return fdsFile


def exampleSaveFile(file=None, outdir=None, outname=None):
    if file is None:
        file = exampleImportFile()
    if outdir is None:
        systemPath = os.path.dirname(os.path.abspath(__file__))
        outdir = os.path.join(systemPath, "generated")
    if outname is None:
        location = os.path.join(outdir, '%s.fds'%(file.head['ID']['CHID']))
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
    if (resultDir is None) and (chid is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
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
    if (resultDir is None) and (chid is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
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
    if (resultDir is None) and (chid is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    if fdsFile is None: fdsFilePath = fds.getFileList(resultDir, chid, 'fds')[0]
    
    datas, times = fds.queryBndf(resultDir, chid, fdsFilePath, fdsQuantities, fdsUnits, axis, value)
    tStartInd = np.argwhere(times >= tStart)[0][0]
    tEndInd = np.argwhere(times <= tEnd)[-1][0]
    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        for mesh in list(datas[qty].keys()):
            data = datas[qty]['DATA']
            x = datas[qty]['X']
            z = datas[qty]['Z']
        meanData = np.mean(data[:, :, tStartInd:tEndInd], axis=2)
        fig = fds.plotSlice(x, z, meanData, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    return datas, fig

def exampleExtractBndfMax(resultDir=None, chid=None, outDir=None,
                          fdsFile=None, smvFile=None, outputNamespace=None,
                          fdsQuantities = ['WALL TEMPERATURE'], fdsUnits = ['C'],
                          tStart=0, tEnd=120, tInt=1, tBand=3, orientations=[0],
                          axis=-2, value=4.4,
                          qnty_mn=20, qnty_mx=100,
                          yticks=[20, 50, 100, 150, 200, 250, 300, 350]):
    if (resultDir is None) and (chid is None) and (outDir is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
        outDir = os.path.join(systemPath, "generated")
    
    if fdsFile is None: fdsFilePath = fds.getFileList(resultDir, chid, 'fds')[0]
    if smvFile is None: smvFilePath = fds.getFileList(resultDir, chid, 'smv')[0]
    if outputNamespace is None: outputNamespace = "%s_max_"%(chid)
    
    outputName = os.path.join(outDir, outputNamespace)
    
    datas = fds.extractMaxBndfValues(fdsFilePath, smvFilePath, resultDir, chid, fdsQuantities,
                    tStart=tStart, tEnd=tEnd, tInt=tInt, tBand=tBand, orientations=orientations)
    
    figs = []
    for qty in fdsQuantities:
        times = datas[qty]['TIMES']
        mPts = datas[qty]['DATA']
        names = datas[qty]['NAMES']
        outName = "%s%s"%(outputName, qty)
        fds.maxValueCSV(times, mPts, names, outName)
        fig = fds.maxValuePlot(times, mPts, names, outName, vName=qty, yticks=yticks)
        figs.append(fig)
    return datas, figs


def example2dSliceToCsv(resultDir=None, outDir=None, chid=None,
                        quantity=None, unit=None, axis=None, value=None,
                        time=None, dt=None):
    if (resultDir is None) and (chid is None) and (outDir is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
        outDir = os.path.join(systemPath, "generated")
    try:
        os.mkdir(outDir)
    except:
        pass
    data = fds.query2dAxisValue(resultDir, chid, quantity, axis, value, time=time, dt=dt)
    fds.renderSliceCsvs(data, chid, outDir)
    fig, ax = fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], axis,
                        clabel="%s (%s)"%(quantity, unit))
    fig.savefig(os.path.join(outDir, '%s_%s_%0.0f_%0.4f_final_frame.png'%(chid, quantity, axis, value)))
    return data, fig


def exampleBndfTimeAverage(resultDir=None, outDir=None, chid=None,
                           quantity=None, dt=None):
    if (resultDir is None) and (chid is None) and (outDir is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
        outDir = os.path.join(systemPath, "generated")
    try:
        os.mkdir(outDir)
    except:
        pass
    
    outFiles, outQty, refFiles, newSmvFile = fds.bndfsTimeAverage(
        resultDir, chid, quantity, dt, outDir=outDir)
    

def exampleWriteToNetCDF4():
    import netCDF4
    
    # Get case information from examples
    systemPath = os.path.dirname(os.path.abspath(__file__))
    chid = "case001"
    resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    qty = 'TEMPERATURE'
    
    # Read data
    grid, data, times = fds.readSLCF3Ddata(chid, resultDir, qty)
    
    # Convert axes to netcdf4 ordering
    data2 = np.moveaxis(data, [3,2,1,0], [0,1,2,3])
    grid2 = np.moveaxis(grid, [2,1,0,3], [0,1,2,3])
    
    # Open netcdf4 file
    rootgrp = netCDF4.Dataset('output.nc','w')
    
    # Establish basic dimension information for netcdf file
    time = rootgrp.createDimension('time', None)
    lon = rootgrp.createDimension('lon', grid.shape[0])
    lat = rootgrp.createDimension('lat', grid.shape[1])
    level = rootgrp.createDimension('level', grid.shape[2])
    
    # Establish netcdf variable
    temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon"))
    temp.units = "K"
    temp[:] = data2 + 273 # Fill temperature variable with values from FDS in K
    
    # Close file
    rootgrp.close()
    

def runExamples():
    systemPath = os.path.dirname(os.path.abspath(__file__))
    exampleInputFdsFile = os.path.join(systemPath, "examples", "case001.fds")
    exampleOutputDir = os.path.join(systemPath, "generated")
    chid = "case001"
    resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    print("Importing model example", flush=True)
    file = exampleImportFile(exampleInputFdsFile)
    
    print("Saving model example", flush=True)
    exampleSaveFile(file=file, outdir=exampleOutputDir)
    
    print("Plotting error example", flush=True)
    exampleErrorCalculation()
    
    print("Importing SLCF2D results example", flush=True)
    datas = exampleReadSlcf2dResults(resultDir=resultDir, chid=chid)
    
    print("Extracting slice from SLCF2D results example", flush=True)
    datas2D, fig = exampleExtract2dFromSlcf2d(datas)
    
    print("Importing SLCF3D results example", flush=True)
    datas = exampleReadSlcf3dResults(resultDir=resultDir, chid=chid)
    
    print("Extracting 2-D slice from SLCF3D results example", flush=True)
    datas2D, fig = exampleExtract2dFromSlcf3d(datas)
    
    print("Importing BNDF results example", flush=True)
    datas, fig = exampleImportBndf(resultDir=resultDir, chid=chid)
    
    print("Extracting max value from BNDF results example", flush=True)
    datas, figs = exampleExtractBndfMax(resultDir=resultDir, chid=chid, outDir=exampleOutputDir)
    
    print("Rendering 2d slice to csv example.", flush=True)
    datas = example2dSliceToCsv(resultDir=resultDir, chid=chid, outDir=exampleOutputDir,
                                axis=1, value=2.45, time=30, dt=60, quantity='TEMPERATURE', unit='C')
    
    print("Time-averaging a boundary file example.", flush=True)
    exampleBndfTimeAverage(dt=30, quantity='WALL TEMPERATURE')
    

if __name__ == '__main__':
    
    runExamples()