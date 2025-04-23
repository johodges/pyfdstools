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
import pandas as pd

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
    
    datas, times = fds.queryBndf(resultDir, chid, fdsFilePath, fdsQuantities, axis, value)
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
    
def stretchedMeshExample(resultDir=None, outDir=None, chid=None,
                           quantity="TEMPERATURE", dt=None, time=None):
    if (resultDir is None) and (chid is None) and (outDir is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "stretched_mesh_example"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
        outDir = os.path.join(systemPath, "generated")
    try:
        os.mkdir(outDir)
    except:
        pass
    
    grid, data, times = fds.readSLCF3Ddata(chid, resultDir, quantity, time=time, dt=dt)
    
    x, z, data_slc = fds.findSliceLocation(grid, data[:, :, :, -1], 1, 0)
    fds.plotSlice(x, z, data_slc, 1, figsize=(10, 4))
    
    x, z, data_slc = fds.findSliceLocation(grid, data[:, :, :, -1], 3, 885)
    fds.plotSlice(x, z, data_slc, 3, figsize=(10, 10))


def exampleAddOccupantFedDevices(indir,outdir):
    name = 3
    height_above_floor = 1.8
    occupants = pd.read_csv(os.path.join(indir,'fed_example_occupants.csv'), header=[0], skiprows=[1])
    t = occupants.loc[occupants['name'] == name,'t'].values
    x = occupants.loc[occupants['name'] == name,'x'].values
    y = occupants.loc[occupants['name'] == name,'y'].values
    z = occupants.loc[occupants['name'] == name,'z'].values + height_above_floor
    
    fdsFile = fds.fdsFileOperations()
    fdsFile.importFile(os.path.join(indir, 'fed_example.fds'))
    fdsFile.addOccupant('Occupant %d'%(name), x, y, z, t, output=True)
    fdsFile.saveModel(1, os.path.join(outdir, "fed_example_out.fds"), allowMeshSplitting=False)
    

def exampleParseS3dFiles(resultDir=None, chid=None):
    if (resultDir is None) and (chid is None):
        systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    values, times = fds.extractS3dValues(resultDir, chid)
    s3dfiles = fds.getFileList(resultDir, chid, 's3d')
    smvFile = fds.getFileList(resultDir, chid, 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    (grids, obsts) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    
    # Render soot file
    meshNum = 0
    quantity = 'SOOT DENSITY'
    data_out = values[meshNum][quantity]
    grid = grids[meshNum]
    dx = np.median(grid[0][1:,1]-grid[0][:-1,1])
    dy = np.median(grid[1][1:,1]-grid[1][:-1,1])
    dz = np.median(grid[2][1:,1]-grid[2][:-1,1])
    dd = (dx*dy*dz)**(1/3)
    file = 'fakeout.s3d'
    print(data_out.shape)
    fds.writeS3dFile(file, times, data_out, quantity, dd)
    
    # To check
    s3dfile = s3dfiles[0]
    f = fds.zopen(s3dfile, 'rb')
    data = f.read()
    f.close()
    
    f1 = fds.zopen(file, 'rb')
    data2 = f1.read()
    f1.close()
    
    print("File %s, "%(s3dfile), data==data2)
    
    
    # Render hrrpuv file
    meshNum = 0
    quantity = 'HRRPUV'
    data_out = values[meshNum][quantity]
    grid = grids[meshNum]
    dx = np.median(grid[0][1:,1]-grid[0][:-1,1])
    dy = np.median(grid[1][1:,1]-grid[1][:-1,1])
    dz = np.median(grid[2][1:,1]-grid[2][:-1,1])
    dd = (dx*dy*dz)**(1/3)
    file = 'fakeout2.s3d'
    
    
    fds.writeS3dFile(file, times, data_out, quantity, dd)

    # To check
    s3dfile = s3dfiles[1]
    f = fds.zopen(s3dfile, 'rb')
    data = f.read()
    f.close()
    
    f1 = fds.zopen(file, 'rb')
    data2 = f1.read()
    f1.close()
    
    print("File %s, "%(s3dfile), data==data2)

    
    # Render temperature file
    meshNum = 0
    quantity = 'TEMPERATURE'
    data_out = values[meshNum][quantity]
    grid = grids[meshNum]
    dx = np.median(grid[0][1:,1]-grid[0][:-1,1])
    dy = np.median(grid[1][1:,1]-grid[1][:-1,1])
    dz = np.median(grid[2][1:,1]-grid[2][:-1,1])
    dd = (dx*dy*dz)**(1/3)
    file = 'fakeout3.s3d'
    fds.writeS3dFile(file, times, data_out, quantity, dd)
    
    # To check
    s3dfile = s3dfiles[2]
    f = fds.zopen(s3dfile, 'rb')
    data = f.read()
    f.close()
    
    f1 = fds.zopen(file, 'rb')
    data2 = f1.read()
    f1.close()
    
    print("File %s, "%(s3dfile), data==data2)

def examplePostProcessVisibility(resultDir=None, chid=None, outDir=None, oldC=3, newC=8):
    if resultDir == None:
        resultDir = 'examples\\visibility_adjustment.zip'
    if chid == None:
        chid = 'visibility_adjustment'
    if outDir == None:
        outDir = 'generated\\'
    qty = 'SOOT VISIBILITY'
    sName = 'rho_C'
    uts = 'kg/m3'
    
    slcfFiles = fds.getFileList(resultDir, chid, 'sf')
    slcfFiles = [x for x in slcfFiles if '_custom.sf' not in x]
    
    endianness = fds.getEndianness(resultDir, chid)
    quantities = []
    for file in slcfFiles:
        f = fds.zopen(file)
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = fds.readSLCFheader(f, endianness)
        quantities.append(qty)
    
    slcfFiles = [x for x,y in zip(slcfFiles, quantities) if qty in y]
    quantities = [y for x,y in zip(slcfFiles, quantities) if qty in y]
    
    smvFile = fds.getFileList(resultDir, chid, 'smv')[0]
    f = fds.zopen(smvFile, 'r')
    smvTxt = f.read()
    f.close()
    if type(smvTxt) is bytes:
        smvTxt = smvTxt.decode('utf-8')
    if '.zip' in smvFile:
        smvFileCustom = smvFile.replace('.smv','_custom.smv').split('.zip')[-1]
        smvFileCustom = os.path.basename(smvFileCustom)
    else:
        smvFileCustom = smvFile.replace('.smv','_custom.smv')
        smvFileCustom = os.path.basename(smvFileCustom)
    smvPathCustom = os.path.join(outDir, smvFileCustom)
    with open(smvPathCustom, 'w') as f:
        f.write(smvTxt)
    
    for slcfFile in slcfFiles:
        lims, data, times = fds.readSingleSlcfFile(slcfFile)
        data = np.squeeze(data)/oldC*newC
        data2 = [np.array(data[:,:,i], dtype=np.float32) for i in range(0, data.shape[2])]
        axis, val = fds.getAxisFromLims(lims)
        if axis < 0:
            print("Warning %s is a 3-D slice, skipping."%(slcfFile))
            continue
        meshnum = int(slcfFile.split('_')[-2])
        if '.zip' in slcfFile:
            outFile = slcfFile.replace('.sf','_custom.sf').split('.zip')[-1]
        else:
            outFile = slcfFile.replace('.sf','_custom.sf')
        outFile = os.path.basename(outFile)
        fds.writeSlice(outFile, outDir, chid, data2, times, axis, val,
                              qty+'C%d'%(newC), sName, uts, meshnum, smvFile=smvFileCustom)

def exampleAddGasPhaseHeatFluxSlice(resultDir=None, chid=None, outDir=None):
    import scipy.interpolate
    if resultDir == None:
        resultDir = 'examples\\hfg_slice.zip'
    if chid == None:
        chid = 'hfg_slice'
    if outDir == None:
        outDir = 'generated\\'
    
    smvFile = fds.getFileList(resultDir, chid, 'smv')[0]
    f = fds.zopen(smvFile, 'r')
    smvTxt = f.read()
    f.close()
    if type(smvTxt) is bytes:
        smvTxt = smvTxt.decode('utf-8')
    if '.zip' in smvFile:
        smvFileCustom = smvFile.replace('.smv','_custom.smv').split('.zip')[-1]
        smvFileCustom = os.path.basename(smvFileCustom)
    else:
        smvFileCustom = smvFile.replace('.smv','_custom.smv')
        smvFileCustom = os.path.basename(smvFileCustom)
    smvPathCustom = os.path.join(outDir, smvFileCustom)
    with open(smvPathCustom, 'w') as f:
        f.write(smvTxt)
    
    ref_slcf_qty = 'TEMPERATURE'
    ref_slcf_axis = 2
    devc_namespace = 'rhf01'
    devc_x0 = -0.95
    devc_dx = 0.1
    devc_y0 = -0.75
    devc_dy = 0.0
    
    quantities, slcfFiles, dimensions, meshes, centers = fds.readSLCFquantities(chid, resultDir)
    
    slcfFiles = [x for x,y in zip(slcfFiles, quantities) if ref_slcf_qty == y]
    
    d = fds.load_csv(resultDir, chid, '_devc')
    c = [x for x in d.columns if devc_namespace in x]
    d_x = [devc_x0+devc_dx*(float(x.split('-')[-1].replace('"',''))-1) for x in c]
    d_y = [devc_y0+devc_dy*(float(x.split('-')[-1].replace('"',''))-1) for x in c]
    d_z = [(float(x.split('-')[-2].replace('p','.'))) for x in c]
    
    d_x = np.round(d_x, decimals=2)
    d_y = np.round(d_y, decimals=2)
    d_z = np.round(d_z, decimals=2)
    
    x_unique = np.unique(d_x)
    y_unique = np.unique(d_y)
    z_unique = np.unique(d_z)
    
    x_grid, y_grid, z_grid = np.meshgrid(x_unique, y_unique, z_unique)
    d_grid = np.zeros((y_grid.shape[0], y_grid.shape[1], y_grid.shape[2], d.values.shape[0]))
    
    for i in range(0, len(d_x)):
        x_ind = np.argwhere(d_x[i] == x_unique)
        y_ind = np.argwhere(d_y[i] == y_unique)
        z_ind = np.argwhere(d_z[i] == z_unique)
        d_grid[y_ind, x_ind, z_ind, :] = d[c[i]].values
    
    d_grid = np.swapaxes(d_grid, 0, 1)
    x_grid = np.swapaxes(x_grid, 0, 1)
    z_grid = np.swapaxes(z_grid, 0, 1)
    
    if ref_slcf_axis == 1:
        xx = y_grid[0,:,:]
        zz = z_grid[0,:,:]
        dd = d_grid[0, :, :, :]
        xxx = y_unique
        zzz = z_unique
    elif ref_slcf_axis == 2:
        xx = x_grid[:,0,:]
        zz = z_grid[:,0,:]
        dd = d_grid[:, 0, :, :]
        xxx = x_unique
        zzz = z_unique
    elif ref_slcf_axis == 3:
        xx = x_grid[:,:,0]
        zz = y_grid[:,:,0]
        dd = d_grid[:, :, 0, :]
        xxx = x_unique
        zzz = y_unique
    #fds.plotSlice(xx, zz, dd[:, :, -1], ref_slcf_axis)
    
    d_time = d['Time'].values
    for slcfFile in slcfFiles:
        lims, _, _ = fds.readSingleSlcfFile(slcfFile)
        x, z, d, times, coords = fds.read2dSliceFile(slcfFile, chid)
        data = np.zeros_like(d)
        for i in range(0, len(d_time)):
            r = scipy.interpolate.RegularGridInterpolator((xxx, zzz), dd[:, :, i], method='linear', bounds_error=False, fill_value=0.0)
            data[:, :, i] = r((x, z))
        
        #fds.plotSlice(x, z, data[:, :, -1], ref_slcf_axis)
        data = np.squeeze(data)
        data2 = [np.array(data[:,:,ii], dtype=np.float32) for ii in range(0, data.shape[2])]
        axis, val = fds.getAxisFromLims(lims)
        meshnum = int(slcfFile.split('_')[-2])
        if '.zip' in slcfFile:
            outFile = slcfFile.replace('.sf','_custom.sf').split('.zip')[-1]
        else:
            outFile = slcfFile.replace('.sf','_custom.sf')
        outFile = os.path.basename(outFile)
        fds.writeSlice(outFile, outDir, chid, data2, times, axis, val,
                              'RADIATIVE HEAT FLUX GAS', 'RHF', 'kW/m2', meshnum, smvFile=smvFileCustom)
        

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
                                axis=1, value=2.55, time=30, dt=60, quantity='TEMPERATURE', unit='C')
    
    print("Time-averaging a boundary file example.", flush=True)
    exampleBndfTimeAverage(dt=30, quantity='WALL TEMPERATURE')
    
    #print("Read SL3D with stretched mesh example", flush=True)
    #stretchedMeshExample()
    
    print("Add Occupant FED calculation example", flush=True)
    exampleAddOccupantFedDevices(os.path.join(systemPath, "examples"),
                                 os.path.join(systemPath, "generated"))
    
    print("Example Post Process Visibility to Change C Factor")
    examplePostProcessVisibility()
    
    #print("Example to parse Smoke 3D files")
    #exampleParseS3dFiles()
    
    print("Example add heat flux slice")
    exampleAddGasPhaseHeatFluxSlice()


if __name__ == '__main__':
    runExamples()
    