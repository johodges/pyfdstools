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
import numpy as np
import matplotlib.pyplot as plt
import glob
import zipfile
import os
import struct
import scipy.interpolate as scpi
import pandas as pd
from collections import defaultdict
from .fdsFileOperations import fdsFileOperations
from .utilities import getDatatypeByEndianness, getEndianness
from .utilities import getFileListFromZip, getFileList, zopen, zreadlines
from .utilities import getFileListFromResultDir
from .utilities import getGridsFromXyzFiles, getAbsoluteGrid, rearrangeGrid
from .utilities import readXYZfile
from .colorSchemes import buildSMVcolormap
from .smokeviewParser import parseSMVFile

def time2str(time, decimals=2):
    """Converts a timestamp to a string
    
    This subroutine converts a float timestamp to a string with the
    period replaced by an underscore.
    
    Parameters
    ----------
    time : float
        Timestamp
    decimals : int, optional
        Number of decimal places to include (default 2)
    
    Returns
    -------
    str
        Timestamp with decimal replaced by underscore
        number of patches
    """
    decStr = '{:.{prec}f}'.format(time - np.floor(time), prec=decimals)
    decStr = decStr.replace('0.','_')
    timeStr = '%08.0f%s'%(np.floor(time), decStr)
    return timeStr

def mesh2str(mesh):
    """Converts a mesh number to string
    
    This subroutine converts an integer mesh number to a string.
    FDS specifies mesh numbers as zero-padded integers with 4 digits.
    
    Parameters
    ----------
    mesh : int
        Integer mesh number
    
    Returns
    -------
    str
        Mesh number as a string zero-padded to 4 digits
    """
    
    meshStr = "%04.0f"%(mesh)
    return meshStr


def readP3Dfile(file):
    """Reads data from plot3D file
    
    This subroutine reads data from a plot3D file.
    TODO: Update this subroutine to use zopen
    
    Parameters
    ----------
    file : str
        String containing the path to a plot3D file
    
    Returns
    -------
    array(NX, NY, NZ, NT)
        Array containing float data in local coordinates for each time
    array()
        Array containing header information from plot3D file
    """
    
    with open(file,'rb') as f:
        header = np.fromfile(f, dtype=np.int32, count=5)
        _ = np.fromfile(f, dtype=np.float32, count=7)
        (nx, ny, nz) = (header[1], header[2], header[3])
        data = np.fromfile(f, dtype=np.float32, count=nx*ny*nz*5)
        data = np.reshape(data, (int(data.shape[0]/5),5), order='F')
    return data, header[1:-1]

def writeP3Dfile(file, data):
    """Writes data to plot3D file
    
    This subroutine writes data to a plot3D file.
    
    Parameters
    ----------
    file : str
        String containing the path to a plot3D file
    
    Returns
    -------
    array(NX, NY, NZ, NT)
        Array containing float data in local coordinates for each time
    array()
        Array containing header information from plot3D file
    """
    
    with open(file,'wb') as f:
        f.write(b'\x0c\x00\x00\x00')
        nx, ny, nz, v = data.shape
        (nx, ny, nz, v) = (int(nx), int(ny), int(nz), int(v))
        f.write(nx.to_bytes(4, 'little'))
        f.write(ny.to_bytes(4, 'little'))
        f.write(nz.to_bytes(4, 'little'))
        f.write(b'\x0c\x00\x00\x00')
        empty = np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.float32)
        empty.tofile(f)
        d1 = data.flatten(order='F')
        d = d1.tobytes()
        f.write(d)


def buildDataFile(meshStr, time):
    """Builds a plot3D datafile name
    
    Parameters
    ----------
    meshStr : mesh number string
    time : float timestamp
    
    Returns
    -------
    str
        String containing the name of the plot3D data file
    """
    
    dataStr = meshStr.replace('.xyz','_%s.q'%(time2str(time)))
    return dataStr
    
def printExtents(grid, data):
    """Prints the extents of plot3D grid and data
    
    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Array containing float global coordinates
    data : array(NX, NY, NZ, 5, NT)
        Array containing float data in local coordinates for each time
    """
    
    print("%0.4f < x < %0.4f"%(np.min(grid[:,0]),np.max(grid[:,0])))
    print("%0.4f < y < %0.4f"%(np.min(grid[:,1]),np.max(grid[:,1])))
    print("%0.4f < z < %0.4f"%(np.min(grid[:,2]),np.max(grid[:,2])))
    print("%0.4f < IBLK < %0.4f"%(np.min(grid[:,3]),np.max(grid[:,3])))
    
    print("")
    
    print("%0.4f < T < %0.4f"%(np.min(data[:,0]),np.max(data[:,0])))
    print("%0.4f < U < %0.4f"%(np.min(data[:,1]),np.max(data[:,1])))
    print("%0.4f < V < %0.4f"%(np.min(data[:,2]),np.max(data[:,2])))
    print("%0.4f < W < %0.4f"%(np.min(data[:,3]),np.max(data[:,3])))
    print("%0.4f < HRRPUV < %0.4f"%(
            np.min(data[:,4]),np.max(data[:,4])))
        

def findSliceLocation(grid, data, axis, value, plot3d=False):
    """Extracts slice location from absolute grid
    
    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Array containing grid coordinates
    data : array(NX, NY, NZ)
        Array containing float data for each grid coordinate
    axis : int
        Axis number to query
    value : float
        Axis value to query
    plot3d : bool
        Flag specifying whether data is plot3d or not
    
    Returns
    -------
    array(NX, NZ)
        Array containing the 2-D slice x coordinates
    array(NX, NZ)
        Array containing the 2-D slice z coordinates
    array(NX, NZ)
        Array containing the 2-D slice data
    """
    
    xGrid = grid[:, :, :, 0]
    yGrid = grid[:, :, :, 1]
    zGrid = grid[:, :, :, 2]
    if axis == 1: 
        nGrid = xGrid
        nValues = nGrid[:, 0, 0]
    if axis == 2: 
        nGrid = yGrid
        nValues = nGrid[0, :, 0]
    if axis == 3: 
        nGrid = zGrid
        nValues = nGrid[0, 0, :]
    ind = np.argmin(abs(nValues - value))
    if axis == 1:
        x = np.squeeze(yGrid[ind, :, :])
        z = np.squeeze(zGrid[ind, :, :])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[ind, :, :, :])
        else:
            data2 = np.squeeze(data[ind, :, :])
    elif axis == 2:
        x = np.squeeze(xGrid[:, ind, :])
        z = np.squeeze(zGrid[:, ind, :])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:, ind, :, :])
        else:
            data2 = np.squeeze(data[:, ind, :])
    elif axis == 3:
        x = np.squeeze(xGrid[:, :, ind])
        z = np.squeeze(yGrid[:, :, ind])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:, :, ind, :])
        else:
            data2 = np.squeeze(data[:, :, ind])
    if plot3d:
        T = np.squeeze(data2[:, :, 0])
        U = np.squeeze(data2[:, :, 1])
        V = np.squeeze(data2[:, :, 2])
        W = np.squeeze(data2[:, :, 3])
        HRR = np.squeeze(data2[:, :, 4])
        return x, z, T, U, V, W, HRR
    else:
        return x, z, data2

def plotSlice(x, z, data_slc, axis, fig=None, ax=None,
              cmap=None, figsize=None, fs=16, figsizeMult=4,
              qnty_mn=None, qnty_mx=None, cbarnumticks=10,
              tickDecimals=None,
              levels=None, cbarticks=None, clabel=None,
              highlightValue=None, highlightWidth=None,
              reverseXY=False, percentile=None,
              xmn=None, xmx=None, zmn=None, zmx=None,
              xlabel=None, zlabel=None,
              addCbar=True, fixXLims=True, fixZLims=True,
              title=None, contour=True, extend='both',
              returnIm=False):
    if qnty_mn == None:
        qnty_mn = np.nanmin(data_slc)
    if qnty_mx == None:
        qnty_mx = np.nanmax(data_slc)
    if qnty_mn == qnty_mx:
        qnty_mn = 0
        qnty_mx = 1
    if highlightValue is not None:
        percentile = (highlightValue - qnty_mn) / (qnty_mx - qnty_mn)
    if (xmn == None): xmn = x.min()
    if (xmx == None): xmx = x.max()
    if (zmn == None): zmn = z.min()
    if (zmx == None): zmx = z.max()
    if cmap == None:
        cmap = buildSMVcolormap(
                percentile=percentile, width=highlightWidth)
    if cmap == 'SMV':
        cmap = buildSMVcolormap(
                percentile=percentile, width=highlightWidth)
    if levels == None:
        levels = 100
    if cbarticks is None:
        cbarticks = np.linspace(qnty_mn, qnty_mx, cbarnumticks)
        if highlightValue is not None:
            cx = cbarticks[1] - cbarticks[0]
            cbarticks = list(cbarticks)
            cbarticks = [x for x in cbarticks if abs(highlightValue-x) > 0.25*cx]
            cbarticks.append(highlightValue)
            cbarticks = sorted(cbarticks)
            
            
    if reverseXY:
        (x1 , z1, d1) = (z, x, data_slc[:, :])
        zrange = xmx-xmn #x.max()-x.min()
        xrange = zmx-zmn #z.max()-z.min()
    else:
        (x1, z1, d1) = (x, z, data_slc[:, :])
        xrange = xmx-xmn #x.max()-x.min()
        zrange = zmx-zmn #z.max()-z.min()
    if figsize == None:
        if zrange > xrange:
            figsize = (figsizeMult, figsizeMult * zrange / xrange)
        else:
            figsize = (figsizeMult * xrange / zrange, figsizeMult)
    if fig == None or ax == None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)
    if contour:
        im = ax.contourf(
            x1, z1, d1, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, 
            levels=np.linspace(qnty_mn, qnty_mx, levels), extend=extend)
    else:
        im = ax.imshow(d1[::-1, :], cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, extent=[x1.min(), x1.max(),z1.min(), z1.max()])
    if addCbar:
        if tickDecimals is not None:
            fmt = "%"+".%d"%(tickDecimals)+"f"
            #fmt = FormatScalarFormatter("%"+fmtStr+"f")
        else:
            fmt = None
        cbar = fig.colorbar(im, cmap=cmap, extend=extend, ticks=cbarticks, format=fmt)
        cbar.ax.tick_params(labelsize=fs)
        if clabel is not None:
            cbar.set_label(clabel, fontsize=fs)
    if xlabel is None:
        if abs(axis) == 1: xlabel = 'y (m)'
        if abs(axis) == 2: xlabel = 'x (m)'
        if abs(axis) == 3: xlabel = 'x (m)'
    if zlabel is None:
        if abs(axis) == 1: zlabel = 'z (m)'
        if abs(axis) == 2: zlabel = 'z (m)'
        if abs(axis) == 3: zlabel = 'y (m)'
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(zlabel, fontsize=fs)
    if fixXLims:
        ax.set_xlim(xmn, xmx)
    if fixZLims:
        ax.set_ylim(zmn, zmx)
    if title is not None:
        ax.set_title(title, fontsize=fs)
    ax.tick_params(labelsize=fs)
    #fig.tight_layout()
    if returnIm:
        return fig, ax, im
    else:
        return fig, ax

def readPlot3Ddata(chid, resultDir, time):
    xyzFiles = glob.glob("%s%s*.xyz"%(resultDir, chid))
    tFiles = glob.glob(xyzFiles[0].replace('.xyz','*.q'))
    
    tNames = [x.split(chid)[-1] for x in tFiles]
    timeNameLength = len(tNames[0].split('_')[-2])
    
    if timeNameLength == 8:
        times = [x.split('.q')[0].split('_') for x in tNames]
        times = [float(x[2])+float(x[3])/100 for x in times]
    elif timeNameLength == 2:
        times = [x.split('.q')[0].split('_') for x in tNames]
        times = [round(float(x[2])+float(x[3])/100,2) for x in times]
    
    grids = []
    datas = []
    
    ind = np.argmin(abs(np.array(times) - time))
    dataFile = tFiles[ind]
    
    for xyzFile in xyzFiles:
        #dataFile = buildDataFile(xyzFile, time)
        
        grid, gridHeader = readXYZfile(xyzFile)
        data, dataHeader = readP3Dfile(dataFile)
        (nx, ny, nz) = (dataHeader[0], dataHeader[1], dataHeader[2])
        
        
        
        printExtents(grid, data)
        
        xGrid, yGrid, zGrid = rearrangeGrid(grid, dataHeader)
        data = np.reshape(data, (nx, ny, nz, 5), order='F')
        grids.append([xGrid, yGrid, zGrid])
        datas.append(data)
        
    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    data_abs = np.zeros((xGrid_abs.shape[0],
                         xGrid_abs.shape[1],
                         xGrid_abs.shape[2],
                         5))
    data_abs[:,:,:,:] = np.nan
    
    for grid, data in zip(grids, datas):
        (xGrid, yGrid, zGrid) = (grid[0], grid[1], grid[2])
        
        xloc = np.where(abs(xGrid_abs-xGrid[0,0,0]) == 0)[0][0]
        yloc = np.where(abs(yGrid_abs-yGrid[0,0,0]) == 0)[1][0]
        zloc = np.where(abs(zGrid_abs-zGrid[0,0,0]) == 0)[2][0]
        
        (NX, NY, NZ) = np.shape(xGrid)
        
        data_abs[xloc:xloc+NX, yloc:yloc+NY, zloc:zloc+NZ,:] = data
    return grid_abs, data_abs

def extractResultDirAndChidFromSlcfName(slcfFile):
    resultDir = os.sep.join(os.path.abspath(slcfFile).split(os.sep)[:-1])
    chid = '_'.join(os.path.abspath(slcfFile).split(os.sep)[-1].split('_')[:-2])
    return resultDir, chid

def readSingleSlcfFile(slcfFile,
                       time=None, dt=None, saveTimesFile=False, endianness="<"):
    if saveTimesFile:
        timesFile = slcfFile.replace('.sf','_times.csv')
    else:
        timesFile = None
    resultDir, chid = extractResultDirAndChidFromSlcfName(slcfFile)
    #endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    timesSLCF = readSLCFtimes(slcfFile, timesFile, endianness=endianness)
    times = []
    f = zopen(slcfFile)
    
    qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
    # Check if slice is 2-dimensional
    threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
    
    if threeDimSlice:
        (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
        # Check if slice is 3-D
        #print(slcfFile, qty, sName, uts, iX, eX, iY, eY, iZ, eZ)
        shape = (NX+1, NY+1, NZ+1)
        if time == None:
            NT = len(timesSLCF)
            datas2 = np.zeros((NX+1, NY+1, NZ+1, NT))
            for i in range(0, NT):
                t, data = readNextTime(f, NX, NY, NZ, datatype)
                data = np.reshape(data, shape, order='F')
                datas2[:, :, :, i] = data
            times = timesSLCF
        elif (time != None) and (dt == None):
            datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
            i = np.argmin(abs(timesSLCF-time))
            f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, 0] = data
            times = [timesSLCF[i]]
        elif (time != None) and (dt != None):
            datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
            i = np.argmin(abs(timesSLCF - (time - dt/2)))
            j = np.argmin(abs(timesSLCF - (time + dt/2)))
            f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
            data = True
            for ii in range(i, j+1):
                if data is not False:
                    t, data = readNextTime(f, NX, NY, NZ, datatype)
                    data = np.reshape(data, shape, order='F')
                    datas2[:, :, :, 0] += data
            if j - i > 0:
                datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
            times = [timesSLCF[i]]
        lims = [iX, eX, iY, eY, iZ, eZ]
    else:
        (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
        #print("2-D slice:", slcfFile)
        shape = (NX+1, NY+1, NZ+1)
        if time == None:
            NT = len(timesSLCF)
            datas2 = np.zeros((NX+1, NY+1, NZ+1, NT))
            for i in range(0, NT):
                t, data = readNextTime(f, NX, NY, NZ, datatype)
                data = np.reshape(data, shape, order='F')
                datas2[:, :, :, i] = data
            times = timesSLCF
        elif (time != None) and (dt == None):
            datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
            i = np.argmin(abs(timesSLCF-time))
            f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, 0] = data
            times = [timesSLCF[i]]
        elif (time != None) and (dt != None):
            datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
            i = np.argmin(abs(timesSLCF - (time - dt/2)))
            j = np.argmin(abs(timesSLCF - (time + dt/2)))
            f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
            for ii in range(i, j+1):
                t, data = readNextTime(f, NX, NY, NZ, datatype)
                data = np.reshape(data, shape, order='F')
                datas2[:, :, :, 0] += data
            if j - i > 0:
                datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
            times = [timesSLCF[i]]
        lims = [iX, eX, iY, eY, iZ, eZ]
    f.close()
    return lims, datas2, times


def readSLCF3Ddata(chid, resultDir, quantityToExport,
                   time=None, dt=None, saveTimesFile=False):
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    if '.zip' in resultDir:
        xyzFiles = getFileListFromZip(resultDir, chid, 'xyz')
    else:
        xyzFiles = glob.glob("%s%s%s*.xyz"%(resultDir,os.sep, chid))
    print(xyzFiles)
    grids = defaultdict(bool)
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid)
        
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s_"%(chid, mesh)
        if '.zip' in resultDir:
            slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s%s*.sf"%(resultDir, os.sep, meshStr))
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        times3D = []
        datas3D = []
        lims3D = []
        print("%s%s_*.sf"%(resultDir, meshStr))
        for slcfFile in slcfFiles:
            if saveTimesFile:
                timesFile = slcfFile.replace('.sf','_times.csv')
            else:
                timesFile = None
            timesSLCF = readSLCFtimes(slcfFile, timesFile, endianness=endianness)
            times = []
            f = zopen(slcfFile)

            qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
            # Check if slice is correct quantity
            correctQuantity = (qty == quantityToExport)
            # Check if slice is 2-dimensional
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            #print(qty, quantityToExport)
            if correctQuantity and threeDimSlice:
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                # Check if slice is 3-D
                #print(slcfFile, qty, sName, uts, iX, eX, iY, eY, iZ, eZ)
                shape = (NX+1, NY+1, NZ+1)
                if time == None:
                    NT = len(timesSLCF)
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, NT))
                    for i in range(0, NT):
                        t, data = readNextTime(f, NX, NY, NZ, datatype)
                        data = np.reshape(data, shape, order='F')
                        datas2[:, :, :, i] = data
                    times = timesSLCF
                elif (time != None) and (dt == None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF-time))
                    f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
                    t, data = readNextTime(f, NX, NY, NZ, datatype)
                    data = np.reshape(data, shape, order='F')
                    datas2[:, :, :, 0] = data
                    times = [timesSLCF[i]]
                elif (time != None) and (dt != None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF - (time - dt/2)))
                    j = np.argmin(abs(timesSLCF - (time + dt/2)))
                    f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
                    data = True
                    for ii in range(i, j+1):
                        if data is not False:
                            t, data = readNextTime(f, NX, NY, NZ, datatype)
                            data = np.reshape(data, shape, order='F')
                            datas2[:, :, :, 0] += data
                    if j - i > 0:
                        datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
                    times = [timesSLCF[i]]
                lims3D.append([iX, eX, iY, eY, iZ, eZ])
                datas3D.append(datas2)
                times3D.append(times)
            f.close()
        grids[meshStr]['datas3D'] = datas3D
        grids[meshStr]['lims3D'] = lims3D
    
    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    if len(grids[list(grids.keys())[0]]['datas3D']) == 0:
        print("No 3D slice data found for qty %s"%(quantityToExport))
        return False, False, False
    tInd = grids[list(grids.keys())[0]]['datas3D'][0].shape[3]
    data_abs = np.zeros((xGrid_abs.shape[0],
                         xGrid_abs.shape[1],
                         xGrid_abs.shape[2],
                         tInd))
    data_abs[:, :, :, :] = np.nan
    
    abs_xs = xGrid_abs[:, 0, 0]
    abs_ys = yGrid_abs[0, :, 0]
    abs_zs = zGrid_abs[0, 0, :]
    
    for key in list(grids.keys()):
        x1 = grids[key]['xGrid'].flatten()
        y1 = grids[key]['yGrid'].flatten()
        z1 = grids[key]['zGrid'].flatten()
        
        xi = np.where(np.logical_or(np.isclose(abs_xs, x1.min()), np.isclose(abs_xs,x1.max())))[0]
        yi = np.where(np.logical_or(np.isclose(abs_ys, y1.min()), np.isclose(abs_ys,y1.max())))[0]
        zi = np.where(np.logical_or(np.isclose(abs_zs, z1.min()), np.isclose(abs_zs,z1.max())))[0]
        
        xg, yg, zg = np.meshgrid(abs_xs[xi[0]:xi[1]+1], abs_ys[yi[0]:yi[1]+1], abs_zs[zi[0]:zi[1]+1])
        xg = np.swapaxes(xg, 0, 1)
        yg = np.swapaxes(yg, 0, 1)
        zg = np.swapaxes(zg, 0, 1)
        
        p = np.zeros((xg.flatten().shape[0],3))
        p[:, 0] = xg.flatten()
        p[:, 1] = yg.flatten()
        p[:, 2] = zg.flatten()
        
        xs = grids[key]['xGrid'][:, 0, 0]
        ys = grids[key]['yGrid'][0, :, 0]
        zs = grids[key]['zGrid'][0, 0, :]        
        d = np.zeros_like(data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, zi[0]:zi[1]+1, 0])
        for t in range(0, tInd):
            d[:, :, :] = np.nan
            if t < grids[key]['datas3D'][0].shape[3]:
                lims3D = grids[key]['lims3D']
                iX, eX, iY, eY, iZ, eZ = lims3D[0]
                #d[iX:eX+1,iY:eY+1,iZ:eZ+1] = grids[key]['datas3D'][0][:, :, :, t]
                d = grids[key]['datas3D'][0][:, :, :, t] #hybrid meshing
                #print(xs.shape, ys.shape, zs.shape, grids[key]['datas3D'][0][:, :, :, t].shape)
                #print(grids[key]['lims3D'])
                #d2 = scpi.interpn((xs, ys, zs), grids[key]['datas3D'][0][:, :, :, t], p, method='linear', fill_value=None, bounds_error=False)
                d2 = scpi.interpn((xs, ys, zs), d, p, method='linear', fill_value=None, bounds_error=False)
                d2r = np.reshape(d2, xg.shape)
                data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, zi[0]:zi[1]+1, t] = d2r
    
    
    
    
    '''
    for key in list(grids.keys()):
        xGrid = grids[key]['xGrid']
        yGrid = grids[key]['yGrid']
        zGrid = grids[key]['zGrid']
        datas3D = grids[key]['datas3D']
        lims3D = grids[key]['lims3D']
        print("Starting Grid %s"%(key))
        for data, lim, times in zip(datas3D, lims3D, times3D):
            xloc_mn = np.where(np.isclose(
                    abs(xGrid_abs - xGrid[lim[0], 0, 0]),
                    0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(
                    abs(xGrid_abs - xGrid[lim[1], 0, 0]),
                    0, atol=1e-04))[0][0]
            yloc_mn = np.where(np.isclose(
                    abs(yGrid_abs - yGrid[0, lim[2], 0]),
                    0, atol=1e-04))[1][0]
            yloc_mx = np.where(np.isclose(
                    abs(yGrid_abs - yGrid[0, lim[3], 0]),
                    0, atol=1e-04))[1][0]
            zloc_mn = np.where(np.isclose(
                    abs(zGrid_abs - zGrid[0, 0, lim[4]]),
                    0, atol=1e-04))[2][0]
            zloc_mx = np.where(np.isclose(
                    abs(zGrid_abs - zGrid[0, 0, lim[5]]),
                    0, atol=1e-04))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            ANX = xloc_mx-xloc_mn + 1
            ANY = yloc_mx-yloc_mn + 1
            ANZ = zloc_mx-zloc_mn + 1
            if (NX != ANX) or (NY != ANY) or (NZ != ANZ):
                x = xGrid[lim[0]:lim[1]+1, 0, 0]
                y = yGrid[0, lim[2]:lim[3]+1, 0]
                z = zGrid[0, 0, lim[4]:lim[5]+1]
                
                xi = grid_abs[xloc_mn:xloc_mx+1, 
                              yloc_mn:yloc_mx+1, 
                              zloc_mn:zloc_mx+1,
                              0].flatten()
                yi = grid_abs[xloc_mn:xloc_mx+1,
                              yloc_mn:yloc_mx+1,
                              zloc_mn:zloc_mx+1,
                              1].flatten()
                zi = grid_abs[xloc_mn:xloc_mx+1,
                              yloc_mn:yloc_mx+1,
                              zloc_mn:zloc_mx+1,
                              2].flatten()
                
                x = np.round(x, decimals=4)
                y = np.round(y, decimals=4)
                z = np.round(z, decimals=4)
                
                xi = np.round(xi, decimals=4)
                yi = np.round(yi, decimals=4)
                zi = np.round(zi, decimals=4)
                
                xi[xi < np.min(x)] = np.min(x)
                xi[xi > np.max(x)] = np.max(x)
                yi[yi < np.min(y)] = np.min(y)
                yi[yi > np.max(y)] = np.max(y)
                zi[zi < np.min(z)] = np.min(z)
                zi[zi > np.max(z)] = np.max(z)
                
                tmpGrid = np.array([xi, yi, zi]).T
                for i in range(0, NT):
                    interpolator = scpi.RegularGridInterpolator(
                            (x, y, z), data[:, :, :, i])
                    data2 = interpolator(tmpGrid)
                    data2 = np.reshape(
                            data2, (ANX, ANY, ANZ), order='C')
                    try:
                        data_abs[xloc_mn:xloc_mx+1,
                                 yloc_mn:yloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading mesh %s at time %0.0f"%(key, i))
            else:
                try:
                    data_abs[xloc_mn:xloc_mx+1,
                             yloc_mn:yloc_mx+1,
                             zloc_mn:zloc_mx+1,
                             :] = data
                except:
                    try:
                        NTT = min([data_abs.shape[3], data.shape[3]])
                        data_abs[xloc_mn:xloc_mx+1,
                                 yloc_mn:yloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 :NTT] = data[:, :, :, :NTT]
                        print("Error loading mesh %s at time %0.0f"%(key, i))
                    except:
                        print("Error loading mesh %s at all times"%(key))
                        print(data.shape, data_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1, :].shape, NTT)
    '''
    return grid_abs, data_abs, times3D[0]





def readSLCF2Ddata(chid, resultDir, quantityToExport,
                   time=None, dt=None):
    fdsFile = fdsFileOperations()
    fdsFile.importFile(resultDir+os.sep+chid+'.fds')
    slcfDir = resultDir
    if fdsFile.dump['ID'] is not False:
        if fdsFile.dump['ID']['RESULTS_DIR'] is not False:
            slcfDir = resultDir + os.sep + fdsFile.dump['ID']['RESULTS_DIR'] + os.sep
    if '.zip' in resultDir:
        xyzFiles = getFileListFromZip(slcfDir, chid, 'xyz')
    else:
        xyzFiles = glob.glob("%s%s%s*.xyz"%(slcfDir, os.sep, chid))
    grids = defaultdict(bool)
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid)
        
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
        if '.zip' in slcfDir:
            slcfFiles = getFileListFromZip(slcfDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s%s_*.sf"%(slcfDir, os.sep, meshStr))
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        datas2D = []
        lims2D = []
        coords2D = []
        times2D = []
        for slcfFile in slcfFiles:
            timesSLCF = readSLCFtimes(slcfFile, None, endianness)
            times = []
            f = zopen(slcfFile)
            
            qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
            # Check if slice is correct quantity
            correctQuantity = (qty == quantityToExport)
            # Check if slice is 2-dimensional
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            if correctQuantity and not threeDimSlice:
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                print("2-D slice:", slcfFile)
                shape = (NX+1, NY+1, NZ+1)
                if time == None:
                    NT = len(timesSLCF)
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, NT), dtype=np.float32)
                    for i in range(0, NT):
                        t, data = readNextTime(f, NX, NY, NZ, datatype)
                        data = np.reshape(data, shape, order='F')
                        datas2[:, :, :, i] = np.array(data, dtype=np.float32)
                    times = timesSLCF
                elif (time != None) and (dt == None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1), dtype=np.float32)
                    i = np.argmin(abs(timesSLCF-time))
                    f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
                    t, data = readNextTime(f, NX, NY, NZ, datatype)
                    data = np.reshape(data, shape, order='F')
                    datas2[:, :, :, 0] = np.array(data, dtype=np.float32)
                    times = [timesSLCF[i]]
                elif (time != None) and (dt != None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1), dtype=np.float32)
                    i = np.argmin(abs(timesSLCF - (time - dt/2)))
                    j = np.argmin(abs(timesSLCF - (time + dt/2)))
                    f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
                    for ii in range(i, j+1):
                        t, data = readNextTime(f, NX, NY, NZ, datatype)
                        data = np.reshape(data, shape, order='F')
                        datas2[:, :, :, 0] += np.array(data, dtype=np.float32)
                    if j - i > 0:
                        datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
                    times = [timesSLCF[i]]
                lims2D.append([iX, eX, iY, eY, iZ, eZ])
                datas2D.append(datas2)
                coords2D.append([xGrid[iX, iY, iZ],
                                 yGrid[iX, iY, iZ],
                                 zGrid[iX, iY, iZ]])
                times2D.append(np.array(times))
                timesOut = np.array(times)
            f.close()
        grids[meshStr]['datas2D'] = datas2D
        grids[meshStr]['lims2D'] = lims2D
        grids[meshStr]['coords2D'] = coords2D
        
    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    tInds = [] 
    for i, key in enumerate(list(grids.keys())):
        try:
            tInd = grids[key]['datas2D'][0].shape[3]
            tInds.append(tInd)
        except:
            pass
    tInd = np.min(tInds)
    data_abs = np.zeros((xGrid_abs.shape[0],
                         xGrid_abs.shape[1],
                         xGrid_abs.shape[2],
                         tInd))
    data_abs[:, :, :, :] = np.nan
    for key in list(grids.keys()):
        xGrid = grids[key]['xGrid']
        yGrid = grids[key]['yGrid']
        zGrid = grids[key]['zGrid']
        datas2D = grids[key]['datas2D']
        lims2D = grids[key]['lims2D']
        coords2D = grids[key]['coords2D']
        
        for data, coord in zip(datas2D, coords2D):
            print(coord)
            xloc = np.where(np.isclose(
                    abs(xGrid_abs - coord[0]), 0, atol=1e-06))[0][0]
            yloc = np.where(np.isclose(
                    abs(yGrid_abs - coord[1]), 0, atol=1e-06))[1][0]
            zloc = np.where(np.isclose(
                    abs(zGrid_abs - coord[2]), 0, atol=1e-06))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            NT = min([NT, data_abs.shape[3]])
            data_abs[xloc:xloc+NX,
                     yloc:yloc+NY,
                     zloc:zloc+NZ,
                     :NT] = data[:, :, :, :NT]
            if (coord[0] == 0.0) and (coord[2] == 0.3514):
                print(xloc, xloc+NX, yloc, yloc+NY, zloc, zloc+NZ )
                plt.imshow(data[:, 0, :, 0])
    return grid_abs, data_abs, timesOut







def extractPoint(point, grid, data):
    ind = np.argmin(np.sum(abs(grid-point),axis=3).flatten())
    ind = np.unravel_index(ind, grid[:,:,:,0].shape)
    x = grid[ind[0],ind[1],ind[2],0]
    y = grid[ind[0],ind[1],ind[2],1]
    z = grid[ind[0],ind[1],ind[2],2]
    d = data[ind[0],ind[1],ind[2],:]
    err = abs(point[0]-x)+abs(point[1]-y)+abs(point[2]-z)
    if err > 0.25:
        errPt = "%0.4f, %0.4f, %0.4f, is %0.4f"%(x, y, z, err)
        print("Warning, error for point %s."%(errPt))
    return d

def readNextTime(f, NX, NY, NZ, datatype):
    _ = np.frombuffer(f.read(8), dtype=datatype)
    time = np.frombuffer(f.read(4), dtype=datatype)
    _ = np.frombuffer(f.read(8), dtype=datatype)
    try:
        data = np.frombuffer(f.read((NX+1)*(NY+1)*(NZ+1)*4), 
                             dtype=datatype)
    except:
        data = False
    return time, data

def readSLCFheader(f, endianness, byteSize=False):
    data = f.read(142)
    header = data[:110]
    size = struct.unpack('%siiiiii'%(endianness), data[118:142])
    tmp = header.split(b'\x1e')
    quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
    shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
    units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
    
    if byteSize:
        return quantity, shortName, units, size
    else:
        iX, eX, iY, eY, iZ, eZ = size
        return quantity, shortName, units, iX, eX, iY, eY, iZ, eZ

def readSLCFtimes(file, timesFile=None, endianness=None):
    if timesFile != None:
        if os.path.exists(timesFile) == True:
            times = np.loadtxt(timesFile, delimiter=',')
            return times
    if endianness == None:
        resultDir, chid = extractResultDirAndChidFromSlcfName(file)
        endianness = getEndianness(resultDir, chid)
    f = zopen(file)
    qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    data = f.read()
    f.close()
    if len(data) % 4 == 0:
        fullFile = np.frombuffer(data, dtype=np.float32)
    else:
        #print(len(data))
        remainder = -1*int(len(data) % 4)
        #print(len(data[:remainder]))
        fullFile = np.frombuffer(data[:remainder], dtype=np.float32)
    times = fullFile[2::(NX+1)*(NY+1)*(NZ+1)+5]
    if timesFile != None:
        np.savetxt(timesFile, times)
    return times

def readSLCFquantities(chid, resultDir):
    try:
        smvFile = os.path.join(resultDir, '%s.smv'%(chid))
        smvData = parseSMVFile(smvFile)
    except:
        smvDir = os.sep.join(os.path.abspath(resultDir).split(os.sep)[:-1])
        smvFile = os.path.join(smvDir, '%s.smv'%(chid))
        smvData = parseSMVFile(smvFile)
        files = smvData['files']
        files2 = defaultdict(bool)
        files2['SLICES'] = defaultdict(bool)
        for key in list(files['SLICES'].keys()):
            n = os.path.basename(key)
            files2['SLICES'][n] = files['SLICES'][key]
        smvData['files'] = files2
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    
    if '.zip' in resultDir:
        slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
        zip = zipfile.ZipFile(resultDir, 'r')
    else:
        slcfFiles = glob.glob("%s%s_*.sf"%(resultDir, chid))
    endianness = getEndianness(resultDir, chid)
    quantities = []
    dimensions = []
    meshes = []
    centers = []
    for file in slcfFiles:
        if '.zip' in resultDir:
            f = zip.open(file.split("%s%s"%('.zip',os.sep))[1])
        else:
            f = open(file, 'rb')
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
        quantities.append(qty)
        dimensions.append([iX, eX, iY, eY, iZ, eZ])
        n = file.split(chid)[-1].split('_')
        meshStr = n[-2]
        meshes.append(meshStr)
        #print(files['SLICES'][file.split(os.sep)[-1]])
        #print(file)
        centers.append(files['SLICES'][file.split(os.sep)[-1]]['CELL_CENTERED'])
        f.close()
    if '.zip' in resultDir:
        zip.close()
    return quantities, slcfFiles, dimensions, meshes, centers

def buildQTYstring(chid, resultDir, qty):
    quantities, slcfFiles = readSLCFquantities(chid, resultDir)
    quantitiesCheck = [True if qty == x else False for x in quantities]
    inds = np.where(quantitiesCheck)[0]
    if len(inds) == 0:
        print("Quantity %s unknown."%(qty))
        print("Known quantities:")
        for qnty in sorted(set(quantities)):
            print(qnty)
    else:
        ind = inds[0]
    quantityStr = slcfFiles[ind].split('.sf')[0].split('_')[-1]
    return quantityStr

def getLimsFromGrid(grid):
    xGrid = grid[:, :, :, 0]
    yGrid = grid[:, :, :, 1]
    zGrid = grid[:, :, :, 2]
    
    (xmn, xmx) = (xGrid[:, 0, 0].min(), xGrid[:, 0, 0].max())
    (ymn, ymx) = (yGrid[0, :, 0].min(), yGrid[0, :, 0].max())
    (zmn, zmx) = (zGrid[0, 0, :].min(), zGrid[0, 0, :].max())
    
    return [xmn, xmx, ymn, ymx, zmn, zmx]

def visualizePlot3D(x, z, T, U, V, W, HRR,
                    qnty_mn=None, qnty_mx=None):
    cmap = buildSMVcolormap()
    
    xrange = x.max()-x.min()
    zrange = z.max()-z.min()
    
    if zrange > xrange:
        plt.figure(figsize=(12*xrange/zrange,12))
    else:
        plt.figure(figsize=(12,12*zrange/xrange))
    if qnty_mn is None:
        qnty_mn = np.nanmin(T)
    if qnty_mx is None:
        qnty_mx = np.nanmax(T)
    levels = np.linspace(qnty_mn, qnty_mx, 100)
    plt.contourf(x, z, T, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx,
                 levels=levels, extend='both')
    plt.colorbar()


def read2dSliceFile(slcfFile, chid, time=None, dt=None, cen=False):
    resultDir = os.sep.join(os.path.abspath(slcfFile).split(os.sep)[:-1])
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    timesSLCF = readSLCFtimes(slcfFile, None, endianness=endianness)
    times = []
    
    xyzFile = '%s%s'%('_'.join(slcfFile.split('_')[:-1]), '.xyz')
    grids = getGridsFromXyzFiles([xyzFile], chid)
    grid = grids[list(grids.keys())[0]]
    xGrid = grid['xGrid']
    yGrid = grid['yGrid']
    zGrid = grid['zGrid']
    
    if cen:
        dx = np.round((grid['xGrid'][-1, 0, 0] - grid['xGrid'][0, 0, 0]) / (grid['xGrid'][:, 0, 0].shape[0]-1), decimals=4)
        dy = np.round((grid['yGrid'][0, -1, 0] - grid['yGrid'][0, 0, 0]) / (grid['yGrid'][0, :, 0].shape[0]-1), decimals=4)
        dz = np.round((grid['zGrid'][0, 0, -1] - grid['zGrid'][0, 0, 0]) / (grid['zGrid'][0, 0, :].shape[0]-1), decimals=4)
    else:
        dx = 0
        dy = 0
        dz = 0
    f = zopen(slcfFile)
    
    qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
    headerSize = 142
    
    (NX, NY, NZ) = (eX - iX, eY - iY, eZ - iZ)
    if (NX == 0):
        slcf_axis = 1
    elif (NY == 0):
        slcf_axis = 2
    elif (NZ == 0):
        slcf_axis = 3
    else:
        slcf_axis = -1
    
    if slcf_axis < 0:
        print("Not a 2-D slice.")
        return None
    shape = (NX+1, NY+1, NZ+1)
    if (time == None) and (dt == None):
        NT = len(timesSLCF)
        datas2 = np.zeros((NX+1, NY+1, NZ+1, NT), dtype=np.float32)
        for i in range(0, NT):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, i] = data
        times = timesSLCF
    elif (time == None) and (dt != None):
        NT = len(timesSLCF)
        datas3 = np.zeros((NX+1, NY+1, NZ+1, NT))
        for i in range(0, NT):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas3[:, :, :, i] = data
        datas2 = datas3.copy()
        slcfDt = np.median(timesSLCF[1:] - timesSLCF[:-1])
        #print(slcfDt, dt)
        for i in range(0, NT):
            t1 = max([0, timesSLCF[i]-dt/2])
            t2 = min([timesSLCF[-1], timesSLCF[i]+dt/2])
            inds = np.where(np.logical_and(timesSLCF >= t1, timesSLCF <= t2))[0]
            datas2[:, :, :, i] = np.nanmean(datas3[:, :, :, inds], axis=3)
            #ind1 = max([0, int(i-(dt/2)/slcfDt)])
            #ind2 = min([NT, int(i+(dt/2)/slcfDt)])
            if i > 55 and i < 130:#(timesSLCF[i] > 55) and (timesSLCF[i] < 65):
                #print(i, timesSLCF[i], dt, slcfDt, ind1, ind2)
                pass
                
            #datas2[:, :, :, i] = np.nanmean(datas3[:, :, :, ind1:ind2], axis=3)
        times = timesSLCF
    elif (time != None) and (dt == None):
        datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
        i = np.argmin(abs(timesSLCF-time))
        f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
        t, data = readNextTime(f, NX, NY, NZ, datatype)
        data = np.reshape(data, shape, order='F')
        datas2[:, :, :, 0] = data
        times = [timesSLCF[i]]
        NT = 1
    elif (time != None) and (dt != None):
        datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
        
        t1 = max([0, time-dt/2])
        t2 = min([timesSLCF[-1], time+dt/2])
        inds = np.where(np.logical_and(timesSLCF >= t1, timesSLCF <= t2))[0]
        
        #print(timesSLCF, t1, t2, timesSLCF[inds])
        ts = []
        for ind in inds:
            f.seek(ind * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)) + headerSize, 0)
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, 0] += data
            ts.append(t)
        datas2[:, :, :, 0] = datas2[:, :, :, 0] / float(len(inds))
        times = [(np.min(ts) + np.max(ts))/2]
        '''
        i = np.argmin(abs(timesSLCF - (time - dt/2)))
        j = np.argmin(abs(timesSLCF - (time + dt/2)))
        f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
        for ii in range(i, j+1):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, 0] += data
        if j - i > 0:
            datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
        times = [(timesSLCF[i] + timesSLCF[j])/2]
        
        print(i, timesSLCF[i], j, timesSLCF[j])
        '''
        NT = 1
    coords = [xGrid[iX, iY, iZ], xGrid[eX, eY, eZ],
              yGrid[iX, iY, iZ], yGrid[eX, eY, eZ],
              zGrid[iX, iY, iZ], zGrid[eX, eY, eZ]]
    f.close()
    if slcf_axis == 1:
        x = yGrid[iX, iY:eY+1, iZ:eZ+1]-dy/2
        z = zGrid[iX, iY:eY+1, iZ:eZ+1]-dz/2
        d = datas2[0, :, :, :]
    elif slcf_axis == 2:
        x = xGrid[iX:eX+1, iY, iZ:eZ+1]-dx/2
        z = zGrid[iX:eX+1, iY, iZ:eZ+1]-dz/2
        d = datas2[:, 0, :, :]
    elif slcf_axis == 3:
        x = xGrid[iX:eX+1, iY:eY+1, iZ]-dx/2
        z = yGrid[iX:eX+1, iY:eY+1, iZ]-dy/2
        d = datas2[:, :, 0, :]
    
    return x, z, d, times, coords

def getAxisAndValueFromXB(XB, grid, cen):
    NX = XB[1] - XB[0]
    NY = XB[3] - XB[2]
    NZ = XB[5] - XB[4]
    
    dx = np.round((grid['xGrid'][-1, 0, 0] - grid['xGrid'][0, 0, 0]) / (grid['xGrid'][:, 0, 0].shape[0]-1), decimals=4)
    dy = np.round((grid['yGrid'][0, -1, 0] - grid['yGrid'][0, 0, 0]) / (grid['yGrid'][0, :, 0].shape[0]-1), decimals=4)
    dz = np.round((grid['zGrid'][0, 0, -1] - grid['zGrid'][0, 0, 0]) / (grid['zGrid'][0, 0, :].shape[0]-1), decimals=4)
    if (NX == 0):
        axis = 1
        if cen:
            value = grid['xGrid'][XB[0], 0, 0] - dx/2
        else:
            value = grid['xGrid'][XB[0], 0, 0]
    elif (NY == 0):
        axis = 2
        if cen:
            value = grid['yGrid'][0, XB[2], 0] - dy/2
        else:
            value = grid['yGrid'][0, XB[2], 0]
    elif (NZ == 0):
        axis = 3
        if cen:
            value = grid['zGrid'][0, 0, XB[4]] - dz/2
        else:
            value = grid['zGrid'][0, 0, XB[4]]
    else:
        axis = -1
        value = -1
    return axis, value

def query2dAxisValue(resultDir, chid, quantity, axis, value, time=None, dt=None, atol=1e-8, printInfo=False):
    xyzFiles = getFileListFromResultDir(resultDir, chid, 'xyz')
    if len(xyzFiles) == 0:
        xyzFiles = getFileListFromResultDir(resultDir+'*'+os.sep, chid, 'xyz')
        resultDir2 = os.path.dirname(xyzFiles[0]) + os.sep 
    else:
        resultDir2 = resultDir
    grids = getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = getAbsoluteGrid(grids)
    if abs(axis) == 1:
        xAbs = grids_abs[0, :, :, 1]
        zAbs = grids_abs[0, :, :, 2]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    elif abs(axis) == 2:
        xAbs = grids_abs[:, 0, :, 0]
        zAbs = grids_abs[:, 0, :, 2]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    elif abs(axis) == 3:
        xAbs = grids_abs[:, :, 0, 0]
        zAbs = grids_abs[:, :, 0, 1]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    quantities, slcfFiles, dimensions, meshes, centers = readSLCFquantities(chid, resultDir2)
    if printInfo:
        print(quantities)
        print(slcfFiles)
    
    datas = defaultdict(bool)
    foundSlice = False
    available_slices = []
    for qty, slcfFile, dim, cen in zip(quantities, slcfFiles, dimensions, centers):
        if qty == quantity:
            n = slcfFile.split(chid)[-1].split('_')
            meshStr = n[-2]
            if meshStr == '': meshStr = chid
            tmp = [len(x) for x in list(grids.keys())]
            if np.max(tmp) < 4:
                meshStr = str(int(meshStr))
            #mesh = slcfFile.split(chid)[-1].split('.sf')[0].split('_')[-2]
            #meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
            
            slcf_axis, slcf_value = getAxisAndValueFromXB(dim, grids[meshStr], cen)
            
            available_slices.append([slcf_axis, slcf_value])
            if np.isclose(axis, abs(slcf_axis)) and np.isclose(slcf_value, value, atol=atol):
                if printInfo:
                    print("Reading %s"%(slcfFile), time, dt)
                x, z, d, times, coords = read2dSliceFile(slcfFile, chid, time=time, dt=dt)
                if cen:
                    x = (x[1:, :]+x[:-1, :])/2
                    z = (z[:,1:]+z[:,:-1])/2
                    x = x[:,:-1]
                    z = z[:-1,:]
                    d = d[1:,1:]
                #print(x.shape, z.shape, d.shape)
                slcfName = '%s_%0.4f_%0.4f_%0.4f_%0.4f_%0.4f_%0.4f'%(qty, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5])
                datas[slcfName] = defaultdict(bool)
                datas[slcfName]['limits'] = coords
                datas[slcfName]['times'] = times
                datas[slcfName]['datas'] = d.copy()
                datas[slcfName]['x'] = x.copy()
                datas[slcfName]['z'] = z.copy()
                datas[slcfName]['center'] = cen
                foundSlice = True
    if not foundSlice:
        print("Warning did not find a 2-D slice of %s on axis %d at value %0.4f"%(quantity, axis, value))
        print("Available slices of qty %s:"%(quantity))
        print('\tAxis\tValue')
        for slcf_axis, slcf_value in available_slices:
            print("\t%d\t%0.4f"%(slcf_axis, slcf_value))
        print("Note an axis of -1 indicates a 3-D slice which is not currently supported in this function.")
        return None        
    data_abs = np.zeros((xAbs.shape[0], xAbs.shape[1], len(times)), dtype=np.float32)
    data_abs_c = np.zeros((xAbs_c.shape[0], xAbs_c.shape[1], len(times)), dtype=np.float32)
    for slcfName in list(datas.keys()):
        x = datas[slcfName]['x']
        z = datas[slcfName]['z']
        d = datas[slcfName]['datas']
        c = datas[slcfName]['center']
        
        if c:
            xloc_mn = np.where(np.isclose(abs(xAbs_c - np.nanmin(x)), 0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xAbs_c - np.nanmax(x)), 0, atol=1e-04))[0][0]
            zloc_mn = np.where(np.isclose(abs(zAbs_c - np.nanmin(z)), 0, atol=1e-04))[1][0]
            zloc_mx = np.where(np.isclose(abs(zAbs_c - np.nanmax(z)), 0, atol=1e-04))[1][0]
        else:
            xloc_mn = np.where(np.isclose(abs(xAbs - np.nanmin(x)), 0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xAbs - np.nanmax(x)), 0, atol=1e-04))[0][0]
            zloc_mn = np.where(np.isclose(abs(zAbs - np.nanmin(z)), 0, atol=1e-04))[1][0]
            zloc_mx = np.where(np.isclose(abs(zAbs - np.nanmax(z)), 0, atol=1e-04))[1][0]
            
        (NX, NZ, NT) = np.shape(d)
        ANX = xloc_mx-xloc_mn + 1
        ANZ = zloc_mx-zloc_mn + 1
        
        if (NX != ANX) or (NZ != ANZ):
            if c:
                xi = xAbs_c[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
                zi = zAbs_c[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
            else:
                xi = xAbs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
                zi = zAbs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
            
            x = np.round(x, decimals=4)
            z = np.round(z, decimals=4)
            
            xi = np.round(xi, decimals=4)
            zi = np.round(zi, decimals=4)
            
            xi[xi < np.min(x)] = np.min(x)
            xi[xi > np.max(x)] = np.max(x)
            zi[zi < np.min(z)] = np.min(z)
            zi[zi > np.max(z)] = np.max(z)
            
            tmpGrid = np.array([xi, zi]).T
            for i in range(0, NT):
                interpolator = scpi.RegularGridInterpolator((x[:, 0], z[0, :]), d[:, :, i])
                data2 = interpolator(tmpGrid)
                data2 = np.reshape(data2, (ANX, ANZ), order='C')
                if c:
                    try:
                        data_abs_c[xloc_mn:xloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading %s at time %0.0f"%(slcfName, i))
                else:
                    try:
                        data_abs[xloc_mn:xloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading %s at time %0.0f"%(slcfName, i))
        else:
            try:
                data_abs[xloc_mn:xloc_mx+1,
                         zloc_mn:zloc_mx+1,
                         :] = d
            except:
                try:
                    NTT = min([data_abs.shape[2], d.shape[2]])
                    data_abs[xloc_mn:xloc_mx+1,
                             zloc_mn:zloc_mx+1,
                             :NTT] = d[:, :, :NTT]
                    print("Error loading %s at time %0.0f"%(slcfName, i))
                except:
                    print("Error loading %s at all times"%(slcfName))
                    print(d.shape, data_abs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1, :].shape, NTT)
                    
    data_abs_out = defaultdict(bool)
    data_abs_out['x'] = xAbs
    data_abs_out['z'] = zAbs
    data_abs_out['datas'] = data_abs
    data_abs_out['times'] = times
    return data_abs_out

def renderSliceCsvs(data, chid, outdir):
    times = data['times']
    xs = data['x'][:, 0]
    zs = data['z'][0, :]
    for i in range(0, len(times)):
        time = times[i]
        outFile = os.path.join(outdir, '%s_%0.4f.csv'%(chid, time))
        print("Rendering time %0.4f to file %s"%(time, outFile))
        d = pd.DataFrame(data['datas'][:, :, i].T, index=zs, columns=xs)
        d.to_csv(outFile)


def writeSLCFheader(f, quantity, shortName, units, size, endianness):
    sz = struct.pack('%s%0.0fi'%(endianness, len(size)), *size)
    qty = str.encode("{:<30}".format(quantity))
    sn = str.encode("{:<30}".format(shortName))
    un = str.encode("{:<30}".format(units))
    f.write(b'\x1e\x00\x00\x00')
    f.write(qty)
    f.write(b'\x1e\x00\x00\x00\x1e\x00\x00\x00')
    f.write(sn)
    f.write(b'\x1e\x00\x00\x00\x1e\x00\x00\x00')
    f.write(un)
    f.write(b'\x1e\x00\x00\x00\x18\x00\x00\x00')
    f.write(sz)
    f.write(b'\x18\x00\x00\x00')

def writeSLCFTime(f, time, data, endianness):
    f.write(b'\x04\x00\x00\x00')
    t = time.tobytes()
    f.write(t)
    f.write(b'\x04\x00\x00\x00')
    f.write(struct.pack('%si'%(endianness), data.shape[0]*4))
    d = data.tobytes()
    f.write(d)
    if data.shape[0]*4 <= 65535:
        f.write(struct.pack('%sHH'%(endianness), data.shape[0]*4,0))
    else:
        f.write(struct.pack('%sI'%(endianness), data.shape[0]*4))

def writeSlice(outFile, resultDir, chid, data, times, axis, val,
                       outQty, sName, uts, meshnum, smvFile=None, endianness="<"):
    outPath = os.path.join(resultDir, outFile)
    smvPath = os.path.join(resultDir, smvFile)
    
    if axis == 1:
        (NX, NY, NZ) = (0, data[0].shape[0]-1, data[0].shape[1]-1)
        size = np.array([val, val, 0, NY, 0, NZ], dtype=np.int32)
        X = [val, val, 0, NY, 0, NZ]
    elif axis == 2:
        (NX, NY, NZ) = (data[0].shape[0]-1, 0, data[0].shape[1]-1)
        size = np.array([0, NX, val, val, 0, NZ], dtype=np.int32)
        X = [0, NX, val, val, 0, NZ]
    elif axis == 3:
        (NX, NY, NZ) = (data[0].shape[0]-1, data[0].shape[1]-1, 0)
        size = np.array([0, NX, 0, NY, val, val], dtype=np.int32)
        X = [0, NX, 0, NY, val, val]
        
    NT = len(times)
    f = zopen(outPath, 'wb')
    writeSLCFheader(f, outQty, sName, uts, size, endianness)
    shape2 = ((NX+1)*(NY+1)*(NZ+1),)
    for i in range(0, NT):
        data_out = np.reshape(data[i], shape2, order='F')
        writeSLCFTime(f, times[i], data_out, endianness)
    f.close()
    
    if smvFile is not None:
        writeSliceToSmv(smvPath, meshnum, X, outQty, outFile, sName, uts)

def writeSliceToSmv(file, meshNum, X, outQty, outFile, sName, uts):
    with open(file, 'a') as f:
        f.write('SLCF     %0.0f # STRUCTURED &     %0.0f    %0.0f     %0.0f    %0.0f     %0.0f    %0.0f !      %0.0f\n'%(meshNum, X[0], X[1], X[2], X[3], X[4], X[5], 1))
        f.write(' %s\n'%(outFile))
        f.write(' %s\n'%(outQty))
        f.write(' %s\n'%(sName))
        f.write(' %s\n\n'%(uts))

def getAxisFromLims(lims):
    iX, eX, iY, eY, iZ, eZ = lims
    (NX, NY, NZ) = (eX - iX+1, eY - iY+1, eZ - iZ+1)
    if (NX == 1):
        slcf_axis = 1
        val = eX
    elif (NY == 1):
        slcf_axis = 2
        val = eY
    elif (NZ == 1):
        slcf_axis = 3
        val = eZ
    else:
        slcf_axis = -1
        val = None
    return slcf_axis, val

def slcfTimeAverage(slcfFile, dt, outFile=None, outQty=None, outdt=None):
    
    # Read the data
    resultDir, chid = extractResultDirAndChidFromSlcfName(slcfFile)
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    timesSLCF = readSLCFtimes(slcfFile, None, endianness=endianness)
    f = zopen(slcfFile)
    qty, sName, uts, size = readSLCFheader(f, endianness, byteSize=True)
    iX, eX, iY, eY, iZ, eZ = size
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    NT = len(timesSLCF)
    datas2 = np.zeros((NX+1, NY+1, NZ+1, NT), dtype=np.float32)
    shape = (NX+1, NY+1, NZ+1)
    for i in range(0, NT):
        t, data = readNextTime(f, NX, NY, NZ, datatype)
        if data is not False:
            data2 = np.reshape(data, shape, order='F')
            datas2[:, :, :, i] = data2
    f.close()
    
    # Average the data
    if outdt is not None:
        outTimes = np.linspace(0, timesSLCF[-1], int((timesSLCF[-1]/outdt)+1), dtype=np.float32)
        
    else:
        outTimes = timesSLCF
    
    NT2 = outTimes.shape[0]
    data_avg = np.zeros((NX+1, NY+1, NZ+1, NT2), dtype=np.float32)
    
    for i in range(0, NT2):
        t = outTimes[i]
        j1 = np.argmin(abs(timesSLCF - t - dt/2))
        j0 = np.argmin(abs(timesSLCF - t + dt/2))
        #i0 = max([int(j - dt/dt2), 0])
        #i1 = min([int(j + dt/dt2), NT])
        #print(i, j0, j1, timesSLCF[j0], timesSLCF[j1])
        data_avg[:, :, :, i] = np.nanmean(datas2[:, :, :, j0:j1], axis=3)
    
    # Write the data
    if outFile is None:
        outFile = slcfFile.replace('.sf', '_avg.sf')
    if outQty is None:
        outQty = sName + ' (%0.0ds Avg)'%(dt)
    f = zopen(outFile, 'wb')
    writeSLCFheader(f, outQty, sName, uts, size)
    shape2 = ((NX+1)*(NY+1)*(NZ+1),)
    for i in range(0, NT2):
        data = np.reshape(data_avg[:,:,:,i], shape2, order='F')
        writeSLCFTime(f, outTimes[i], data)
    f.close()
    if outdt is None:
        print("Wrote %s with %0.1f average window"%(outFile, dt))
    else:
        print("Wrote %s with %0.1f average window at interval %0.1f"%(outFile, dt, outdt))
    header = qty, sName, uts, iX, eX, iY, eY, iZ, eZ
    return outFile, header
    

def slcfsTimeAverage(resultDir, chid, fdsQuantity, dt, outDir=None, outQty=None, outdt=None):
    slcfFiles = getFileList(resultDir, chid, 'sf')
    filesWithQueriedQuantity = []
    endianness = getEndianness(resultDir, chid)
    for sliceFile in slcfFiles:
        f = zopen(sliceFile)
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
        f.close()
        # Check if slice is correct quantity
        correctQuantity = (qty == fdsQuantity)
        if correctQuantity:
            filesWithQueriedQuantity.append(sliceFile)
    
    outFiles = []
    meshInfo = []
    for sliceFile in slcfFiles:
        print("Starting to avg file %s"%(sliceFile))
        outFile, header = slcfTimeAverage(sliceFile, dt, outFile=None, outQty=outQty, outdt=outdt)
        outFiles.append(outFile)
        print("wrote %s"%(outFile))
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = header
        (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
        meshNum = int(sliceFile.split('_')[-2])
        meshInfo.append([meshNum, NX, NY, NZ])
    
    # Add new slices to smokeview file
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    linesSMV = zreadlines(smvFile)
    smvData = parseSMVFile(smvFile)
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    slices = files['SLICES']
    if outQty is None:
        outQty = sName + ' (%0.0ds Avg)'%(dt)
    
    for outFile, info in zip(outFiles, meshInfo):
        fname = outFile.replace('_avg.sf','.sf').split(os.sep)[-1]
        lineText = slices[fname]['LINETEXT']
        lineText = lineText.replace('.sf', '_avg.sf')
        lineText = lineText[:-4] + ('%s'%(int(lineText.replace('\n','').split()[-1])+1000)).ljust(4) + '\n'
        linesSMV.append(lineText.replace('.sf', '_avg.sf'))
        meshNum, NX, NY, NZ = info
        #linesSMV.append('SLCF     %0.0f # STRUCTURED &     0    %0.0f     0    %0.0f     0    %0.0f !      %0.0f\n'%(meshNum, NX, NY, NZ, meshNum))
        linesSMV.append(' %s\n'%(outFile.split(os.sep)[-1]))
        linesSMV.append(' %s\n'%(outQty))
        linesSMV.append(' %s\n'%(sName))
        linesSMV.append(' %s\n'%(uts))
    
    smvText = '\n'.join(linesSMV)
    smvText = smvText + '\n'
    smvText = smvText.replace('\n\n','\n')
    
    if outDir is None:
        newSmvFile = outFile.split(os.sep)[:-1]
        newSmvFile.append(smvFile.split(os.sep)[-1].replace('.smv','_avg.smv'))
        newSmvFile = os.sep.join(newSmvFile)
    else:
        newSmvFile = outDir.split(os.sep)
        newSmvFile.append(smvFile.split(os.sep)[-1].replace('.smv','_avg.smv'))
        newSmvFile = os.sep.join(newSmvFile)
        
    with open(newSmvFile, 'w') as f:
        f.write(smvText)
    print("wrote to %s"%(newSmvFile))
    
    return outFiles, outQty, filesWithQueriedQuantity, newSmvFile


def writeXYZfile(file, grid):
    """Writes grid to an xyz file
    
    This subroutine writes a grid to an xyz file.
    
    Parameters
    ----------
    file : str
        String containing the path to an xyz file
    
    grid : array(NX, NY, NZ, 3)
        Array containing float global coordinates
    """
    nx = np.unique(grid[:, 0]).shape[0]
    ny = np.unique(grid[:, 1]).shape[0]
    nz = np.unique(grid[:, 2]).shape[0]
    v = 4
    with open(file,'wb') as f:
        f.write(b'\x0c\x00\x00\x00')
        (nx, ny, nz, v) = (int(nx), int(ny), int(nz), int(v))
        f.write(nx.to_bytes(4, 'little'))
        f.write(ny.to_bytes(4, 'little'))
        f.write(nz.to_bytes(4, 'little'))
        f.write(b'\x0c\x00\x00\x00')
        #empty = np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.float32)
        empty = np.array([0], dtype=np.float32)
        empty.tofile(f)
        d1 = grid.flatten(order='F')
        d = d1.tobytes()
        f.write(d)
        byteStr = struct.pack('I', nx*ny*nz*v*4)
        f.write(byteStr)