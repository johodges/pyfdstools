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
from collections import defaultdict
from .utilities import getFileListFromZip, zopen
from .colorSchemes import buildSMVcolormap

def time2str(time):
    timeStr = "%08.0f_%02.0f"%(np.floor(time), (time-np.floor(time))*100)
    return timeStr

def mesh2str(mesh):
    meshStr = "%04.0f"%(mesh)
    return meshStr

def readXYZfile(file):
    f = zopen(file)
    header = struct.unpack('<iiiiif', f.read(24))
    (nx, ny, nz) = (header[1], header[2], header[3])
    data = np.frombuffer(f.read(nx*ny*nz*4*4), dtype=np.float32)
    grid = np.reshape(data, (int(data.shape[0]/4),4),order='F')
    f.close()
    return grid, header[1:-1]

def readP3Dfile(file):
    with open(file,'rb') as f:
        header = np.fromfile(f,dtype=np.int32,count=5)
        _ = np.fromfile(f,dtype=np.float32,count=7)
        (nx, ny, nz) = (header[1], header[2], header[3])
        data = np.fromfile(f,dtype=np.float32,count=nx*ny*nz*5)
        data = np.reshape(data, (int(data.shape[0]/5),5),order='F')
    return data, header[1:-1]

def rearrangeGrid(grid, header):
    (xs, ys, zs) = (np.unique(grid[:,0]), np.unique(grid[:,1]), np.unique(grid[:,2]))
    xGrid, yGrid, zGrid = np.meshgrid(xs, ys, zs)
    xGrid = np.swapaxes(xGrid, 0, 1)
    yGrid = np.swapaxes(yGrid, 0, 1)
    zGrid = np.swapaxes(zGrid, 0, 1)
    return xGrid, yGrid, zGrid

def buildDataFile(meshStr, time):
    dataStr = meshStr.replace('.xyz','_%s.q'%(time2str(time)))
    return dataStr
    
def printExtents(grid, data):
    print("%0.4f < x < %0.4f"%(np.min(grid[:,0]),np.max(grid[:,0])))
    print("%0.4f < y < %0.4f"%(np.min(grid[:,1]),np.max(grid[:,1])))
    print("%0.4f < z < %0.4f"%(np.min(grid[:,2]),np.max(grid[:,2])))
    print("%0.4f < IBLK < %0.4f"%(np.min(grid[:,3]),np.max(grid[:,3])))
    
    print("")
    
    print("%0.4f < T < %0.4f"%(np.min(data[:,0]),np.max(data[:,0])))
    print("%0.4f < U < %0.4f"%(np.min(data[:,1]),np.max(data[:,1])))
    print("%0.4f < V < %0.4f"%(np.min(data[:,2]),np.max(data[:,2])))
    print("%0.4f < W < %0.4f"%(np.min(data[:,3]),np.max(data[:,3])))
    print("%0.4f < HRRPUV < %0.4f"%(np.min(data[:,4]),np.max(data[:,4])))
        
def getAbsoluteGrid(grids):
    mins = []
    maxs = []
    deltas = []
    for key in list(grids.keys()):
        xGrid = grids[key]['xGrid']
        yGrid = grids[key]['yGrid']
        zGrid = grids[key]['zGrid']
        mins.append([xGrid.min(), yGrid.min(), zGrid.min()])
        maxs.append([xGrid.max(), yGrid.max(), zGrid.max()])
        dx = np.round(xGrid[1,0,0]-xGrid[0,0,0],decimals=4)
        dy = np.round(yGrid[0,1,0]-yGrid[0,0,0],decimals=4)
        dz = np.round(zGrid[0,0,1]-zGrid[0,0,0],decimals=4)
        deltas.append([dx, dy, dz])
    absMins = np.min(mins,axis=0)
    absMaxs = np.max(maxs,axis=0)
    absDeltas = np.min(deltas, axis=0)
    
    Nx = int(np.round((absMaxs[0]-absMins[0])/absDeltas[0]) + 1)
    Ny = int(np.round((absMaxs[1]-absMins[1])/absDeltas[1]) + 1)
    Nz = int(np.round((absMaxs[2]-absMins[2])/absDeltas[2]) + 1)
    
    xs = np.linspace(absMins[0],absMaxs[0], Nx)
    ys = np.linspace(absMins[1],absMaxs[1], Ny)
    zs = np.linspace(absMins[2],absMaxs[2], Nz)
    
    xGrid_abs, yGrid_abs, zGrid_abs = np.meshgrid(xs, ys, zs)
    xGrid_abs = np.swapaxes(xGrid_abs, 0, 1)
    yGrid_abs = np.swapaxes(yGrid_abs, 0, 1)
    zGrid_abs = np.swapaxes(zGrid_abs, 0, 1)
    
    grid_abs = np.zeros((xGrid_abs.shape[0], xGrid_abs.shape[1], xGrid_abs.shape[2], 3))
    grid_abs[:,:,:,0] = xGrid_abs
    grid_abs[:,:,:,1] = yGrid_abs
    grid_abs[:,:,:,2] = zGrid_abs
    
    return grid_abs

def findSliceLocation(grid, data, axis, value, plot3d=False):
    (xGrid, yGrid, zGrid) = (grid[:,:,:,0], grid[:,:,:,1], grid[:,:,:,2])
    if axis == 1: 
        nGrid = xGrid
        nValues = nGrid[:, 0, 0]
    if axis == 2: 
        nGrid = yGrid
        nValues = nGrid[0, :, 0]
    if axis == 3: 
        nGrid = zGrid
        nValues = nGrid[0, 0, :]
    ind = np.argmin(abs(nValues-value))
    if axis == 1:
        x = np.squeeze(yGrid[ind,:,:])
        z = np.squeeze(zGrid[ind,:,:])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[ind,:,:,:])
        else:
            data2 = np.squeeze(data[ind,:,:])
    elif axis == 2:
        x = np.squeeze(xGrid[:,ind,:])
        z = np.squeeze(zGrid[:,ind,:])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:,ind,:,:])
        else:
            data2 = np.squeeze(data[:,ind,:])
    elif axis == 3:
        x = np.squeeze(xGrid[:,:,ind])
        z = np.squeeze(yGrid[:,:,ind])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:,:,ind,:])
        else:
            data2 = np.squeeze(data[:,:,ind])
    if plot3d:
        T = np.squeeze(data2[:,:,0])
        U = np.squeeze(data2[:,:,1])
        V = np.squeeze(data2[:,:,2])
        W = np.squeeze(data2[:,:,3])
        HRR = np.squeeze(data2[:,:,4])
        return x, z, T, U, V, W, HRR
    else:
        return x, z, data2

def plotSlice(x, z, data_slc, axis,
              cmap=None, figsize=None, fs=16,
              qnty_mn=None, qnty_mx=None,
              levels=None, cbarticks=None, clabel=None):
    if cmap == None: cmap = buildSMVcolormap()
    if cmap == 'SMV': cmap = buildSMVcolormap()
    if levels == None: levels = 100
    if cbarticks is None: cbarticks = np.linspace(qnty_mn, qnty_mx, 10)
    if figsize == None:
        xrange = x.max()-x.min()
        zrange = z.max()-z.min()
        if zrange > xrange:
            fig, ax = plt.subplots(1, 1, figsize=(12*xrange/zrange,12))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(12,12*zrange/xrange))
    else:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    if qnty_mn == None: qnty_mn = np.nanmin(data_slc)
    if qnty_mx == None: qnty_mx = np.nanmax(data_slc)
    im = ax.contourf(x, z, data_slc[:,:], cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, levels=np.linspace(qnty_mn,qnty_mx,levels), extend='both')
    cbar = fig.colorbar(im, cmap=cmap, extend='both', ticks=cbarticks)
    cbar.ax.tick_params(labelsize=fs)
    if clabel is not None: cbar.set_label(clabel, fontsize=fs)
    if abs(axis) == 1: ax.set_xlabel('y (m)', fontsize=fs)
    if (abs(axis) == 2) or (abs(axis) == 3): ax.set_xlabel('x (m)', fontsize=fs)
    if (abs(axis) == 1) or (abs(axis) == 2): ax.set_ylabel('z (m)', fontsize=fs)
    if (abs(axis) == 3): ax.set_ylabel('y (m)', fontsize=fs)
    ax.tick_params(labelsize=fs)
    return fig, ax

def readPlot3Ddata(chid, resultDir, time):
    xyzFiles = glob.glob("%s%s*.xyz"%(resultDir, chid))
    timeFiles = glob.glob(xyzFiles[0].replace('.xyz','*.q'))
    
    times = [x.split(chid)[1].split('.q')[0].split('_') for x in timeFiles]
    times = [float(x[2])+float(x[3])/100 for x in times]
    
    grids = []
    datas = []
    
    for xyzFile in xyzFiles:
        dataFile = buildDataFile(xyzFile, time)
        
        grid, gridHeader = readXYZfile(xyzFile)
        data, dataHeader = readP3Dfile(dataFile)
        (nx, ny, nz) = (dataHeader[0], dataHeader[1], dataHeader[2])
        
        printExtents(grid, data)
        
        xGrid, yGrid, zGrid = rearrangeGrid(grid, dataHeader)
        data = np.reshape(data, (nx, ny, nz, 5), order='F')
        grids.append([xGrid, yGrid, zGrid])
        datas.append(data)

    grid_abs = getAbsoluteGrid(grids)
    (xGrid_abs, yGrid_abs, zGrid_abs) = (grid_abs[:,:,:,0], grid_abs[:,:,:,1], grid_abs[:,:,:,2])
    data_abs = np.zeros((xGrid_abs.shape[0], xGrid_abs.shape[1], xGrid_abs.shape[2], 5))
    data_abs[:,:,:,:] = np.nan
    
    for grid, data in zip(grids, datas):
        (xGrid, yGrid, zGrid) = (grid[0], grid[1], grid[2])
        
        xloc = np.where(abs(xGrid_abs-xGrid[0,0,0]) == 0)[0][0]
        yloc = np.where(abs(yGrid_abs-yGrid[0,0,0]) == 0)[1][0]
        zloc = np.where(abs(zGrid_abs-zGrid[0,0,0]) == 0)[2][0]
        
        (NX, NY, NZ) = np.shape(xGrid)
        
        data_abs[xloc:xloc+NX, yloc:yloc+NY, zloc:zloc+NZ,:] = data
    return grid_abs, data_abs

def readSLCF3Ddata(chid, resultDir, quantityToExport, time=None, dt=None):
    if '.zip' in resultDir:
        xyzFiles = getFileListFromZip(resultDir, chid, 'xyz')
    else:
        xyzFiles = glob.glob("%s%s*.xyz"%(resultDir, chid))
    grids = defaultdict(bool)
    
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid, gridHeader)
        
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
        if '.zip' in resultDir:
            slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s_*.sf"%(resultDir, meshStr))
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        times3D = []
        datas3D = []
        lims3D = []
        for slcfFile in slcfFiles:
            timesSLCF = readSLCFtimes(slcfFile)
            times = []
            f = zopen(slcfFile)

            #with open(slcfFile, 'rb') as f:
            quantity, shortName, units, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f)
            # Check if slice is correct quantity
            correctQuantity = (quantity == quantityToExport)
            # Check if slice is 2-dimensional
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            
            if correctQuantity and threeDimSlice:
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                # Check if slice is 3-D
                print(slcfFile, quantity, shortName, units, iX, eX, iY, eY, iZ, eZ)
                if time == None:
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, len(timesSLCF)))
                    for i in range(0, len(timesSLCF)):
                        t, data = readNextTime(f, NX, NY, NZ)
                        data = np.reshape(data, (NX+1, NY+1, NZ+1),order='F')
                        datas2[:,:,:,i] = data
                    times = timesSLCF
                elif (time != None) and (dt == None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF-time))
                    f.seek(i*4*(5+(NX+1)*(NY+1)*(NZ+1)),1)
                    t, data = readNextTime(f, NX, NY, NZ)
                    data = np.reshape(data, (NX+1, NY+1, NZ+1),order='F')
                    datas2[:,:,:,0] = data
                    times = [timesSLCF[i]]
                elif (time != None) and (dt != None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF-(time-dt/2)))
                    j = np.argmin(abs(timesSLCF-(time+dt/2)))
                    f.seek(i*4*(5+(NX+1)*(NY+1)*(NZ+1)),1)
                    data = True
                    for ii in range(i,j+1):
                        if data is not False:
                            t, data = readNextTime(f, NX, NY, NZ)
                            data = np.reshape(data, (NX+1, NY+1, NZ+1), order='F')
                            datas2[:,:,:,0] = datas2[:,:,:,0] + data
                    if j-i > 0: datas2[:,:,:,0] = datas2[:,:,:,0] / (j-i)
                    times = [timesSLCF[i]]
                lims3D.append([iX, eX, iY, eY, iZ, eZ])
                datas3D.append(datas2)
                times3D.append(times)
            f.close()
        grids[meshStr]['datas3D'] = datas3D
        grids[meshStr]['lims3D'] = lims3D
        
    grid_abs = getAbsoluteGrid(grids)
    (xGrid_abs, yGrid_abs, zGrid_abs) = (grid_abs[:,:,:,0], grid_abs[:,:,:,1], grid_abs[:,:,:,2])
    tInd = grids[list(grids.keys())[0]]['datas3D'][0].shape[3]
    data_abs = np.zeros((xGrid_abs.shape[0], xGrid_abs.shape[1], xGrid_abs.shape[2], tInd))
    data_abs[:,:,:,:] = np.nan
    for key in list(grids.keys()):
        (xGrid, yGrid, zGrid) = (grids[key]['xGrid'], grids[key]['yGrid'], grids[key]['zGrid'])
        (datas3D, lims3D) = (grids[key]['datas3D'], grids[key]['lims3D'])
        for data, lim, times in zip(datas3D, lims3D, times3D):
            xloc_mn = np.where(np.isclose(abs(xGrid_abs-xGrid[lim[0],0,0]),0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xGrid_abs-xGrid[lim[1],0,0]),0, atol=1e-04))[0][0]
            yloc_mn = np.where(np.isclose(abs(yGrid_abs-yGrid[0,lim[2],0]),0, atol=1e-04))[1][0]
            yloc_mx = np.where(np.isclose(abs(yGrid_abs-yGrid[0,lim[3],0]),0, atol=1e-04))[1][0]
            zloc_mn = np.where(np.isclose(abs(zGrid_abs-zGrid[0,0,lim[4]]),0, atol=1e-04))[2][0]
            zloc_mx = np.where(np.isclose(abs(zGrid_abs-zGrid[0,0,lim[5]]),0, atol=1e-04))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            (ANX, ANY, ANZ) = (xloc_mx-xloc_mn+1, yloc_mx-yloc_mn+1, zloc_mx-zloc_mn+1)
            if (NX != ANX) or (NY != ANY) or (NZ != ANZ):
                (x, y, z) = (xGrid[lim[0]:lim[1]+1, 0, 0], yGrid[0, lim[2]:lim[3]+1, 0], zGrid[0, 0, lim[4]:lim[5]+1])
                
                xi = grid_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,0].flatten()
                yi = grid_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,1].flatten()
                zi = grid_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,2].flatten()
                
                (x, y, z) = (np.round(x, decimals=4), np.round(y, decimals=4), np.round(z, decimals=4))
                (xi, yi, zi) = (np.round(xi, decimals=4), np.round(yi, decimals=4), np.round(zi, decimals=4))
                
                (xi[xi < np.min(x)], xi[xi > np.max(x)]) = (np.min(x), np.max(x))
                (yi[yi < np.min(y)], yi[yi > np.max(y)]) = (np.min(y), np.max(y))
                (zi[zi < np.min(z)], zi[zi > np.max(z)]) = (np.min(z), np.max(z))
                
                tmpGrid = np.array([xi, yi, zi]).T
                for i in range(0, NT):
                    interpolator = scpi.RegularGridInterpolator((x, y, z), data[:, :, :, i])
                    data2 = interpolator(tmpGrid)
                    data2 = np.reshape(data2, (ANX, ANY, ANZ), order='C')
                    data_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,i] = data2
            else:
                data_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,:] = data
        
    return grid_abs, data_abs, times3D[0]





def readSLCF2Ddata(chid, resultDir, quantityToExport, time=None, dt=None):
    if '.zip' in resultDir:
        xyzFiles = getFileListFromZip(resultDir, chid, 'xyz')
    else:
        xyzFiles = glob.glob("%s%s*.xyz"%(resultDir, chid))
    grids = defaultdict(bool)
    
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid, gridHeader)
        
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
        if '.zip' in resultDir:
            slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s_*.sf"%(resultDir, meshStr))
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        datas2D = []
        lims2D = []
        coords2D = []
        times2D = []
        for slcfFile in slcfFiles:
            timesSLCF = readSLCFtimes(slcfFile)
            times = []
            f = zopen(slcfFile)
            
            quantity, shortName, units, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f)
            # Check if slice is correct quantity
            correctQuantity = (quantity == quantityToExport)
            # Check if slice is 2-dimensional
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            
            if correctQuantity and not threeDimSlice:
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                print("2-D slice:", slcfFile)
                if time == None:
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, len(timesSLCF)))
                    for i in range(0, len(timesSLCF)):
                        t, data = readNextTime(f, NX, NY, NZ)
                        data = np.reshape(data, (NX+1, NY+1, NZ+1), order='F')
                        datas2[:, :, :, i] = data
                    times = timesSLCF
                elif (time != None) and (dt == None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF-time))
                    f.seek(i*4*(5+(NX+1)*(NY+1)*(NZ+1)),1)
                    t, data = readNextTime(f, NX, NY, NZ)
                    data = np.reshape(data, (NX+1, NY+1, NZ+1), order='F')
                    datas2[:, :, :, 0] = data
                    times = [timesSLCF[i]]
                elif (time != None) and (dt != None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF-(time-dt/2)))
                    j = np.argmin(abs(timesSLCF-(time+dt/2)))
                    f.seek(i*4*(5+(NX+1)*(NY+1)*(NZ+1)),1)
                    for ii in range(i,j+1):
                        t, data = readNextTime(f, NX, NY, NZ)
                        data = np.reshape(data, (NX+1, NY+1, NZ+1),order='F')
                        datas2[:,:,:,0] = datas2[:,:,:,0] + data
                        
                    if j-i > 0: datas2[:,:,:,0] = datas2[:,:,:,0] / (j-i)
                    times = [timesSLCF[i]]
                lims2D.append([iX, eX, iY, eY, iZ, eZ])
                datas2D.append(datas2)
                coords2D.append([xGrid[iX, iY, iZ], yGrid[iX, iY, iZ], zGrid[iX, iY, iZ]])
                times2D.append(np.array(times))
            f.close()
        grids[meshStr]['datas2D'] = datas2D
        grids[meshStr]['lims2D'] = lims2D
        grids[meshStr]['coords2D'] = coords2D
        
    grid_abs = getAbsoluteGrid(grids)
    (xGrid_abs, yGrid_abs, zGrid_abs) = (grid_abs[:,:,:,0], grid_abs[:,:,:,1], grid_abs[:,:,:,2])
    tInd = grids[list(grids.keys())[0]]['datas2D'][0].shape[3]
    data_abs = np.zeros((xGrid_abs.shape[0], xGrid_abs.shape[1], xGrid_abs.shape[2], tInd))
    data_abs[:,:,:,:] = np.nan
    for key in list(grids.keys()):
        (xGrid, yGrid, zGrid) = (grids[key]['xGrid'], grids[key]['yGrid'], grids[key]['zGrid'])
        (datas2D, lims2D, coords2D) = (grids[key]['datas2D'], grids[key]['lims2D'], grids[key]['coords2D'])
        
        for data, coord in zip(datas2D, coords2D):
            xloc = np.where(np.isclose(abs(xGrid_abs-coord[0]),0, atol=1e-04))[0][0]
            yloc = np.where(np.isclose(abs(yGrid_abs-coord[1]),0, atol=1e-04))[1][0]
            zloc = np.where(np.isclose(abs(zGrid_abs-coord[2]),0, atol=1e-04))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            data_abs[xloc:xloc+NX, yloc:yloc+NY, zloc:zloc+NZ,:NT] = data
        
    return grid_abs, data_abs, times2D[0]







def extractPoint(point, grid, data):
    ind = np.argmin(np.sum(abs(grid-point),axis=3).flatten())
    ind = np.unravel_index(ind, grid[:,:,:,0].shape)
    x = grid[ind[0],ind[1],ind[2],0]
    y = grid[ind[0],ind[1],ind[2],1]
    z = grid[ind[0],ind[1],ind[2],2]
    d = data[ind[0],ind[1],ind[2],:]
    err = abs(point[0]-x)+abs(point[1]-y)+abs(point[2]-z)
    if err > 0.25:
        print("Warning, error for point %0.4f, %0.4f, %0.4f, is %0.4f."%(x, y, z, err))
    return d

def readNextTime(f, NX, NY, NZ):
    _ = np.frombuffer(f.read(8), dtype=np.float32)
    time = np.frombuffer(f.read(4), dtype=np.float32)
    _ = np.frombuffer(f.read(8), dtype=np.float32)
    try:
        data = np.frombuffer(f.read((NX+1)*(NY+1)*(NZ+1)*4), dtype=np.float32)
    except:
        data = False
    return time, data

def readSLCFheader(f):
    data = f.read(142)
    header = data[:110]
    size = struct.unpack('>iiiiii', data[115:139])
    tmp = header.split(b'\x1e')
    quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
    shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
    units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
    
    iX, eX, iY, eY, iZ, eZ = size
    
    return quantity, shortName, units, iX, eX, iY, eY, iZ, eZ

def readSLCFtimes(file):
    f = zopen(file)
    quantity, shortName, units, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f)
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    data = f.read()
    f.close()
    if len(data) % 4 == 0:
        fullFile = np.frombuffer(data, dtype=np.float32)
    else:
        print(len(data))
        remainder = -1*int(len(data) % 4)
        print(len(data[:remainder]))
        fullFile = np.frombuffer(data[:remainder], dtype=np.float32)
    times = fullFile[2::(NX+1)*(NY+1)*(NZ+1)+5]
    return times

def readSLCFquantities(chid, resultDir):
    if '.zip' in resultDir:
        slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
        zip = zipfile.ZipFile(resultDir, 'r')
    else:
        slcfFiles = glob.glob("%s%s_*.sf"%(resultDir, chid))
    quantities = []
    dimensions = []
    for file in slcfFiles:
        if '.zip' in resultDir:
            f = zip.open(file.split("%s%s"%('.zip',os.sep))[1])
        else:
            f = open(file, 'rb')
        quantity, shortName, units, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f)
        quantities.append(quantity)
        dimensions.append([iX, eX, iY, eY, iZ, eZ])
        f.close()
    if '.zip' in resultDir:
        zip.close()
    return quantities, slcfFiles, dimensions

def buildQTYstring(chid, resultDir, quantity):
    quantities, slcfFiles = readSLCFquantities(chid, resultDir)
    quantitiesCheck = [True if quantity == x else False for x in quantities]
    inds = np.where(quantitiesCheck)[0]
    if len(inds) == 0:
        print("Quantity %s unknown."%(quantity))
        print("Known quantities:")
        for qnty in sorted(set(quantities)):
            print(qnty)
    else:
        ind = inds[0]
    quantityStr = slcfFiles[ind].split('.sf')[0].split('_')[-1]
    return quantityStr

def getLimsFromGrid(grid):
    (xGrid, yGrid, zGrid) = (grid[:,:,:,0], grid[:,:,:,1], grid[:,:,:,2])
    
    (xmn, xmx) = (xGrid[:,0,0].min(), xGrid[:,0,0].max())
    (ymn, ymx) = (yGrid[0,:,0].min(), yGrid[0,:,0].max())
    (zmn, zmx) = (zGrid[0,0,:].min(), zGrid[0,0,:].max())
    
    return [xmn, xmx, ymn, ymx, zmn, zmx]

def visualizePlot3D(x, z, T, U, V, W, HRR):
    cmap = buildSMVcolormap()
    
    xrange = x.max()-x.min()
    zrange = z.max()-z.min()
    
    if zrange > xrange:
        plt.figure(figsize=(12*xrange/zrange,12))
    else:
        plt.figure(figsize=(12,12*zrange/xrange))
    (qnty_mn, qnty_mx) = (np.nanmin(T), np.nanmax(T))
    qnty_mn = 20
    qnty_mx = 370
    plt.contourf(x, z, T, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, levels=np.linspace(qnty_mn, qnty_mx,100), extend='both')
    plt.colorbar()
