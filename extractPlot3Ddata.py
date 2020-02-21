# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 09:35:22 2019

@author: JHodges
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import zipfile
import os
import struct
import scipy.interpolate as scpi
from collections import defaultdict

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

def buildSMVcolormap():
    from matplotlib.colors import ListedColormap
    newcmp = np.zeros((256,4))
    
    colors = np.array([[0,0,1,1],
              [0,1,1,1],
              [0,1,0,1],
              [1,1,0,1],
              [1,0,0,1],])
    colorInds = np.array([0, 64, 128, 192, 255])
    
    j = 0
    for i in range(0,len(newcmp)):
        if i == colorInds[j]:
            newcmp[i,:] = colors[j,:]
            j = j + 1
        else:
            m = (colors[j,:]-colors[j-1,:])/(colorInds[j]-colorInds[j-1])
            newcmp[i] = colors[j-1,:]+m*(i-colorInds[j-1])
    cmap = ListedColormap(newcmp)
    
    return cmap

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
    if axis == 0: nGrid = xGrid
    if axis == 1: nGrid = yGrid
    if axis == 2: nGrid = zGrid
    inds = np.where(np.isclose(nGrid, value, atol=1e-04))
    indDiff = [np.max(x)-np.min(x) for x in inds]
    indCheck = [True if x == 0 else False for x in indDiff]
    ind = inds[np.where(indCheck)[0][0]][0]
    if axis == 0:
        x = np.squeeze(yGrid[ind,:,:])
        z = np.squeeze(zGrid[ind,:,:])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[ind,:,:,:])
        else:
            data2 = np.squeeze(data[ind,:,:])
    elif axis == 1:
        x = np.squeeze(xGrid[:,ind,:])
        z = np.squeeze(zGrid[:,ind,:])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:,ind,:,:])
        else:
            data2 = np.squeeze(data[:,ind,:])
    elif axis == 2:
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
    if axis == 0: ax.set_xlabel('y (m)', fontsize=fs)
    if (axis == 1) or (axis == 2): ax.set_xlabel('x (m)', fontsize=fs)
    if (axis == 0) or (axis == 1): ax.set_ylabel('z (m)', fontsize=fs)
    if (axis == 2): ax.set_ylabel('y (m)', fontsize=fs)
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

def getFileListFromZip(filename, chid, extension):
    filelist = []
    with zipfile.ZipFile(filename, 'r') as zip:
        for info in zip.infolist():
            if info.filename.split('.')[-1] == extension:
                if chid in info.filename:
                    filelist.append("%s%s%s"%(filename, os.sep, info.filename))
    return filelist

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
        
        datas3D = []
        datas2D = []
        lims2D = []
        lims3D = []
        coords2D = []
        for slcfFile in slcfFiles:
            timesSLCF = readSLCFtimes(slcfFile)
            times = []
            f = zopen(slcfFile)

            #with open(slcfFile, 'rb') as f:
            quantity, shortName, units, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f)
            # Check if slice is correct quantity
            if (quantity == quantityToExport):
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                # Check if slice is 3-D
                if (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0):
                    print(slcfFile, quantity, shortName, units, iX, eX, iY, eY, iZ, eZ)
                    if time == None:
                        datas2 = np.zeros((NX+1, NY+1, NZ+1, len(timesSLCF)))
                        for i in range(0, len(timesSLCF)):
                            t, data = readNextTime(f, NX, NY, NZ)
                            data = np.reshape(data, (NX+1, NY+1, NZ+1),order='F')
                            datas2[:,:,:,i] = data
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
                else:
                    '''
                    print("2-D slice:", slcfFile)
                    if time == None:
                        datas2 = np.zeros((NX+1, NY+1, NZ+1, len(timesSLCF)))
                        for i in range(0, len(timesSLCF)):
                            t, data = readNextTime(f, NX, NY, NZ)
                            data = np.reshape(data, (NX+1, NY+1, NZ+1),order='F')
                            datas2[:,:,:,i] = data
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
                        for ii in range(i,j+1):
                            t, data = readNextTime(f, NX, NY, NZ)
                            data = np.reshape(data, (NX+1, NY+1, NZ+1),order='F')
                            datas2[:,:,:,0] = datas2[:,:,:,0] + data
                            
                        if j-i > 0: datas2[:,:,:,0] = datas2[:,:,:,0] / (j-i)
                        times = [timesSLCF[i]]
                    lims2D.append([iX, eX, iY, eY, iZ, eZ])
                    datas2D.append(datas2)
                    coords2D.append([xGrid[iX, iY, iZ], yGrid[iX, iY, iZ], zGrid[iX, iY, iZ]])
                    '''
            f.close()
        grids[meshStr]['datas3D'] = datas3D
        grids[meshStr]['lims3D'] = lims3D
        grids[meshStr]['datas2D'] = datas2D
        grids[meshStr]['lims2D'] = lims2D
        grids[meshStr]['coords2D'] = coords2D
        
    grid_abs = getAbsoluteGrid(grids)
    (xGrid_abs, yGrid_abs, zGrid_abs) = (grid_abs[:,:,:,0], grid_abs[:,:,:,1], grid_abs[:,:,:,2])
    #print(xGrid_abs.shape)
    #print(len(datas))
    #print(datas[0].shape)
    #print(np.where(np.isnan(datas[0])))
    #print(lims)
    tInd = grids[list(grids.keys())[0]]['datas3D'][0].shape[3]
    data_abs = np.zeros((xGrid_abs.shape[0], xGrid_abs.shape[1], xGrid_abs.shape[2], tInd))
    data_abs[:,:,:,:] = np.nan
    for key in list(grids.keys()):
        (xGrid, yGrid, zGrid) = (grids[key]['xGrid'], grids[key]['yGrid'], grids[key]['zGrid'])
        (datas3D, lims3D) = (grids[key]['datas3D'], grids[key]['lims3D'])
        (datas2D, lims2D, coords2D) = (grids[key]['datas2D'], grids[key]['lims2D'], grids[key]['coords2D'])
        for data, lim in zip(datas3D, lims3D):
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
                    #print(x.shape, y.shape, z.shape, data.shape)
                    #print(xGrid[lim[0],0,0], xGrid[lim[1],0,0], yGrid[0,lim[2],0], yGrid[0,lim[3],0], zGrid[0,0,lim[4]], zGrid[0,0,lim[5]])
                    #print(xGrid_abs[xloc_mn,0,0], xGrid_abs[xloc_mx,0,0], yGrid_abs[0,yloc_mn,0], yGrid_abs[0,yloc_mx,0], zGrid_abs[0,0,zloc_mn], zGrid_abs[0,0,zloc_mx])
                    interpolator = scpi.RegularGridInterpolator((x, y, z), data[:, :, :, i])
                    data2 = interpolator(tmpGrid)
                    data2 = np.reshape(data2, (ANX, ANY, ANZ), order='C')
                    data_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,i] = data2
            else:
                data_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1,:] = data
            #data_abs[xloc:xloc+NX, yloc:yloc+NY, zloc:zloc+NZ,:] = data
        
        for data, coord in zip(datas2D, coords2D):
            xloc = np.where(np.isclose(abs(xGrid_abs-coord[0]),0, atol=1e-04))[0][0]
            yloc = np.where(np.isclose(abs(yGrid_abs-coord[1]),0, atol=1e-04))[1][0]
            zloc = np.where(np.isclose(abs(zGrid_abs-coord[2]),0, atol=1e-04))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            data_abs[xloc:xloc+NX, yloc:yloc+NY, zloc:zloc+NZ,:NT] = data
        
    return grid_abs, data_abs, times

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
    #data = np.fromfile(f,dtype=np.float32,count=(NX+1)*(NY+1)*(NZ+1))
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

def zopen(file):
    if '.zip' in file:
        zname = '%s.zip'%(file.split('.zip')[0])
        fname = file.split('.zip%s'%(os.sep))[1]
        zip = zipfile.ZipFile(zname, 'r')
        f = zip.open(fname)
    else:
        f = open(file, 'rb')
    return f

def readSLCFtimes(file):
    f = zopen(file)
    quantity, shortName, units, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f)
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    data = f.read()
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
    (qnty_mn, qnty_mx) = (np.nanmin(data_slc), np.nanmax(data_slc))
    qnty_mn = 20
    qnty_mx = 370
    plt.contourf(x, z, T, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, levels=np.linspace(qnty_mn,qnty_mx,100), extend='both')
    plt.colorbar()

if __name__ == "__main__":
    smvFile = "E:\\projects\\mineFireLearning\\mineValidation\\scenario06_v2\\Mine_fire_s6.smv"
    fdsFile = "E:\\projects\\mineFireLearning\\mineValidation\\scenario06_v2\\Mine_fire_s6.fds"
    
    meshNumber = 2
    time = 90
    dt = 10
    chid = 'Mine_fire_s6'
    quantity = 'TEMPERATURE'
    resultDir = "E:\\projects\\mineFireLearning\\mineValidation\\scenario06_v2\\"
    
    (axis, value) = (1, 15.3)
    #(axis, value) = (0, 3.0)
    (axis, value) = (2, 2.0)
    
    quantities = ['TEMPERATURE',
                  'U-VELOCITY',
                  'V-VELOCITY',
                  'W-VELOCITY',
                  'HRRPUV',
                  'OXYGEN VOLUME FRACTION',
                  'REAC_FUEL VOLUME FRACTION',
                  'SOOT VOLUME FRACTION',
                  'PRESSURE']
    
    quantities = ['TEMPERATURE']
    
    grids_abs, data_abs = readPlot3Ddata(chid, resultDir, time)
    
    x, z, T, U, V, W, HRR = findSliceLocation(grids_abs, data_abs, axis, value)
    assert False, "Stopped"
    
    #grid, T, times = readSLCF3Ddata(chid, resultDir, 'FAKENEWS')
    for time in np.linspace(0,100,11):
        for i in range(0,len(quantities)):
            quantity = quantities[i]
            grid, data, times = readSLCF3Ddata(chid, resultDir, quantity, time=time, dt=dt)
            if i == 0:
                datas = np.expand_dims(data, axis=4)
            else:
                data = np.expand_dims(data, axis=4)
                datas = np.concatenate((datas, data), axis=4)
    
        x, z, data_slc = findSliceLocation(grid, datas[:,:,:,:,0], axis, value, plot3d=False)
        #x2, z2, data2 = findSliceLocation(grid2, T2, axis, value, plot3d=False)
        
        size = 5
        
        cx = np.random.randint((size-1)/2,grid.shape[0]-(size-1)/2)
        cy = np.random.randint((size-1)/2,grid.shape[1]-(size-1)/2)
        cz = np.random.randint((size-1)/2,grid.shape[2]-(size-1)/2)
        
        (xmn, xmx) = (int(cx-(size-1)/2), int(cx+(size-1)/2+1))
        (ymn, ymx) = (int(cy-(size-1)/2), int(cy+(size-1)/2+1))
        (zmn, zmx) = (int(cz-(size-1)/2), int(cz+(size-1)/2+1))
        
        #s_data = data_abs[xmn:xmx,ymn:ymx,zmn:zmx,:]
    
        print(cx, cy, cz)
        
        cmap = buildSMVcolormap()
        
        xrange = x.max()-x.min()
        zrange = z.max()-z.min()
        
        if zrange > xrange:
            plt.figure(figsize=(12*xrange/zrange,12))
        else:
            plt.figure(figsize=(12,12*zrange/xrange))
        (qnty_mn, qnty_mx) = (np.nanmin(data_slc), np.nanmax(data_slc))
        qnty_mn = 20
        qnty_mx = 370
        plt.contourf(x, z, data_slc[:,:, 500], cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, levels=np.linspace(qnty_mn,qnty_mx,100), extend='both')
        plt.colorbar()
    