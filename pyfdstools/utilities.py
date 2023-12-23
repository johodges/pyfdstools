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
# This script has utilities used to throughout the package for tasks
# associated with visualization and filtering.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import numpy as np
import matplotlib.colors as pltc
import scipy.spatial as scsp
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
import os
import zipfile
import glob
import struct
from collections import defaultdict
from .colorSchemes import getVTcolors

def timeAverage2(data, times, window):
    # Reshape data into array for time averaging
    sz = data.shape
    pts = np.product(sz[:-1])
    data2 = np.reshape(data, (pts,sz[-1]))
    
    # Time average the array
    data3 = np.zeros_like(data2) #data2.copy()
    dt = window/2
    tmax = len(times)
    tind = np.where(times-window > 0)[0][0]
    data3[:, 0] = np.nanmean(data2[:, :tind], axis=1)
    for i in range(1, tmax):
        data_dt = times[i]-times[i-1]
        data3[:,i] = (data3[:, i-1]*(dt-data_dt)+ data2[:, i]*data_dt)/dt
    data4 = np.reshape(data3, sz)
    return data4

def timeAverage(data, times, window, outdt=-1, smoothEnds=False, queryTime=-1):
    tmax = np.nanmax(times)
    tmin = np.nanmin(times)
    if window > (tmax-tmin):
        print("Warning, windowSize > time interval. Not time averaging.")
        return data, times
    
    dts = times[1:] - times[:-1]
    dt = np.nanmin(dts[dts > 0])
    
    if queryTime > 0:
        tmin = np.floor((queryTime-window/2)/dt)*dt
        tmax = np.ceil((queryTime+window/2)/dt)*dt
        data2 = np.zeros((data.shape[0], data.shape[1], 1))
        t1 = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
        for i in range(0, data.shape[0]):
            for j in range(0, data.shape[1]):
                v1 = np.interp(t1, times, data[i, j, :])
                data2[i, j, 0] = np.nanmean(v1)
        return data2, queryTime
    
    tmin = np.floor(tmin/dt)*dt
    tmax = np.ceil(tmax/dt)*dt
    
    t1 = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
    N = int(np.round(window/dt))
    N2 = int(N/2)
    f = np.zeros((N)) + 1
    f = f/f.sum()
    
    if data.shape[-1] < times.shape[-1]:
        print("Warning, data shape is less than timesteps, truncating.")
        times = times[:data.shape[-1]]
    
    data2 = np.zeros((data.shape[0], data.shape[1], int(len(t1)-N+1)))
    for i in range(0, data.shape[0]):
        for j in range(0, data.shape[1]):
            v1 = np.interp(t1, times, data[i, j, :])
            data2[i, j, :] = np.convolve(v1, f, mode='valid')
            if smoothEnds:
                for k in range(0, N2):
                    data2[i,j,k] = np.nanmean(v1[k:k+N2])
                for k in range(data2.shape[2]-N, data2.shape[2]):
                    data2[i,j,k] = np.nanmean(v1[k:k+N2])
            else:
                data2[i, j, :N2] = v1[:N2]
                data2[i ,j,-N2:] = v1[-N2:]
    
    if outdt > 0:
        tmin = np.floor(tmin/outdt)*outdt
        tmax = np.ceil(tmax/outdt)*outdt
        t2 = np.linspace(tmin, tmax, int((tmax-tmin)/outdt + 1))
        data3 = np.zeros((data.shape[0], data.shape[1], t2.shape[0]))
        for i in range(0, data.shape[0]):
            for j in range(0, data.shape[1]):
                data3[i, j, :] = np.interp(t2, t1[:data2.shape[2]], data2[i, j, :])
        return data3, t2
    else:
        return data2, t1[:data2.shape[2]]
        
    
    return data2, t1[:data2.shape[2]]

def kalmanFilter(z, Q=1e-5, R=0.5**2):
    # This subroutine applies a kalman filter to an input set of data.
    #
    #Inputs:
    #    z: series of data to be filtered
    #    Q: process variance
    #    R: measurement variance
    #Outputs:
    #    xhat: filtered series of data
    
    # intial parameters
    sz = z.shape[0] # size of array
    
    # allocate space for arrays
    xhat=np.zeros(sz)      # a posteri estimate of x
    P=np.zeros(sz)         # a posteri error estimate
    xhatminus=np.zeros(sz) # a priori estimate of x
    Pminus=np.zeros(sz)    # a priori error estimate
    K=np.zeros(sz)         # gain or blending factor
    
    # intial guesses
    xhat[0] = z[0]
    P[0] = 1.0
    
    for k in range(1,sz):
        # time update
        xhatminus[k] = xhat[k-1]
        Pminus[k] = P[k-1]+Q
    
        # measurement update
        K[k] = Pminus[k]/( Pminus[k]+R )
        xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
        P[k] = (1-K[k])*Pminus[k]
    
    return xhat

def smvVisual(obstructions,surfaces,namespace,fs=16,fig=None,ax=None,
              limits=[0,15,0,8,0,5]):
    if fig is None: fig = plt.figure(figsize=(12,12))
    if ax is None: ax = a3.Axes3D(fig)
    
    for obst in obstructions:
        pts, colors = getPtsFromObst(obst,surfaces)
        print(pts)
        print(colors)
        for pt, color in zip(pts,colors):
            f = a3.art3d.Poly3DCollection(pt)
            f.set_color(color)
            f.set_edgecolor('k')
            #f.set_alpha(1)
            ax.add_collection3d(f)
    plt.xlim(limits[0],limits[1])
    plt.ylim(limits[2],limits[3])
    ax.set_zlim(limits[4],limits[5])
    plt.xlabel('x (m)',fontsize=fs)
    plt.ylabel('y (m)',fontsize=fs)
    ax.set_zlabel('z (m)',fontsize=fs)
    plt.tick_params(labelsize=fs)
    plt.savefig('%s_smvvisual.png'%(namespace),dpi=300)
    
    return fig, ax

def buildSMVgeometry(file):
    with open(file,'r') as f:
        lines = f.readlines()
    inds = []
    for i in range(0,len(lines)):
        if lines[i][0] != ' ' and lines[i][0] != '-':
            inds.append(i)
    surfaces = []
    obstructions = []
    for ind in inds:
        if 'SURFACE' in lines[ind]:
            sname = ' '.join(lines[ind+1].split())
            (Tign,eps) = (lines[ind+2].split()[0],lines[ind+2].split()[1])
            (stype,t_width,t_height) = (lines[ind+3].split()[0],lines[ind+3].split()[1],lines[ind+3].split()[2])
            (c1,c2,c3,c4) = (lines[ind+3].split()[3],lines[ind+3].split()[4],lines[ind+3].split()[5],lines[ind+3].split()[6])
            surfaces.append([sname,Tign,eps,stype,t_width,t_height,c1,c2,c3,c4])
        if 'OBST' in lines[ind] and '_OBST' not in lines[ind]:
            nObst = int(lines[ind+1].split()[0])
            for i in range(0,nObst):
                obst = [float(x) for x in lines[ind+i+2].split()]
                obstructions.append(obst)
    return surfaces, obstructions

def getPtsFromObst(obst,surfaces):
    pts = []
    colors = []
    pts = np.array([[obst[0],obst[2],obst[4]],
                   [obst[0],obst[2],obst[5]],
                   [obst[1],obst[2],obst[5]],
                   [obst[1],obst[2],obst[4]],
                   [obst[0],obst[3],obst[5]],
                   [obst[1],obst[3],obst[5]],
                   [obst[1],obst[3],obst[4]],
                   [obst[0],obst[3],obst[4]]])
    
    # y-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[1],obst[2],obst[4]),
    #            (obst[1],obst[2],obst[5]),(obst[0],obst[2],obst[5])])
    surf = surfaces[int(obst[7])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # y-positive surface
    #pts.append([(obst[0],obst[3],obst[4]),(obst[1],obst[3],obst[4]),
    #            (obst[1],obst[3],obst[5]),(obst[0],obst[3],obst[5])])
    surf = surfaces[int(obst[8])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # x-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[0],obst[2],obst[5]),
    #            (obst[0],obst[3],obst[5]),(obst[0],obst[3],obst[4])])
    surf = surfaces[int(obst[9])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # x-positive surface
    #pts.append([(obst[1],obst[2],obst[4]),(obst[1],obst[2],obst[5]),
    #            (obst[1],obst[3],obst[5]),(obst[1],obst[3],obst[4])])
    surf = surfaces[int(obst[10])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # z-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[1],obst[2],obst[4]),
    #            (obst[1],obst[3],obst[4]),(obst[0],obst[3],obst[4])])
    surf = surfaces[int(obst[11])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # z-positive surface
    #pts.append([(obst[0],obst[2],obst[5]),(obst[1],obst[2],obst[5]),
    #            (obst[1],obst[3],obst[5]),(obst[0],obst[3],obst[5])])
    surf = surfaces[int(obst[12])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    
    return pts, colors

def maxValueCSV(times, mPts, names, namespace):
    '''  mPts rows correlated to times, columns correlated to different groups. '''
    numberOfGroups = mPts.shape[1]
    header = 'Time,'
    for i in range(0,numberOfGroups):
        name = names[i].replace(',','_')
        header = header+name+','
    header = header[:-1]+'\n'
    data = np.append(np.reshape(times,(times.shape[0],1)),mPts,axis=1)
    csvName = '%s.csv'%(namespace)
    print("Saving max value csv to %s"%(csvName))
    np.savetxt(csvName, data, delimiter=',', header=header)
    return '%s.csv'%(namespace)   

def maxValuePlot(times, mPts, names, namespace, fs=16, lw=3, pcs=None, vName='',
                 yticks=None, xticks=None):
    '''  mPts rows correlated to times, columns correlated to different groups. '''
    numberOfGroups = mPts.shape[1]
    if pcs is None:
        pcs = getVTcolors()
        if len(pcs) < numberOfGroups: pcs = getPlotColors(numberOfGroups)
    fig = plt.figure(figsize=(12,8))
    for i in range(0,numberOfGroups):
        plt.plot(times,mPts[:,i],color=pcs[i],label=names[i],linewidth=lw)
    if yticks is not None: plt.yticks(yticks)
    if xticks is not None: plt.xticks(xticks)
    plt.legend(fontsize=fs)
    plt.xlabel('time (s)',fontsize=fs)
    plt.ylabel('%s'%(vName),fontsize=fs)
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    figName = '%s_maxTPlot.png'%(namespace)
    print("Saving max value figure to %s"%(figName))
    plt.savefig(figName, dpi=300)
    plt.show()
    return fig

def getPlotColors(numberOfGroups):
    pcs = []
    print(numberOfGroups)
    for i in range(0,numberOfGroups):
        v = np.random.rand(3)
        su = np.sum(v)
        while (su > 2.7) or (su < 0.3):
            v = np.random.rand(3)
            su = np.sum(v)
        tmp = pltc.rgb2hex(v)
        pcs.append(tmp)
    return pcs

def pointsFromXB(XB,extend=[0,0,0]):
    ''' This routine builds a list of XYZ points from an obstruction XB
    
    Inputs:
        XB: Septuplet containing [xmin, xmax, ymin, ymax, zmin, zmax]
        extend: Float array containing amount to extend polygon along each axis
    Outputs:
        pts: List of corner points
    '''
    pts = [[XB[0]-extend[0],XB[2]-extend[1],XB[4]-extend[2]],
           [XB[0]-extend[0],XB[2]-extend[1],XB[5]+extend[2]],
           [XB[0]-extend[0],XB[3]+extend[1],XB[4]-extend[2]],
           [XB[0]-extend[0],XB[3]+extend[1],XB[5]+extend[2]],
           [XB[1]+extend[0],XB[2]-extend[1],XB[4]-extend[2]],
           [XB[1]+extend[0],XB[2]-extend[1],XB[5]+extend[2]],
           [XB[1]+extend[0],XB[3]+extend[1],XB[4]-extend[2]],
           [XB[1]+extend[0],XB[3]+extend[1],XB[5]+extend[2]]]
    return pts

def in_hull(p, hull):
    if not isinstance(hull,scsp.Delaunay):
        hull = scsp.Delaunay(hull)
    return hull.find_simplex(p)>=0

def pts2polygons(groups):
    '''
    Build polygons from series of points.
    '''
    polygons = []
    
    for group in groups:
        linkedPolygons = []
        for pts in group:
            try:
                linkedPolygons.append(scsp.ConvexHull(pts))
            except:
                print("Failed points:")
                print(pts)
        polygons.append(linkedPolygons)
        
    return polygons, len(polygons)

def getFileList(resultDir, chid, extension):
    if '.zip' in resultDir:
        files = getFileListFromZip(resultDir, chid, extension)
    else:
        path = os.path.join(resultDir, '%s*.%s'%(chid, extension))
        files = glob.glob(path)
    return files

def getFileListFromZip(filename, chid, extension):
    filelist = []
    with zipfile.ZipFile(filename, 'r') as zip:
        for info in zip.infolist():
            if info.filename.split('.')[-1] == extension:
                if chid in info.filename:
                    filelist.append("%s%s%s"%(filename, os.sep, info.filename))
    return filelist

def zreadlines(file):
    f = zopen(file, readtype='r')
    lines = f.readlines()
    if '.zip' in file:
        lines = [line.decode("utf-8").replace('\r','').replace('\n','') for line in lines]
    f.close()
    return lines

def zopen(file, readtype='rb'):
    if '.zip' in file:
        zname = '%s.zip'%(file.split('.zip')[0])
        fname = file.split('.zip%s'%(os.sep))[1]
        zip = zipfile.ZipFile(zname, 'r')
        f = zip.open(fname)
    else:
        f = open(file, readtype)
    return f


def getEndianness(resultDir, chid):
    endFiles = getFileList(resultDir, chid, 'end')
    if len(endFiles) == 0:
        #print('Unable to find endianness file, %s.'%(os.path.join(resultDir, chid+'.end')))
        #print("Assuming little-endian.")
        return "<"
        
    f = zopen(endFiles[0])
    data = f.read()
    f.close()
    if struct.unpack('<i',data[4:8])[0] == 1:
        return '<'
    elif struct.unpack('>i',data[4:8])[0] == 1:
        return ">"
    else:
        print("Unable to determine endianness from file, %s."%(os.path.join(resultDir, chid+'.end')))
        print("Assuming little-endian.")
        return "<"

def getDatatypeByEndianness(datatype1, endianness):
    if endianness == '>':
        datatype2 = np.dtype(datatype1).newbyteorder('>')
    elif endianness == '<':
        datatype2 = np.dtype(datatype1).newbyteorder('<')
    else:
        print("Endianness unknown, could not determine datatype")
    return datatype2


def getFileListFromResultDir(resultDir, chid, ext):
    if '.zip' in resultDir:
        fileList = getFileListFromZip(resultDir, chid, ext)
    else:
        fileList = glob.glob(os.path.join(resultDir, '%s*.%s'%(chid, ext)))
    return fileList
    

def getGridsFromXyzFiles(xyzFiles, chid):
    grids = defaultdict(bool)
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid)
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else mesh
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
    return grids


def rearrangeGrid(grid, iblock=False):
    """Builds a meshgrid based on grid array
    
    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Array containing float global coordinates
    
    Returns
    -------
    array(NX, NY, NZ)
        Array containing float global coordinates for x-axis
    array(NX, NY, NZ)
        Array containing float global coordinates for y-axis
    array(NX, NY, NZ)
        Array containing float global coordinates for z-axis
    """
    
    xs = np.unique(grid[:,0])
    ys = np.unique(grid[:,1])
    zs = np.unique(grid[:,2])
    xGrid, yGrid, zGrid = np.meshgrid(xs, ys, zs)
    xGrid = np.swapaxes(xGrid, 0, 1)
    yGrid = np.swapaxes(yGrid, 0, 1)
    zGrid = np.swapaxes(zGrid, 0, 1)
    
    if iblock:
        iGrid = np.zeros_like(xGrid)
        for i in range(0, xGrid.shape[0]):
            for j in range(0, xGrid.shape[1]):
                for k in range(0, xGrid.shape[2]):
                    dist = abs(xs-xGrid[i, j, k]) + abs(ys-yGrid[i, j, k]) + abs(zs-zGrid[i, j, k])
                    ind = np.argmin(dist)
                    iGrid[i, j, k] = grid[ind, 3]
        return xGrid, yGrid, zGrid, iGrid
    else:
        return xGrid, yGrid, zGrid


def readXYZfile(file):
    """Reads points from an xyz file
    
    This subroutine reads grid coordinates from an xyz file. Note,
    xyz file can be generated in FDS by adding WRITE_XYZ=.TRUE. in the
    &DUMP namelist.
    
    Parameters
    ----------
    file : str
        String containing the path to an xyz file or xyz file in an
        archive
    
    Returns
    -------
    array(NX, NY, NZ, 3)
        Array containing float global coordinates
    array()
        Array containing header information from xyz file
    """
    
    try:
        f = zopen(file)
    except FileNotFoundError:
        tmp = file[:-4].split('_')
        meshStr = str(int(tmp[-1]))
        file2 = '%s_%s.xyz'%('_'.join(tmp[:-1]), meshStr)
        f = zopen(file2)
    header = struct.unpack('<iiiiif', f.read(24))
    (nx, ny, nz) = (header[1], header[2], header[3])
    data = np.frombuffer(f.read(nx*ny*nz*4*4), dtype=np.float32)
    grid = np.reshape(data, (int(data.shape[0]/4), 4),order='F')
    f.close()
    return grid, header[1:-1]
    

def getAbsoluteGrid(grids, makeUniform=False):
    """Builds absolute grid from defaultdict of local grids
    
    Parameters
    ----------
    grids : defaultdict
        Dictionary containing arrays(NX, NY, NZ) for each grid
    makeUniform : bool
        Boolean flag indicating whether the absolute grid should be
        interpolated to be rectangular at the most resolved mesh.
    Returns
    -------
    array(NX, NY, NZ, 3)
        Array containing the absolute grid coordinates
    """
    
    if makeUniform:
        mins = []
        maxs = []
        deltas = []
        for key in list(grids.keys()):
            xGrid = grids[key]['xGrid']
            yGrid = grids[key]['yGrid']
            zGrid = grids[key]['zGrid']
            mins.append([xGrid.min(), yGrid.min(), zGrid.min()])
            maxs.append([xGrid.max(), yGrid.max(), zGrid.max()])
            dx = np.round(xGrid[1, 0, 0] - xGrid[0, 0, 0], decimals=4)
            dy = np.round(yGrid[0, 1, 0] - yGrid[0, 0, 0], decimals=4)
            dz = np.round(zGrid[0, 0, 1] - zGrid[0, 0, 0], decimals=4)
            deltas.append([dx, dy, dz])
        absMins = np.min(mins, axis=0)
        absMaxs = np.max(maxs, axis=0)
        absDeltas = np.min(deltas, axis=0)
        
        Nx = int(np.round((absMaxs[0] - absMins[0]) / absDeltas[0]) + 1)
        Ny = int(np.round((absMaxs[1] - absMins[1]) / absDeltas[1]) + 1)
        Nz = int(np.round((absMaxs[2] - absMins[2]) / absDeltas[2]) + 1)
        
        xs = np.linspace(absMins[0], absMaxs[0], Nx)
        ys = np.linspace(absMins[1], absMaxs[1], Ny)
        zs = np.linspace(absMins[2], absMaxs[2], Nz)
        
        xGrid_abs, yGrid_abs, zGrid_abs = np.meshgrid(xs, ys, zs)
        xGrid_abs = np.swapaxes(xGrid_abs, 0, 1)
        yGrid_abs = np.swapaxes(yGrid_abs, 0, 1)
        zGrid_abs = np.swapaxes(zGrid_abs, 0, 1)
        
        grid_abs = np.zeros((xGrid_abs.shape[0],
                             xGrid_abs.shape[1],
                             xGrid_abs.shape[2],
                             3))
        grid_abs[:, :, :, 0] = xGrid_abs
        grid_abs[:, :, :, 1] = yGrid_abs
        grid_abs[:, :, :, 2] = zGrid_abs
    else:
        all_xs = []
        all_ys = []
        all_zs = []
        for key in list(grids.keys()):
            xs = grids[key]['xGrid'][:, 0, 0]
            ys = grids[key]['yGrid'][0, :, 0]
            zs = grids[key]['zGrid'][0, 0, :]
            all_xs.extend(xs)
            all_ys.extend(ys)
            all_zs.extend(zs)
        
        abs_xs = np.unique(all_xs)
        abs_ys = np.unique(all_ys)
        abs_zs = np.unique(all_zs)
        xGrid_abs, yGrid_abs, zGrid_abs = np.meshgrid(abs_xs, abs_ys, abs_zs)
        
        xGrid_abs = np.swapaxes(xGrid_abs, 0, 1)
        yGrid_abs = np.swapaxes(yGrid_abs, 0, 1)
        zGrid_abs = np.swapaxes(zGrid_abs, 0, 1)
        
        grid_abs = np.zeros((xGrid_abs.shape[0],
                             xGrid_abs.shape[1],
                             xGrid_abs.shape[2],
                             3))
        grid_abs[:, :, :, 0] = xGrid_abs
        grid_abs[:, :, :, 1] = yGrid_abs
        grid_abs[:, :, :, 2] = zGrid_abs
    return grid_abs


def getTwoZone(z, val, lowtohigh=True):
    if lowtohigh:
        z = z[::-1]
        val = val[::-1]
        val_low = val[-1]
    else:
        val_low = val[0]
    H = z.max()
    H0 = z.min()
    tmpZ = np.linspace(0, H, 101)
    tmpV = np.interp(tmpZ, z, val)
    if np.isclose(np.nanmin(tmpV), np.nanmax(tmpV), atol=1e-1):
        zInt = H
        return val_low, val_low, zInt
    
    I1 = np.trapz(tmpV, tmpZ)
    I2 = np.trapz(1/tmpV, tmpZ)
    zInt = val_low*(I1*I2-H**2)/(I1+I2*val_low**2-2*val_low*H)
    
    zU = np.linspace(zInt, H, num=50)
    
    val_high_tmp = np.interp(zU, z, val)
    val_high = np.trapz(val_high_tmp, zU)/(H-zInt)
    
    zL = np.linspace(0, zInt, num=50)
    
    val_low_tmp = np.interp(zL, z, val)
    val_low = np.trapz(val_low_tmp, zL)/(zInt-H0)
    
    return val_low, val_high, zInt
