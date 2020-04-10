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
import os
import numpy as np
import yaml
from collections import defaultdict
import glob
from .utilities import zopen, in_hull, zreadlines, getFileList, pts2polygons
from .fdsFileOperations import fdsFileOperations
from .smokeviewParser import parseSMVFile


class fdspatch(object):
    """
    A class used to represent a patch from an FDS boundary file

    ...

    Attributes
    ----------
    data : int array(NX, NY, NT)
        array containing data from the patch
    lims : float array(6)
        Six component array containing X_min, X_max, Y_min, Y_max,
        Z_min, Z_max coordinates
    orientation : int
        Integer specifying the orientation of the patch
    x : array(NX, NY)
        array containing x-coordinates of the patch
    y : array(NX, NY)
        array containing y-coordinates of the patch
    z : array(NX, NY)
        array containing z-coordinates of the patch


    Methods
    -------
    append(data, iT)
        Adds data from one time stamp to the data attribute
    average(inds)
        Averages data from multiple time stamps
    buildSpace()
        Builds x, y, z coordinates
    extractPoints()
        Extracts points from data array
    """
    
    def __init__(self, NX, NY, NT, DS, OR):
        """
        Parameters
        ----------
        NX : int
            Number of cells in the x-direction
        NY : int
            Number of cells in the y-direction
        NT : int
            Number of time stamps
        DS : array(6)
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates
        OR : int
            Integer specifying the orientation of the patch
        """
        
        self.data = np.zeros((NX, NY, NT))
        self.lims = DS
        self.orientation = OR
        
        
    def append(self, data, iT):
        """Adds data from one time stamp to the data attribute
        
        Parameters
        ----------
        data : int array(NX, NY, NT)
            array containing data from the patch
        iT : int
            time index
        """
        
        if len(data.shape) == 3:
            self.data[:, :, iT] = data[:, :, 0]
        else:
            self.data[:, :, iT] = data
        
        
    def average(self, inds):
        """Averages data from multiple time stamps
        
        Parameters
        ----------
        inds : list
            List of indices to include in average
        
        Returns
        -------
        array(NX, NY)
            Array containing time-averaged data
        """
        
        return np.mean(self.data[:, :, inds], axis=2)
    
    
    def buildSpace(self):
        """Builds x, y, z coordinates
        
        This function builds an array of the x, y, and z coordinates and
        stores them in the x, y, and z attributes
        """
        
        (NX, NY) = (self.data.shape[0], self.data.shape[1])
        if self.lims[0] == self.lims[1]:
            xGrid = np.zeros((NX, NY)) + self.lims[0]
            y = np.linspace(self.lims[2], self.lims[3], NX)
            z = np.linspace(self.lims[4], self.lims[5], NY)
            yGrid, zGrid = np.meshgrid(y, z)
        elif self.lims[2] == self.lims[3]:
            yGrid = np.zeros((NX, NY)) + self.lims[2]
            x = np.linspace(self.lims[0], self.lims[1], NX)
            z = np.linspace(self.lims[4], self.lims[5], NY)
            xGrid, zGrid = np.meshgrid(x, z)
        elif self.lims[4] == self.lims[5]:
            zGrid = np.zeros((NX, NY)) + self.lims[4]
            x = np.linspace(self.lims[0], self.lims[1], NX)
            y = np.linspace(self.lims[2], self.lims[3], NY)
            xGrid, yGrid = np.meshgrid(x, y)
        (self.x, self.y, self.z) = (xGrid, yGrid, zGrid)
        
        
    def extractPoints(self):
        """Extracts points from patches
        
        Returns
        -------
        array(NX*NY, 3)
            Array containing x, y, z coordinates of points from patch
        array(NX*NY, NT)
            Array containing data for each point for each time stamp
        array(NX*NY, 1)
            Array containing orientation for each point
        """
        
        NX, NY, NT = self.data.shape
        pts = np.zeros((NX*NY, NT))
        coords = np.zeros((NX*NY, 3))
        for iT in range(0, NT):
            pts[:, iT] = self.data[:, :, iT].flatten()
        coords[:, 0] = self.x.flatten()
        coords[:, 1] = self.y.flatten()
        coords[:, 2] = self.z.flatten()
        orients = np.zeros((NX*NY,)) + self.orientation
        return coords, pts, orients
    
    
def getPatches(bndfFile, smvFile, axis, value, meshNum,
               xmin=999, xmax=-999, ymin=999, ymax=-999,
               zmin=999, zmax=-999, dx=999, dz=999):
    """Read patches from boundary file
    
    Parameters
    ----------
    bndfFile : str
        String containing the path to an archive or boundary file
    smvFile : str
        String containing the path to an archive or boundary file
    axis : int
        Integer representing the axis of the patch
    value : float
        Location along axis of the patch
    meshNum : int
        Integer number of the mesh
    xmin : float, optional
        Initialized minimum x value of the patches (default 999)
    xmax : float, optional
        Initialized maximum x value of the patches (default -999)
    ymin : float, optional
        Initialized minimum y value of the patches (default 999)
    ymax : float, optional
        Initialized maximum y value of the patches (default -999)
    zmin : float, optional
        Initialized minimum z value of the patches (default 999)
    zmax : float, optional
        Initialized maximum z value of the patches (default -999)
    dx : float, optional
        Initialized x-delta of the patch (default 999)
    dz : float, optional
        Initialized z-delta of the patch (default 999)
    
    Returns
    -------
    list
        List of floats containing the time stamps
    list
        List of patches
    float
        Minimum x value of the patches
    float
        Maximum x value of the patches
    float
        Minimum y value of the paches
    float
        Maximum y value of the patches
    float
        Minimum z value of the patches
    float
        Maximum z value of the patches
    float
        X-axis delta of the patches
    float
        Z-axis delta of the patches
    """
    
    f = zopen(bndfFile)
    quantity, shortName, units, npatch = parseBndfHeader(f)
    patchInfo, data = parseBndfPatches(f, npatch)
    f.close()
    (patchPts,patchDs,patchIors) = patchInfo
    
    times, patches = importBoundaryFile(
            bndfFile, smvFile, gridNum=meshNum)
    allPatches = []
    for i in range(0, len(patches)):
        patches[i].data[patches[i].data < -1e10] = np.nan
        lims = patches[i].lims
        if lims[0] == lims[1]:
            pass
        if patches[i].data.shape[0]*patches[i].data.shape[1] > -1:
            check = False
            if (abs(axis) == 1):
                if (lims[0] == value) and (lims[1] == value):
                    check = True
            if (abs(axis) == 2):
                if (lims[2] == value) and (lims[3] == value):
                    check = True
            if (abs(axis) == 3):
                if (lims[4] == value) and (lims[5] == value):
                    check = True
            if check and (axis == patchIors[i]):
                allPatches.append(patches[i])
                xmin = min([xmin, lims[0]])
                xmax = max([xmax, lims[1]])
                ymin = min([ymin, lims[2]])
                ymax = max([ymax, lims[3]])
                zmin = min([zmin, lims[4]])
                zmax = max([zmax, lims[5]])
                NX, NZ, NT = patches[i].data.shape
                if abs(axis) == 1:
                    dx = np.round((lims[3]-lims[2])/NX, decimals=4)
                    dz = np.round((lims[5]-lims[4])/NZ, decimals=4)
                elif abs(axis) == 2:
                    dx = np.round((lims[1]-lims[0])/NX, decimals=4)
                    dz = np.round((lims[5]-lims[4])/NZ, decimals=4)
                elif abs(axis) == 3:
                    dx = np.round((lims[1]-lims[0])/NX, decimals=4)
                    dz = np.round((lims[3]-lims[2])/NZ, decimals=4)
                    
                allPatches.append(patches[i])
    return times, allPatches, xmin, xmax, ymin, ymax, zmin, zmax, dx, dz


def buildAbsPatch(patches, xmin, xmax, ymin, ymax, zmin, zmax,
                  dx, dz, axis):
    """Converts local patches to absolute patches
    
    Parameters
    ----------
    patches : list
        List of patches
    xmin : float
        Minimum x value of the patches
    xmax : float
        Maximum x value of the patches
    ymin : float
        Minimum y value of the patches
    ymax : float
        Maximum y value of the patches
    zmin : float
        Minimum z value of the patches
    zmax : float, optional
        Maximum z value of the patches
    dx : float
        X-axis delta of the patches
    dz : float
        Z-axis delta of the patches
    
    Returns
    -------
    array(NXA, NZA)
        Array containing global x-axis coordinates
    array(NXA, NZA)
        Array containing global z-axis coordinates
    array(NXZ, NZA, NT)
        Array containing data in global coordinates for each timestamp
    """
    
    if abs(axis) == 1:
        x_abs = np.linspace(ymin, ymax, int(np.round((ymax-ymin)/dx)+1))
        z_abs = np.linspace(zmin, zmax, int(np.round((zmax-zmin)/dz)+1))
    if abs(axis) == 2:
        x_abs = np.linspace(xmin, xmax, int(np.round((xmax-xmin)/dx)+1))
        z_abs = np.linspace(zmin, zmax, int(np.round((zmax-zmin)/dz)+1))
    if abs(axis) == 3:
        x_abs = np.linspace(xmin, xmax, int(np.round((xmax-xmin)/dx)+1))
        z_abs = np.linspace(ymin, ymax, int(np.round((ymax-ymin)/dz)+1))
    
    x_grid_abs, z_grid_abs = np.meshgrid(x_abs, z_abs)
    NXA, NZA = x_grid_abs.shape
    NX, NZ, NT = patches[0].data.shape
    data_abs = np.zeros((NXA, NZA, NT))
    data_abs[:, :, :] = np.nan
    for patch in patches:
        lims = patch.lims
        if abs(axis) == 1:
            (xMin, xMax) = (lims[2], lims[3])
            (zMin, zMax) = (lims[4], lims[5])
        elif abs(axis) == 2:
            (xMin, xMax) = (lims[0], lims[1])
            (zMin, zMax) = (lims[4], lims[5])
        elif abs(axis) == 3:
            (xMin, xMax) = (lims[0], lims[1])
            (zMin, zMax) = (lims[2], lims[3])
        else:
            estr = "Axis %0.0f not in [-1, -2, -3, 1, 2, 3]."%(axis)
            assert False, estr
        xInd1 = np.argwhere(np.isclose(x_grid_abs, xMin))[0][1]
        xInd2 = np.argwhere(np.isclose(x_grid_abs, xMax))[0][1]
        zInd1 = np.argwhere(np.isclose(z_grid_abs, zMin))[1][0]
        zInd2 = np.argwhere(np.isclose(z_grid_abs, zMax))[1][0]
        for t in range(0, patch.data.shape[2]):
            pdata = patch.data[:, :, t].T
            data_abs[zInd1:zInd2, xInd1:xInd2, t] = pdata
    return x_grid_abs, z_grid_abs, data_abs


def parseBndfHeader(f):
    """Parse header from boundary file
    
    Parameters
    ----------
    f : file
        Binary file open for reading
    
    Returns
    -------
    str
        String containing the boundary quantity
    str
        String containing the quantity short name
    str
        String containing the quantity units
    int
        Number of patches
    """
    
    data = f.read(130)
    header = data[:110]
    patchInfo = data[110:]
    tmp = header.split(b'\x1e')
    
    quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
    shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
    units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
    
    data = np.frombuffer(patchInfo, dtype=np.int32, count=5)
    npatch = data[2]
    return quantity, shortName, units, npatch


def parseBndfPatches(f, npatch):
    pts = 0
    patchPts = []
    patchDs = []
    patchIors = []
    for k in range(0, npatch):
        data = np.frombuffer(f.read(44), dtype=np.int32, count=11)
        #data = np.fromfile(f,dtype=np.int32,count=11)
        #print(k, data)
        (dx,dy,dz) = (data[1]-data[0],data[3]-data[2],data[5]-data[4])
        if abs(data[6]) == 1:
            dy = dy+1
            dz = dz+1
        elif abs(data[6]) == 2:
            dx = dx+1
            dz = dz+1
        elif abs(data[6]) == 3:
            dx = dx+1
            dy = dy+1
        nPts = np.max([dx,1])*np.max([dy,1])*np.max([dz,1])+2
        patchPts.append([pts,pts+nPts])
        patchDs.append([dx,dy,dz,data[0],data[1],data[2],data[3],data[4],data[5]])
        patchIors.append(data[6])
        pts = pts+nPts
    #data = np.fromfile(f,dtype=np.float32,count=-1)
    data = np.frombuffer(f.read(), dtype=np.float32)
    patchInfo = (patchPts, patchDs, patchIors)
    return patchInfo, data

def readBoundaryHeader(file):
    f = zopen(file)
    quantity, shortName, units, npatch = parseBndfHeader(f)
    f.close()
    return quantity, shortName, units, npatch

def readBoundaryQuantities(dataDir, chid):
    files = glob.glob("%s%s*.bf"%(dataDir, chid))
    quantities = defaultdict(bool)
    for file in files:
        quantity, shortName, units, npatch = readBoundaryHeader(file)
        if quantities[quantity]:
            quantities[quantity].append(file)
        else:
            quantities[quantity] = [file]
    return quantities

def importBoundaryFile(fname, smvFile=None, gridNum=0, grid=None):
    try:
        f = zopen(fname)
    except:
        print("Unable to open file %s."%(fname))
        return None, None
    quantity, shortName, units, npatch = parseBndfHeader(f)
    patchInfo, data = parseBndfPatches(f, npatch)
    f.close()
    
    (patchPts,patchDs,patchIors) = patchInfo
    if (grid is None) and (smvFile is not None):
        grid, obst, bndfs, surfs = parseSMVFile(smvFile)
    elif (grid is None) and (smvFile is None):
        print("Either smokeview file or grid must be provided.")
        return None, None
    times, patches = buildPatches(patchPts, patchDs, patchIors, data, grid[gridNum])
    return times, patches

def buildPatches(patchPts, patchDs, patchIors, data, grid):
    pts = patchPts[-1][1]
        
    timeSteps = int((data.shape[0]+1)/(pts+3))
    patches = []
    times = np.zeros((timeSteps-1,))
    for i in range(1,timeSteps):
        ind1 = i*(pts+3)
        ind2 = (i+1)*(pts+3)
        thisData = data[ind1:ind2]
        t = thisData[0]
        times[i-1] = t
        for k in range(0,len(patchPts)):
            (pind1,pind2) = (patchPts[k][0]+3,patchPts[k][1]+1)
            (pdx,pdy,pdz) = (patchDs[k][0],patchDs[k][1],patchDs[k][2])
            patchData = thisData[pind1:pind2].copy()
            if abs(patchIors[k]) == 1:
                patchData = patchData.reshape((pdy,pdz),order='F')
            elif abs(patchIors[k]) == 2:
                patchData = patchData.reshape((pdx,pdz),order='F')
            elif abs(patchIors[k]) == 3:
                patchData = patchData.reshape((pdx,pdy),order='F')
            patchData = patchData[:-1,:-1]
            if i == 1:
                #print(k)
                lims = getLimsFromGrid(patchDs[k][3:],grid)
                patches.append(fdspatch(patchData.shape[0],patchData.shape[1],timeSteps-1,lims,patchIors[k]))
                patches[k].append(patchData,i-1)
            else:
                patches[k].append(patchData,i-1)
    return times, patches

def getLimsFromGrid(data,grid):
    try:
        ind1 = np.argwhere(grid[0][:,0] == data[0])[0][0]
        ind2 = np.argwhere(grid[0][:,0] == data[1])[0][0]
        ind3 = np.argwhere(grid[1][:,0] == data[2])[0][0]
        ind4 = np.argwhere(grid[1][:,0] == data[3])[0][0]
        ind5 = np.argwhere(grid[2][:,0] == data[4])[0][0]
        ind6 = np.argwhere(grid[2][:,0] == data[5])[0][0]
        
        lims = [grid[0][ind1,1],grid[0][ind2,1],
                grid[1][ind3,1],grid[1][ind4,1],
                grid[2][ind5,1],grid[2][ind6,1]]
    except:
        print(grid[0][:,0])
        print(grid[1][:,0])
        print(grid[2][:,0])
        print("Failed to make lims, returning null")
        lims = [0,0,0,0,0,0]
        assert False, "Stopped"
    return lims


def readInputFile(file):
    ''' Load input file '''
    params = defaultdict(bool,yaml.load(open(file,'r')))
    return params

def parseFDSforPts(fdsObsts, smvObsts, names, extend=[0,0,0],fileSMV=None):
    ''' This routine parses an FDS file looking for a list of names and
    stores a list of points defining polygons for each name which is
    found.
    
    The extend argument allows the polygon to be extended by a number of grid
    cells. This is useful since FDS snaps obstructions to the grid which means
    the actual coordinate location of the data from FDS may not align exactly
    with the input file.
    '''

    polygons = []
    for name in names:
        linkedPolygons = []
        for key in list(fdsObsts.keys()):
            if key == name:
                coord = fdsObsts[name]['XB']
                snapInd = np.argmin(np.sum(abs(smvObsts[:,:6]-coord)**2,axis=1)**0.5)
                snapPts = smvObsts[snapInd,13:19].copy()
                pts = [[snapPts[0]-extend[0],snapPts[2]-extend[1],snapPts[4]-extend[2]],
                       [snapPts[0]-extend[0],snapPts[2]-extend[1],snapPts[5]+extend[2]],
                       [snapPts[0]-extend[0],snapPts[3]+extend[1],snapPts[4]-extend[2]],
                       [snapPts[0]-extend[0],snapPts[3]+extend[1],snapPts[5]+extend[2]],
                       [snapPts[1]+extend[0],snapPts[2]-extend[1],snapPts[4]-extend[2]],
                       [snapPts[1]+extend[0],snapPts[2]-extend[1],snapPts[5]+extend[2]],
                       [snapPts[1]+extend[0],snapPts[3]+extend[1],snapPts[4]-extend[2]],
                       [snapPts[1]+extend[0],snapPts[3]+extend[1],snapPts[5]+extend[2]]]
                linkedPolygons.append(pts)
        polygons.append(linkedPolygons)
    return polygons

def loadBNDFdata_lessParams(tStart, tEnd, tInt, tBand, bndfs, smvGrids, smvObsts, orientations, polygons):
    coords2, pts2, times, orients2 = getPointsFromFiles(bndfs, smvGrids, smvObsts, tStart, tEnd, tBand, tInt)
    if orientations[0] == 0:
        orientations = [-3, -2, -1, 1, 2, 3]
    orientInds = np.where([True if x in orientations else False for x in orients2])[0]
    
    coords = coords2[orientInds,:]
    pts = pts2[orientInds,:]
    orients = orients2[orientInds]
    
    # Generate point mask for polygons
    masks = getCoordinateMasks(coords, polygons)
    
    mPts = np.zeros((pts.shape[1],masks.shape[1]))
    for i in range(0,masks.shape[1]):
        if np.where(masks[:,i] == 1)[0].shape[0] > 1:
            mPts[:,i] = np.nanmax(pts[masks[:,i] == 1],axis=0)
        else:
            mPts[:,i] = -1
    
    # Remove last time if it is zero
    if times[-1] == 0:
        mPts = np.delete(mPts,mPts.shape[0]-1,axis=0)
        times = np.delete(times,times.shape[0]-1,axis=0)
        
    return times, mPts, orients

def readBoundaryFile(fname):
    
    try:
        f = zopen(fname)
    except:
        return None, None, None

    
    with open(fname,'rb') as f:
        bytes_read = f.read(110)
        tmp = bytes_read.split(b'\x1e')

        quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
        shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
        units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
        
        data = np.fromfile(f,dtype=np.int32,count=5)
        npatch = data[2]
        pts = 0
        patchPts = []
        patchDs = []
        patchIors = []
        for k in range(0,npatch):
            data = np.fromfile(f,dtype=np.int32,count=11)
            (dx,dy,dz) = (data[1]-data[0],data[3]-data[2],data[5]-data[4])
            
            if abs(data[6]) == 1:
                dy = dy+1
                dz = dz+1
            elif abs(data[6]) == 2:
                dx = dx+1
                dz = dz+1
            elif abs(data[6]) == 3:
                dx = dx+1
                dy = dy+1
            nPts = np.max([dx,1])*np.max([dy,1])*np.max([dz,1])+2
            patchPts.append([pts,pts+nPts])
            patchDs.append([dx,dy,dz,data[0],data[1],data[2],data[3],data[4],data[5]])
            patchIors.append(data[6])
            pts = pts+nPts
        data = np.fromfile(f,dtype=np.float32,count=-1)
    bndfInfo = (quantity,shortName,units,npatch)
    patchInfo = (patchPts,patchDs,patchIors)
    return bndfInfo, patchInfo, data

def extractTime(tStart, tEnd, tBand, tInt, patches, times):
    NT = int((tEnd-tStart)/tInt+1)
    newPatches = []
    for patch in patches:
        newPatches.append(fdspatch(patch.data.shape[0],patch.data.shape[1],NT,patch.lims,patch.orientation))
    tCurrent = tStart
    iT = 0
    newTimes = np.zeros((NT,))
    while tCurrent < tEnd:
        t1 = max([tCurrent-tBand/2,0])
        t2 = min([tCurrent+tBand/2,tEnd])
        newTimes[iT] = (t1+t2)/2
        #print(t1,t2)
        tCurrent = tCurrent+tInt
        tMask = np.argwhere((times >= t1) & (times <= t2))
        for i in range(0,len(patches)):
            newPatches[i].append(patches[i].average(tMask),iT)
        iT = iT+1
    return newTimes, newPatches

def getPatchesFromMesh(grid,obst,bndfFile):
    bndfInfo, patchInfo, data = readBoundaryFile(bndfFile)
    if bndfInfo == None:
        return [None], [None]
    (quantity,shortName,units,npatch) = bndfInfo
    (patchPts,patchDs,patchIors) = patchInfo
    times, patches = buildPatches(patchPts, patchDs, patchIors, data, grid)
    return times, patches

def buildSpace(patches):
    for i in range(0,len(patches)):
        patches[i].buildSpace()

def extractPoints(patches):
    allCoords = []
    allPoints = []
    allOrients = []
    for patch in patches:
        coords, pts, orients = patch.extractPoints()
        if len(allCoords) == 0:
            allCoords = coords
            allPoints = pts
            allOrients = orients
        else:
            allCoords = np.append(allCoords,coords,axis=0)
            allPoints = np.append(allPoints,pts,axis=0)
            allOrients = np.append(allOrients,orients,axis=0)
    return allCoords, allPoints, allOrients

def getPointsFromFiles(bndfs, grids, obsts, tStart, tEnd, tBand, tInt):
    allCoords = []
    allPts = []
    allOrients = []
    newTimes = []
    for file, mesh in bndfs:
        mesh = int(mesh)
        #times, patches = getPatchesFromMesh(grids[mesh-1], obsts[mesh-1], file)
        times, patches = importBoundaryFile(file, gridNum=mesh, grid=grids)
        if len(times) > 1:
            newTimes, newPatches = extractTime(tStart, tEnd, tBand, tInt, patches, times)
            buildSpace(newPatches)
            #savePatches(dataDir,chid,times,patches)
            coords, pts, orients = extractPoints(newPatches)
            
            if len(allCoords) == 0:
                allCoords = coords
                allPts = pts
                allOrients = orients
            else:
                allCoords = np.append(allCoords,coords,axis=0)
                allPts = np.append(allPts,pts,axis=0)
                allOrients = np.append(allOrients,orients,axis=0)
                
    if len(newTimes) == 0:
        print(bndfs)
        assert False, "No valid patches found."
    return allCoords, allPts, newTimes, allOrients

def extractTimeParams(params):
    ''' Extract time parameters from parameter dictionary '''
    initialTime = params['time']['initial']
    finalTime = params['time']['final']
    intervalTime = params['time']['interval']
    bandTime = params['time']['band']
    return initialTime,finalTime,intervalTime,bandTime

def getCoordinateMasks(coords,polygons):
    masks = np.zeros((coords.shape[0],len(polygons)))
    for i in range(0,len(polygons)):
        linkedpolygons = polygons[i]
        for p in linkedpolygons:
            masks[np.where(in_hull(coords,p.points)),i] = 1
    return masks

def queryBndf(resultDir, chid, fdsFilePath, fdsQuantities, fdsUnits, axis, value):
    datas = defaultdict(bool)
    
    bndfFiles = getFileList(resultDir, chid, 'bf')
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsFilePath)
    meshes = list(fdsFile.meshes.keys())
    if 'unknownCounter' in meshes: meshes.remove('unknownCounter')
    numberOfMeshes = len(meshes)
    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        allPatches = []
        (xmin, xmax, ymin, ymax, zmin, zmax, dx, dz) = (999, -999, 999, -999, 999, -999, 999, 999)
        for file in bndfFiles:
            bndfFile = os.path.abspath(file)
            quantity, shortName, units, npatch = readBoundaryHeader(file)    
            meshNumber = 0 if numberOfMeshes == 1 else int(file.split('_')[-2])-1
            if quantity == qty:
                times, patches, xmin1, xmax1, ymin1, ymax1, zmin1, zmax1, dx1, dz1 = getPatches(bndfFile, smvFile, axis, value, meshNumber)
                (xmin, xmax) = (min([xmin, xmin1]), max([xmax, xmax1]))
                (ymin, ymax) = (min([ymin, ymin1]), max([ymax, ymax1]))
                (zmin, zmax) = (min([zmin, zmin1]), max([zmax, zmax1]))
                (dx, dz) = (min([dx, dx1]), min([dz, dz1]))
                for patch in patches: allPatches.append(patch)
        x_grid_abs, z_grid_abs, data_abs = buildAbsPatch(allPatches, xmin, xmax, ymin, ymax, zmin, zmax, dx, dz, axis)
        datas[qty] = defaultdict(bool)
        datas[qty]["MESH-%04.0f"%(meshNumber)] = defaultdict(bool)
        datas[qty]["MESH-%04.0f"%(meshNumber)]['X'] = x_grid_abs
        datas[qty]["MESH-%04.0f"%(meshNumber)]['Z'] = z_grid_abs
        datas[qty]["MESH-%04.0f"%(meshNumber)]['DATA'] = data_abs
    return datas, times

def linkBndfFileToMesh(meshes, bndfs, fdsQuantities):
    if 'unknownCounter' in meshes: meshes.remove('unknownCounter')
    numberOfMeshes = len(meshes)
    
    bndf_dic = defaultdict(bool)
    for qty in fdsQuantities:
        bndf_qty = []
        for bndf in bndfs:
            quantity, shortName, units, npatch = readBoundaryHeader(bndf)
            if quantity == qty:
                meshNumber = 0 if numberOfMeshes == 1 else int(bndf.split('_')[-2])-1
                bndf_qty.append([bndf, meshNumber])
        bndf_dic[qty] = bndf_qty
    return bndf_dic

def extractMaxBndfValues(fdsFilePath, smvFilePath, resultDir, chid, fdsQuantities,
                         tStart=0, tEnd=120, tInt=1, tBand=3, orientations=[0]):
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsFilePath)
    meshes = list(fdsFile.meshes.keys())
    names = fdsFile.getPolygonNamesFromFdsFile()
    smvGrids, smvObsts, smvBndfs, smvSurfs = parseSMVFile(smvFilePath)
    #linesSMV = zreadlines(smvFilePath)
    points = parseFDSforPts(fdsFile.obsts, smvObsts, names, extend=[0,0,0])
    polygons, numberOfGroups = pts2polygons(points)
    
    bndfs = getFileList(resultDir, chid, 'bf')
    bndf_dic = linkBndfFileToMesh(meshes, bndfs, fdsQuantities)
    
    smvGrids, smvObsts, smvBndfs, smvSurfs = parseSMVFile(smvFilePath)
    
    datas = defaultdict(bool)
    for qty in fdsQuantities:
        datas[qty] = defaultdict(bool)
        times, mPts, orients = loadBNDFdata_lessParams(tStart, tEnd, tInt, tBand, bndf_dic[qty], smvGrids, smvObsts, orientations, polygons)
        datas[qty]['TIMES'] = times
        datas[qty]['NAMES'] = names
        datas[qty]['DATA'] = mPts
    return datas