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
from collections import defaultdict
import glob
import cv2
import pandas as pd
from .utilities import zopen, in_hull, getFileList, pts2polygons, zreadlines
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
    times : float array(NT)
        array containing the timestamps of the patch
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
        if abs(OR) == 1:
            self.dx = (DS[3]-DS[2])/(NX-1)
            self.dz = (DS[5]-DS[4])/(NY-1)
        elif abs(OR) == 2:
            self.dx = (DS[1]-DS[0])/(NX-1)
            self.dz = (DS[5]-DS[4])/(NY-1)
        elif abs(OR) == 3:
            self.dx = (DS[1]-DS[0])/(NX-1)
            self.dz = (DS[3]-DS[2])/(NY-1)
        
        
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
        
        (NX, NY) = (self.data.shape[0]+1, self.data.shape[1]+1)
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
        #coords = np.zeros((NX*NY, 3))
        coords = np.zeros(((NX+1)*(NY+1), 3))
        for iT in range(0, NT):
            pts[:, iT] = self.data[:, :, iT].flatten()
        coords[:, 0] = self.x.flatten() #((self.x[:-1,:-1]+self.x[1:,1:])/2).flatten() #self.x.flatten()
        coords[:, 1] = self.y.flatten() #((self.y[:-1,:-1]+self.y[1:,1:])/2).flatten() #self.y.flatten()
        coords[:, 2] = self.z.flatten() #((self.z[:-1,:-1]+self.z[1:,1:])/2).flatten() #self.z.flatten()
        orients = np.zeros((NX*NY,)) + self.orientation
        return coords, pts, orients
    
    
def getPatchOptions(bndfFile, smvFile, meshNum, smvData=None):
    """Read patches from boundary file and return options for each axis
    """
    f = zopen(bndfFile)
    quantity, shortName, units, npatch = parseBndfHeader(f)
    patchInfo, data = parseBndfPatches(f, npatch)
    f.close()
    (patchPts, patchDs, patchIors, patchNBs, patchNMs) = patchInfo
    
    times, patches, units = importBoundaryFile(
            bndfFile, smvFile, gridNum=meshNum, smvData=smvData)
    xoptions = []
    yoptions = []
    zoptions = []
    for i in range(0, len(patches)):
        patches[i].data[patches[i].data < -1e10] = np.nan
        lims = patches[i].lims
        if lims[0] == lims[1]: xoptions.append([lims, patchIors[i]])
        if lims[2] == lims[3]: yoptions.append([lims, patchIors[i]])
        if lims[4] == lims[5]: zoptions.append([lims, patchIors[i]])
    return xoptions, yoptions, zoptions

        
def getPatches(bndfFile, smvFile, axis, value, meshNum,
               xmin=999, xmax=-999, ymin=999, ymax=-999,
               zmin=999, zmax=-999, dx=999, dz=999,
               decimals=4, smvData=None,
               showOptions=False):
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
    decimals : int, optional
        Number of decimals to use in rounding (default 4)
    
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
    (patchPts, patchDs, patchIors, patchNBs, patchNMs) = patchInfo
    
    times, patches, units = importBoundaryFile(
            bndfFile, smvFile, gridNum=meshNum, smvData=smvData)
    allPatches = []
    if patches is None:
        return None, None, None, None, None, None, None, None, None, None
    for i in range(0, len(patches)):
        patches[i].data[patches[i].data < -1e10] = np.nan
        #print(patches[i].lims, patchIors[i])
        lims = patches[i].lims
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
                '''
                if abs(axis) == 1:
                    dx = np.round((lims[3]-lims[2])/(NX-1), decimals=decimals+2)
                    dz = np.round((lims[5]-lims[4])/(NZ-1), decimals=decimals+2)
                elif abs(axis) == 2:
                    dx = np.round((lims[1]-lims[0])/(NX-1), decimals=decimals+2)
                    dz = np.round((lims[5]-lims[4])/(NZ-1), decimals=decimals+2)
                elif abs(axis) == 3:
                    dx = np.round((lims[1]-lims[0])/(NX-1), decimals=decimals+2)
                    dz = np.round((lims[3]-lims[2])/(NZ-1), decimals=decimals+2)
                '''
                dx = patches[i].dx
                dz = patches[i].dz
                allPatches.append(patches[i])
    return times, allPatches, xmin, xmax, ymin, ymax, zmin, zmax, dx, dz


def buildAbsPatch(patches, xmin, xmax, ymin, ymax, zmin, zmax,
                  dx, dz, axis, decimals=4):
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
    decimals : int, optional
        Number of decimals to round (default 4)
    
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
        if patches[0].cell_centered == True:
            #x_abs = np.linspace(ymin, ymax-dx, int(np.round((ymax-ymin)/dx)))
            #z_abs = np.linspace(zmin, zmax-dz, int(np.round((zmax-zmin)/dz)))
            x_abs = np.linspace(ymin, ymax, int(np.round((ymax-ymin)/dx)+1))
            z_abs = np.linspace(zmin, zmax, int(np.round((zmax-zmin)/dz)+1))
        else:
            x_abs = np.linspace(ymin, ymax, int(np.round((ymax-ymin)/dx)+1))
            z_abs = np.linspace(zmin, zmax, int(np.round((zmax-zmin)/dz)+1))
    if abs(axis) == 2:
        if patches[0].cell_centered == True:
            x_abs = np.linspace(xmin, xmax, int(np.round((xmax-xmin)/dx)+1))
            z_abs = np.linspace(zmin, zmax, int(np.round((zmax-zmin)/dz)+1))
        else:
            x_abs = np.linspace(xmin, xmax, int(np.round((xmax-xmin)/dx)+1))
            z_abs = np.linspace(zmin, zmax, int(np.round((zmax-zmin)/dz)+1))
    if abs(axis) == 3:
        if patches[0].cell_centered == True:
            #x_abs = np.linspace(xmin, xmax-dx, int(np.round((xmax-xmin)/dx)))
            #z_abs = np.linspace(ymin, ymax-dz, int(np.round((ymax-ymin)/dz)))
            x_abs = np.linspace(xmin, xmax, int(np.round((xmax-xmin)/dx)+1))
            z_abs = np.linspace(ymin, ymax, int(np.round((ymax-ymin)/dz)+1))
        else:
            x_abs = np.linspace(xmin, xmax, int(np.round((xmax-xmin)/dx)+1))
            z_abs = np.linspace(ymin, ymax, int(np.round((ymax-ymin)/dz)+1))
    x_abs = np.round(x_abs, decimals=decimals)
    z_abs = np.round(z_abs, decimals=decimals)
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
        xInd1 = np.argwhere(np.isclose(x_grid_abs, xMin, atol=10**(-1*decimals)))[0][1]
        if xInd1 == 0:
            xInd2 = np.argwhere(np.isclose(x_grid_abs, xMax, atol=10**(-1*decimals)))[0][1]
        else:
            #print(patch.cell_centered)
            #print(x_grid_abs)
            #print(xMax)
            xInd2 = np.argwhere(np.isclose(x_grid_abs, xMax, atol=10**(-1*decimals)))[0][1]
        zInd1 = np.argwhere(np.isclose(z_grid_abs, zMin, atol=10**(-1*decimals)))[1][0]
        if zInd1 == 0:
            zInd2 = np.argwhere(np.isclose(z_grid_abs, zMax, atol=10**(-1*decimals)))[1][0]
        else:
            zInd2 = np.argwhere(np.isclose(z_grid_abs, zMax, atol=10**(-1*decimals)))[1][0]
        NT = min([patch.data.shape[2], data_abs[zInd1:zInd2, xInd1:xInd2, :].shape[2]])
        for t in range(0, NT):
            pdata = patch.data[:, :, t].T
            if pdata.shape[0] != data_abs[zInd1:zInd2, xInd1:xInd2, t].shape[0]:
                tmp = data_abs[zInd1:zInd2, xInd1:xInd2, t].shape
                tmp2 = (tmp[1], tmp[0])
                pdata2 = cv2.resize(pdata, dsize=tmp2, interpolation=cv2.INTER_LINEAR)
                data_abs[zInd1:zInd2, xInd1:xInd2, t] = pdata2
            else:
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
    """Parse patches from boundary file
    
    Parameters
    ----------
    f : file
        Binary file open for reading
    npatch : int
        Number of patches in boundary file
    
    Returns
    -------
    tuple
        Tuple containing points, limits, orientations, boundary number, and
        mesh number
    float array(NX*NZ*NT)
        Array containing data from patches
    """
    
    pts = 0
    patchPts = []
    patchDs = []
    patchIors = []
    patchNBs = []
    patchNMs = []
    for k in range(0, npatch):
        data = np.frombuffer(f.read(44), dtype=np.int32, count=11)
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
        patchPts.append([pts, pts+nPts])
        patchDs.append([dx,dy,dz,
                        data[0],data[1],data[2],
                        data[3],data[4],data[5]])
        patchIors.append(data[6])
        patchNBs.append(data[7])
        patchNMs.append(data[8])
        pts = pts+nPts
    dat = f.read()
    if len(dat) % 4 != 0:
        data = np.frombuffer(dat[:int(len(dat)/4)*4], dtype=np.float32)
    else:
        data = np.frombuffer(dat, dtype=np.float32)
    patchInfo = (patchPts, patchDs, patchIors, patchNBs, patchNMs)
    return patchInfo, data


def readBoundaryHeader(file):
    """Read boundary header from boundary file
    
    Parameters
    ----------
    file : str
        Path to boundary file to read
    
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
    
    f = zopen(file)
    quantity, shortName, units, npatch = parseBndfHeader(f)
    f.close()
    return quantity, shortName, units, npatch


def readBoundaryQuantities(dataDir, chid):
    """Returns quantities available in boundary files
    
    This function searches through a data directory or archive and
    determines which QUANTITIES are available in all files.
    
    Parameters
    ----------
    dataDir : str
        Path to data file
    chid : str
        CHID of simulation to check
    
    Returns
    -------
    defaultdict
        Dictionary containing the quantities available in the files
    """
    
    files = getFileList(dataDir, chid, 'bf')
    quantities = defaultdict(bool)
    for file in files:
        quantity, shortName, units, npatch = readBoundaryHeader(file)
        if quantities[quantity]:
            quantities[quantity].append(file)
        else:
            quantities[quantity] = [file]
    return quantities


def importBoundaryFile(fname, smvFile=None, gridNum=0, smvData=None):
    """Import patches from a single boundary file
    
    Parameters
    ----------
    fname : str
        Path to boundary file to import
    smvFile : str
        Path to smokeview file
    gridNum : int, optional
        Grid corresponding to the patch. Must be specified if multiple
        meshes are used.
    smvData : dict, optional
        Previously imported smvData
    
    Returns
    -------
    list
        List of float time stamps
    list
        List of patches
    """
    
    try:
        f = zopen(fname)
    except:
        print("Unable to open file %s."%(fname))
        return None, None, None
    quantity, shortName, units, npatch = parseBndfHeader(f)
    patchInfo, data = parseBndfPatches(f, npatch)
    f.close()
    (patchPts, patchDs, patchIors, patchNBs, patchNMs) = patchInfo
    if (smvData is None) and (smvFile is not None):
        smvData = parseSMVFile(smvFile)
    elif (smvData is None) and (smvFile is None):
        print("Either smokeview file or smokeview data must be provided.")
        return None, None, None
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    times, patches = buildPatches(
            patchPts, patchDs, patchIors, data, grid[gridNum])
    fname2 = fname.split('.zip')[-1].split(os.sep)[-1]
    #print(bndfs)
    bndf_extract = [[x[1], x[4]] for x in bndfs]
    #print(bndf_extract)
    bnd_types = [x[1] for x in bndf_extract if x[0] == fname2]
    #print('bnd_types')
    #print(fname2)
    #print(bnd_types)
    if len(bnd_types) == 0: return None, None, None
    bnd_type = bnd_types[0]
    
    for i in range(0, len(patches)):
        if bnd_type == 'BNDC':
            patches[i].data = patches[i].data[:-1,:-1, :]
            patches[i].cell_centered = True
            XB = patches[i].lims
            IOR = patches[i].orientation
            (dx, dz) = (patches[i].dx, patches[i].dz)
            #print("PRE", IOR, XB, dx, dy)
            if abs(IOR) == 1:
                (XB[2], XB[3]) = (XB[2]+dx/2, XB[3]+dx/2)
                (XB[4], XB[5]) = (XB[4]+dz/2, XB[5]+dz/2)
            elif abs(IOR) == 2:
                (XB[0], XB[1]) = (XB[0]+dx/2, XB[1]+dx/2)
                (XB[4], XB[5]) = (XB[4]+dz/2, XB[5]+dz/2)
            elif abs(IOR) == 3:
                (XB[0], XB[1]) = (XB[0]+dx/2, XB[1]+dx/2)
                (XB[2], XB[3]) = (XB[2]+dz/2, XB[3]+dz/2)
            #print("POST", IOR, XB)
            patches[i].lims = XB
        else:
            patches[i].cell_centered = False
    return times, patches, units


def buildPatches(patchPts, patchDs, patchIors, data, grid):
    """Build spatial patches from patch points
    
    Parameters
    ----------
    patchPts : list
        Point indicies
    patchDs : list
        Patch deltas and limits
    patchIors : list
        Patch orientations
    data : float array(NX*NZ*NT)
        Array containing data from the boundary file
    grid : list
        List of arrays containing the grid of the mesh
        
    Returns
    -------
    float array(NT)
        Array containing time stamps
    list
        List containing each patch
    """
    
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
            #patchData = patchData[:-1,:-1]
            if i == 1:
                lims = getLimsFromGrid(patchDs[k][3:],grid)
                patch = fdspatch(patchData.shape[0],
                                 patchData.shape[1],
                                 timeSteps-1, lims, patchIors[k])
                patches.append(patch)
                patches[k].append(patchData,i-1)
            else:
                patches[k].append(patchData,i-1)
    for k in range(0, len(patches)):
        patches[k].times = times
    return times, patches


def getLimsFromGrid(data, grid):
    """Extracts limits from mesh grid and patch data
    
    Parameters
    ----------
    data : list
        List containing limits from patch
    grid : list
        List of arrays containing the grid of the mesh
        
    Returns
    -------
    float array(6)
        Array containing limit extents
    """
    
    try:
        ind1 = np.argwhere(grid[0][:,0] == data[0])[0][0]
        #print("ind1 " + str(ind1))
        ind2 = np.argwhere(grid[0][:,0] == data[1])[0][0]
        #print("ind2 " + str(ind2))
        ind3 = np.argwhere(grid[1][:,0] == data[2])[0][0]
        #print("ind3 " + str(ind3))
        ind4 = np.argwhere(grid[1][:,0] == data[3])[0][0]
        #print("ind4 " + str(ind4))
        ind5 = np.argwhere(grid[2][:,0] == data[4])[0][0]
        #print("ind5 " + str(ind5))
        ind6 = np.argwhere(grid[2][:,0] == data[5])[0][0]
        #print("ind6 " + str(ind6))
        
        lims = [grid[0][ind1,1], grid[0][ind2,1],
                grid[1][ind3,1], grid[1][ind4,1],
                grid[2][ind5,1], grid[2][ind6,1]]
    except:
        print(grid[0][:,0])
        print(grid[1][:,0])
        print(grid[2][:,0])
        print(data)
        print("Failed to make lims, returning null")
        lims = [0,0,0,0,0,0]
        assert False, "Stopped"
    return lims


def parseFDSforPts(fdsObsts, smvObsts, names, extend=[0,0,0]):
    """ Builds a list of polygons from an fds input file
    
    This routine parses an FDS file looking for a list of names and
    stores a list of points defining polygons for each name which is
    found.
    
    The extend argument allows the polygon to be extended by a number of
    grid cells. This is useful since FDS snaps obstructions to the grid
    which means the actual coordinate location of the data from FDS may
    not align exactly with the input file.
    
    Parameters
    ----------
    fdsObsts : dict
        Dictionary containing obstruction definitions from fds input
    smvObsts : dict
        Dictionary containing obstruction definitions from smv input
    names : list
        List of obstruction names of interest
    extend : float array(3), optional
        Distance to extend polygon along each axis
    
    Returns
    -------
    list
        List of polygons corresponding to each name
    """
    
    polygons = []
    for name in names:
        linkedPolygons = []
        obstNames = list(fdsObsts.keys())
        if 'unknownCounter' in obstNames: obstNames.remove('unknownCounter')
        for key in obstNames:
            if fdsObsts[key]['ID'] == name:
                coord = fdsObsts[key]['XB']
                obstdif = abs(smvObsts[:,:6]-coord)
                obstsum = np.sum((obstdif)**2, axis=1)**0.5
                snapInd = np.argmin(obstsum)
                snapPts = smvObsts[snapInd,13:19].copy()
                p1 = snapPts[0]-extend[0]
                p2 = snapPts[2]-extend[1]
                p3 = snapPts[4]-extend[2]
                p4 = snapPts[5]+extend[2]
                p5 = snapPts[3]+extend[1]
                p6 = snapPts[1]+extend[0]
                pts = [[p1, p2, p3], [p1, p2, p4], [p1, p5, p3],
                       [p1, p5, p4], [p6, p2, p3], [p6, p2, p4],
                       [p6, p5, p3], [p6, p5, p4]]
                linkedPolygons.append(pts)
        polygons.append(linkedPolygons)
    return polygons


def loadBNDFdata(tStart, tEnd, tInt, tBand, bndfs, smvData,
                 orientations, polygons):
    """Loads maximum boundary data versus time from boundary files
    
    This routine will read all boundary data from boundary files for
    a simulation.
    
    Parameters
    ----------
    tStart : float
        Timestamp to start extracting data
    tEnd : float
        Timestamp to finish extracting data
    tInt : float
        Timestamp interval to include in extracted data
    tBand : float
        Time averaging window to include in extracted data
    bndfs : list
        List containing string names of boundary files to be loaded
    smvData : dict
        Dictionary containing data imported from smokeview
    orientations : list
        List containing plane orientations to include in loading
    polygons : list
        List of polygons corresponding to each name
    
    Returns
    -------
    float array(NT)
        Array containing time stamps
    float array(NT)
        Array containing max temperature from patch
    float array(NT)
        Array containing orientations
    """
    
    coords2, pts2, times, orients2, units = getPointsFromFiles(
            bndfs, smvData, tStart, tEnd, tBand, tInt)
    #import pandas as pd
    #pd.DataFrame(coords2).to_csv('E:\\projects\\kansas_city_fire_modeling\\fromcluster\\coords.csv')
    #pd.DataFrame(pts2).to_csv('E:\\projects\\kansas_city_fire_modeling\\fromcluster\\pts.csv')
    #pd.DataFrame(times).to_csv('E:\\projects\\kansas_city_fire_modeling\\fromcluster\\times.csv')
    #pd.DataFrame(orients2).to_csv('E:\\projects\\kansas_city_fire_modeling\\fromcluster\\orients.csv')
    #print(coords2.shape, pts2.shape, times.shape, orients2.shape)
    if orientations[0] == 0:
        orientations = [-3, -2, -1, 1, 2, 3]
    mask = [True if x in orientations else False for x in orients2]
    orientInds = np.where(mask)[0]
    
    coords = coords2[orientInds,:]
    pts = pts2[orientInds,:]
    orients = orients2[orientInds]
    
    # Generate point mask for polygons
    masks = getCoordinateMasks(coords, polygons)
    
    #pd.DataFrame(coords).to_csv('test_points.csv')
    #pd.DataFrame(orients).to_csv('test_orients.csv')
    #pd.DataFrame(masks).to_csv('test_masks.csv')
    #pd.DataFrame(pts).to_csv('test_pts.csv')
    
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
        
    return times, mPts, orients, units


def readBoundaryFile(fname):
    """Loads data from a boundary file
    
    This routine will read all boundary data from a boundary file
    
    Parameters
    ----------
    fname : str
        Path to boundary file to import
    
    Returns
    -------
    tuple
        Tuple containing boundary quantity, short name, nunits, and
        number of patches
    tuple
        Tuple containing patchPts, patchDs, and patchIors
    float array(NX, NZ, NT)
        Array containing arranged data from boundary file
    """
    
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
        
        data = np.fromfile(f, dtype=np.int32, count=5)
        npatch = data[2]
        pts = 0
        patchPts = []
        patchDs = []
        patchIors = []
        patchNBs = []
        patchNMs = []
        for k in range(0, npatch):
            data = np.fromfile(f, dtype=np.int32, count=11)
            dx = data[1] - data[0]
            dy = data[3] - data[2]
            dz = data[5] - data[4]
            
            if abs(data[6]) == 1:
                dy = dy+1
                dz = dz+1
            elif abs(data[6]) == 2:
                dx = dx+1
                dz = dz+1
            elif abs(data[6]) == 3:
                dx = dx+1
                dy = dy+1
            nPts = np.max([dx, 1])*np.max([dy, 1])*np.max([dz, 1]) + 2
            patchPts.append([pts, pts+nPts])
            patchDs.append([dx, dy, dz, data[0], data[1], data[2],
                            data[3], data[4], data[5]])
            patchIors.append(data[6])
            patchNBs.append(data[7])
            patchNMs.append(data[8])
            pts = pts+nPts
        data = np.fromfile(f, dtype=np.float32, count=-1)
    bndfInfo = (quantity, shortName, units, npatch)
    patchInfo = (patchPts, patchDs, patchIors, patchNBs, patchNMs)
    return bndfInfo, patchInfo, data


def extractTime(tStart, tEnd, tBand, tInt, patches, times):
    """Extracts a timestamp interval from a patch time series
    
    Parameters
    ----------
    tStart : float
        Timestamp to start extracting data
    tEnd : float
        Timestamp to finish extracting data
    tInt : float
        Timestamp interval to include in extracted data
    tBand : float
        Time averaging window to include in extracted data
    patches : list
        List of patches
    times : list
        List of timestamps
    
    Returns
    -------
    list
        List containing new time stamps
    list
        List containing new patches
    """
    
    NT = int((tEnd-tStart)/tInt + 1)
    newPatches = []
    for patch in patches:
        newPatch = fdspatch(patch.data.shape[0],
                            patch.data.shape[1],
                            NT, patch.lims, patch.orientation)
        newPatches.append(newPatch)
    tCurrent = tStart
    iT = 0
    newTimes = np.zeros((NT,))
    while tCurrent < tEnd:
        t1 = max([tCurrent - tBand/2, 0])
        t2 = min([tCurrent + tBand/2, tEnd])
        newTimes[iT] = (t1 + t2)/2
        tCurrent = tCurrent + tInt
        tMask = np.argwhere((times >= t1) & (times <= t2))
        for i in range(0, len(patches)):
            newPatches[i].append(patches[i].average(tMask), iT)
        iT = iT + 1
    return newTimes, newPatches


def getPatchesFromMesh(grid, obst, bndfFile):
    """Extracts patches from a mesh
    
    Parameters
    ----------
    grid : list
        List of arrays containing grid coordinates
    obst : list
        List of obstructions from smokeview file
    bndfFile : str
        String path to boundary file
    
    Returns
    -------
    list
        List containing time stamps
    list
        List containing patches
    """
    
    bndfInfo, patchInfo, data = readBoundaryFile(bndfFile)
    if bndfInfo == None:
        return [None], [None]
    (quantity, shortName, units, npatch) = bndfInfo
    (patchPts, patchDs, patchIors, patchNBs, patchNMs) = patchInfo
    times, patches = buildPatches(
            patchPts, patchDs, patchIors, data, grid)
    return times, patches


def buildSpace(patches):
    """Builds spatial grid for each patch in a list
    
    Parameters
    ----------
    patches : list
        List of patches
    """
    
    for i in range(0,len(patches)):
        patches[i].buildSpace()
    
    
def extractPoints(patches):
    """Extracts points from a list of patches
    
    Parameters
    ----------
    patches : list
        List of patches
    
    Returns
    -------
    list
        List containing all coordinates
    list
        List containing all points
    list
        List containing all orientations
    """
    
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
            allCoords = np.append(allCoords, coords, axis=0)
            allPoints = np.append(allPoints, pts, axis=0)
            allOrients = np.append(allOrients, orients, axis=0)
    return allCoords, allPoints, allOrients


def getPointsFromFiles(bndfs, smvData, tStart, tEnd, tBand, tInt):
    """Extracts a timestamp interval from a patch time series
    
    Parameters
    ----------
    bndfs : list
        List of filenames for boundary files
    smvData : dict
        Dictionary containing data from smokeview file
    tStart : float
        Timestamp to start extracting data
    tEnd : float
        Timestamp to finish extracting data
    tInt : float
        Timestamp interval to include in extracted data
    tBand : float
        Time averaging window to include in extracted data
    
    Returns
    -------
    list
        List containing all coordinates
    list
        List containing all points
    list
        List containing all times
    list
        List containing all orientations
    """
    allCoords = []
    allPts = []
    allOrients = []
    newTimes = []
    for file, mesh in bndfs:
        mesh = int(mesh)
        times, patches, units = importBoundaryFile(
                file, gridNum=mesh, smvData=smvData)
        if len(times) > 1:
            if (tInt != None) and (tBand != None):
                newTimes, newPatches = extractTime(
                        tStart, tEnd, tBand, tInt, patches, times)
                buildSpace(newPatches)
            else:
                newTimes, newPatches = times, patches
                buildSpace(newPatches)
            coords, pts, orients = extractPoints(newPatches)
            
            if len(allCoords) == 0:
                allCoords = coords
                allPts = pts
                allOrients = orients
            else:
                allCoords = np.append(allCoords, coords, axis=0)
                allPts = np.append(allPts, pts, axis=0)
                allOrients = np.append(allOrients, orients, axis=0)
                
    if len(newTimes) == 0:
        print(bndfs)
        assert False, "No valid patches found."
    return allCoords, allPts, newTimes, allOrients, units


def getCoordinateMasks(coords, polygons):
    """Returns a mask of coordinates contained in a polygon
    
    Parameters
    ----------
    coords : list
        List of coordinates
    polygons : list
        List of polygons
    
    Returns
    -------
    list
        List of coordinate masks
    """
    
    masks = np.zeros((coords.shape[0], len(polygons)))
    for i in range(0,len(polygons)):
        linkedpolygons = polygons[i]
        for p in linkedpolygons:
            masks[np.where(in_hull(coords, p.points)), i] = 1
    return masks


def timeAverageBndfs(resultDir, chid, fdsFilePath, fdsQuantity, dt):
    """Timev average boundary files
    
    Parameters
    ----------
    resultDir : str
        String containing path to results directory
    chid : str
        String containing CHID from simulation
    fdsFilePath : str
        String containing path to fds file or achive
    fdsQuantity : str
        String containing the quantity to time average
    """
    bndfFiles = getFileList(resultDir, chid, 'bf')
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsFilePath)
    meshes = list(fdsFile.meshes.keys())
    if 'unknownCounter' in meshes: meshes.remove('unknownCounter')
    numberOfMeshes = len(meshes)
    
    for file in bndfFiles:
        quantity, shortName, units, npatch = readBoundaryHeader(file)
        if numberOfMeshes == 1:
            meshNumber = 0
        else:
            meshNumber = int(file.split('_')[-2]) - 1
            
        if quantity == fdsQuantity:
            times, patches, units = importBoundaryFile(
                file, smvFile, gridNum=meshNumber)
            
            
            
    

def queryBndf(resultDir, chid, fdsFilePath, fdsQuantities, 
              axis, value, decimals=4):
    """Query boundary files
    
    Parameters
    ----------
    resultDir : str
        String containing path to results directory
    chid : str
        String containing CHID from simulation
    fdsFilePath : str
        String containing path to fds file or achive
    fdsQuantities : list
        List containing quantities to query from boundary file
    axis : int
        Integer specifying axis to query
    value : float
        Axis value to query
    decimals : int
        Number of decimals to include in rounding
    
    Returns
    -------
    dict
        Dictionary containty data for each quantity
    list
        List of timestamps
    """
    
    datas = defaultdict(bool)
    
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsFilePath)
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    smvData = parseSMVFile(smvFile)
    bndfDir = resultDir
    if fdsFile.dump['ID'] is not False:
        if fdsFile.dump['ID']['RESULTS_DIR'] is not False:
            bndfDir = resultDir + os.sep + fdsFile.dump['ID']['RESULTS_DIR'] + os.sep
            for i in range(0, len(smvData['bndfs'])):
                f = smvData['bndfs'][i][1]
                f = os.path.abspath(f)
                f = f.split(os.sep)[-1]
                smvData['bndfs'][i][1] = f
    #print(bndfDir)
    bndfFiles = getFileList(bndfDir, chid, 'bf')
    #print('bndfFiles',bndfFiles)
    #print(resultDir, chid)
    
    meshes = list(fdsFile.meshes.keys())
    if 'unknownCounter' in meshes: meshes.remove('unknownCounter')
    numberOfMeshes = 0
    for mesh in meshes:
        if fdsFile.meshes[mesh]['MULT_ID'] is not False:
            mult_id = fdsFile.meshes[mesh]['MULT_ID']
            mult = fdsFile.mult[mult_id]
            (I_UPPER, J_UPPER, K_UPPER) = (0, 0, 0)
            if mult['I_UPPER'] is not False: I_UPPER = mult['I_UPPER']
            if mult['J_UPPER'] is not False: J_UPPER = mult['J_UPPER']
            if mult['K_UPPER'] is not False: K_UPPER = mult['K_UPPER']
            numberOfMeshes += (I_UPPER+1)*(J_UPPER+1)*(K_UPPER+1)
        else:
            numberOfMeshes += 1
    
    #numberOfMeshes = len(meshes)
    #print(numberOfMeshes)
    ts_out = None
    for qty in fdsQuantities:
        allPatches = []
        (xmin, xmax) = (999, -999)
        (ymin, ymax) = (999, -999)
        (zmin, zmax) = (999, -999)
        (dx, dz) = (999, 999)
        qtyFound = False
        for file in bndfFiles:
            quantity, shortName, units, npatch = readBoundaryHeader(
                    file)
            if numberOfMeshes == 1:
                meshNumber = 0
            else:
                meshNumber = int(file.split('_')[-2]) - 1
            #print(quantity, qty, shortName, units, npatch)
            #print(file, meshNumber)
            if quantity == qty:
                #print(axis, value)
                ts, ps, x1, x2, y1, y2, z1, z2, dx1, dz1 = getPatches(
                        file, smvFile, axis, value, meshNumber, decimals=decimals, smvData=smvData)
                #print(file, dx1, dz1)
                if ts is None: continue
                ts_out = ts  
                #print(file, meshNumber, x1, x2, y1, y2, z1, z2)
                #print(file.split('_')[-2], numberOfMeshes)
                #print(quantity, qty, shortName, units, npatch)
                (xmin, xmax) = (min([xmin, x1]), max([xmax, x2]))
                (ymin, ymax) = (min([ymin, y1]), max([ymax, y2]))
                (zmin, zmax) = (min([zmin, z1]), max([zmax, z2]))
                (dx, dz) = (min([dx, dx1]), min([dz, dz1]))
                for patch in ps:
                    allPatches.append(patch)
                    qtyFound = True
                #print(len(allPatches))
                outUnits = units
        if len(allPatches) == 0:
            for file in bndfFiles:
                quantity, shortName, units, npatch = readBoundaryHeader(
                        file)
                if numberOfMeshes == 1:
                    meshNumber = 0
                else:
                    meshNumber = int(file.split('_')[-2]) - 1
                if quantity == qty:
                    xoptions, yoptions, zoptions = getPatchOptions(
                            file, smvFile, meshNumber, smvData=smvData)
                    print("Queried axis and value not found.")
                    print("Options for IOR %0.0f in %s"%(axis, file))
                    for option in xoptions:
                        if option[1] == axis:
                            print(option[0])
                    for option in yoptions:
                        if option[1] == axis:
                            print(option[0])
                    for option in zoptions:
                        if option[1] == axis:
                            print(option[0])
                    
        #print("ABS INFO: ")
        #print(xmin, xmax, ymin, ymax, zmin, zmax, dx, dz)
        if qtyFound:
            x_grid_abs, z_grid_abs, data_abs = buildAbsPatch(
                    allPatches, xmin, xmax, ymin, ymax, zmin, zmax,
                    dx, dz, axis, decimals=decimals)
            datas[qty] = defaultdict(bool)
            datas[qty]['X'] = x_grid_abs
            datas[qty]['Z'] = z_grid_abs
            datas[qty]['DATA'] = data_abs
            datas[qty]['UNITS'] = outUnits
        else:
            print("Quantity %s not found in boundary data"%(qty))
            quantities = readBoundaryQuantities(resultDir, chid)
            print("Available quantites are:")
            for q in quantities:
                print("\t%s"%(q))
        #datas[qty]["MESH-%04.0f"%(meshNumber)] = defaultdict(bool)
        #datas[qty]["MESH-%04.0f"%(meshNumber)]['X'] = x_grid_abs
        #datas[qty]["MESH-%04.0f"%(meshNumber)]['Z'] = z_grid_abs
        #datas[qty]["MESH-%04.0f"%(meshNumber)]['DATA'] = data_abs
    tmax = len(ts)
    for qty in list(datas.keys()):
        tmax = min([datas[qty]['DATA'].shape[2], tmax])
    for qty in list(datas.keys()):
        datas[qty]['DATA'] = datas[qty]['DATA'][:, :, :tmax]
    ts = ts[:tmax]
    return datas, ts


def linkBndfFileToMesh(meshes, bndfs, fdsQuantities):
    """Links boundary file back to mesh
    
    Parameters
    ----------
    meshes : list
        List containing mesh definitions
    bndfs : list
        List containing path to each boundary file to query
    fdsQuantities : list
        List containing quantities to query from boundary file
    
    Returns
    -------
    dict
        Dictionary containty mesh for each boundary file
    """
    
    if 'unknownCounter' in meshes: meshes.remove('unknownCounter')
    numberOfMeshes = len(meshes)
    
    bndf_dic = defaultdict(bool)
    for qty in fdsQuantities:
        bndf_qty = []
        for bndf in bndfs:
            quantity, shortName, units, npatch = readBoundaryHeader(
                    bndf)
            if quantity == qty:
                if numberOfMeshes == 1:
                    meshNumber = 0
                else:
                    meshNumber = int(bndf.split('_')[-2]) - 1
                bndf_qty.append([bndf, meshNumber])
        bndf_dic[qty] = bndf_qty
    return bndf_dic


def extractMaxBndfValues(fdsF, smvF, resultDir, chid, quantities,
                         tStart=0, tEnd=120, tInt=1, tBand=3,
                         orientations=[0], extend=[0, 0, 0],
                         names=None):
    """Extract maximum value from boundary file
    
    Parameters
    ----------
    fdsF : str
        String containing path to fds file or archive
    smvF : str
        String containing path to smv file or archive
    resultDir : str
        String containing path to results directory
    chid : str
        String containing CHID from simulation
    quantities : list
        List containing quantities to query from boundary file
    tStart : float
        Timestamp to start extracting data
    tEnd : float
        Timestamp to finish extracting data
    tInt : float
        Timestamp interval to include in extracted data
    tBand : float
        Time averaging window to include in extracted data
    orientations : list
        List of integers specifying which orientations to consider
    names : list
        List of names specifiynig which obstructions to extract
    
    Returns
    -------
    dict
        Dictionary containty data for each quantity
    """
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsF)
    meshes = list(fdsFile.meshes.keys())
    if names == None:
        #print(fdsFile.obsts)
        names = fdsFile.getPolygonNamesFromFdsFile()
    
    #print(names)
    
    #names = ['CE212319', 'CE210006', 'CE212298']

    smvData = parseSMVFile(smvF)
    (smvGrids, smvObsts) = (smvData['grids'], smvData['obsts'])
    (smvBndfs, smvSurfs) = (smvData['bndfs'], smvData['surfs'])
    (smvFiles, bndes) = (smvData['files'], smvData['bndes'])
    
    fdsObsts = fdsFile.obsts
    points = parseFDSforPts(fdsObsts, smvObsts, names, extend=extend)
    polygons, numberOfGroups = pts2polygons(points)
    
    bndfs = getFileList(resultDir, chid, 'bf')
    bndf_dic = linkBndfFileToMesh(meshes, bndfs, quantities)
    #pd.DataFrame(names).to_csv('test_names.csv')
    datas = defaultdict(bool)
    for qty in quantities:
        datas[qty] = defaultdict(bool)
        times, mPts, orients, units = loadBNDFdata(
                tStart, tEnd, tInt, tBand, bndf_dic[qty], 
                smvData, orientations, polygons)
        datas[qty]['TIMES'] = times
        datas[qty]['NAMES'] = names
        datas[qty]['DATA'] = mPts
        datas[qty]['UNITS'] = units
    return datas


def bndfsTimeAverage(resultDir, chid, fdsQuantity, dt, outDir=None, outQty=None):
    bndfFiles = getFileList(resultDir, chid, 'bf')
    filesWithQueriedQuantity = []
    for boundaryFile in bndfFiles:
        f = zopen(boundaryFile)
        quantity, shortName, units, npatch = parseBndfHeader(f)
        f.close()
        if quantity == fdsQuantity:
            filesWithQueriedQuantity.append(boundaryFile)
    
    outFiles = []
    for boundaryFile in filesWithQueriedQuantity:
        print("Starting to avg file %s"%(boundaryFile))
        if outDir is None:
            outFile = None
        else:
            tmp = boundaryFile.replace('.bf','_avg.bf')
            tmp = tmp.split(os.sep)[-1]
            outFile = os.sep.join([outDir, tmp])
            tmp = outFile.split(os.sep)
            outFile = os.sep.join([x for x in tmp if '.zip' not in x])
        outFile, outQty = bndfTimeAverage(boundaryFile, dt, outFile, outQty)
        print("wrote %s to %s"%(outQty, outFile))
        outFiles.append(outFile)
    
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    linesSMV = zreadlines(smvFile)
    smvData = parseSMVFile(smvFile)
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    bndf_ref_files = [x[1] for x in bndfs]
    meshes = [int(x[0]) for x in bndfs]
    vnums = [int(x[3]) for x in bndfs]
    
    bndfLines = []
    for refFile, outFile in zip(filesWithQueriedQuantity, outFiles):
        rfile = refFile.split(os.sep)[-1]
        mesh, vnum = [(y, z) for x, y, z in zip(bndf_ref_files, meshes, vnums) if rfile in x][0]
        outFile = outFile.split(os.sep)[-1]
        bndfLine = buildBndfSmvLine(mesh, vnum, outQty, outFile, shortName, units)
        bndfLines.extend(bndfLine)
    
    if outDir is None:
        newSmvFile = outFile.split(os.sep)[:-1]
        newSmvFile.append(smvFile.split(os.sep)[-1].replace('.smv','_avg.smv'))
        newSmvFile = os.sep.join(newSmvFile)
    else:
        newSmvFile = outDir.split(os.sep)
        newSmvFile.append(smvFile.split(os.sep)[-1].replace('.smv','_avg.smv'))
        newSmvFile = os.sep.join(newSmvFile)
    
    linesSMV.extend(bndfLines)
    smvText = '\n'.join(linesSMV)
    smvText = smvText + '\n'
    smvText = smvText.replace('\n\n','\n')
    with open(newSmvFile, 'w') as f:
        f.write(smvText)
    print("wrote to %s"%(newSmvFile))
    
    return outFiles, outQty, filesWithQueriedQuantity, newSmvFile
        

def buildBndfSmvLine(mesh, vnum, qty, file, shortName, units):
    bndfLine = "BNDF" + str(mesh).rjust(6) + str(vnum).rjust(6) + "\n"
    bndfLine = bndfLine + " " + file + "\n"
    bndfLine = bndfLine + " " + qty + "\n"
    bndfLine = bndfLine + " " + shortName + "\n"
    bndfLine = bndfLine + " " + units
    bndfLine = bndfLine.split('\n')
    return bndfLine
    

def bndfTimeAverage(boundaryFile, dt, outFile=None, outQty=None):
    f = zopen(boundaryFile)
    quantity, shortName, units, npatch = parseBndfHeader(f)
    patchInfo, data = parseBndfPatches(f, npatch)
    f.close()
    pts = patchInfo[0][-1][1]
    
    times = data[::(pts+3)]
    
    # Reshape data into array for time averaging
    data2 = np.zeros((pts-3, len(times)))
    for i, time in enumerate(times):
        ind1 = i*(pts+3)
        ind2 = (i+1)*(pts+3)
        data2[:, i] = data[ind1+3:ind2-3]
    
    # Time average the array
    data3 = data2.copy()
    for i, time in enumerate(times):
        t1 = max([(time-dt/2), times[0]])
        t2 = min([(time+dt/2), times[-1]])
        tinds = np.logical_and(times>=t1, times<=t2)
        
        data3[:, i] = np.nanmean(data2[:, tinds], axis=1)
    
    # Reshape data back to 1-D array for binary write
    data4 = np.array(data, dtype=np.float32)
    for i, time in enumerate(times):
        ind1 = i*(pts+3)
        ind2 = (i+1)*(pts+3)
        
        data4[ind1+3:ind2-3] = data3[:, i]
    
    data5 = data4.tobytes()
    if outFile is None:
        outFile = boundaryFile.replace('.bf','_avg.bf')
        if '.zip' in outFile:
            tmp = outFile.split(os.sep)
            outFile = os.sep.join([x for x in tmp if '.zip' not in x])
    if outQty is None:
        outQty = shortName + ' (%0.0ds Avg)'%(dt)
    
    with open(outFile, 'wb') as f:
        header = writeBndfHeader(outQty, shortName, units, npatch)
        darray = writeBndfPatchInfo(patchInfo)
        
        f.write(header)
        f.write(darray)
        f.write(data5)
    return outFile, outQty


def writeBndfHeader(quantity, shortName, units, npatch):
    
    header = b'\x1e\x00\x00\x00'
    header = header + (quantity.ljust(30, ' ')).encode('utf-8')
    header = header + b'\x1e\x00\x00\x00' + b'\x1e\x00\x00\x00'
    header = header + (shortName.ljust(30, ' ')).encode('utf-8')
    header = header + b'\x1e\x00\x00\x00' + b'\x1e\x00\x00\x00'
    header = header + (units.ljust(30, ' ')).encode('utf-8')
    #header = header + b'\x1e\x00\x00\x00\x04\x00\x00\x00'
    #header = header + b'\xbf\x0f\x00\x00\x04\x00\x00\x00$\x00\x00\x00'
    
    npatch_array = np.array([30, 4, npatch, 4, 36], dtype=np.int32)
    
    npatch_binary = npatch_array.tobytes()
    
    header = header + npatch_binary
    
    return header


def writeBndfPatchInfo(patchInfo):
    patchDs = patchInfo[1]
    patchIors = patchInfo[2]
    patchNBs = patchInfo[3]
    patchNMs = patchInfo[4]
    
    patchHeader = b''
    for i in range(0, len(patchDs)):
        patchD = patchDs[i]
        darray = np.array(patchD[3:], dtype=np.int32)
        if i < len(patchDs) - 1:
            darray = np.append(darray, [patchIors[i], patchNBs[i], patchNMs[i], 36, 36])
        else:
            darray = np.append(darray, [patchIors[i], patchNBs[i], patchNMs[i], 36, 4])
        patchHeader = patchHeader + darray.tobytes()
    return patchHeader
