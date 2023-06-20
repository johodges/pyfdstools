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
# This script parses the smokeview file generated by FDS.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import numpy as np
from .utilities import zreadlines
from collections import defaultdict


def parseGRID(lines, i):
    """This function parses grid data from smokeview lines.
    
    Parameters
    ----------
    lines : list
        List of strings corresponding to lines from a smokeview file
    i : int
        Index of grid line
    
    Returns
    -------
    array(I)
        Array containing x-coordinates of grid
    array(J)
        Array containing y-coordinates of grid
    array(K)
        Array containing z-coordinates of grid
    """
    
    gridPts = [int(x) for x in lines[i+1].replace('\n','').split()]
    (gridTRNX, gridTRNY, gridTRNZ) = ([], [], [])
    (xind0, xind1) = (i+8, i+9+gridPts[0])
    (yind0, yind1) = (i+12+gridPts[0], i+13+gridPts[0]+gridPts[1])
    zind0 = i+16+gridPts[0]+gridPts[1]
    zind1 = i+17+gridPts[0]+gridPts[1]+gridPts[2]
    xlines = lines[xind0:xind1]
    ylines = lines[yind0:yind1]
    zlines = lines[zind0:zind1]
    for x in xlines:
        gridTRNX.append([float(y) for y in x.replace('\n','').split()])
    for x in ylines:
        gridTRNY.append([float(y) for y in x.replace('\n','').split()])
    for x in zlines:
        gridTRNZ.append([float(y) for y in x.replace('\n','').split()])
    gridTRNX = np.array(gridTRNX)
    gridTRNY = np.array(gridTRNY)
    gridTRNZ = np.array(gridTRNZ)
    return gridTRNX, gridTRNY, gridTRNZ


def calculateDeltas(gridTRNX, gridTRNY, gridTRNZ):
    """This function calculates the deltas for each axis
    
    Parameters
    ----------
    gridTRNX : array(I)
        Array containing x-coordinates of grid
    gridTRNY : array(J)
        Array containing y-coordinates of grid
    gridTRNZ : array(K)
        Array containing z-coordinates of grid
    
    Returns
    -------
    float
        Delta for x-axis
    float
        Delta for y-axis
    float
        Delta for z-axis
    """
    
    dx = (gridTRNX.max()-gridTRNX.min())/(gridTRNX.shape[0]-1)
    dy = (gridTRNY.max()-gridTRNY.min())/(gridTRNY.shape[0]-1)
    dz = (gridTRNZ.max()-gridTRNZ.min())/(gridTRNZ.shape[0]-1)
    return dx, dy, dz


def parseOBST(lines, i, gridTRNX, gridTRNY, gridTRNZ):
    """This function parses and obstruction from a someview line
    
    Parameters
    ----------
    lines : list
        List of strings corresponding to lines from a smokeview file
    i : int
        Index of obstruction line
    gridTRNX : array(I)
        Array containing x-coordinates of grid
    gridTRNY : array(J)
        Array containing y-coordinates of grid
    gridTRNZ : array(K)
        Array containing z-coordinates of grid
    
    Returns
    -------
    array(M, 20)
        Array containing obstruction information from smokeview file
    """
    
    dx, dy, dz = calculateDeltas(gridTRNX, gridTRNY, gridTRNZ)
    numOBST = int(lines[i+1].replace(' ',''))
    tmp1 = lines[i+2:i+2+numOBST]
    tmp2 = lines[i+2+numOBST:i+2+numOBST+numOBST]
    tmp1 = [x.replace('\n','').split('!')[0] for x in tmp1]
    tmp2 = [x.replace('\n','').split('!')[0] for x in tmp2]
    tmp1 = [[float(y) for y in x.split()] for x in tmp1]
    tmp2 = [[float(y) for y in x.split()] for x in tmp2]
    for j in range(0, len(tmp2)):
        if len(tmp2[j]) > 8:
            tmp2[j] = tmp2[j][:8]
    smvObj = np.array([x1+x2 for x1, x2 in zip(tmp1,tmp2)])
    for j in range(0, smvObj.shape[0]):
        pts = smvObj[j,13:19]
        x1 = gridTRNX[np.where(gridTRNX[:,0] == pts[0])[0][0],1]
        x2 = gridTRNX[np.where(gridTRNX[:,0] == pts[1])[0][0],1]
        y1 = gridTRNY[np.where(gridTRNY[:,0] == pts[2])[0][0],1]
        y2 = gridTRNY[np.where(gridTRNY[:,0] == pts[3])[0][0],1]
        z1 = gridTRNZ[np.where(gridTRNZ[:,0] == pts[4])[0][0],1]
        z2 = gridTRNZ[np.where(gridTRNZ[:,0] == pts[5])[0][0],1]
        newPts = np.array([x1,x2,y1,y2,z1,z2])
        if newPts[0] == newPts[1]: newPts[1] = newPts[1] + dx
        if newPts[2] == newPts[3]: newPts[3] = newPts[3] + dy
        if newPts[4] == newPts[5]: newPts[5] = newPts[5] + dz
        smvObj[j,13:19] = newPts
    return smvObj


def parseBNDF(lines, i):
    """This function parses the smokeview file to find bndf files
    
    Parameters
    ----------
    lines : list
        List of strings corresponding to lines from a smokeview file
    i : int
        Index of bndf line
    
    Returns
    -------
    float
        Mesh number for bndf file
    string
        Boundary file name
    string
        String variable number
    float
        Variable number
    """
    
    (_,mesh,vNum) = lines[i-1].split()
    bndfName = lines[i].split(' ')[1].replace('\n','')
    vID = ' '.join(lines[i+1].split(' ')[1:]).replace('\n','')
    (mesh, vNum) = (float(mesh), float(vNum))
    return mesh, bndfName, vID, vNum


def parseBNDE(lines, i):
    """This function parses the smokeview file to find bnde files
    
    Parameters
    ----------
    lines : list
        List of strings corresponding to lines from a smokeview file
    i : int
        Index of bndf line
    
    Returns
    -------
    float
        Mesh number for bndf file
    string
        Boundary file name
    string
        String variable number
    float
        Variable number
    """
    
    (_,mesh,vNum) = lines[i].split()
    bndeName = lines[i+1].split(' ')[1].replace('\n','')
    gbfName = lines[i+2].split(' ')[1].replace('\n','')
    vID = ' '.join(lines[i+3].split(' ')[1:]).replace('\n','')
    (mesh, vNum) = (float(mesh), float(vNum))
    return mesh, bndeName, gbfName, vID, vNum


def parseSMVFile(smvFile):
    """This function parses a smokeview file
    
    Parameters
    ----------
    smvFile : str
        String containing path to archive or smokeview file
    
    Returns
    -------
    list
        List of grids
    list
        List of obstructions
    list
        List of boundary files
    list
        List of surfaces
    """
    
    linesSMV = zreadlines(smvFile)
    grids = []
    obsts = []
    bndfs = []
    surfs = []
    bndes = []
    files = defaultdict(bool)
    files['SLICES'] = defaultdict(bool)
    for i in range(0, len(linesSMV)):
        line2 = linesSMV[i]
        
        if ("GRID" in line2):
            gridTRNX, gridTRNY, gridTRNZ = parseGRID(linesSMV, i)
            grids.append([gridTRNX.copy(),
                          gridTRNY.copy(),
                          gridTRNZ.copy()])
            dx, dy, dz = calculateDeltas(gridTRNX, gridTRNY, gridTRNZ)
        if ("OBST" in line2[:6]) and ("HIDE_OBST" not in line2):
            smvObj = parseOBST(
                    linesSMV, i, gridTRNX, gridTRNY, gridTRNZ)
            if len(smvObj) > 0:
                if len(obsts) == 0:
                    obsts = smvObj
                else:
                    obsts = np.append(obsts, smvObj, axis=0)
        if (".bf" in line2):
            mesh, bndfName, vID, vNum = parseBNDF(linesSMV, i)
            bndfs.append([mesh, bndfName, vID, vNum])
        if ("BNDE" in line2[:6]):
            mesh, bndeName, gbfName, vID, vNum = parseBNDE(linesSMV, i)
            bndes.append([mesh, bndeName, gbfName, vID, vNum])
        if 'SURFACE\n' in linesSMV[i]:
            sname = ' '.join(linesSMV[i+1].split())
            Tign = linesSMV[i+2].split()[0]
            eps = linesSMV[i+2].split()[1]
            stype = linesSMV[i+3].split()[0]
            t_width = linesSMV[i+3].split()[1]
            t_height = linesSMV[i+3].split()[2]
            c1 = linesSMV[i+3].split()[3]
            c2 = linesSMV[i+3].split()[4]
            c3 = linesSMV[i+3].split()[5]
            c4 = linesSMV[i+3].split()[6]
            surfs.append([sname, Tign, eps, stype, t_width, t_height, 
                          c1, c2, c3, c4])
        if 'SLCF' in linesSMV[i]:
            file = '%s.sf'%(linesSMV[i+1][1:].split('.sf')[0])
            files['SLICES'][file] = defaultdict(bool)
            files['SLICES'][file]['CELL_CENTERED'] = False
            files['SLICES'][file]['QUANTITY'] = linesSMV[i+2].strip()
            files['SLICES'][file]['SHORTNAME'] = linesSMV[i+3].strip()
            files['SLICES'][file]['UNITS'] = linesSMV[i+4].strip()
            files['SLICES'][file]['LINETEXT'] = linesSMV[i]
        if 'SLCC' in linesSMV[i]:
            file = '%s.sf'%(linesSMV[i+1][1:].split('.sf')[0])
            files['SLICES'][file] = defaultdict(bool)
            files['SLICES'][file]['CELL_CENTERED'] = True
            files['SLICES'][file]['QUANTITY'] = linesSMV[i+2].strip()
            files['SLICES'][file]['SHORTNAME'] = linesSMV[i+3].strip()
            files['SLICES'][file]['UNITS'] = linesSMV[i+4].strip()
            files['SLICES'][file]['LINETEXT'] = linesSMV[i]
    return grids, obsts, bndfs, surfs, files, bndes

