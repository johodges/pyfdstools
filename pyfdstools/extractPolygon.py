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
# This script extracts boundary data from defined polygons.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import numpy as np
import os
from collections import defaultdict
#from . import utilities as ut

from .fdsFileOperations import fdsFileOperations
from .utilities import in_hull, zreadlines, getFileList, pts2polygons
from .extractBoundaryData import linkBndfFileToMesh, loadBNDFdata_lessParams
from .smokeviewParser import parseSMVFile

def extractMaxBndfValues(fdsFilePath, smvFilePath, resultDir, chid, fdsQuantities,
                         tStart=0, tEnd=120, tInt=1, tBand=3, orientations=[0]):
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsFilePath)
    meshes = list(fdsFile.meshes.keys())
    names = getPolygonNamesFromFdsFile(fdsFile)
    linesSMV = zreadlines(smvFilePath)
    points = parseFDSforPts(fdsFile, linesSMV, names, extend=[0,0,0])
    polygons, numberOfGroups = pts2polygons(points)
    
    bndfs = getFileList(resultDir, chid, 'bf')
    bndf_dic = linkBndfFileToMesh(meshes, bndfs, fdsQuantities)
    
    smvGrids, smvObsts, smvBndfs, smvSurfs, bndes = parseSMVFile(smvFilePath)
    
    datas = defaultdict(bool)
    for qty in fdsQuantities:
        datas[qty] = defaultdict(bool)
        times, mPts, orients = loadBNDFdata_lessParams(tStart, tEnd, tInt, tBand, bndf_dic[qty], smvGrids, smvObsts, orientations, polygons)
        datas[qty]['TIMES'] = times
        datas[qty]['NAMES'] = names
        datas[qty]['DATA'] = mPts
    return datas

def getPolygonNamesFromFdsFile(file):
    names = []
    obstList = list(file.obsts.keys())
    if 'unknownCounter' in obstList: obstList.remove('unknownCounter')
    for key in obstList:
        if file.obsts[key]['BNDF_OBST']:
            names.append(file.obsts[key]["ID"])
    names = list(set(names))
    return names

def parseFDSforPts(fileFDS, linesSMV, names, extend=[0,0,0]):
    ''' This routine parses an FDS file looking for a list of names and
    stores a list of points defining polygons for each name which is
    found.
    
    The extend argument allows the polygon to be extended by a number of grid
    cells. This is useful since FDS snaps obstructions to the grid which means
    the actual coordinate location of the data from FDS may not align exactly
    with the input file.
    '''
    smvObjs = []
    for i in range(0,len(linesSMV)):
        line2 = linesSMV[i]
        if "GRID" in line2:
            gridPts = [int(x) for x in linesSMV[i+1].replace('\n','').split()]
            gridTRNX = np.array([[float(y) for y in x.replace('\n','').split()] for x in linesSMV[i+8:i+9+gridPts[0]]])
            gridTRNY = np.array([[float(y) for y in x.replace('\n','').split()] for x in linesSMV[i+12+gridPts[0]:i+13+gridPts[0]+gridPts[1]]])
            gridTRNZ = np.array([[float(y) for y in x.replace('\n','').split()] for x in linesSMV[i+16+gridPts[0]+gridPts[1]:i+17+gridPts[0]+gridPts[1]+gridPts[2]]])
            dx = (gridTRNX.max()-gridTRNX.min())/(gridTRNX.shape[0]-1)
            dy = (gridTRNY.max()-gridTRNY.min())/(gridTRNY.shape[0]-1)
            dz = (gridTRNZ.max()-gridTRNZ.min())/(gridTRNZ.shape[0]-1)
        if "OBST" in line2 and "HIDE_OBST" not in line2:
            try:
                numOBST = int(linesSMV[i+1].replace(' ',''))
            except:
                print(linesSMV[i-2:i+2])
                assert False, "Stopped"
            tmp1 = linesSMV[i+2:i+2+numOBST]
            tmp2 = linesSMV[i+2+numOBST:i+2+numOBST+numOBST]
            tmp1 = [x.replace('\n','') for x in tmp1]
            tmp2 = [x.replace('\n','') for x in tmp2]
            tmp1 = [[float(y) for y in x.split()] for x in tmp1]
            tmp2 = [[float(y) for y in x.split()] for x in tmp2]
            for i in range(0, len(tmp1)):
                if len(tmp2[i]) > 8:
                    tmp2[i] = tmp2[i][:8]
                    
            smvObj = np.array([x1+x2 for x1, x2 in zip(tmp1,tmp2)])
            for j in range(0,smvObj.shape[0]):
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
                #print("Pre-snap:",smvObj[j,:6])
                #print("Post-snap:",newPts)
                smvObj[j,13:19] = newPts
            if len(smvObjs) == 0:
                smvObjs = smvObj
            else:
                smvObjs = np.append(smvObjs,smvObj,axis=0)
    
    obstKeys = list(fileFDS.obsts.keys())
    if 'unknownCounter' in obstKeys: obstKeys.remove('unknownCounter')
    polygons = []
    for name in names:
        linkedPolygons = []
        for key in obstKeys:
            if name in fileFDS.obsts[key]['ID']:
                coord = fileFDS.obsts[key]['XB']
                snapInd = np.argmin(np.sum(abs(smvObjs[:,:6]-coord)**2,axis=1)**0.5)
                snapPts = smvObjs[snapInd,13:19].copy()
                
                pts = [[snapPts[0]-extend[0],snapPts[2]-extend[1],snapPts[4]-extend[2]],
                       [snapPts[0]-extend[0],snapPts[2]-extend[1],snapPts[5]+extend[2]],
                       [snapPts[0]-extend[0],snapPts[3]+extend[1],snapPts[4]-extend[2]],
                       [snapPts[0]-extend[0],snapPts[3]+extend[1],snapPts[5]+extend[2]],
                       [snapPts[1]+extend[0],snapPts[2]-extend[1],snapPts[4]-extend[2]],
                       [snapPts[1]+extend[0],snapPts[2]-extend[1],snapPts[5]+extend[2]],
                       [snapPts[1]+extend[0],snapPts[3]+extend[1],snapPts[4]-extend[2]],
                       [snapPts[1]+extend[0],snapPts[3]+extend[1],snapPts[5]+extend[2]]]
                #print("Before snapping:",coord)
                #print("After snapping:",newPts)
                linkedPolygons.append(pts)
        polygons.append(linkedPolygons)
    return polygons

def parseFDSforVID(file,vName):
    '''
    '''
    with open(file,'r') as f:
        lines = f.readlines()
    vIDCounter = 0
    for line2 in lines:
        line = line2.replace('/\n','')
        if '&BNDF' in line:
            vIDCounter = vIDCounter + 1
            if vName in line:
                vID = vIDCounter
    return vID

def getCoordinateMasks(coords,polygons):
    masks = np.zeros((coords.shape[0],len(polygons)))
    for i in range(0,len(polygons)):
        linkedpolygons = polygons[i]
        for p in linkedpolygons:
            masks[np.where(in_hull(coords,p.points)),i] = 1
    return masks

def getCoordinateMasks2(coords,polygons):
    masks = np.zeros((coords.shape[0],len(polygons)))
    for i in range(0,len(polygons)):
        linkedpolygons = polygons[i]
        for p in linkedpolygons:
            for j in range(0,coords.shape[0]):
                if ut.pnt_in_cvex_hull(p, coords[j,:]):
                    masks[j,i] = 1
    return masks


