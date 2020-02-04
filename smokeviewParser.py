# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 08:17:40 2019

@author: JHodges
"""

import numpy as np

def parseGRID(lines,i):
    gridPts = [int(x) for x in lines[i+1].replace('\n','').split()]
    gridTRNX = np.array([[float(y) for y in x.replace('\n','').split()] for x in lines[i+8:i+9+gridPts[0]]])
    gridTRNY = np.array([[float(y) for y in x.replace('\n','').split()] for x in lines[i+12+gridPts[0]:i+13+gridPts[0]+gridPts[1]]])
    gridTRNZ = np.array([[float(y) for y in x.replace('\n','').split()] for x in lines[i+16+gridPts[0]+gridPts[1]:i+17+gridPts[0]+gridPts[1]+gridPts[2]]])
    return gridTRNX, gridTRNY, gridTRNZ

def calculateDeltas(gridTRNX, gridTRNY, gridTRNZ):
    dx = (gridTRNX.max()-gridTRNX.min())/(gridTRNX.shape[0]-1)
    dy = (gridTRNY.max()-gridTRNY.min())/(gridTRNY.shape[0]-1)
    dz = (gridTRNZ.max()-gridTRNZ.min())/(gridTRNZ.shape[0]-1)
    return dx, dy, dz

def parseOBST(lines, i, gridTRNX, gridTRNY, gridTRNZ):
    dx, dy, dz = calculateDeltas(gridTRNX, gridTRNY, gridTRNZ)
    numOBST = int(lines[i+1].replace(' ',''))
    tmp1 = lines[i+2:i+2+numOBST]
    tmp2 = lines[i+2+numOBST:i+2+numOBST+numOBST]
    tmp1 = [x.replace('\n','') for x in tmp1]
    tmp2 = [x.replace('\n','') for x in tmp2]
    tmp1 = [[float(y) for y in x.split()] for x in tmp1]
    tmp2 = [[float(y) for y in x.split()] for x in tmp2]
    for i in range(0, len(tmp2)):
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
        smvObj[j,13:19] = newPts
    return smvObj

def parseBNDF(lines, i):
    (_,mesh,vNum) = lines[i-1].split()
    bndfName = lines[i].split(' ')[1].replace('\n','')
    vID = ' '.join(lines[i+1].split(' ')[1:]).replace('\n','')
    (mesh, vNum) = (float(mesh), float(vNum))
    return mesh, bndfName, vID, vNum

def parseSMVFile(smvFile):
    with open(smvFile,'r') as f:
        linesSMV = f.readlines()
    grids = []
    obsts = []
    bndfs = []
    surfs = []
    for i in range(0,len(linesSMV)):
        line2 = linesSMV[i]
        if ("GRID" in line2):
            gridTRNX, gridTRNY, gridTRNZ = parseGRID(linesSMV, i)
            grids.append([gridTRNX.copy(), gridTRNY.copy(), gridTRNZ.copy()])
            dx, dy, dz = calculateDeltas(gridTRNX, gridTRNY, gridTRNZ)
        if ("OBST" in line2) and ("HIDE_OBST" not in line2):
            smvObj = parseOBST(linesSMV, i, gridTRNX, gridTRNY, gridTRNZ)
            if len(obsts) == 0:
                obsts = smvObj
            else:
                obsts = np.append(obsts,smvObj,axis=0)
            #obsts.append(smvObj)
        if (".bf" in line2):
            mesh, bndfName, vID, vNum = parseBNDF(linesSMV, i)
            bndfs.append([mesh, bndfName, vID, vNum])
        if 'SURFACE\n' in linesSMV[i]:
            sname = ' '.join(linesSMV[i+1].split())
            (Tign,eps) = (linesSMV[i+2].split()[0],linesSMV[i+2].split()[1])
            (stype,t_width,t_height) = (linesSMV[i+3].split()[0],linesSMV[i+3].split()[1],linesSMV[i+3].split()[2])
            (c1,c2,c3,c4) = (linesSMV[i+3].split()[3],linesSMV[i+3].split()[4],linesSMV[i+3].split()[5],linesSMV[i+3].split()[6])
            surfs.append([sname,Tign,eps,stype,t_width,t_height,c1,c2,c3,c4])
    return grids, obsts, bndfs, surfs
