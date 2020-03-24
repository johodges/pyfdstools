#----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#----------------------------------------------------------------------
#======================================================================
# 
# DESCRIPTION:
# This software is part of a python library to assist in developing and
# analyzing simulation results from Fire Dynamics Simulator (FDS).
# FDS is an open source software package developed by NIST. The source
# code is available at: https://github.com/firemodels/fds
# 
# This script extracts boundary data from defined polygons.
#
#=========================================================================
# # IMPORTS
#=========================================================================
import numpy as np
import yaml
import matplotlib.colors as pltc
import os
from collections import defaultdict
import sys
from . import utilities as ut

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl

import scipy.optimize as scop

def setupQT():
    app = QtGui.QApplication([])
    w = gl.GLViewWidget()
    w.show()
    w.setWindowTitle('pyqtgraph example: GLMeshItem')
    w.setCameraPosition(distance=40)
    
    g = gl.GLGridItem()
    g.scale(2,2,1)
    w.addItem(g)
    return app, w

def getRunInfo(params):
    ''' Return run information from parameter dictionary '''
    chid = params['chid']
    orientations = params['orientations']
    outputs = defaultdict(bool,params['outputs'])
    
    if params['space']:
        dLimited = params['space']['limitDomain']
        if dLimited:
            dLimits = params['space']['domainLimits']
        else:
            dLimits=[0,1,0,1,0,1]
    else:
        dLimited = False
        dLimits = [0,1,0,1,0,1]
    return chid, orientations, outputs, dLimited, dLimits

def getPolysFromInput(params):
    ''' Get points from entries in parameter dictionary '''
    allPoints = []
    for name in params['polygons']:
        linkedPolygons = params[name]['points']
        linkedPts = []
        for polygon in linkedPolygons:
            pts = polygon
            linkedPts.append(pts)
        allPoints.append(linkedPts)
    return allPoints

def getNamesFromInput(params,file):
    ''' Get names from entries in parameter dictionary '''
    if params['polygons']:
        names = []
        for p in params['polygons']:
            if params[p]:
                pDict = defaultdict(bool,params[p])
                if pDict['name']:
                    names.append(pDict['name'])
                else:
                    names.append(p)
            else:
                names.append(p)
                params[p] = defaultdict(bool,{'name':p,'threshold':params['defaultThreshold']})
                print(params[p])
    else:
        with open(file,'r') as f:
            lines = f.readlines()
        comments = []
        for x in lines:
            sLine = x.replace('\n','').split('/')
            if len(sLine) <= 1:
                comments.append('')
            elif len(sLine) == 2:
                comments.append(sLine[1])
            else:
                comments.append('/'.join(sLine[1:]))
        #comments = [x.replace('\n','').split('/')[1] for x in lines]
        names = []
        for line, cmt in zip(lines,comments):
            if 'BNDF_OBST=.TRUE.' in line:
                names.append(cmt)
        names = list(set(names))
        params['polygons'] = names
        for name in names:
            params[name] = defaultdict(bool,{'threshold':params['defaultThreshold']})
    return names, params

def parseFDSforPts_old(file,names,extend=[0,0,0]):
    ''' This routine parses an FDS file looking for a list of names and
    stores a list of points defining polygons for each name which is
    found.
    
    The extend argument allows the polygon to be extended by a number of grid
    cells. This is useful since FDS snaps obstructions to the grid which means
    the actual coordinate location of the data from FDS may not align exactly
    with the input file.
    '''
    with open(file,'r') as f:
        lines = f.readlines()
    (xmn,xmx) = (100000,-100000)
    (ymn,ymx) = (100000,-100000)
    (zmn,zmx) = (100000,-100000)
    
    for line2 in lines:
        line = line2.replace('/\n','')
        if '&MESH' in line:
            inds = [float(x) for x in line.split('IJK=')[1].split(',')[:3]]
            limsStr = line.split('XB=')[1].split(',')[:6]
            try:
                float(limsStr[5])
            except:
                limsStr[5] = limsStr[5].split('/')[0]
            lims = [float(x) for x in limsStr]
            dx = (lims[1]-lims[0])/inds[0]
            dy = (lims[3]-lims[2])/inds[1]
            dz = (lims[5]-lims[4])/inds[2]
            (xmn,xmx) = (min([xmn,lims[0]]),max([xmx,lims[1]]))
            (ymn,ymx) = (min([ymn,lims[2]]),max([ymx,lims[3]]))
            (zmn,zmx) = (min([zmn,lims[4]]),max([zmx,lims[5]]))
    xs = np.linspace(xmn,xmx,(xmx-xmn)/dx+1)
    ys = np.linspace(ymn,ymx,(ymx-ymn)/dy+1)
    zs = np.linspace(zmn,zmx,(zmx-zmn)/dz+1)
    print(xs)
    print(ys)
    print(zs)
    #print(xs)
    comments = []
    for x in lines:
        sLine = x.replace('\n','').split('/')
        if len(sLine) <= 1:
            comments.append('')
        elif len(sLine) == 2:
            comments.append(sLine[1])
        else:
            comments.append('/'.join(sLine[1:]))
    polygons = []
    for name in names:
        linkedPolygons = []
        for line, cmt in zip(lines,comments):
            #print(name,cmt)
            if name in cmt and "&OBST" in line:
                print(name,cmt,line)
                coord = [float(x) for x in line.split('XB=')[1].split(',')[:6]]
                print("Before snapping:",coord)
                #print(abs(xs-coord[0]))
                #print(abs(xs-coord[1]))
                
                coord[0] = xs[np.where(np.round(abs(xs-coord[0]),3) == np.min(np.round(abs(xs-coord[0]),3)))[0][-1]]
                coord[1] = xs[np.where(np.round(abs(xs-coord[1]),3) == np.min(np.round(abs(xs-coord[1]),3)))[0][0]]
                coord[2] = ys[np.where(np.round(abs(ys-coord[2]),3) == np.min(np.round(abs(ys-coord[2]),3)))[0][-1]]
                coord[3] = ys[np.where(np.round(abs(ys-coord[3]),3) == np.min(np.round(abs(ys-coord[3]),3)))[0][0]]
                coord[4] = zs[np.where(np.round(abs(zs-coord[4]),3) == np.min(np.round(abs(zs-coord[4]),3)))[0][-1]]
                coord[5] = zs[np.where(np.round(abs(zs-coord[5]),3) == np.min(np.round(abs(zs-coord[5]),3)))[0][0]]
                
                #print(np.argmin(np.round(abs(xs-coord[0]),3)))
                #print(np.argmin(np.round(abs(xs-coord[1]),3)))
                '''
                coord[0] = xs[np.argmin(abs(xs-coord[0]))]
                coord[1] = xs[np.argmin(abs(xs-coord[1]))]
                coord[2] = ys[np.argmin(abs(ys-coord[2]))]
                coord[3] = ys[np.argmin(abs(ys-coord[3]))]
                coord[4] = zs[np.argmin(abs(zs-coord[4]))]
                coord[5] = zs[np.argmin(abs(zs-coord[5]))]
                '''
                print("After snapping:",coord)
                pts = [[coord[0]-dx*extend[0],coord[2]-dy*extend[1],coord[4]-dz*extend[2]],
                       [coord[0]-dx*extend[0],coord[2]-dy*extend[1],coord[5]+dz*extend[2]],
                       [coord[0]-dx*extend[0],coord[3]+dy*extend[1],coord[4]-dz*extend[2]],
                       [coord[0]-dx*extend[0],coord[3]+dy*extend[1],coord[5]+dz*extend[2]],
                       [coord[1]+dx*extend[0],coord[2]-dy*extend[1],coord[4]-dz*extend[2]],
                       [coord[1]+dx*extend[0],coord[2]-dy*extend[1],coord[5]+dz*extend[2]],
                       [coord[1]+dx*extend[0],coord[3]+dy*extend[1],coord[4]-dz*extend[2]],
                       [coord[1]+dx*extend[0],coord[3]+dy*extend[1],coord[5]+dz*extend[2]]]
                linkedPolygons.append(pts)
        polygons.append(linkedPolygons)
    return polygons

def parseFDSforPts(fileFDS,names,extend=[0,0,0],fileSMV=None):
    ''' This routine parses an FDS file looking for a list of names and
    stores a list of points defining polygons for each name which is
    found.
    
    The extend argument allows the polygon to be extended by a number of grid
    cells. This is useful since FDS snaps obstructions to the grid which means
    the actual coordinate location of the data from FDS may not align exactly
    with the input file.
    '''
    with open(fileFDS,'r') as f:
        lines = f.readlines()
    if fileSMV == None:
        for line in lines:
            if "&HEAD" in line:
                tmp = line.split("CHID")[1]
                tmp = tmp.split(',')[0]
                if "'" in tmp:
                    chid = tmp.split("'")[1]
                elif '"' in tmp:
                    chid = tmp.split('"')[1]
                else:
                    assert False, "Could not determine CHID from fds input file."
                tmp = fileFDS.split(os.sep)
                tmp[-1] = "%s.smv"%(chid)
                fileSMV = os.sep.join(tmp)
        print(fileSMV)

    with open(fileSMV,'r') as f:
        linesSMV = f.readlines()
    
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
    print(smvObjs.shape)
    comments = []
    for x in lines:
        sLine = x.replace('\n','').split('/')
        if len(sLine) <= 1:
            comments.append('')
        elif len(sLine) == 2:
            comments.append(sLine[1])
        else:
            comments.append('/'.join(sLine[1:]))
            
    polygons = []
    print(names)
    for name in names:
        linkedPolygons = []
        for line, cmt in zip(lines,comments):
            #print(name,cmt)
            if name in cmt and "&OBST" in line:
                #print(name,cmt,line)
                coord = np.array([float(x) for x in line.split('XB=')[1].split(',')[:6]])
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

def getMaxValueVersusTime(timeParameters,otherParameters,polygons):
    (initialTime,finalTime) = (timeParameters[0],timeParameters[1])
    (intervalTime,bandTime) = (timeParameters[2],timeParameters[3])
    (dataDir,chid,vIDs) = (otherParameters[0],otherParameters[1],otherParameters[2])
    (outName,orientations) = (otherParameters[3],otherParameters[4])
    (dLimited,dLimits) = (otherParameters[5],otherParameters[6])
    (saveOutputs,names) = (otherParameters[7],otherParameters[8])
    (outDir,vName) = (otherParameters[9],otherParameters[10])
    (chid) = (otherParameters[11])
    times = [initialTime]
    mPts = []
    while times[-1] + bandTime/2-intervalTime < finalTime:
        timeStart = max([times[-1]-bandTime/2,initialTime])
        timeEnd = min([times[-1]+bandTime/2,finalTime])
        oPts = []
        allPts = [[] for x in polygons]
        for orientation in orientations:
            pts = []
            for vID in vIDs:
                inName = buildFDS2asciiInput(dataDir,chid,vID,outName,
                                             timeStart,timeEnd,orientation,
                                             domainLimited=dLimited,
                                             domainLimits=dLimits)
                logfile = runFDS2ascii(dataDir,inName,outName)
                vpts = loadFDS2asciiPoints(dataDir+os.sep+outName)
                deleteFDS2asciiOutput(dataDir,inName,outName,logfile)
                if len(pts) == 0:
                    pts = vpts
                else:
                    pts.extend(vpts)
            pPts = []
            for i in range(0,len(polygons)):
                linkedpolygons = polygons[i]
                vpts = []
                for p in linkedpolygons:
                    for pt in pts:
                        if ut.pnt_in_cvex_hull(p, pt[:-1]):
                            vpts.append(pt)
                #vpts = np.array(vpts)
                if len(vpts) > 0:
                    mPt = np.max(np.array(vpts)[:,3])
                else:
                    mPt = -1
                pPts.append(mPt)
                allPts[i].append(vpts)
            oPts.append(pPts)
        oPts = np.array(oPts)
        if saveOutputs:
            for i in range(0,len(polygons)):
                thesePoints = []
                for pt in allPts[i]:
                    if len(thesePoints) == 0:
                        thesePoints = pt
                    else:
                        thesePoints.extend(pt)
                np.savetxt(outDir+os.sep+chid+'_'+names[i]+'_t_%.0f_%.0f.csv'%(timeStart,timeEnd),
                           np.squeeze(thesePoints),delimiter=',',header='X (m),Y (m), Z (m), %s\n'%(vName))
        
        mPts.append(np.max(oPts,axis=0))
        times.append(times[-1]+intervalTime)
    mPts = np.array(mPts)
    times = np.array(times)[:-1]
    return mPts, times

def getMeshes(params,fdsInputFile):
    meshes = params['meshes']
    print("INPUT FILE MESHES:",meshes)
    if meshes[0] == -1:
        with open(fdsInputFile,'r') as f:
            lines = f.readlines()
        meshes = []
        meshCounter = 0
        for line in lines:
            if '&MESH' in line:
                meshCounter = meshCounter+1
                meshes.append(int(meshCounter))
    print("OUTPUT MESHES:",meshes)
    return meshes

def processInputFileCustom(argsFile):
    params = readInputFile(argsFile)
    
    # Load run information
    outDir = getOutputDirectory(argsFile,params)
    dataDir = getDataDirectory(argsFile,params)
    chid, orientations, outputs, dLimited, dLimits = getRunInfo(params)
    fdsPath = getFdsPath(params)
    vName = getVariableName(params)
    outName = getOutName(params)
    
    # Load polygon information
    names, params = getNamesFromInput(params)
    thresholds = getThresholdsFromInput(params)
    if params['fdsInputFile']:
        fdsInputFile = os.path.abspath(os.sep.join(os.path.abspath(argsFile).split(os.sep)[:-1])+os.sep+params['fdsInputFile'])
        points = parseFDSforPts(fdsInputFile,names,extend=[0.0,0.0,0.0])
    else:
        print("FDS file not included in input file. Using points defined in input file.")
        points = getPolysFromInput(params)
    
    # Load meshes
    meshes = getMeshes(params,fdsInputFile)
    print(meshes)
    
    # Generate polygons
    polygons, numberOfGroups = ut.pts2polygons(points)
    
    # Load time information
    timeParameters = extractTimeParams(params)    
    tStart,tEnd,tBand,tInt = timeParameters
    # Run fds2ascii the first time to generate logfile
    #setupEnvironment(path=fdsPath)
    #inName = buildFDS2asciiInput(dataDir,chid,1,outName,
    #                             timeParameters[0],timeParameters[1],orientations[0],
    #                             domainLimited=dLimited,domainLimits=dLimits)
    #logfile = runFDS2ascii(dataDir,inName,outName)
    
    # Find variable ID
    vID = parseFDSforVID(fdsInputFile,vName)
    
    # Load boundary files
    coords, pts, times = fa.getPointsFromMeshes(dataDir,chid,meshes,vID,tStart,tEnd,tBand,tInt)
    
    # Generate point mask for polygons
    masks = getCoordinateMasks(coords,polygons)
    
    mPts = np.zeros((pts.shape[1],masks.shape[1]))
    for i in range(0,masks.shape[1]):
        mPts[:,i] = np.nanmax(pts[masks[:,i] == 1],axis=0)
    
    # Remove last time if it is zero
    if times[-1] == 0:
        mPts = np.delete(mPts,mPts.shape[0]-1,axis=0)
        times = np.delete(times,times.shape[0]-1,axis=0)
    
    # Make consistent plot colors
    pcs = []
    for i in range(0,numberOfGroups): pcs.append(pltc.rgb2hex(np.random.rand(3)))
    
    # Generate outputs
    namespace = "%s%s"%(outDir+os.sep,chid)
    
    return times, mPts, names, thresholds, polygons, namespace, params




def loadBNDFdata(argsFile,polygons):
    params = readInputFile(argsFile)
    
    # Load run information
    outDir = getOutputDirectory(argsFile,params)
    dataDir = getDataDirectory(argsFile,params)
    chid, orientations, outputs, dLimited, dLimits = getRunInfo(params)
    print(params)
    print(orientations)
    vName = getVariableName(params)
    
    # Load fdsInputFile Name
    fdsInputFile = os.path.abspath(os.sep.join(os.path.abspath(argsFile).split(os.sep)[:-1])+os.sep+params['fdsInputFile'])
    
    # Find polygon names
    names, params = getNamesFromInput(params,fdsInputFile)
    
    # Load polygon information
    thresholds = getThresholdsFromInput(params)
    
    # Load meshes
    meshes = getMeshes(params,fdsInputFile)
    print(meshes)
    
    # Load time information
    tStart,tEnd,tInt,tBand = extractTimeParams(params)    
    
    # Find variable ID
    vID = parseFDSforVID(fdsInputFile,vName)
    
    # Load boundary files
    coords2, pts2, times, orients2 = fa.getPointsFromMeshes(dataDir,chid,meshes,vID,tStart,tEnd,tBand,tInt)
    
    print(coords2.shape, pts2.shape, times.shape, orients2.shape)
    
    orientInds = np.where([True if x in orientations else False for x in orients2])[0]
    
    coords = coords2[orientInds,:]
    pts = pts2[orientInds,:]
    orients = orients2[orientInds]
    print(np.unique(orients))
    print(coords.shape, pts.shape, times.shape, orients.shape)
    
    # Generate point mask for polygons
    masks = getCoordinateMasks(coords,polygons)
    
    # Output specific polygon mask
    #nameCheck = ['Riser A-TTCY' in x for x in names]
    #ind = np.where(nameCheck)[0][0]
    #print("SAVING: E:\\projects\\khnp\\run0002\\POINT_COORDINATES.csv")
    #np.savetxt("E:\\projects\\khnp\\run0002\\POINT_COORDINATES.csv",coords[masks[:,ind] == 1],delimiter=',')
    #print("SAVING: E:\\projects\\khnp\\run0002\\POINT_VALUES.csv")
    #np.savetxt("E:\\projects\\khnp\\run0002\\POINT_VALUES.csv",pts[masks[:,ind] == 1],delimiter=',')
    
    mPts = np.zeros((pts.shape[1],masks.shape[1]))
    for i in range(0,masks.shape[1]):
        if np.where(masks[:,i] == 1)[0].shape[0] > 1:
            mPts[:,i] = np.nanmax(pts[masks[:,i] == 1],axis=0)
            np.savetxt('polygon%02.0f.csv'%(i),coords[masks[:,i] == 1,:],delimiter=',',header='x,y,z')
        else:
            mPts[:,i] = -1
    
    # Remove last time if it is zero
    if times[-1] == 0:
        mPts = np.delete(mPts,mPts.shape[0]-1,axis=0)
        times = np.delete(times,times.shape[0]-1,axis=0)
    
    # Generate outputs
    namespace = "%s%s"%(outDir+os.sep,chid)
    
    return times, mPts, names, thresholds, namespace, params




def getArgFilePolygons(argsFile):
    params = readInputFile(argsFile)
    
    # Load run information
    chid, orientations, outputs, dLimited, dLimits = getRunInfo(params)
    
    # Load fdsInputFile Name
    fdsInputFile = os.path.abspath(os.sep.join(os.path.abspath(argsFile).split(os.sep)[:-1])+os.sep+params['fdsInputFile'])
    
    # Find polygon names
    names, params = getNamesFromInput(params,fdsInputFile)
    
    # Load polygon information
    points = parseFDSforPts(fdsInputFile,names,extend=[0.0,0.0,0.0])
    
    # Generate polygons
    polygons, numberOfGroups = ut.pts2polygons(points)
    
    # Generate color scheme for each polygon
    pcs = []
    for i in range(0,numberOfGroups): pcs.append(pltc.rgb2hex(np.random.rand(3)))
    
    return polygons, pcs


def getSmvFile(argsFile):
    params = readInputFile(argsFile)
    dataDir = getDataDirectory(argsFile,params)
    chid, orientations, outputs, dLimited, dLimits = getRunInfo(params)
    smvFile = os.path.abspath(dataDir+os.sep+chid+'.smv')
    return smvFile

def getCoordinateMasks(coords,polygons):
    masks = np.zeros((coords.shape[0],len(polygons)))
    for i in range(0,len(polygons)):
        linkedpolygons = polygons[i]
        for p in linkedpolygons:
            masks[np.where(in_hull(coords,p.points)),i] = 1
    return masks

def in_hull(p,hull):
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0

def getCoordinateMasks2(coords,polygons):
    masks = np.zeros((coords.shape[0],len(polygons)))
    for i in range(0,len(polygons)):
        linkedpolygons = polygons[i]
        for p in linkedpolygons:
            for j in range(0,coords.shape[0]):
                if ut.pnt_in_cvex_hull(p, coords[j,:]):
                    masks[j,i] = 1
    return masks

def processInputFileReturn(argsFile):
    params = readInputFile(argsFile)
    
    # Load run information
    outDir = getOutputDirectory(argsFile,params)
    dataDir = getDataDirectory(argsFile,params)
    chid, orientations, outputs, dLimited, dLimits = getRunInfo(params)
    fdsPath = getFdsPath(params)
    vName = getVariableName(params)
    outName = getOutName(params)
    
    # Load polygon information
    names, params = getNamesFromInput(params)
    thresholds = getThresholdsFromInput(params)
    if params['fdsInputFile']:
        fdsInputFile = os.path.abspath(os.sep.join(os.path.abspath(argsFile).split(os.sep)[:-1])+os.sep+params['fdsInputFile'])
        points = parseFDSforPts(fdsInputFile,names,extend=[0.5,0.5,0.5])
    else:
        print("FDS file not included in input file. Using points defined in input file.")
        points = getPolysFromInput(params)
    
    # Generate polygons
    polygons, numberOfGroups = ut.pts2polygons(points)
    
    # Load time information
    timeParameters = extractTimeParams(params)    
    
    # Run fds2ascii the first time to generate logfile
    setupEnvironment(path=fdsPath)
    inName = buildFDS2asciiInput(dataDir,chid,1,outName,
                                 timeParameters[0],timeParameters[1],orientations[0],
                                 domainLimited=dLimited,domainLimits=dLimits)
    logfile = runFDS2ascii(dataDir,inName,outName)
    
    # Find variable ID in logfile
    checkForLogfile(logfile)
    vIDs = findVariableInLog(logfile,vName)
    vID = checkVIDs(vIDs,logfile,vName)
    
    # Call fds2ascii for each time and orientation configuration
    saveOutputs = True if outputs['allPoints'] else False
    otherParameters = [dataDir,chid,vID,outName,orientations,dLimited,dLimits,
                       saveOutputs,names,outDir,vName,chid]
    mPts,times = getMaxValueVersusTime(timeParameters,otherParameters,polygons)
    
    # Make consistent plot colors
    pcs = []
    for i in range(0,numberOfGroups): pcs.append(pltc.rgb2hex(np.random.rand(3)))
    
    # Generate outputs
    namespace = "%s%s"%(outDir+os.sep,chid)
    
    smvFile = os.path.abspath(dataDir+os.sep+chid+'.smv')
    surfaces, obstructions = ut.buildSMVgeometry(smvFile)
    
    return times, mPts, names, thresholds, polygons, surfaces, obstructions, namespace, params

def processInputFile(argsFile):
    params = readInputFile(argsFile)
    
    # Load run information
    outDir = getOutputDirectory(argsFile,params)
    dataDir = getDataDirectory(argsFile,params)
    chid, orientations, outputs, dLimited, dLimits = getRunInfo(params)
    fdsPath = getFdsPath(params)
    vName = getVariableName(params)
    outName = getOutName(params)
    
    # Load polygon information
    names = getNamesFromInput(params)
    thresholds = getThresholdsFromInput(params)
    if params['fdsInputFile']:
        fdsInputFile = os.path.abspath(os.sep.join(os.path.abspath(argsFile).split(os.sep)[:-1])+os.sep+params['fdsInputFile'])
        points = parseFDSforPts(fdsInputFile,names,extend=[0.5,0.5,0.5])
    else:
        print("FDS file not included in input file. Using points defined in input file.")
        points = getPolysFromInput(params)
    
    # Generate polygons
    polygons, numberOfGroups = ut.pts2polygons(points)
    
    # Load time information
    timeParameters = extractTimeParams(params)    
    
    # Run fds2ascii the first time to generate logfile
    setupEnvironment(path=fdsPath)
    inName = buildFDS2asciiInput(dataDir,chid,1,outName,
                                 timeParameters[0],timeParameters[1],orientations[0],
                                 domainLimited=dLimited,domainLimits=dLimits)
    logfile = runFDS2ascii(dataDir,inName,outName)
    
    # Find variable ID in logfile
    checkForLogfile(logfile)
    vIDs = findVariableInLog(logfile,vName)
    vID = checkVIDs(vIDs,logfile)
    
    # Call fds2ascii for each time and orientation configuration for each mesh
    otherParameters = [dataDir,chid,vID,outName,orientations,dLimited,dLimits]
    mPts,times = getMaxValueVersusTime(timeParameters,otherParameters,polygons)
    
    # Make consistent plot colors
    pcs = []
    for i in range(0,numberOfGroups): pcs.append(pltc.rgb2hex(np.random.rand(3)))
    
    # Generate outputs
    namespace = "%s%s"%(outDir+os.sep,chid)
    if outputs['maxTCSV']: ut.maxValueCSV(times,mPts,names,namespace)
    if outputs['failureTime']: ut.failureTimesCSV(times,mPts,names,thresholds,namespace)
    if outputs['maxTPlot']: ut.maxValuePlot(times,mPts,names,namespace,pcs=pcs,vName=vName)
    if outputs['polygonVisual']:
        app, w = setupQT()
        for i in range(0,len(polygons)):
            for j in range(0,len(polygons[i])):
                p = polygons[i][j]
                c = np.array([pltc.to_rgba(pcs[i]) for x in p.simplices])
                m1 = gl.GLMeshItem(vertexes=p.points, faces=p.simplices, faceColors=c, smooth=False)
                m1.setGLOptions('opaque')
                w.addItem(m1)
        
        smvFile = os.path.abspath(dataDir+os.sep+chid+'.smv')
        surfaces, obstructions = ut.buildSMVgeometry(smvFile)
        
        for i in range(0,len(obstructions)):#obst in obstructions:
            obst = obstructions[i]
            pts, colors = ut.getPtsFromObst(obst,surfaces)
            for pt, color in zip(pts,colors):
                c = np.array([[color[0],color[1],color[2],color[3]] for x in color])
                p = np.array([[x[0],x[1],x[2]] for x in pt])
                faces = np.array([[0,1,2],[0,2,3]])
                m1 = gl.GLMeshItem(vertexes=p, faces=faces, faceColors=c, smooth=False)
                if c[0][3] < 0.1:
                    m1.setVisible(False)
                else:
                    m1.setGLOptions('opaque')
                w.addItem(m1)
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()
    
    
if __name__ == '__main__':
    
    app = QApplication(sys.argv)

    w = QWidget()
    w.resize(250, 150)
    w.move(300, 300)
    w.setWindowTitle('Simple')
    w.show()
    
    sys.exit(app.exec_())
    
    # Read input file
    args = sys.argv
    systemDir, defaultParamFile = determineSystemDir(args)
    argsFile = checkCommandLineArgs(args)
    
    # Process input file
    print(argsFile)
    processInputFile(argsFile)

        
    '''
    app = QApplication(sys.argv)
    window = QDialog()
    ui = Ui_MainWindow()
    ui.setupUi(window)
    
    window.show()
    sys.exit(app.exec_())
    '''
    
    
    #fig, ax = ut.polygonVisual(polygons,namespace,pcs=pcs)
    #fig, ax = ut.smvVisual(obstructions,surfaces,namespace,fig=fig,ax=ax)
    #plt.show()