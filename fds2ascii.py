# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 10:04:29 2018

@author: jhodges
"""

#import numpy as np
#import glob
#import os
#import matplotlib.pyplot as plt
#import subprocess as sp

class fdspatch(object):
    def __init__(self,NX,NY,NT,DS,OR):
        self.data = np.zeros((NX,NY,NT))
        self.lims = DS
        self.orientation = OR
    def append(self,data,iT):
        self.data[:,:,iT] = data
    def average(self,inds):
        return np.mean(self.data[:,:,inds],axis=2)[:,:,0]
    def buildSpace(self):
        (NX,NY) = (self.data.shape[0],self.data.shape[1])
        if self.lims[0] == self.lims[1]:
            xGrid = np.zeros((NX,NY))+self.lims[0]
            y = np.linspace(self.lims[2],self.lims[3],NX)
            z = np.linspace(self.lims[4],self.lims[5],NY)
            yGrid,zGrid = np.meshgrid(y,z)
        elif self.lims[2] == self.lims[3]:
            yGrid = np.zeros((NX,NY))+self.lims[2]
            x = np.linspace(self.lims[0],self.lims[1],NX)
            z = np.linspace(self.lims[4],self.lims[5],NY)
            xGrid,zGrid = np.meshgrid(x,z)
        elif self.lims[4] == self.lims[5]:
            zGrid = np.zeros((NX,NY))+self.lims[4]
            x = np.linspace(self.lims[0],self.lims[1],NX)
            y = np.linspace(self.lims[2],self.lims[3],NY)
            xGrid,yGrid = np.meshgrid(x,y)
        self.x = xGrid
        self.y = yGrid
        self.z = zGrid
    def extractPoints(self):
        pts = np.zeros((self.data.shape[0]*self.data.shape[1],self.data.shape[2]))
        for iT in range(0,self.data.shape[2]):
            pts[:,iT] = self.data[:,:,iT].flatten()
        coords = np.zeros((self.data.shape[0]*self.data.shape[1],3))
        coords[:,0] = self.x.flatten()
        coords[:,1] = self.y.flatten()
        coords[:,2] = self.z.flatten()
        orients = np.zeros((self.data.shape[0]*self.data.shape[1],))+self.orientation
        return coords, pts, orients
        
def bytesToNum(byte):
    byte = byte.decode('utf-8')
    byte = byte.replace('\x00','0')
    byte = byte.replace('\x01','1')
    byte = byte.replace('\x02','2')
    byte = byte.replace('\x03','3')
    byte = byte.replace('\x04','4')
    byte = byte.replace('\x05','5')
    byte = byte.replace('\x06','6')
    byte = byte.replace('\x07','7')
    byte = byte.replace('\x08','8')
    byte = byte.replace('\x09','9')
    return byte

def readBoundaryFile(fname):
    if not os.path.isfile(fname):
        return None, None, None
    
    with open(fname,'rb') as f:
        bytes_read = f.read(110)
        tmp = bytes_read.split(b'\x1e')

        quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
        shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
        units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
        
        data = np.fromfile(f,dtype=np.int32,count=5)
        npatch = data[2]
        #print(quantity, shortName, units, npatch)
        pts = 0
        patchPts = []
        patchDs = []
        patchIors = []
        for k in range(0,npatch):
            data = np.fromfile(f,dtype=np.int32,count=11)
            #data2 = np.fromfile(f,dtype=np.int32,count=2)
            #data3 = np.fromfile(f,dtype=np.int8,count=8)
            #data2 = np.fromfile(f,dtype=np.int32,count=4)
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
            #print(data,data2,data3,nPts)
            patchPts.append([pts,pts+nPts])
            patchDs.append([dx,dy,dz,data[0],data[1],data[2],data[3],data[4],data[5]])
            patchIors.append(data[6])
            pts = pts+nPts
        data = np.fromfile(f,dtype=np.float32,count=-1)
    bndfInfo = (quantity,shortName,units,npatch)
    patchInfo = (patchPts,patchDs,patchIors)
    return bndfInfo, patchInfo, data

def getLimsFromGrid(data,grid):
    #if len(data[0].shape) > 0:
    try:
        #print(grid)
        #print(data)
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
        print("Failed to make lims, returning null")
        '''
        print(data)
        print("Grid0:")
        print(grid[0])
        print("Grid1:")
        print(grid[1])
        print("Grid2:")
        print(grid[2])
        print(len(grid))
        print( np.argwhere(grid[0][:,0] == data[0]))
        print(ind1)
        print(ind2)
        print(ind3)
        print(ind4)
        print(ind5)
        print(ind6)
        assert False, "Stopped"
        '''
        lims = [0,0,0,0,0,0]
    #else:
    #    print(data)
    #    assert False, "Stopped"
    #    lims = [0,0.01,0,0.01,0,0.0]
    return lims

def getPatchObstId(data,obst,ior):
    xmn,xmx,ymn,ymx,zmn,zmx = data
    data = np.array(data)
    #print(obst.shape,data.shape)
    #print("OBST:")
    #print(obst[:,:6])
    print("DATA:")
    print(data)
    diff = abs(obst[:,:6]-data)
    #print("DIFF:")
    #print(diff)
    #print(ior)
    if ior == 1:
        dist = diff[:,2]+diff[:,3]+diff[:,4]+diff[:,5]+diff[:,1]
    elif ior == -1:
        dist = diff[:,2]+diff[:,3]+diff[:,4]+diff[:,5]+diff[:,2]
    elif ior == 2:
        dist = diff[:,0]+diff[:,1]+diff[:,4]+diff[:,5]+diff[:,3]
    elif ior == -2:
        dist = diff[:,0]+diff[:,1]+diff[:,4]+diff[:,5]+diff[:,2]
    elif ior == 3:
        dist = diff[:,0]+diff[:,1]+diff[:,2]+diff[:,3]+diff[:,5]
    elif ior == -3:
        dist = diff[:,0]+diff[:,1]+diff[:,2]+diff[:,3]+diff[:,4]
    ind = np.argmin(abs(dist))
    print("Best OBST:")
    print(obst[ind,:6])
    print("Best DIST:",dist[ind])
    print("OBST ID:",obst[ind,6])
    return obst[ind,6]

def buildPatches(patchPts,patchDs,patchIors,data,grid,obst):
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
                lims = getLimsFromGrid(patchDs[k][3:],grid)
                #patchId = getPatchObstId(lims,obst,patchIors[k])
                patches.append(fdspatch(patchData.shape[0],patchData.shape[1],timeSteps-1,lims,patchIors[k]))
                patches[k].append(patchData,i-1)
            else:
                patches[k].append(patchData,i-1)
    return times, patches

def extractTime(tStart,tEnd,tBand,tInt,patches,times):
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

def buildSpace(patches):
    for i in range(0,len(patches)):
        patches[i].buildSpace()

def buildMesh_old(smvFile,mesh):
    with open(smvFile,'r') as f:
        file = f.read()
    gFile = file.split('GRID')
    for i in range(1,len(gFile)):
        meshStr = gFile[i].split('\n')[0].split()[0]
        passed = False
        while not passed:
            try:
                float(meshStr)
                passed = True
            except:
                meshStr=meshStr[1:]
                if len(meshStr) == 0:
                    passed = True
                    meshStr = mesh
        if float(meshStr) == float(mesh):
            lines = gFile[i].split('\n')
    inds = []
    for i in range(0,len(lines)):
        if len(lines[i]) > 0:
            if lines[i][0] != ' ' and lines[i][0] != '-':
                inds.append(i)
    for i in range(0,len(inds)):
        ind = inds[i]
        if 'TRNX' in lines[ind]:
            xGrid = np.array([[float(y) for y in x.split()] for x in lines[inds[i]+2:inds[i+1]-1]])
        if 'TRNY' in lines[ind]:
            yGrid = np.array([[float(y) for y in x.split()] for x in lines[inds[i]+2:inds[i+1]-1]])
        if 'TRNZ' in lines[ind]:
            zGrid = np.array([[float(y) for y in x.split()] for x in lines[inds[i]+2:inds[i+1]-1]])
        if 'OBST' in lines[ind] and '_OBST' not in lines[ind]:
            num = int(lines[inds[i]+1].split()[0])
            tmp = lines[inds[i]+2:inds[i]+2+num]
            obst = np.array([[float(y) for y in x.split()] for x in tmp])
    grids = [xGrid,yGrid,zGrid]
    return grids, obst

def buildMesh(smvFile):
    with open(smvFile,'r') as f:
        linesSMV = f.readlines()
    
    grids = []
    obsts = []
    bndfs = []
    for i in range(0,len(linesSMV)):
        line2 = linesSMV[i]
        if "GRID" in line2:
            gridPts = [int(x) for x in linesSMV[i+1].replace('\n','').split()]
            gridTRNX = np.array([[float(y) for y in x.replace('\n','').split()] for x in linesSMV[i+8:i+9+gridPts[0]]])
            gridTRNY = np.array([[float(y) for y in x.replace('\n','').split()] for x in linesSMV[i+12+gridPts[0]:i+13+gridPts[0]+gridPts[1]]])
            gridTRNZ = np.array([[float(y) for y in x.replace('\n','').split()] for x in linesSMV[i+16+gridPts[0]+gridPts[1]:i+17+gridPts[0]+gridPts[1]+gridPts[2]]])
            grids.append([gridTRNX.copy(),gridTRNY.copy(),gridTRNZ.copy()])
            dx = (gridTRNX.max()-gridTRNX.min())/(gridTRNX.shape[0]-1)
            dy = (gridTRNY.max()-gridTRNY.min())/(gridTRNY.shape[0]-1)
            dz = (gridTRNZ.max()-gridTRNZ.min())/(gridTRNZ.shape[0]-1)
        if "OBST" in line2 and "HIDE_OBST" not in line2:
            numOBST = int(linesSMV[i+1].replace(' ',''))
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
            obsts.append(smvObj)
        if ".bf" in line2:
            (_,mesh,_) = linesSMV[i-1].split()
            bndfName = linesSMV[i].split(' ')[1].replace('\n','')
            vNum = bndfName.split('_')[-1].replace('.bf','')
            vID = ' '.join(linesSMV[i+1].split(' ')[1:]).replace('\n','')
            bndfs.append([float(mesh),bndfName,vID,float(vNum)])
    #grid = grids[mesh]
    #obst = obsts[mesh]
    return grids, obsts, bndfs

def buildMesh2(smvFile):
    with open(smvFile,'r') as f:
        lines = f.readlines()
    inds = []
    for i in range(0,len(lines)):
        if lines[i][0] != ' ' and lines[i][0] != '-':
            inds.append(i)
    for i in range(0,len(inds)):
        ind = inds[i]
        if 'TRNX' in lines[ind]:
            xGrid = np.array([[float(y) for y in x.split()] for x in lines[inds[i]+2:inds[i+1]]])
        if 'TRNY' in lines[ind]:
            yGrid = np.array([[float(y) for y in x.split()] for x in lines[inds[i]+2:inds[i+1]]])
        if 'TRNZ' in lines[ind]:
            zGrid = np.array([[float(y) for y in x.split()] for x in lines[inds[i]+2:inds[i+1]]])
    grids = [xGrid,yGrid,zGrid]
    return grids

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

def getPatchesFromMesh(grid,obst,bndfFile):
    bndfInfo, patchInfo, data = readBoundaryFile(bndfFile)
    if bndfInfo == None:
        return [None], [None]
    (quantity,shortName,units,npatch) = bndfInfo
    (patchPts,patchDs,patchIors) = patchInfo
    times, patches = buildPatches(patchPts,patchDs,patchIors,data,grid,obst)
    return times, patches

def savePatches(dataDir,chid,times,patches):
    for j in range(len(patches)):
        p = patches[j]
        for i in range(0,len(times)):
            t = times[i]
            nome = "%s%s_p%02.0f_t_%04.4f"%(dataDir,chid,j,t)
            nome = nome.replace('.','p')
            np.savetxt("%s.csv"%(nome),p.data[:,:,i],delimiter=',')

def getPointsFromMeshes(dataDir,chid,meshes,vID,tStart,tEnd,tBand,tInt):
    allCoords = []
    allPts = []
    allOrients = []
    newTimes = []
    
    smvFile = os.path.abspath("%s%s%s.smv"%(dataDir,os.sep,chid))
    grids, obsts, bndfs = buildMesh(smvFile)
    
    for mesh in meshes:
        for bndf in bndfs:
            if (float(bndf[0]) == float(mesh)) and (float(vID) == float(bndf[3])):
                bndfFile = "%s%s%s"%(dataDir,os.sep,bndf[1])
                times, patches = getPatchesFromMesh(grids[mesh-1],obsts[mesh-1],bndfFile)
                if len(times) > 1:
                    newTimes, newPatches = extractTime(tStart,tEnd,tBand,tInt,patches,times)
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
        print(meshes)
        print(bndfs)
        assert False, "No valid patches found."
    return allCoords, allPts, newTimes, allOrients

def getPointsFromFiles(bndfs, grids, obsts, tStart, tEnd, tBand, tInt):
    allCoords = []
    allPts = []
    allOrients = []
    newTimes = []
    
    for file, mesh in bndfs:
        mesh = int(mesh)
        times, patches = getPatchesFromMesh(grids[mesh-1], obsts[mesh-1], file)
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
        print(meshes)
        print(bndfs)
        assert False, "No valid patches found."
    return allCoords, allPts, newTimes, allOrients

def buildFDS2asciiInput(outDir,chid,variableID,outName,
                        timeStart,timeEnd,orientation,
                        fileType=3,downsampleRate=1,
                        domainLimited=False,domainLimits=[0,1,0,1,0,1],
                        fname='fds2asciiInputFile.txt'):
    '''
    This subroutine builds the fds2ascii input file and returns the filename.
    '''
    front = "%s\n%.0f\n%.0f\n"%(chid,fileType,downsampleRate)
    if domainLimited:
        (xmn,xmx) = (domainLimits[0],domainLimits[1])
        (ymn,ymx) = (domainLimits[2],domainLimits[3])
        (zmn,zmx) = (domainLimits[4],domainLimits[5])
        dls = "%s\n%.2f %.2f %.2f %.2f %.2f %.2f\n"%('y',
                                                     xmn,xmx,ymn,ymx,zmn,zmx)
    else:
        dls = "%s\n"%('n')
    ts = "%.6f\n%.6f\n"%(timeStart,timeEnd)
    back = "%.0f\n%.0f\n%.0f\n%s"%(orientation,1,variableID,outName)
    string = "%s%s%s%s"%(front,dls,ts,back)
    with open(outDir+os.sep+fname,'w') as f:
        f.write(string)
    return fname

def runFDS2ascii(dataDir,inputFile,outputFile,logname='fds2ascii.log'):
    '''
    This subroutine runs fds2ascii and returns the log filename.
    '''
    p = sp.Popen(['fds2ascii.exe','<',inputFile,'>',logname],
                 cwd=dataDir,shell=True)
    p.wait()
    return dataDir+os.sep+logname

def loadFDS2asciiPoints(fName):
    '''
    Reads and returns all points and the value from the fds2ascii output file.
    '''
    with open(fName,'r') as f:
        lines = f.readlines()
    pts = []
    for line in lines:
        s = line.split(',')
        try:
            (x,y,z,v) = (float(s[0]),float(s[1]),float(s[2]),float(s[3]))
            pts.append([x,y,z,v])
        except:
            pass
    return pts

def deleteFDS2asciiOutput(dataDir,inName,outName,logfile):
    '''
    Delete fds2ascii input and output after loading. This is important so that
    an old file is not read in for a future time step.
    '''
    os.remove(dataDir+os.sep+inName)
    os.remove(dataDir+os.sep+outName)
    os.remove(logfile)

def findVariableInLog(logfile,vName):
    '''
    Searches fds2ascii logfile for the string variable name and returns the
    variable IDs associated with it.
    '''
    with open(logfile,'r') as f:
        lines = f.readlines()
    vIDs = []
    for line in lines:
        if vName in line:
            vIDs.append(int(line.split()[0]))
        else:
            pass
    return vIDs

def setupEnvironment(path=r"C:\Program Files\FDS\FDS6\bin"):
    '''
    This subroutine changes the path to use the specified FDS directory.
    '''
    oldpath = os.environ['PATH']
    os.environ['PATH'] = path + ";" + oldpath
    os.environ['OMP_NUM_THREADS'] = '1'

def getFdsPath(params):
    ''' Return fds path from parameter dictionary '''
    if params['fdsPath']:
        fdsPath = params['fdsPath']
    else:
        fdsPath = r"C:\Program Files\FDS\FDS6\bin"
    return fdsPath

def getVariableName(params):
    ''' Return fds variable name from parameter dictionary '''
    vName = params['variableName']
    return vName

def checkForLogfile(file):
    if not os.path.isfile(file):
        assert False, "Could not find %s, did FDS2ascii run correctly?"%(file)

def getOutName(params):
    if params['outName']:
        outName = params['outName']
    else:
        outName = 'f2aoutput.csv'
    return outName

def checkVIDs(vIDs,logfile,vName):
    ''' Make sure variable ID is found in the log file. If it is found
    multiply times error since multiple meshes have not been implemented yet.
    '''
    if len(vIDs) > 1:
        print("Warning: Multiple meshes not supported at this time.")
        vID = vIDs
    elif len(vIDs) == 0:
        assert False, "Could not find %s in %s"%(vName,logfile)
    else:
        vID = vIDs[0]
    return vID

if __name__ == "__main__":
    case = 19
    GHF_version = 2
    if case == 0:
        dataDir = "../../SFPE/simulations/"
        chid = "preliminaryStudy_nominal"
        vID = 7
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
    elif case == 1:
        dataDir = "F:\\railcarScaling\\nfpa286\\"
        chid = "nfpa286_w00_fs_inert"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 45
    elif case == 2:
        dataDir = "F:\\railcarScaling\\nfpa286\\"
        chid = "nfpa286_w00_p5s_inert"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 46
    elif case == 3:
        dataDir = "F:\\railcarScaling\\nfpa286\\"
        chid = "nfpa286_w00_p2s_inert"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 45
    elif case == 4:
        dataDir = "F:\\railcarScaling\\nfpa286\\"
        chid = "nfpa286_w02_fs_inert"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 46
    elif case == 5:
        dataDir = "F:\\railcarScaling\\nfpa286\\"
        chid = "nfpa286_w02_p5s_inert"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 47
    elif case == 6:
        dataDir = "F:\\railcarScaling\\nfpa286\\"
        chid = "nfpa286_w02_p2s_inert"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 45
    elif case == 7:
        dataDir = "F:\\railcarScaling\\nfpa286\\FullScale\\"
        chid = "nfpa286_w00_fs_inert_q1000"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 45
    elif case == 8:
        dataDir = "F:\\railcarScaling\\nfpa286\\FullScale\\"
        chid = "nfpa286_w02_fs_inert_q1000"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 46
    elif case == 9:
        dataDir = "F:\\railcarScaling\\nfpa286\\FullScale\\"
        chid = "nfpa286_w00_fs_inert_q3000"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 32
    elif case == 10:
        dataDir = "F:\\railcarScaling\\nfpa286\\FullScale\\"
        chid = "nfpa286_w02_fs_inert_q3000"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 32
        
    elif case == 11:
        dataDir = "F:\\railcarScaling\\nfpa286\\HalfScale\\"
        chid = "nfpa286_w00_p5s_inert_q250"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 46
    elif case == 12:
        dataDir = "F:\\railcarScaling\\nfpa286\\HalfScale\\"
        chid = "nfpa286_w02_p5s_inert_q250"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 47
    elif case == 13:
        dataDir = "F:\\railcarScaling\\nfpa286\\HalfScale\\"
        chid = "nfpa286_w00_p5s_inert_q750"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 37
    elif case == 14:
        dataDir = "F:\\railcarScaling\\nfpa286\\HalfScale\\"
        chid = "nfpa286_w02_p5s_inert_q750"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 39
        
    elif case == 15:
        dataDir = "F:\\railcarScaling\\nfpa286\\FifthScale\\outputw00\\"
        chid = "nfpa286_w00_p2s_inert_q40"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 45
    elif case == 16:
        dataDir = "F:\\railcarScaling\\nfpa286\\FifthScale\\outputw02\\"
        chid = "nfpa286_w02_p2s_inert_q40"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 45
    elif case == 17:
        dataDir = "F:\\railcarScaling\\nfpa286\\FifthScale\\outputw00\\"
        chid = "nfpa286_w00_p2s_inert_q120"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 32
    elif case == 18:
        dataDir = "F:\\railcarScaling\\nfpa286\\FifthScale\\outputw02\\"
        chid = "nfpa286_w02_p2s_inert_q120"
        vID = 4
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 32
    elif case == 19:
        dataDir = "E:\\projects\\khnp\\run0008\\"
        chid = "khnp_run0008"
        vID = 1
    
        meshes = [1]
        (tStart,tEnd,tBand,tInt) = (0,10,3,1)
        patchID = 32
        
    allCoords, allPts, newTimes = getPointsFromMeshes(dataDir,chid,meshes,vID,tStart,tEnd,tBand,tInt)
    
    mesh = meshes[0]
    vIDstr = "%02.0f"%(vID)
    
    if mesh == 1:
        files = glob.glob(dataDir+chid+"_*_*.bf")
        if len(files) == 0:
            meshStr = ""
        else:
            meshStr = "_%04.0f"%(mesh)
    else:
        meshStr = "_%04.0f"%(mesh)
    bndfFile = "%s%s%s%s_%s.bf"%(dataDir,os.sep,chid,meshStr,vIDstr)
    bndfFile = os.path.abspath(bndfFile)
    smvFile = os.path.abspath("%s%s%s.smv"%(dataDir,os.sep,chid))
    bndfInfo, patchInfo, data = readBoundaryFile(bndfFile)
    (quantity,shortName,units,npatch) = bndfInfo
    (patchPts,patchDs,patchIors) = patchInfo
    grid, obst, bndfs = buildMesh(smvFile)
    times, patches = buildPatches(patchPts,patchDs,patchIors,data,grid[1],obst[1])
    
    if GHF_version == 1:
        ghf = np.max(np.max(patches[patchID].data,axis=0),axis=0)
        ghfName = "_FLOOR_GHF"
    elif GHF_version == 2:
        sz = patches[patchID].data.shape
        xDir = int(np.round(sz[1]/4,decimals=0))
        yDir = int(np.round(sz[0]/2,decimals=0))
        ghf = patches[patchID].data[yDir,xDir,:]
        ghfName = "_FLOOR_GHF_v2"
    #saveData = np.array([times,ghf]).T
    #np.savetxt(dataDir+chid+ghfName+".csv",saveData,delimiter=',')
    #print("MADE: %s"%(dataDir+chid+ghfName+".csv"))
    
    #times, patches = getPatchesFromMesh(dataDir,chid,mesh,vID)
    #newTimes, newPatches = extractTime(tStart,tEnd,tBand,tInt,patches,times)
    #buildSpace(newPatches)
    
    #coords, pts = extractPoints(newPatches)
    
