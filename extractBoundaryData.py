import sys
import os
import matplotlib.colors as pltc
import matplotlib.pyplot as plt
import numpy as np
import yaml
from collections import defaultdict
import glob

class fdspatch(object):
    def __init__(self,NX,NY,NT,DS,OR):
        self.data = np.zeros((NX,NY,NT))
        self.lims = DS
        self.orientation = OR
    def append(self,data,iT):
        self.data[:,:,iT] = data
    def average(self,inds):
        return np.mean(self.data[:,:,inds],axis=2)
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

def parseBndfHeader(f):
    bytes_read = f.read(110)
    tmp = bytes_read.split(b'\x1e')

    quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
    shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
    units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
    data = np.fromfile(f,dtype=np.int32,count=5)
    npatch = data[2]
    return quantity, shortName, units, npatch

def parseBndfPatches(f, npatch):
    pts = 0
    patchPts = []
    patchDs = []
    patchIors = []
    for k in range(0, npatch):
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
    patchInfo = (patchPts, patchDs, patchIors)
    return patchInfo, data

def readBoundaryHeader(file):
    with open(file, 'rb') as f:
        quantity, shortName, units, npatch = parseBndfHeader(f)
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

def importBoundaryFile(fname, smvFile):
    if not os.path.isfile(fname):
        return None, None
    
    with open(fname,'rb') as f:
        quantity, shortName, units, npatch = parseBndfHeader(f)
        patchInfo, data = parseBndfPatches(f, npatch)
        
    (patchPts,patchDs,patchIors) = patchInfo
    grid, obst, bndfs = buildMesh(smvFile)
    times, patches = buildPatches(patchPts,patchDs,patchIors,data,grid[0],obst[0])
    return times, patches

def buildPatches(patchPts, patchDs, patchIors, data, grid, obst):
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
        print("Failed to make lims, returning null")
        lims = [0,0,0,0,0,0]
    return lims

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

def readInputFile(file):
    ''' Load input file '''
    params = defaultdict(bool,yaml.load(open(file,'r')))
    return params

def getPolygonNames(params, obsts):
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
    else:
        names = []
        for key in list(fdsObsts.keys()):
            if fdsObsts[key]['BNDF_OBST']:
                names.append(fdsObsts[key]['ID'])
        params['polygons'] = names
    return names, params

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

def loadBNDFdata(params, bndfs, smvGrids, smvObsts):
    tStart, tEnd, tBand, tInt = extractTimeParams(params)
    coords2, pts2, times, orients2 = f2a.getPointsFromFiles(bndfs, smvGrids, smvObsts, tStart, tEnd, tBand, tInt)
    orientInds = np.where([True if x in params['orientations'] else False for x in orients2])[0]
    
    coords = coords2[orientInds,:]
    pts = pts2[orientInds,:]
    orients = orients2[orientInds]
    
    # Generate point mask for polygons
    masks = getCoordinateMasks(coords,polygons)
    
    mPts = np.zeros((pts.shape[1],masks.shape[1]))
    for i in range(0,masks.shape[1]):
        if np.where(masks[:,i] == 1)[0].shape[0] > 1:
            mPts[:,i] = np.nanmax(pts[masks[:,i] == 1],axis=0)
            if params['savePolygons']:
                np.savetxt('polygon%02.0f.csv'%(i),coords[masks[:,i] == 1,:],delimiter=',',header='x,y,z')
        else:
            mPts[:,i] = -1
    
    # Remove last time if it is zero
    if times[-1] == 0:
        mPts = np.delete(mPts,mPts.shape[0]-1,axis=0)
        times = np.delete(times,times.shape[0]-1,axis=0)
        
    return times, mPts, orients

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
            masks[np.where(in_hull2(coords,p.points)),i] = 1
    return masks

def in_hull2(p,hull):
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0

if __name__ == '__main__':
    systemDir = os.path.abspath(os.path.abspath(sys.argv[0])).replace(os.sep+'extractBoundaryData.py','').replace('extractBoundaryData.exe','')
    #defaultParamFile, serverParamFile, configPath, databasePath = determineLocalPaths(systemDir)
    inputFile = "E:\\projects\\1MJP00013.000.001\\fromcluster\\Bay708B_run0001.yaml"
    inputFile = os.path.abspath(inputFile)
    configDir = os.sep.join(inputFile.split(os.sep)[:-1])
    params = readInputFile(inputFile)
    params['configDir'] = configDir
    params['smvfile'] = os.path.abspath(configDir+os.sep+params['chid']+'.smv')
    params['fdsInputFile'] = os.path.abspath(configDir+os.sep+params['fdsInputFile'])
    params['dataDir'] = os.path.abspath(configDir+os.sep+params['dataDir'])
    params['orientations'] = [-3, -2, -1, 1, 2, 3] if not params['orientations'] else params['orientations']
    params['orientations'] = [-3, -2, -1, 1, 2, 3] if (0 in params['orientations']) else params['orientations']
    
    tStart,tEnd,tInt,tBand = extractTimeParams(params)
    smvGrids, smvObsts, smvBndfs, smvSurfs = smv.parseSMVFile(params['smvfile'])
    fdsChid, fdsDevcs, fdsObsts, fdsVents, fdsSurfs, fdsRamps, fdsCtrls = fds.parseFDSFile(params['fdsInputFile'])
    
    polygonNames, params = getPolygonNames(params, fdsObsts)
    points = parseFDSforPts(fdsObsts, smvObsts, polygonNames)
    
    polygons, numberOfGroups = ut.pts2polygons(points)
    
    # Generate color scheme for each polygon
    pcs = []
    for i in range(0,numberOfGroups): pcs.append(pltc.rgb2hex(np.random.rand(3)))
    
    # Get boundary file names
    bndfs = [[os.sep.join([params['dataDir'],x[1]]), x[0]] for x in smvBndfs if (params['variableName'] in x[2])]
    
    # Read boundary data
    times, mPts, orients = loadBNDFdata(params, bndfs, smvGrids, smvObsts)
    
    # Output data
    namespace = params['fdsInputFile'].replace('.fds','')
    ut.maxValueCSV(times, mPts, polygonNames, "%s_%s_max"%(namespace, params['variableName']))
    
    # Plot data
    if not params['plotStyle']:
        params['plotStyle'] = defaultdict(bool)
        params['plotStyle']['fontSize'] = 16
        params['plotStyle']['figureSize'] = [12, 12]
        params['plotStyle']['lineWidth'] = 3
    plt.figure(figsize=(params['plotStyle']['figureSize'][0],
                        params['plotStyle']['figureSize'][1]))
    for i in range(0, mPts.shape[1]):
        c = pltc.to_rgba(pcs[i])
        #c = (int(255*c[0]),int(255*c[1]),int(255*c[2]))
        plt.plot(times, mPts[:,i], c=c, label=polygonNames[i], linewidth=params['plotStyle']['lineWidth'])
    plt.legend(fontsize=params['plotStyle']['fontSize'])
    plt.xlabel('Time (s)', fontsize=params['plotStyle']['fontSize'])
    plt.ylabel(params['variableName'], fontsize=params['plotStyle']['fontSize'])
    plt.tick_params(labelsize=params['plotStyle']['fontSize'])
    plt.tight_layout()
        
