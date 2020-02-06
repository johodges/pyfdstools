# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 08:17:40 2019

@author: JHodges
"""

import numpy as np
import pandas as pd
from collections import defaultdict
import datetime
import re
import scipy.spatial as scsp
import os
import zipfile

def load_devc(file):
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data

def load_hrr(file):
    with open(file, 'r') as f:
        line = f.readline()
        line = f.readline()
    header = line.split(',')
    data = pd.read_csv(file, delimiter=',', names=header, skiprows=2)
    return data

def dictMerge(a, b, path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                dictMerge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass
            else:
                a[key] = b[key]
        else:
            a[key] = b[key]
    return a

class fdsFileOperations(object):
    def __init__(self):
        self.head = defaultdict(bool)
        self.devcs = defaultdict(bool)
        self.inits = defaultdict(bool)
        self.obsts = defaultdict(bool)
        self.holes = defaultdict(bool)
        self.vents = defaultdict(bool)
        self.surfs = defaultdict(bool)
        self.ramps = defaultdict(bool)
        self.ctrls = defaultdict(bool)
        self.meshes = defaultdict(bool)
        self.slcfs = defaultdict(bool)
        self.bndfs = defaultdict(bool)
        self.time = defaultdict(bool)
        self.dump = defaultdict(bool)
        self.misc = defaultdict(bool)
        self.zones = defaultdict(bool)
        self.reacs = defaultdict(bool)
        self.matls = defaultdict(bool)
        self.radis = defaultdict(bool)
        self.pres = defaultdict(bool)
        
        self.devcUnknownCounter = 0
        self.obstUnknownCounter = 0
        self.ventUnknownCounter = 0
        self.meshUnknownCounter = 0
        self.holeUnknownCounter= 0
        
        self.slcfCounter = 0
        self.bndfCounter = 0
        
        self.meshOrder = False
    
    def zopen(self, file):
        if '.zip' in file:
            zname = '%s.zip'%(file.split('.zip')[0])
            fname = file.split('.zip%s'%(os.sep))[1]
            zip = zipfile.ZipFile(zname, 'r')
            f = zip.open(fname)
        else:
            f = open(file, 'rb')
        return f
    
    def importFile(self, file=None, text=None, textList=None):
        if file != None:
            f = self.zopen(file)
            textFDS = f.read()
            textFDS = textFDS.decode("utf-8")
        elif text != None:
            textFDS = text
        elif textList != None:
            textFDS = '\n'.join(textList)
        lines = self.makeFDSLines(textFDS)
        self.parseFDSLines(lines)
    
    def makeFDSLines(self, textFDS):
        linesFDS = [x.split('/')[0] for x in textFDS.split("&")[1:]]
        for i in range(0, len(linesFDS)):
            line2 = linesFDS[i]
            line2 = "%s,"%(line2) if line2[-1] != ',' else line2
            linesFDS[i] = line2
        return linesFDS
    
    def addOBST(self, ID, XB, SURF_IDS=None, SURF_ID=None, SURF_ID6=None,
                BNDF_OBST=True, THICKEN=None, TRANSPARENCY=None,
                COLOR=None):
        obst = defaultdict(bool)
        obst['XB'] = XB
        obst['ID'] = ID
        if SURF_IDS != None: obst['SURF_IDS'] = SURF_IDS
        if SURF_ID != None: obst['SURF_ID'] = SURF_ID
        if SURF_ID6 != None: obst['SURF_ID6'] = SURF_ID6
        if BNDF_OBST: obst['BNDF_OBST'] = True
        if THICKEN != None: obst['THICKEN'] = THICKEN
        if TRANSPARENCY != None: obst['TRANSPARENCY'] = TRANSPARENCY
        if COLOR != None: obst['COLOR'] = COLOR
        if self.obsts[ID]:
            counter = self.obsts[ID]['counter']
            counter += 1
            self.obsts["%s-%0.0f"%(ID, counter)] = obst
            self.obsts[ID]['counter'] = counter
        else:
            obst['counter'] = 0
            self.obsts[ID] = obst
    
    def addHOLE(self, ID, XB):
        hole = defaultdict(bool)
        hole['XB'] = XB
        hole['ID'] = ID
        self.holes[ID] = hole
    
    def addHEAD(self, chid, title=None):
        head = defaultdict(bool)
        head['CHID'] = chid
        if title != None:
            head['TITLE'] = title
        else:
            head['TITLE'] = chid
        self.head['ID'] = head
    
    def addTIME(self, T_END=0.0, T_BEGIN=0.0):
        time = defaultdict(bool)
        time['T_BEGIN'] = T_BEGIN
        time['T_END'] = T_END
        self.time['ID'] = time
    
    def addMISC(self, BNDF_DEFAULT=None, TMPA=None):
        misc = defaultdict(bool)
        if BNDF_DEFAULT != None: misc['BNDF_DEFAULT'] = BNDF_DEFAULT
        if TMPA != None: misc['TMPA'] = TMPA
        self.misc['ID'] = misc
    
    def addDUMP(self, RENDER_FILE=None, COLUMN_DUMP_LIMIT=False,
                WRITE_XYZ=False, DT_PL3D=None, DT_SL3D=None, DT_SLCF=None,
                DT_BNDF=None, DT_DEVC=None, DT_CTRL=None, DT_HRR=None,
                DT_RESTART=None):
        dump = defaultdict(bool)
        if RENDER_FILE != None: dump['RENDER_FILE'] = RENDER_FILE
        if COLUMN_DUMP_LIMIT: dump['COLUMN_DUMP_LIMIT'] = COLUMN_DUMP_LIMIT
        if WRITE_XYZ: dump['WRITE_XYZ'] = WRITE_XYZ
        if DT_PL3D != None: dump['DT_PL3D'] = DT_PL3D
        if DT_SL3D != None: dump['DT_SL3D'] = DT_SL3D
        if DT_SLCF != None: dump['DT_SLCF'] = DT_SLCF
        if DT_BNDF != None: dump['DT_BNDF'] = DT_BNDF
        if DT_DEVC != None: dump['DT_DEVC'] = DT_DEVC
        if DT_CTRL != None: dump['DT_CTRL'] = DT_CTRL
        if DT_HRR != None: dump['DT_HRR'] = DT_HRR
        if DT_RESTART != None: dump['DT_RESTART'] = DT_RESTART
        self.dump['ID'] = dump
    
    def addCTRL(self, ID, FUNCTION_TYPE, INPUT_ID, DELAY=None):
        ctrl = defaultdict(bool)
        ctrl['ID'] = ID
        ctrl['FUNCTION_TYPE'] = FUNCTION_TYPE
        ctrl['INPUT_ID'] = INPUT_ID
        if DELAY != None: ctrl['DELAY'] = DELAY
        self.ctrls[ID] = ctrl
    
    def addPRES(self, VELOCITY_TOLERANCE=None, MAX_PRESSURE_ITERATIONS=None):
        pres = defaultdict(bool)
        if VELOCITY_TOLERANCE != None: pres['VELOCITY_TOLERANCE'] = VELOCITY_TOLERANCE
        if MAX_PRESSURE_ITERATIONS != None: pres['MAX_PRESSURE_ITERATIONS'] = MAX_PRESSURE_ITERATIONS
        self.pres['ID'] = pres
    
    def addDEVC(self, ID, QUANTITY, XYZ=None, XB=None, IOR=None, SPEC_ID=None,
                TIME_AVERAGED=None, SPATIAL_STATISTIC=None,
                INITIAL_STATE=None, INIT_ID=None, SETPOINT=None,
                DUCT_ID=None):
        devc = defaultdict(bool)
        devc['ID'] = ID
        devc['QUANTITY'] = QUANTITY
        if XYZ != None: devc['XYZ'] = XYZ
        if XB != None: devc['XB'] = XB
        if INITIAL_STATE != None: devc['INITIAL_STATE'] = INITIAL_STATE
        if INIT_ID != None: devc['INIT_ID'] = INIT_ID
        if SETPOINT != None: devc['SETPOINT'] = SETPOINT
        if IOR != None: devc['IOR'] = IOR
        if TIME_AVERAGED != None: devc['TIME_AVERAGED'] = TIME_AVERAGED
        if SPATIAL_STATISTIC != None: devc['SPATIAL_STATISTIC'] = SPATIAL_STATISTIC
        if DUCT_ID != None: devc['DUCT_ID'] = DUCT_ID
        if SPEC_ID != None: devc['SPEC_ID'] = SPEC_ID
        self.devcs[ID] = devc
    
    def addMESHfromCOMP(self, ID, x0, y0, z0, x1, y1, z1, doors,
                        extSURF_IDS=['OPEN','OPEN','OPEN','OPEN','ADIABATIC','ADIABATIC'],
                        wallExtension=defaultdict(bool), dx=0.1, dy=0.1, dz=0.1):
        mesh = defaultdict(bool)
        IJK = [(x1-x0)/dx, (y1-y0)/dy, (z1-z0)/dz]
        XB = [x0, x1, y0, y1, z0, z1]
        mesh['ID'] = ID
        mesh['IJK'] = IJK
        mesh['XB'] = XB
        self.meshes[ID] = mesh
        
        for key in list(wallExtension.keys()):
            (xmn, xmx, ymn, ymx, zmn, zmx) = (x0, x1, y0, y1, z0, z1)
            if (wallExtension[key] > 0) and (doors[key]['vtype'] != 'MECHANICAL'):
                mesh = defaultdict(bool)
                meshID = "%s-%s"%(ID, key)
                mesh['ID'] = meshID
                
                door = doors[key]
                xs = door['xs']
                ys = door['ys']
                
                if (key == 'west'):
                    xmx = xmn
                    xmn = xmn-wallExtension[key]
                    ymn = max(ys[0]-3*dx, y0)
                    ymx = min(ys[1]+3*dx, y1)
                    self.addVENT("%s-SOUTH"%(meshID), extSURF_IDS[2], [xmn, xmx, ymn, ymn, zmn, zmx])
                    self.addVENT("%s-NORTH"%(meshID), extSURF_IDS[3], [xmn, xmx, ymx, ymx, zmn, zmx])
                    self.addVENT("%s-WEST"%(meshID), extSURF_IDS[0], [xmn, xmn, ymn, ymx, zmn, zmx])
                elif (key == 'east'):
                    xmn = xmx
                    xmx = xmx+wallExtension[key]
                    ymn = max(ys[0]-3*dx, y0)
                    ymx = min(ys[1]+3*dx, y1)
                    self.addVENT("%s-SOUTH"%(meshID), extSURF_IDS[2], [xmn, xmx, ymn, ymn, zmn, zmx])
                    self.addVENT("%s-NORTH"%(meshID), extSURF_IDS[3], [xmn, xmx, ymx, ymx, zmn, zmx])
                    self.addVENT("%s-EAST"%(meshID), extSURF_IDS[0], [xmx, xmx, ymn, ymx, zmn, zmx])
                elif (key == 'north'):
                    ymn = ymx
                    ymx = ymx+wallExtension[key]
                    xmn = max(xs[0]-3*dx, x0)
                    xmx = min(xs[1]+3*dx, x1)
                    self.addVENT("%s-EAST"%(meshID), extSURF_IDS[0], [xmx, xmx, ymn, ymx, zmn, zmx])
                    self.addVENT("%s-WEST"%(meshID), extSURF_IDS[0], [xmn, xmn, ymn, ymx, zmn, zmx])
                    self.addVENT("%s-NORTH"%(meshID), extSURF_IDS[3], [xmn, xmx, ymx, ymx, zmn, zmx])
                elif (key == 'south'):
                    ymx = ymn
                    ymn = ymn-wallExtension[key]
                    xmn = max(xs[0]-3*dx, x0)
                    xmx = min(xs[1]+3*dx, x1)
                    self.addVENT("%s-EAST"%(meshID), extSURF_IDS[0], [xmx, xmx, ymn, ymx, zmn, zmx])
                    self.addVENT("%s-WEST"%(meshID), extSURF_IDS[0], [xmn, xmn, ymn, ymx, zmn, zmx])
                    self.addVENT("%s-SOUTH"%(meshID), extSURF_IDS[2], [xmn, xmx, ymn, ymn, zmn, zmx])
                elif (key == 'floor'):
                    zmx = zmn
                    zmn = zmn-wallExtension[key]
                elif (key == 'ceiling'):
                    zmn = zmx
                    zmx = zmx+wallExtension[key]
                IJK = [(xmx-xmn)/dx, (ymx-ymn)/dy, (zmx-zmn)/dz]
                XB = [xmn, xmx, ymn, ymx, zmn, zmx]
                mesh['IJK'] = IJK
                mesh['XB'] = XB
                self.meshes[meshID] = mesh
                self.addVENT("%s-FLOOR"%(meshID), extSURF_IDS[4], [xmn, xmx, ymn, ymx, zmn, zmn])
                self.addVENT("%s-CEILING"%(meshID), extSURF_IDS[4], [xmn, xmx, ymn, ymx, zmx, zmx])
    
    def addOBSTfromCOMP(self, ID, x0, y0, z0, x1, y1, z1, doors, dx=0.1, dy=0.1, dz=0.1, SURF_ID=None, SURF_IDS=None, SURF_ID6=None, BNDF_OBST=True):
        self.addOBST('%s-WEST'%(ID), [x0, x0+dx, y0, y1, z0, z1], SURF_IDS=SURF_IDS, SURF_ID=SURF_ID, SURF_ID6=SURF_ID6, BNDF_OBST=BNDF_OBST)
        self.addOBST('%s-EAST'%(ID), [x1-dx, x1, y0, y1, z0, z1], SURF_IDS=SURF_IDS, SURF_ID=SURF_ID, SURF_ID6=SURF_ID6, BNDF_OBST=BNDF_OBST)
        self.addOBST('%s-NORTH'%(ID), [x0, x1, y1-dy, y1, z0, z1], SURF_IDS=SURF_IDS, SURF_ID=SURF_ID, SURF_ID6=SURF_ID6, BNDF_OBST=BNDF_OBST)
        self.addOBST('%s-SOUTH'%(ID), [x0, x1, y0, y0+dx, z0, z1], SURF_IDS=SURF_IDS, SURF_ID=SURF_ID, SURF_ID6=SURF_ID6, BNDF_OBST=BNDF_OBST)
        
        self.addVENT('%s-FLOOR'%(ID), SURF_ID, XB=[x0, x1, y0, y1, z0, z0])
        self.addVENT('%s-CEILING'%(ID), SURF_ID, XB=[x0, x1, y0, y1, z1, z1])
        #self.addOBST('%s-FLOOR'%(ID), [x0, x1, y0, y1, z0-dz, z0], SURF_IDS=SURF_IDS, SURF_ID=SURF_ID, SURF_ID6=SURF_ID6, BNDF_OBST=BNDF_OBST)
        #self.addOBST('%s-CEILING'%(ID), [x0, x1, y0, y1, z1, z1+dz], SURF_IDS=SURF_IDS, SURF_ID=SURF_ID, SURF_ID6=SURF_ID6, BNDF_OBST=BNDF_OBST)
        
        for key in list(doors.keys()):
            door = doors[key]
            if door['exist'] == 1:
                xs = door['xs']
                ys = door['ys']
                doorX = abs(xs[1]-xs[0])
                doorY = abs(ys[1]-ys[0])
                compX = abs(x1-x0)
                compY = abs(y1-y0)
                zs = [0, 2.5]
                if door['vtype'] != 'MECHANICAL':
                    if doorX < dx: xs = [xs[0]-2*dx, xs[1]+2*dx]
                    if doorY < dy: ys = [ys[0]-2*dy, ys[1]+2*dy]
                    if abs(doorX-compX) < dx: xs = [xs[0]+dx, xs[1]-dx]
                    elif (xs[0] == x0): xs = [xs[0] + dx, xs[1]]
                    elif (xs[1] == x1): xs = [xs[0], xs[1] - dx]
                    if abs(doorY-compY) < dy: ys = [ys[0]+dy, ys[1]-dy]
                    elif (ys[0] == y0): ys = [ys[0] + dy, ys[1]]
                    elif (ys[1] == y1): ys = [ys[0], ys[1] - dy]
                    XB = [xs[0], xs[1], ys[0], ys[1], zs[0], zs[1]]
                    self.addHOLE('%s-DOOR-%s'%(ID, key),XB)
                else:
                    if (key == 'south'): XB = [xs[0], xs[1], ys[0]+dx, ys[1]+dx, zs[0], zs[1]]
                    if (key == 'north'): XB = [xs[0], xs[1], ys[0]-dx, ys[1]-dx, zs[0], zs[1]]
                    if (key == 'west'): XB = [xs[0]+dx, xs[1]+dx, ys[0], ys[1], zs[0], zs[1]]
                    if (key == 'east'): XB = [xs[0]-dx, xs[1]-dx, ys[0], ys[1], zs[0], zs[1]]
                    self.addVENT('%s-DOOR-%s'%(ID, key), 'MECH-VENT', XB)
    
    def checkOverlappingMESH(self):
        def in_hull(p,hull):
            if not isinstance(hull,scsp.Delaunay):
                hull = scsp.Delaunay(hull)
            return hull.find_simplex(p)>=0
        def pointsFromXB(XB,extend=[0.05, -0.05, 0.05, -0.05, 0, 0]):
            pts = [[XB[0]+extend[0],XB[2]+extend[2],XB[4]+extend[4]],
                   [XB[0]+extend[0],XB[2]+extend[2],XB[5]+extend[5]],
                   [XB[0]+extend[0],XB[3]+extend[3],XB[4]+extend[4]],
                   [XB[0]+extend[0],XB[3]+extend[3],XB[5]+extend[5]],
                   [XB[1]+extend[1],XB[2]+extend[2],XB[4]+extend[4]],
                   [XB[1]+extend[1],XB[2]+extend[2],XB[5]+extend[5]],
                   [XB[1]+extend[1],XB[3]+extend[3],XB[4]+extend[4]],
                   [XB[1]+extend[1],XB[3]+extend[3],XB[5]+extend[5]]]
            return pts
        meshHulls = defaultdict(bool)
        for key in list(self.meshes.keys()):
            pts = pointsFromXB(self.meshes[key]['XB'])
            meshHull = scsp.Delaunay(pts)
            meshHulls[key] = meshHull
        overlap = False
        for key1 in list(self.meshes.keys()):
            for key2 in list(self.meshes.keys()):
                if (key1 != key2):
                    extend = [0.05, -0.05, 0.05, -0.05, 0, 0]
                    if ('east' in key2): extend = [0.05, 0.1, 0.05, -0.05, 0, 0]
                    if ('west' in key2): extend = [-0.1, -0.05, 0.05, -0.05, 0, 0]
                    if ('north' in key2): extend = [0.05, -0.05, 0.05, 0.1, 0, 0]
                    if ('south' in key2): extend = [0.05, -0.05, -0.1, -0.05, 0, 0]
                    pts = pointsFromXB(self.meshes[key2]['XB'], extend=extend)
                    for p in pts:
                        if in_hull(p, meshHulls[key1]):
                            overlap = True
        return overlap
    
    def calculateMeshCells(self):
        meshes = []
        numCells = []
        for key in list(self.meshes.keys()):
            IJK = self.meshes[key]['IJK']
            numCells.append(IJK[0]*IJK[1]*IJK[2])
            meshes.append(key)
        return meshes, numCells
    
    def addMPIprocesses(self, numberOfProcesses):
        meshes, numCells = self.calculateMeshCells()
        cellsPerProcess = np.sum(numCells)/numberOfProcesses
        mpiConverged = False
        splitConverged = False
        splitMultiplier = 1.20
        while not mpiConverged:
            mpiConverged = True
            while not splitConverged:
                splitConverged = True
                meshes, numCells = self.calculateMeshCells()
                for mesh, numCell in zip(meshes, numCells):
                    if numCell > cellsPerProcess*splitMultiplier:
                        self.splitMESHonce(self.meshes[mesh])
                        splitConverged = False
            
            meshes, numCells = self.calculateMeshCells()
            mpiProcessInds = np.zeros((len(numCells),))-1
            mpiProcess = np.zeros((numberOfProcesses,))
            while np.argwhere(mpiProcessInds == -1).shape[0] > 0:
                ind = np.argmax(numCells)
                ind2 = np.argmin(mpiProcess)
                mpiProcessInds[ind] = ind2
                mpiProcess[ind2] += numCells[ind]
                numCells[ind] = 0
            if np.max(mpiProcess) > cellsPerProcess*1.20:
                mpiConverged = False
                splitConverged = False
                splitMultiplier = splitMultiplier*0.9
        for key, mp in zip(meshes, mpiProcessInds):
            self.meshes[key]['MPI_PROCESS'] = mp
        self.mpiProcesses = numberOfProcesses
        self.meshOrder = np.argsort(mpiProcessInds)
        
    def splitMESHonce(self, mesh):
        IJK = np.round(mesh['IJK'])
        XB = mesh['XB']
        dxs = [(XB[1]-XB[0])/float(IJK[0]), (XB[3]-XB[2])/float(IJK[1]), (XB[5]-XB[4])/float(IJK[2])]
        ind = np.argmax(IJK)
        if ind == 2:
            IJK_temp = list(IJK)
            IJK_temp[2] = 0
            ind = np.argmax(IJK_temp)
        
        IJK2 = list(IJK)
        XB2 = list(XB)
        IJK2[ind] = int(IJK[ind]/2)
        if IJK2[ind] % 2 > 0: IJK2[ind] = IJK2[ind]-1
        XB2[int(2*ind+1)] = XB2[int(2*ind)] + dxs[ind]*float(IJK2[ind])
        
        IJK3 = list(IJK)
        XB3 = list(XB)
        IJK3[ind] = IJK[ind] - IJK2[ind]
        XB3[int(2*ind)] = XB2[int(2*ind+1)]
        
        mesh2 = defaultdict(bool)
        mesh2['ID'] = "%s-A"%(mesh["ID"])
        mesh2['XB'] = XB2
        mesh2['IJK'] = IJK2
        
        mesh3 = defaultdict(bool)
        mesh3['ID'] = "%s-B"%(mesh["ID"])
        mesh3['XB'] = XB3
        mesh3['IJK'] = IJK3
        
        self.meshes.pop(mesh['ID'], False)
        self.meshes[mesh2['ID']] = mesh2
        self.meshes[mesh3['ID']] = mesh3
    
    def addMESH(self, ID, IJK, XB):
        mesh = defaultdict(bool)
        mesh['ID'] = ID
        mesh['IJK'] = IJK
        mesh['XB'] = XB
        self.meshes[ID] = mesh
        
    def addREAC(self, ID, FUEL=None, FORMULA=None, AIT=None, SY=None, COY=None, HOC=None,
                C=None, H=None, O=None, N=None, FYI=None, RF=None):
        reac = defaultdict(bool)
        reac['ID'] = ID
        if FUEL != None: reac['FUEL'] = FUEL
        if FORMULA != None: reac['FORMULA'] = FORMULA
        if AIT != None: reac['AUTO_IGNITION_TEMPERATURE'] = AIT
        if SY != None: reac['SOOT_YIELD'] = SY
        if COY != None: reac['CO_YIELD'] = COY
        if HOC != None: reac['HEAT_OF_COMBUSTION'] = HOC
        if C != None: reac['C'] = C
        if H != None: reac['H'] = H
        if O != None: reac['O'] = O
        if N != None: reac['N'] = N
        if FYI != None: reac['FYI'] = FYI
        if RF != None: reac['RADIATIVE_FRACTION'] = RF
        self.reacs[ID] = reac       
        
    def addMATL(self, ID, Emi=None, Den=None, Con=None, Spe=None, kramp=None, cpramp=None, fyi=None):
        matl = defaultdict(bool)
        matl['ID'] = ID
        if Emi != None: matl['EMISSIVITY'] = Emi
        if Den != None: matl['DENSITY'] = Den
        if Con != None: matl['CONDUCTIVITY'] = Con
        if Spe != None: matl['SPECIFIC_HEAT'] = Spe
        if kramp != None: matl['CONDUCTIVITY_RAMP'] = kramp
        if cpramp != None: matl['SPECIFIC_HEAT_RAMP'] = cpramp
        if fyi != None: matl['FYI'] = fyi
        self.matls[ID] = matl
    
    def addRAMP(self, ID, T, F):
        if self.ramps[ID]:
            Ts = self.ramps[ID]['T']
            Fs = self.ramps[ID]['F']
            for t, f in zip(T, F):
                Ts.append(t)
                Fs.append(f)
            self.ramps[ID]['T'] = Ts
            self.ramps[ID]['F'] = Fs
        else:
            self.ramps[ID] = defaultdict(bool)
            self.ramps[ID]['T'] = T
            self.ramps[ID]['F'] = F
            self.ramps[ID]['ID'] = ID
    
    def addSURF(self, ID, Mid=None, Col=None, Thi=None, Bac=None, Geo=None,
                Fyi=None, Len=None, LeaPat=None, Hrrpua=None, qramp=None,
                Rgb=None, adiabatic=False, VOLUME_FLOW=None,
                VEL_T=None):
        surf = defaultdict(bool)
        surf['ID'] = ID
        if Mid != None: surf['MATL_ID'] = Mid
        if Col != None: surf['COLOR'] = Col
        if Thi != None: surf['THICKNESS'] = Thi
        if Bac != None: surf['BACKING'] = Bac
        if Geo != None: surf['GEOMETRY'] = Geo
        if Fyi != None: surf['FYI'] = Fyi
        if Len != None: surf['LENGTH'] = Len
        if LeaPat != None: surf['LEAK_PATH'] = LeaPat
        if Hrrpua != None: surf['HRRPUA'] = Hrrpua
        if qramp != None: surf['RAMP_Q'] = qramp
        if Rgb != None: surf['RGB'] = Rgb
        if adiabatic: surf['ADIABATIC'] = True
        if VOLUME_FLOW != None: surf['VOLUME_FLOW'] = VOLUME_FLOW
        if VEL_T != None: surf['VEL_T'] = VEL_T
        self.surfs[ID] = surf
        
    def addVENT(self, ID, SURF_ID, XB=None, CTRL_ID=None, MB=None, IOR=None):
        vent = defaultdict(bool)
        vent['ID'] = ID
        vent['SURF_ID'] = SURF_ID
        if XB != None: vent['XB'] = XB
        if CTRL_ID != None: vent['CTRL_ID'] = CTRL_ID
        if MB != None: vent['MB'] = MB
        if IOR != None: vent['IOR'] = IOR
        if self.vents[ID]:
            counter = self.vents[ID]['counter']
            counter += 1
            self.vents["%s-%0.0f"%(ID, counter)] = vent
            self.vents[ID]['counter'] = counter
        else:
            vent['counter'] = 0
            self.vents[ID] = vent
    
    def addZONE(self, ID, XB, LEAK_AREA=None):
        zone = defaultdict(bool)
        zone['ID'] = ID
        zone['XB'] = XB
        if LEAK_AREA != None: zone['LEAK_AREA'] = LEAK_AREA
        self.zones[ID] = zone
    
    def addSLCF(self, Qty, PBX=None, PBY=None, PBZ=None, Vec=False, XB=None, SPEC_ID=None):
        slcf = defaultdict(bool)
        slcf['ID'] = "SLCF-%05.0f"%(self.slcfCounter)
        slcf['QUANTITY'] = Qty
        if PBX != None: slcf['PBX'] = PBX
        if PBY != None: slcf['PBY'] = PBY
        if PBZ != None: slcf['PBZ'] = PBZ
        if SPEC_ID != None: slcf['SPEC_ID'] = SPEC_ID
        if Vec: slcf['VECTOR'] = 'TRUE'
        if XB != None: slcf['XB'] = XB
        self.slcfCounter += 1
        self.slcfs[slcf['ID']] = slcf
    
    def addBNDF(self, Qty, CELL_CENTERED=None):
        bndf = defaultdict(bool)
        bndf['ID'] = "BNDF-%05.0f"%(self.bndfCounter)
        bndf['QUANTITY'] = Qty
        if CELL_CENTERED != None: bndf['CELL_CENTERED'] = CELL_CENTERED
        self.bndfCounter += 1
        self.bndfs[bndf['ID']] = bndf
    
    def getLineType(self, line):
        lineType = line[:4]
        return lineType
    
    def parseHEAD(self, line):
        head = parseHEAD(line)
        if self.head['ID']:
            self.head['ID'] = dictMerge(self.head['ID'], head)
        else:
            self.head['ID'] = head
    
    def parseDEVC(self, line):
        devc = parseDEVC(line)
        if self.devcs[devc['ID']]:
            self.devcs['NULL%04.0f'%(self.devcUnknownCounter)] = devc
            self.devcUnknownCounter += 1
        else:
            self.devcs[devc['ID']] = devc
    
    def parseINIT(self, line):
        init = parseINIT(line)
        self.init[init['ID']] = init
    
    def parseOBST(self, line):
        obst = parseOBST(line)
        if not obst['ID']:
            obst['ID'] = 'NULL%04.0f'%(self.obstUnknownCounter)
            self.obsts['NULL%04.0f'%(self.obstUnknownCounter)] = obst
            self.obstUnknownCounter += 1
        else:
            if not self.obsts[obst['ID']]:
                obst['number'] = 0
                self.obsts[obst['ID']] = obst
            else:
                obstNum = self.obsts[obst['ID']]['number']+1
                self.obsts["%s-%04.0f"%(obst['ID'],obstNum)] = obst
                self.obsts[obst['ID']]['number'] = obstNum
                
    def parseHOLE(self, line):
        hole = parseHOLE(line)
        if not hole['ID']:
            hole['ID'] = 'NULL%04.0f'%(self.holeUnknownCounter)
            self.holes['NULL%04.0f'%(self.holeUnknownCounter)] = hole
            self.holeUnknownCounter += 1
        else:
            if not self.holes[hole['ID']]:
                hole['number'] = 0
                self.holes[hole['ID']] = hole
            else:
                holeNum = self.holes[hole['ID']]['number']+1
                self.holes["%s-%04.0f"%(hole['ID'],holeNum)] = hole
                self.holes[hole['ID']]['number'] = holeNum
                
    def parseVENT(self, line):
        vent = parseVENT(line)
        if not vent['ID']:
            vent['ID'] = 'NULL%04.0f'%(self.ventUnknownCounter)
            self.vents['NULL%04.0f'%(self.ventUnknownCounter)] = vent
            self.ventUnknownCounter += 1
        else:
            if not self.vents[vent['ID']]:
                vent['number'] = 0
                self.vents[vent['ID']] = vent
            else:
                ventNum = self.vents[vent['ID']]['number']+1
                self.vents["%s-%04.0f"%(vent['ID'],ventNum)] = vent
                self.vents[vent['ID']]['number'] = ventNum

    def parseSURF(self, line):
        surf = parseSURF(line)
        self.surfs[surf['ID']] = surf

    def parseRAMP(self, line):
        ramp = parseRAMP(line)
        if self.ramps[ramp['ID']]:
            Ts = self.ramps[ramp['ID']]['T']
            Ts.append(ramp['T'])
            Fs = self.ramps[ramp['ID']]['F']
            Fs.append(ramp['F'])
            self.ramps[ramp['ID']]['T'] = Ts
            self.ramps[ramp['ID']]['F'] = Fs
        else:
            ramp['T'] = [ramp['T']]
            ramp['F'] = [ramp['F']]
            self.ramps[ramp['ID']] = ramp

    def parseCTRL(self, line):
        ctrl = parseCTRL(line)
        self.ctrls[ctrl["ID"]] = ctrl
    
    def parseMESH(self, line):
        mesh = parseMESH(line)
        if self.meshes[mesh["ID"]]:
            self.meshes["%s-%04.0f"%(mesh["ID"], self.meshUnknownCounter)] = mesh
            self.meshUnknownCounter += 1
        else:
            self.meshes[mesh["ID"]] = mesh
    
    def parsePRES(self, line):
        pres = parsePRES(line)
        if self.pres['ID']:
            self.pres['ID'] = dictMerge(self.pres['ID'], pres)
        else:
            self.pres['ID'] = pres
    
    def parseSLCF(self, line):
        slcf = parseSLCF(line)
        slcf['ID'] = "SLCF-%05.0f"%(self.slcfCounter)
        self.slcfs[slcf['ID']] = slcf
        self.slcfCounter += 1
    
    def parseBNDF(self, line):
        bndf = parseBNDF(line)
        bndf['ID'] = "BNDF-%05.0f"%(self.bndfCounter)
        self.bndfs[bndf['ID']] = bndf
        self.bndfCounter += 1
    
    def parseTIME(self, line):
        time = parseTIME(line)
        if self.time['ID']:
            self.time['ID'] = dictMerge(self.time['ID'], time)
        else:
            self.time['ID'] = time
        
    def parseDUMP(self, line):
        dump = parseDUMP(line)
        if self.dump['ID']:
            self.dump['ID'] = dictMerge(self.dump["ID"], dump)
        else:
            self.dump["ID"] = dump
        
    def parseMISC(self, line):
        misc = parseMISC(line)
        if self.misc["ID"]:
            self.misc["ID"] = dictMerge(self.misc["ID"], misc)
        else:
            self.misc["ID"] = misc
    
    def parseZONE(self, line):
        zone = parseZONE(line)
        self.zones[zone['ID']] = zone
    
    def parseREAC(self, line):
        reac = parseREAC(line)
        self.reacs[reac['ID']] = reac
    
    def parseMATL(self, line):
        matl = parseMATL(line)
        self.matls[matl['ID']] = matl
    
    def parseRADI(self, line):
        radi = parseRADI(line)
        self.radis[radi['ID']] = radi
    
    def parseFDSLines(self, lines):
        for line in lines:
            lineType = self.getLineType(line)
            if lineType == 'HEAD': self.parseHEAD(line)
            if lineType == 'DEVC': self.parseDEVC(line)
            if lineType == 'INIT': self.parseINIT(line)
            if lineType == 'OBST': self.parseOBST(line)
            if lineType == 'VENT': self.parseVENT(line)
            if lineType == 'SURF': self.parseSURF(line)
            if lineType == 'RAMP': self.parseRAMP(line)
            if lineType == 'CTRL': self.parseCTRL(line)
            if lineType == 'MESH': self.parseMESH(line)
            if lineType == 'SLCF': self.parseSLCF(line)
            if lineType == 'BNDF': self.parseBNDF(line)
            if lineType == 'TIME': self.parseTIME(line)
            if lineType == 'DUMP': self.parseDUMP(line)
            if lineType == 'MISC': self.parseMISC(line)
            if lineType == 'ZONE': self.parseZONE(line)
            if lineType == 'REAC': self.parseREAC(line)
            if lineType == 'MATL': self.parseMATL(line)
            if lineType == 'RADI': self.parseRADI(line)
            if lineType == 'PRES': self.parsePRES(line)
            if lineType == 'HOLE': self.parseHOLE(line)
        for key in list(self.devcs.keys()):
            if self.devcs[key]['INIT_ID']:
                self.devcs[key]['XYZ'] = self.inits[self.devcs[key]['INIT_ID']]['XYZ']
    
    def getMeshLimits(self):
        meshLimits = defaultdict(bool)
        limitingXB = [100000, -100000, 100000, -100000, 100000, -100000]
        for key in list(self.meshes.keys()):
            mesh = self.meshes[key]
            XB = mesh['XB']
            limitingXB[0] = min([limitingXB[0], XB[0]])
            limitingXB[1] = max([limitingXB[1], XB[1]])
            limitingXB[2] = min([limitingXB[2], XB[2]])
            limitingXB[3] = max([limitingXB[3], XB[3]])
            limitingXB[4] = min([limitingXB[4], XB[4]])
            limitingXB[5] = max([limitingXB[5], XB[5]])
            meshLimits[key] = mesh
        meshLimits['Overall'] = defaultdict(bool)
        meshLimits['Overall']['XB'] = limitingXB
        return meshLimits
    
    def generateFDStext(self):
        date = datetime.date.today()
        text = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        text = "%s!!!!! Input file generated with fdsTools v1 on %04.0f-%02.0f-%02.0f          !!!!!\n"%(text, date.year, date.month, date.day)
        text = "%s!!!!! Copyright JENSEN HUGHES %04.0f All right reserved.             !!!!!\n"%(text, date.year)
        text = "%s!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"%(text)
        
        text = "%s%s"%(text, makeLinesFromDict(self.head, getHEADtypes(), '&HEAD', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.time, getTIMEtypes(), '&TIME', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.misc, getMISCtypes(), '&MISC', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.inits, getINITtypes(), '&INIT', newline=True))
        text = "%s%s"%(text, makeLinesFromDict(self.dump, getDUMPtypes(), '&DUMP', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.zones, getZONEtypes(), '&ZONE', newline=True))
        text = "%s%s"%(text, makeLinesFromDict(self.pres, getPREStypes(), '&PRES', newline=True))
        
        if not self.meshOrder:
            self.addMPIprocesses(1)
        text = "%s%s"%(text, makeMESH(self.meshes, meshOrder=self.meshOrder))
        text = "%s%s"%(text, makeLinesFromDict(self.reacs, getREACtypes(), '&REAC', newline=True))
        text = "%s%s"%(text, makeLinesFromDict(self.radis, getRADItypes(), '&RADI', newline=True))
        text = "%s%s"%(text, makeLinesFromDict(self.matls, getMATLtypes(), '&MATL', newline=True))
        text = "%s%s"%(text, makeLinesFromDict(self.surfs, getSURFtypes(), '&SURF', newline=True))
        text = "%s%s"%(text, makeRAMP(self.ramps))
        
        text = "%s%s"%(text, makeLinesFromDict(self.obsts, getOBSTtypes(), '&OBST', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.holes, getHOLEtypes(), '&HOLE', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.vents, getVENTtypes(), '&VENT', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.devcs, getDEVCtypes(), '&DEVC', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.ctrls, getCTRLtypes(), '&CTRL', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.bndfs, getBNDFtypes(), '&BNDF', newline=False))
        text = "%s%s"%(text, makeLinesFromDict(self.slcfs, getSLCFtypes(), '&SLCF', newline=False))
        
        return text

def keyAssist(text, types, dic, internalKeys=['counter'], newline=False):
    keys = list(dic.keys())
    keys.sort()
    if 'ID' in keys: keys.insert(0, keys.pop(keys.index('ID')))
    for key in internalKeys:
        if key in keys:
            keys.remove(key)
    for key2 in keys:
        if (types[key2] == 'ignore'):
            pass
        elif (types[key2] == 'string'):
            text = "%s%s='%s', "%(text, key2, dic[key2])
        elif (types[key2] == 'float'):
            text = "%s%s=%0.4f, "%(text, key2, dic[key2])
        elif (types[key2] == 'int'):
            text = "%s%s=%0.0f, "%(text, key2, dic[key2])
        elif (types[key2] == 'bool'):
            if (dic[key2] is True) or (dic[key2] == 'TRUE') or (dic[key2] == '.TRUE.'):
                text = "%s%s=.TRUE., "%(text, key2)
            else:
                text = "%s%s=.FALSE., "%(text, key2)
        elif ('listind' in types[key2]):
            temp = dic[key2]
            for t in temp:
                tempTxt = "%s(%0.0f)="%(key2, t[0])
                if ('string' in types[key2]): tempTxt = "%s '%s',"%(tempTxt, t[1])
                if ('float' in types[key2]): tempTxt = "%s %0.4f,"%(tempTxt, t[1])
                if ('int' in types[key2]): tempTxt = "%s %0.0f,"%(tempTxt, t[1])
                text = "%s%s "%(text, tempTxt)
        elif ('list' in types[key2]):
            temp = dic[key2]
            tempTxt = "%s="%(key2)
            for t in temp:
                if ('string' in types[key2]): tempTxt = "%s '%s',"%(tempTxt, t)
                if ('float' in types[key2]): tempTxt = "%s %0.4f,"%(tempTxt, t)
                if ('int' in types[key2]): tempTxt = "%s %0.0f,"%(tempTxt, t)
            text = "%s%s "%(text, tempTxt)
        else:
            print(key2, types[key2])
        if newline and (types[key2] != 'ignore'):
            text = "%s\n      "%(text)
    return text

def makeLinesFromDict(items, types, prefix, newline=False):
    text = ''
    keys = list(items.keys())
    keys.sort()
    for key in keys:
        text = "%s%s "%(text, prefix)
        text = keyAssist(text, types, items[key], newline=newline)
        text = "%s /\n"%(text)
    return text

def makeMESH(meshes, meshOrder=False):
    text = ''
    meshList = list(meshes.keys())
    if (meshOrder is not False): meshList = [meshList[x] for x in meshOrder]
    for key in meshList:
        meshTypes = getMESHtypes()
        text = "%s&MESH "%(text)
        text = keyAssist(text, meshTypes, meshes[key])
        text = "%s /\n"%(text)
    return text

def makeRAMP(ramps):
    text = ''
    for key in list(ramps.keys()):
        ID = ramps[key]['ID']
        for F, T in zip(ramps[key]['F'], ramps[key]['T']):
            text = "%s&RAMP ID='%s', T = %0.4f, F = %0.4f, /\n"%(text, ID, T, F)
    return text

def getMATLtypes():
    matlTypes = defaultdict(bool)
    matlTypes['ID'] = 'string'
    matlTypes['FYI'] = 'string'
    matlTypes['SPECIFIC_HEAT'] = 'float'
    matlTypes['CONDUCTIVITY'] = 'float'
    matlTypes['DENSITY'] = 'float'
    matlTypes['EMISSIVITY'] = 'float'
    matlTypes['SPECIFIC_HEAT_RAMP'] = 'string'
    matlTypes['CONDUCTIVITY_RAMP'] = 'string'
    return matlTypes

def getREACtypes():
    reacTypes = defaultdict(bool)
    reacTypes['FUEL'] = 'string'
    reacTypes['ID'] = 'string'
    reacTypes['FORMULA'] = 'string'
    reacTypes['AUTO_IGNITION_TEMPERATURE'] = 'float'
    reacTypes['SOOT_YIELD'] = 'float'
    reacTypes['CO_YIELD'] = 'float'
    reacTypes['HEAT_OF_COMBUSTION'] = 'float'
    reacTypes['C'] = 'float'
    reacTypes['H'] = 'float'
    reacTypes['O'] = 'float'
    reacTypes['N'] = 'float'
    reacTypes['FYI'] = 'string'
    reacTypes['RADIATIVE_FRACTION'] = 'float'
    reacTypes['IDEAL'] = 'bool'
    return reacTypes

def getRADItypes():
    radiTypes = defaultdict(bool)
    radiTypes['RADIATIVE_FRACTION'] = 'float'
    return radiTypes

def getMESHtypes():
    meshTypes = defaultdict(bool)
    meshTypes['ID'] = 'string'
    meshTypes['IJK'] = 'listint'
    meshTypes['XB'] = 'listfloat'
    meshTypes['MPI_PROCESS'] = 'int'
    return meshTypes

def getSURFtypes():
    surfTypes = defaultdict(bool)
    surfTypes['ID'] = 'string'
    surfTypes['MATL_ID'] = 'liststring'
    surfTypes['THICKNESS'] = 'listfloat'
    surfTypes['RGB'] = 'listfloat'
    surfTypes['COLOR'] = 'string'
    surfTypes['BACKING'] = 'string'
    surfTypes['GEOMETRY'] = 'string'
    surfTypes['FYI'] = 'string'
    surfTypes['LENGTH'] = 'float'
    surfTypes['C_FORCED_RE'] = 'float'
    surfTypes['C_FORCED_CONSTANT'] = 'float'
    surfTypes['C_FORCED_RE_EXP'] = 'float'
    surfTypes['C_FORCED_PR_EXP'] = 'float'
    surfTypes['CONVECTION_LENGTH_SCALE'] = 'float'
    surfTypes['C_HORIZONTAL'] = 'float'
    surfTypes['C_VERTICAL'] = 'float'
    surfTypes['EMISSIVITY'] = 'float'
    surfTypes['LEAK_PATH'] = 'listint'
    surfTypes['HRRPUA'] = 'float'
    surfTypes['RAMP_Q'] = 'string'
    surfTypes['ADIABATIC'] = 'bool'
    surfTypes['DEFAULT'] = 'bool'
    surfTypes['VOLUME_FLOW'] = 'float'
    surfTypes['VEL_T'] = 'listfloat'
    return surfTypes

def getVENTtypes():
    surfTypes = defaultdict(bool)
    surfTypes['ID'] = 'string'
    surfTypes['XB'] = 'listfloat'
    surfTypes['MB'] = 'string'
    surfTypes['CTRL_ID'] = 'string'
    surfTypes['SURF_ID'] = 'string'
    surfTypes['IOR'] = 'int'
    surfTypes['number'] = 'ignore'
    return surfTypes

def getPREStypes():
    presTypes = defaultdict(bool)
    presTypes['SOLVER'] = 'string'
    presTypes['VELOCITY_TOLERANCE'] = 'float'
    presTypes['MAX_PRESSURE_ITERATIONS'] = 'int'
    return presTypes

def getOBSTtypes():
    obstTypes = defaultdict(bool)
    obstTypes['ID'] = 'string'
    obstTypes['XB'] = 'listfloat'
    obstTypes['CTRL_ID'] = 'string'
    obstTypes['DEVC_ID'] = 'string'
    obstTypes['SURF_ID'] = 'string'
    obstTypes['SURF_ID6'] = 'liststring'
    obstTypes['SURF_IDS'] = 'liststring'
    obstTypes['BNDF_OBST'] = 'bool'
    obstTypes['PERMIT_HOLE'] = 'bool'
    obstTypes['THICKEN'] = 'bool'
    obstTypes['TRANSPARENCY'] = 'float'
    obstTypes['COLOR'] = 'string'
    obstTypes['number'] = 'ignore'
    return obstTypes

def getHOLEtypes():
    holeTypes = defaultdict(bool)
    holeTypes['ID'] = 'string'
    holeTypes['XB'] = 'listfloat'
    holeTypes['CTRL_ID'] = 'string'
    holeTypes['number'] = 'ignore'
    return holeTypes

def getINITtypes():
    initTypes = defaultdict(bool)
    initTypes['ID'] = 'string'
    initTypes['XB'] = 'listfloat'
    initTypes['XYZ'] = 'listfloat'
    initTypes['PART_ID'] = 'string'
    initTypes['N_PARTICLES'] = 'int'
    return initTypes

def getDEVCtypes():
    devcTypes = defaultdict(bool)
    devcTypes['ID'] = 'string'
    devcTypes['XB'] = 'listfloat'
    devcTypes['XYZ'] = 'listfloat'
    devcTypes['QUANTITY'] = 'string'
    devcTypes['INITIAL_STATE'] = 'bool'
    devcTypes['INIT_ID'] = 'string'
    devcTypes['SETPOINT'] = 'float'
    devcTypes['IOR'] = 'int'
    devcTypes['TIME_AVERAGED'] = 'bool'
    devcTypes['SPATIAL_STATISTIC'] = 'string'
    devcTypes['DUCT_ID'] = 'string'
    devcTypes['SPEC_ID'] = 'string'
    return devcTypes

def getBNDFtypes():
    bndfTypes = defaultdict(bool)
    bndfTypes['QUANTITY'] = 'string'
    bndfTypes['CELL_CENTERED'] = 'bool'
    bndfTypes['ID'] = 'ignore'
    return bndfTypes

def getSLCFtypes():
    slcfTypes = defaultdict(bool)
    slcfTypes['QUANTITY'] = 'string'
    slcfTypes['PBX'] = 'float'
    slcfTypes['PBY'] = 'float'
    slcfTypes['PBZ'] = 'float'
    slcfTypes['XB'] = 'listfloat'
    slcfTypes['VECTOR'] = 'bool'
    slcfTypes['SPEC_ID'] = 'string'
    slcfTypes['ID'] = 'ignore'
    return slcfTypes

def getCTRLtypes():
    ctrlTypes = defaultdict(bool)
    ctrlTypes['ID'] = 'string'
    ctrlTypes['FUNCTION_TYPE'] = 'string'
    ctrlTypes['INPUT_ID'] = 'liststring'
    ctrlTypes['DELAY'] = 'float'
    return ctrlTypes

def getZONEtypes():
    zoneTypes = defaultdict(bool)
    zoneTypes['ID'] = 'string'
    zoneTypes['XB'] = 'listfloat'
    zoneTypes['LEAK_AREA'] = 'listindfloat'
    return zoneTypes

def getDUMPtypes():
    dumpTypes = defaultdict(bool)
    dumpTypes['DT_CTRL'] = 'float'
    dumpTypes['DT_HRR'] = 'float'
    dumpTypes['DT_DEVC'] = 'float'
    dumpTypes['DT_BNDF'] = 'float'
    dumpTypes['DT_SLCF'] = 'float'
    dumpTypes['DT_SL3D'] = 'float'
    dumpTypes['DT_RESTART'] = 'float'
    dumpTypes['WRITE_XYZ'] = 'bool'
    dumpTypes['ID'] = 'ignore'
    return dumpTypes

def getTIMEtypes():
    timeTypes = defaultdict(bool)
    timeTypes['T_BEGIN'] = 'float'
    timeTypes['T_END'] = 'float'
    timeTypes['ID'] = 'ignore'
    return timeTypes

def getMISCtypes():
    miscTypes = defaultdict(bool)
    miscTypes['TMPA'] = 'float'
    miscTypes['HUMIDITY'] = 'float'
    miscTypes['SUPPRESSION'] = 'bool'
    miscTypes['BNDF_DEFAULT'] = 'bool'
    miscTypes['ID'] = 'ignore'
    return miscTypes

def getHEADtypes():
    headTypes = defaultdict(bool)
    headTypes['CHID'] = 'string'
    headTypes['TITLE'] = 'string'
    headTypes['ID'] = 'ignore'
    return headTypes
















def parseHEAD(line):
    head = defaultdict(bool)
    if "CHID" in line: head['CHID'] = line.split("CHID")[1].split("'")[1]
    if "TITLE" in line: head['TITLE'] = line.split("TITLE")[1].split("'")[1]
    return head

def parseDEVC(line):
    devc = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        keyTxt = key.split('=')[0].replace(' ','')
        if keyTxt == 'ID':
            devc['ID'] = key.split("'")[1]
        elif keyTxt == 'XYZ':
            devc['XYZ'] = [float(x) for x in key.split('=')[1].replace(' ','').split(',')[:3]]
        elif keyTxt == 'XB':
            devc['XB'] = [float(x) for x in key.split('=')[1].replace(' ','').split(',')[:6]]
        elif keyTxt == "INITIAL_STATE":
            devc['INITIAL_STATE'] = ".%s."%(key.split("=")[1].split(".")[1])
        elif keyTxt == 'QUANTITY':
            devc['QUANTITY'] = key.split('=')[1].split("'")[1]
        elif keyTxt == 'INIT_ID':
            devc['INIT_ID'] = key.split('=')[1].split("'")[1]
        elif keyTxt == 'SETPOINT':
            devc['SETPOINT'] = float(key.split('=')[1].split(',')[0].replace(" ",""))
        elif keyTxt == 'IOR':
            devc['IOR'] = float(key.split('=')[1].split(',')[0].replace(" ",""))
        elif keyTxt == 'TIME_AVERAGED':
            devc['TIME_AVERAGED'] = ".%s."%(key.split("=")[1].split(".")[1])
    return devc

def parseINIT(line):
    init = defaultdict(bool)
    if "XYZ" in line:
        tmp = line.split("XYZ")[1]
        tmp = tmp.replace("=",'').split(',')[:3]
        coords = [float(x) for x in tmp]
        init['XYZ'] = coords
    if "XB" in line:
        tmp = line.split("XB")[1]
        tmp = tmp.replace("=",'').split(',')[:6]
        coords = [float(x) for x in tmp]
        init['XB'] = coords
    if "ID" in line:
        tmp = line.split("ID")[1]
        devID = tmp.split("'")[1]
        init['ID'] = devID
    if "PART_ID" in line:
        tmp = line.split("PART_ID")[1]
        partID = tmp.split("'")[1]
        init['PART_ID'] = partID
    if "N_PARTICLES" in line:
        tmp = line.split("N_PARTICLES")[1]
        tmp = tmp.replace("=",'')
        nParticles = int(tmp.split(",")[0])
        init['N_PARTICLES'] = nParticles
    return init

def parseOBST(line):
    obst = defaultdict(bool)
    if "SURF_ID6" in line:
        obst['SURF_ID6'] = line.split("SURF_ID6")[1].split("'")[1]
    elif "SURF_IDS" in line:
        tmp = line.split('SURF_IDS')[1].split("'")
        tmp = [tmp[1], tmp[3], tmp[5]]
        obst['SURF_IDS'] = tmp
    elif "SURF_ID" in line:
        obst['SURF_ID'] = line.split("SURF_ID")[1].split("'")[1]
    if "ID" in line.replace("SURF_IDS","").replace("SURF_ID6","").replace("SURF_ID","").replace("CTRL_ID",""):
        obst['ID'] = line.replace("SURF_ID","").split("ID")[1].split("'")[1]
    if "XB" in line:
        tmp = np.array([float(x) for x in line.split("XB")[1].replace("=","").split(',')[:6]])
        obst['XB'] =  tmp
    if "BNDF_OBST" in line:
        txt = line.split('BNDF_OBST=')[1].split('.')[1]
        if "TRUE" in txt:
            obst['BNDF_OBST'] = True
        if "FALSE" in txt:
            obst['BNDF_OBST'] = False
    if "CTRL_ID" in line:
        txt = line.split('CTRL_ID')[1].split("'")[1]
        obst['CTRL_ID'] = txt
    if "PERMIT_HOLE" in line:
        txt = line.split('PERMIT_HOLE=')[1].split('.')[1]
        if "TRUE" in txt:
            obst['PERMIT_HOLE'] = True
        if "FALSE" in txt:
            obst['PERMIT_HOLE'] = False
    return obst

def parseHOLE(line):
    hole = defaultdict(bool)
    if "ID" in line.replace('CTRL_ID',''):
        hole['ID'] = line.split("ID")[1].split("'")[1]
    if "XB" in line:
        tmp = np.array([float(x) for x in line.split("XB")[1].replace("=","").split(',')[:6]])
        hole['XB'] =  tmp
    if "CTRL_ID" in line:
        txt = line.split('CTRL_ID')[1].split("'")[1]
        hole['CTRL_ID'] = txt
    return hole

def parseVENT(line):
    vent = defaultdict(bool)
    if "SURF_ID" in line:
        vent['SURF_ID'] = line.split("SURF_ID")[1].split("'")[1]
    if "ID" in line.replace("SURF_ID",""):
        vent['ID'] = line.replace("SURF_ID","").split("ID")[1].split("'")[1]
    if "XB" in line:
        tmp = np.array([float(x) for x in line.split("XB")[1].replace("=","").split(',')[:6]])
        vent['XB'] =  tmp
    if "CTRL_ID" in line:
        tmp = line.split('CTRL_ID')[1].split('=')[1].split(',')[0].split("'")[1]
        vent['CTRL_ID'] = tmp
    if 'MB' in line:
        vent['MB'] = line.split('MB')[1].split("'")[1]
    return vent

def parseSURF(line2):
    line = line2.replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',')
    line = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
    
    regex1 = r"\(\s*\d+\s*,\s*\d+\s*\)"
    regex2 = r"\(\s*\d+\s*:\s*\d+\s*\)"
    try:
        line = re.sub(regex1, regex2, line)
    except:
        pass
    
    keys = line.split(',')
    keys[0] = keys[0].replace('SURF','')
    updatedKeys = []
    txt = ''
    
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    surf = defaultdict(bool)
    for key in updatedKeys:
        if ("ID" in key) and ('MATL_ID' not in key) and ('VOID' not in key) and ('WIDTH' not in key):
            surf['ID'] = key.split('ID')[1].split("'")[1]
        if "MATL_ID" in key:
            #print(key)
            tmp = key.split("=")[1].split(",")
            if '' in tmp: tmp.remove('')
            if ' ' in tmp: tmp.remove(' ')
            tmp = [x.split("'")[1] for x in tmp]
            surf['MATL_ID'] = tmp
        if "COLOR" in key:
            surf['COLOR'] = key.split('COLOR')[1].split("'")[1]
        if "THICKNESS" in key:
            tmp = key.split('=')[1].replace(' ','').split(',')
            for i in range(0, 5):
                if '' in tmp: tmp.remove('')
            surf['THICKNESS'] = [float(x) for x in tmp]                
        if "BACKING" in key:
            surf['BACKING'] = key.split('BACKING')[1].split("'")[1]
        if "GEOMETRY" in key:
            surf['GEOMETRY'] = key.split('GEOMETRY')[1].split("'")[1]
        if "FYI" in key:
            surf['FYI'] = key.split('FYI')[1].split("'")[1]
        if ("LENGTH" in key) and ("CONVECTION_LENGTH_SCALE" not in key):
            surf['LENGTH'] = float(key.split('LENGTH')[1].replace(' ','').replace('=','').split(',')[0])
        if ("CONVECTION_LENGTH_SCALE" in key):
            surf['CONVECTION_LENGTH_SCALE'] = float(key.split('CONVECTION_LENGTH_SCALE')[1].replace(' ','').replace('=','').split(',')[0])
        if ("C_FORCED_RE" in key and "C_FORCED_RE_EXP" not in key):
            surf['C_FORCED_RE'] = float(key.split('=')[1].split(',')[0].replace(' ',''))
        if ("C_FORCED_RE_EXP" in key):
            surf['C_FORCED_RE_EXP'] = float(key.split('=')[1].split(',')[0].replace(' ',''))
        if ("C_FORCED_PR_EXP" in key):
            surf['C_FORCED_PR_EXP'] = float(key.split('=')[1].split(',')[0].replace(' ',''))
        if ("C_HORIZONTAL" in key):
            surf['C_HORIZONTAL'] = float(key.split('=')[1].split(',')[0].replace(' ',''))
        if ("C_VERTICAL" in key):
            surf['C_VERTICAL'] = float(key.split('=')[1].split(',')[0].replace(' ',''))
        if ("EMISSIVITY" in key):
            surf['EMISSIVITY'] = float(key.split('=')[1].split(',')[0].replace(' ',''))
        if "LEAK_PATH" in key:
            tmp = key.split('LEAK_PATH')[1].replace(' ','').replace('=','').split(',')
            if '' in tmp: tmp.remove('')
            tmp = [float(x) for x in tmp]
            surf['LEAK_PATH'] = tmp
        if "HRRPUA" in key:
            tmp = key.split('=')[1].split(',')[0].replace(' ','')
            surf['HRRPUA'] = float(tmp)
        if "RAMP_Q" in key:
            tmp = key.split('=')[1].split("'")[1]
            surf['RAMP_Q'] = tmp
        if "RGB" in key:
            surf['RGB'] = [float(x) for x in key.split('=')[1].split(',')[:3]]
        if ("ADIABATIC" in key) and ("ID" not in key):
            surf['ADIABATIC'] = key.split('.')[1]
        if "VOLUME_FLOW" in key:
            surf['VOLUME_FLOW'] = float(key.split('=')[1].split(',')[0])
            
    return surf

def parseRAMP(line):
    keys = line.split(',')
    updatedKeys = []
    txt = ''
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    ramp = defaultdict(bool)
    for key in updatedKeys:
        if ("ID" in key) and ("DEVC_ID" not in key):
            ramp['ID'] = key.split("'")[1]
        elif ("DEVC_ID" in key):
            ramp['DEVC_ID'] = key.split("'")[1]
        elif "T" in key:
            ramp['T'] = float(key.split("T")[1].split('=')[1].split(',')[0].replace(' ',''))
        elif "F" in key:
            ramp['F'] = float(key.split("F")[1].split('=')[1].split(',')[0].replace(' ',''))
    return ramp

def parseCTRL(line2):
    line = line2.replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',')
    line = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
    
    keys = line.split(',')
    keys[0] = keys[0].replace('CTRL ','')
    updatedKeys = []
    txt = ''
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    ctrl = defaultdict(bool)
    for key in updatedKeys:
        if ("ID" in key) and ('INPUT_ID' not in key):
            ctrl['ID'] = key.split('ID')[1].split("'")[1]
        if "INPUT_ID" in key:
            tmp = key.split("=")[1].split("'")[1:]
            if ',' in tmp: tmp.remove(',')
            if '' in tmp: tmp.remove('')
            if ' ' in tmp: tmp.remove(' ')
            ctrl['INPUT_ID'] = tmp
        if "FUNCTION_TYPE" in key:
            tmp = key.split('=')[1].split("'")[1]
            ctrl['FUNCTION_TYPE'] = tmp
        if ("DELAY" in key) and ("'TIME_DELAY'" not in key):
            tmp = float(key.split('=')[1].split(',')[0].replace(" ",""))
            ctrl['DELAY'] = tmp
    return ctrl

def splitLineIntoKeys(line2):
    line = line2.replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',')
    line = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
    
    keys = line.split(',')
    keys[0] = keys[0][4:]
    updatedKeys = []
    txt = ''
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    return updatedKeys

def parseSLCF(line):
    slcf = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'QUANTITY' in key: slcf['QUANTITY'] = key.split('=')[1].split("'")[1]
        if 'PBX' in key: slcf['PBX'] = float(key.split('=')[1].split(',')[0])
        if 'PBY' in key: slcf['PBY'] = float(key.split('=')[1].split(',')[0])
        if 'PBZ' in key: slcf['PBZ'] = float(key.split('=')[1].split(',')[0])
        if 'VECTOR' in key: slcf['VECTOR'] = key.split('=')[1].split('.')[1]
        if "XB" in key:
            slcf['XB'] = [float(x) for x in key.split('XB')[1].replace('=','').split(',')[:6]]
        if 'SPEC_ID' in key:
            slcf['SPEC_ID'] = key.split("'")[1]
    return slcf

def parseBNDF(line):
    bndf = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'QUANTITY' in key: bndf['QUANTITY'] = key.split('=')[1].split("'")[1]
    return bndf

def parseTIME(line):
    time = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'T_BEGIN' in key: time['T_BEGIN'] = float(key.split('=')[1].split(',')[0])
        if 'T_END' in key: time['T_END'] = float(key.split('=')[1].split(',')[0])
    return time

def parseDUMP(line):
    dump = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'RENDER_FILE' in key: dump['RENDER_FILE'] = key.split("'")[1]
        if 'COLUMN_DUMP_LIMIT' in key: dump['COLUMN_DUMP_LIMIT'] = key.split('.')[1]
        if 'DT_CTRL' in key: dump['DT_CTRL'] = float(key.split('=')[1].split(',')[0])
        if 'DT_DEVC' in key: dump['DT_DEVC'] = float(key.split('=')[1].split(',')[0])
        if 'DT_BNDF' in key: dump['DT_BNDF'] = float(key.split('=')[1].split(',')[0])
        if 'DT_SLCF' in key: dump['DT_SLCF'] = float(key.split('=')[1].split(',')[0])
        if 'DT_SL3D' in key: dump['DT_SL3D'] = float(key.split('=')[1].split(',')[0])
        if 'DT_RESTART' in key: dump['DT_RESTART'] = float(key.split('=')[1].split(',')[0])
        if 'WRITE_XYZ' in key: dump['WRITE_XYZ'] = key.split('.')[1]
    return dump

def parseMISC(line):
    misc = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'TMPA' in key: misc['TMPA'] = float(key.split('=')[1].split(',')[0])
        if 'HUMIDITY' in key: misc['HUMIDITY'] = float(key.split('=')[1].split(',')[0])
        if 'SUPPRESSION' in key: misc['SUPPRESSION'] = key.split('.')[1]
        if 'BNDF_DEFAULT' in key: misc['BNDF_DEFAULT'] = key.split('.')[1]
    return misc

def parseZONE(line):
    zone = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'ID' in key: zone['ID'] = line.split("'")[1]
        if 'XB' in key: zone['XB'] = [float(x) for x in key.split('XB')[1].replace('=','').split(',')[:6]]
    return zone

def parseMESH(line):
    mesh = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if "ID" in key:
            mesh['ID'] = key.split('ID')[1].replace('=','').split("'")[1]
        if 'IJK' in key:
            mesh['IJK'] = [float(x) for x in key.split('IJK')[1].replace('=','').split(',')]
        if 'XB' in key:
            mesh['XB'] = [float(x) for x in key.split('XB')[1].replace('=','').split(',')[:6]]
        if 'MPI_PROCESS' in key:
            mesh['MPI_PROCESS'] = key.split('MPI_PROCESS')[1].replace('=','').replace(',','')
    return mesh

def parseREAC(line):
    reac = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        keyTxt = key.split('=')[0].replace(' ','')
        if (keyTxt == 'C'): reac['C'] = float(key.split('=')[1].split(',')[0])
        if (keyTxt == 'H'): reac['H'] = float(key.split('=')[1].split(',')[0])
        if (keyTxt == 'O'): reac['O'] = float(key.split('=')[1].split(',')[0])
        if (keyTxt == 'N'): reac['N'] = float(key.split('=')[1].split(',')[0])
        if (keyTxt == 'FYI'): reac['FYI'] = key.split('=')[1].split("'")[1]
        if ("ID" in key) and ("IDEAL" not in key): reac['ID'] = key.split('ID')[1].replace('=','').split("'")[1]
        if "FUEL" in key: reac['FUEL'] = key.split("'")[1]
        if "FORMULA" in key: reac['FORMULA'] = key.split("'")[1]
        if "CO_YIELD" in key: reac['CO_YIELD'] = float(key.split('=')[1].split(',')[0])
        if "SOOT_YIELD" in key: reac['SOOT_YIELD'] = float(key.split('=')[1].split(',')[0])
        if "IDEAL" in key: reac['IDEAL'] = key.split('.')
        if "HEAT_OF_COMBUSTION" in key: reac['HEAT_OF_COMBUSTION'] = float(key.split('=')[1].split(',')[0])
    return reac

def parseRADI(line):
    radi = defaultdict(bool)
    radi['ID'] = 'RADIATION'
    keys = splitLineIntoKeys(line)
    for key in keys:
        if ('RADIATIVE_FRACTION' in key): radi['RADIATIVE_FRACTION'] = float(key.split('=')[1].split(',')[0])
    return radi

def parseMATL(line):
    matl = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if "ID" in key: matl['ID'] = key.split("'")[1]
        if "FYI" in key: matl['FYI'] = key.split("'")[1]
        if ("SPECIFIC_HEAT_RAMP" in key):
            matl['SPECIFIC_HEAT_RAMP'] = key.split('=')[1].split(',')[0]
        elif ("SPECIFIC_HEAT" in key):
            matl['SPECIFIC_HEAT'] = float(key.split('=')[1].split(',')[0])
        if ("CONDUCTIVITY_RAMP" in key):
            matl['CONDUCTIVITY_RAMP'] = key.split('=')[1].split(',')[0]
        elif ("CONDUCTIVITY" in key):
            matl['CONDUCTIVITY'] = float(key.split('=')[1].split(',')[0])
        if "DENSITY" in key: matl['DENSITY'] = float(key.split('=')[1].split(',')[0])
        if "EMISSIVITY" in key: matl['EMISSIVITY'] = float(key.split('=')[1].split(',')[0])
    return matl

def parsePRES(line):
    pres = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        keyTxt = key.split('=')[0].replace(' ','')
        if (keyTxt == 'SOLVER'): pres['SOLVER'] = key.split('=')[1].split("'")[1]
    return pres
    
    