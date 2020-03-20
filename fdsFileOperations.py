# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 08:17:40 2019

@author: JHodges
"""

import numpy as np
from collections import defaultdict
import datetime
import re
import scipy.spatial as scsp
import os
import zipfile
from .fdsTypes import fdsLineTypes

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
        self.customLines = []
        
        '''
        self.head['ID'] = defaultdict(bool)
        self.time['ID'] = defaultdict(bool)
        self.dump['ID'] = defaultdict(bool)
        self.misc['ID'] = defaultdict(bool)
        self.pres['ID'] = defaultdict(bool)
        '''
        
        self.devcs['unknownCounter'] = 0
        self.obsts['unknownCounter'] = 0
        self.holes['unknownCounter'] = 0
        self.vents['unknownCounter'] = 0
        self.meshes['unknownCounter'] = 0
        self.slcfs['unknownCounter'] = 0
        self.bndfs['unknownCounter'] = 0
        
        self.meshOrder = False
        self.version = "6.7.4"
    
    def dictMerge(self, a, b, path=None):
        "merges b into a"
        if path is None: path = []
        for key in b:
            if key in a:
                if isinstance(a[key], dict) and isinstance(b[key], dict):
                    self.dictMerge(a[key], b[key], path + [str(key)])
                elif a[key] == b[key]:
                    pass
                else:
                    a[key] = b[key]
            else:
                a[key] = b[key]
        return a
    
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
        #linesFDS = [x.split('/')[0] for x in textFDS.split("&")[1:]]
        linesFDS = [x for x in textFDS.split("&")[1:]]
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
                TIME_AVERAGED=None, SPATIAL_STATISTIC=None, STATISTICS=None,
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
        if STATISTICS != None: devc["STATISTICS"] = STATISTICS
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
        meshKeys = list(self.meshes.keys())
        try:
            meshKeys.remove('unknownCounter')
        except:
            pass
        for key in meshKeys:
            IJK = self.meshes[key]['IJK']
            numCells.append(IJK[0]*IJK[1]*IJK[2])
            meshes.append(key)
        return meshes, numCells
    
    def addMPIprocesses(self, numberOfProcesses, allowMeshSplitting=True):
        meshes, numCells = self.calculateMeshCells()
        cellsPerProcess = np.sum(numCells)/numberOfProcesses
        mpiConverged = False
        splitConverged = False
        splitMultiplier = 1.20
        while not mpiConverged:
            mpiConverged = True
            while not splitConverged and allowMeshSplitting:
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
            if np.max(mpiProcess) > cellsPerProcess*1.20 and allowMeshSplitting:
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
        mesh2['ID'] = "%s-00"%(mesh["ID"])
        mesh2['XB'] = XB2
        mesh2['IJK'] = IJK2
        
        mesh3 = defaultdict(bool)
        mesh3['ID'] = "%s-01"%(mesh["ID"])
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
    
    def keyFromLineType(self, lineType):
        if lineType == 'HEAD': key = 'head'
        if lineType == 'DEVC': key = 'devcs'
        if lineType == 'INIT': key = 'init'
        if lineType == 'OBST': key = 'obsts'
        if lineType == 'VENT': key = 'vents'
        if lineType == 'SURF': key = 'surfs'
        if lineType == 'RAMP': key = 'ramps'
        if lineType == 'CTRL': key = 'ctrls'
        if lineType == 'MESH': key = 'meshes'
        if lineType == 'SLCF': key = 'slcfs'
        if lineType == 'BNDF': key = 'bndfs'
        if lineType == 'TIME': key = 'time'
        if lineType == 'DUMP': key = 'dump'
        if lineType == 'MISC': key = 'misc'
        if lineType == 'ZONE': key = 'zones'
        if lineType == 'REAC': key = 'reacs'
        if lineType == 'MATL': key = 'matls'
        if lineType == 'RADI': key = 'radis'
        if lineType == 'PRES': key = 'pres'
        if lineType == 'HOLE': key = 'holes'
        return key
    
    def mergeTypeFromLineType(self, lineType):
        key = 'unknown'
        if lineType == 'HEAD': key = 'merge'
        if lineType == 'DEVC': key = 'enumerate'
        if lineType == 'INIT': key = 'enumerate'
        if lineType == 'OBST': key = 'enumerate'
        if lineType == 'VENT': key = 'enumerate'
        if lineType == 'SURF': key = 'enumerate'
        if lineType == 'RAMP': key = 'append'
        if lineType == 'CTRL': key = 'enumerate'
        if lineType == 'MESH': key = 'enumerate'
        if lineType == 'SLCF': key = 'enumerate'
        if lineType == 'BNDF': key = 'enumerate'
        if lineType == 'TIME': key = 'merge'
        if lineType == 'DUMP': key = 'merge'
        if lineType == 'MISC': key = 'merge'
        if lineType == 'ZONE': key = 'enumerate'
        if lineType == 'REAC': key = 'enumerate'
        if lineType == 'MATL': key = 'enumerate'
        if lineType == 'RADI': key = 'merge'
        if lineType == 'PRES': key = 'merge'
        if lineType == 'HOLE': key = 'enumerate'
        return key
    
    def parseLine(self, line, lineType, types, key):
        lineDict = self.dictFromLine(line, lineType, types)
        tmp = getattr(self, key)
        mergeType = self.mergeTypeFromLineType(lineType)
        if mergeType == 'merge':
            if not tmp['ID']: tmp['ID'] = defaultdict(bool)
            tmp['ID'] = self.dictMerge(tmp['ID'], lineDict)
            setattr(self, key, tmp)
        elif mergeType == 'append':
            ID = lineDict['ID']
            if tmp[ID]:
                for keyID2 in list(lineDict.keys()):
                    keyType = getattr(types, lineType.lower())[keyID2]
                    keyValue = lineDict[keyID2]
                    if (keyType == 'listrowfloat'):
                        for v in keyValue:
                            tmp[ID][keyID2].append(v)
            else:
                tmp[ID] = lineDict
        elif mergeType == 'enumerate':
            ID = lineDict['ID']
            if ID is False: ID = "ID"
            if tmp[ID]:
                counter = tmp[ID]['counter']
                tmp["%s-%04.0f"%(ID, counter)] = lineDict
                tmp[ID]['counter'] += 1
                pass
            else:
                tmp[ID] = lineDict
                tmp[ID]['counter'] = 0
        else:
            assert False, "Stopped"
    
    def interpretKey(self, key, lineType, types):
        keyID = key.split('=')[0]
        keyValue = '='.join(key.split('=')[1:])
        #keyValue = key.replace(keyID, '')[1:]
        regex1 = r"\(\s*.*\)"
        regex2 = r""
        try:
            keyID2 = re.sub(regex1, regex2, keyID)
        except:
            keyID2 = keyID
        while keyID2[-1] == ' ':
            keyID2 = keyID2[:-1]
        keyType = getattr(types, lineType.lower())[keyID2]
        return keyID, keyID2, keyType, keyValue
    
    def dictFromLine(self, line, lineType, types):
        lineDict = defaultdict(bool)
        keys = self.splitLineIntoKeys(line)
        
        for key in keys:
            keyID, keyID2, keyType, keyValue = self.interpretKey(key, lineType, types)
            if keyType == 'string':
                keyValue = keyValue.split("'")[1]
            elif keyType == 'float':
                keyValue = float(keyValue.replace(' ', '').replace(',','').replace('/',''))
            elif keyType == 'int':
                keyValue = int(keyValue.replace(' ', '').replace(',','').replace('/',''))
            elif keyType == 'bool':
                keyValue = keyValue.split(".")[1]
            elif (keyType == 'liststring') or (keyType == 'listfloat') or (keyType == 'listint'):
                vals = []
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                for t in keyValues:
                    if 'string' in keyType: preprocess = t.split("'")[1]
                    if 'float' in keyType: preprocess = float(t.replace(' ', '').replace(',','').replace('/',''))
                    if 'int' in keyType: preprocess = int(t.replace(' ', '').replace(',','').replace('/',''))
                    vals.append(preprocess)
                keyValue = vals
            elif (keyType == 'listindfloat'):
                vals = []
                ind = int(keyID.replace(keyID2, '').replace('(','').replace(')',''))
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                for t in keyValues:
                    if 'string' in keyType: preprocess = t
                    if 'float' in keyType: preprocess = float(t.replace(' ', '').replace(',','').replace('/',''))
                    if 'int' in keyType: preprocess = int(t.replace(' ', '').replace(',','').replace('/',''))
                    vals.append([ind, preprocess])
                keyValue = vals
            elif (keyType == 'listrowfloat'):
                vals = []
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                for t in keyValues:
                    if 'string' in keyType: preprocess = t
                    if 'float' in keyType: preprocess = float(t.replace(' ', '').replace(',','').replace('/',''))
                    if 'int' in keyType: preprocess = int(t.replace(' ', '').replace(',','').replace('/',''))
                    vals.append(preprocess)
                keyValue = vals
            else:
                print(lineType.lower(), keyID, keyID2, keyType)
                print(len(keyID))
                assert False, "Stopped"
            lineDict[keyID2] = keyValue
        return lineDict
        
    def parseFDSLines(self, lines):
        for line in lines:
            lineType = self.getLineType(line)
            key = self.keyFromLineType(lineType)
            types = fdsLineTypes(version=self.version)
            self.parseLine(line, lineType, types, key)
        devcKeys = list(self.devcs.keys())
        devcKeys.remove('unknownCounter')
        for key in devcKeys:
            if self.devcs[key]['INIT_ID']:
                self.devcs[key]['XYZ'] = self.inits[self.devcs[key]['INIT_ID']]['XYZ']
            else:
                self.devcs[key].pop('INIT_ID')
    
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
        types = fdsLineTypes(version=self.version)
        text = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        text = "%s!!!!! Input file generated with fdsTools v1 on %04.0f-%02.0f-%02.0f          !!!!!\n"%(text, date.year, date.month, date.day)
        text = "%s!!!!! Copyright JENSEN HUGHES %04.0f All right reserved.             !!!!!\n"%(text, date.year)
        text = "%s!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"%(text)
        
        text = "%s%s"%(text, self.makeLinesFromDict(self.head, types.head, '&HEAD', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.time, types.time, '&TIME', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.misc, types.misc, '&MISC', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.inits, types.init, '&INIT', newline=True))
        text = "%s%s"%(text, self.makeLinesFromDict(self.dump, types.dump, '&DUMP', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.zones, types.zone, '&ZONE', newline=True))
        text = "%s%s"%(text, self.makeLinesFromDict(self.pres, types.pres, '&PRES', newline=True))
        
        if self.meshOrder is False:
            self.addMPIprocesses(1)
        text = "%s%s"%(text, self.makeMESH(self.meshes, types.mesh, meshOrder=self.meshOrder))
        text = "%s%s"%(text, self.makeLinesFromDict(self.reacs, types.reac, '&REAC', newline=True))
        text = "%s%s"%(text, self.makeLinesFromDict(self.radis, types.radi, '&RADI', newline=True))
        text = "%s%s"%(text, self.makeLinesFromDict(self.matls, types.matl, '&MATL', newline=True))
        text = "%s%s"%(text, self.makeLinesFromDict(self.surfs, types.surf, '&SURF', newline=True))
        text = "%s%s"%(text, self.makeRAMP(self.ramps))
        
        text = "%s%s"%(text, self.makeLinesFromDict(self.obsts, types.obst, '&OBST', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.holes, types.hole, '&HOLE', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.vents, types.vent, '&VENT', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.devcs, types.devc, '&DEVC', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.ctrls, types.ctrl, '&CTRL', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.bndfs, types.bndf, '&BNDF', newline=False))
        text = "%s%s"%(text, self.makeLinesFromDict(self.slcfs, types.slcf, '&SLCF', newline=False))
        
        for line in self.customLines:
            text = "%s%s\n"%(text, line)
        
        return text

    def keyAssist(self, text, types, dic, internalKeys=['counter'], newline=False):
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
    
    def makeLinesFromDict(self, items, types, prefix, newline=False):
        text = ''
        keys = list(items.keys())
        keys.sort()
        if 'unknownCounter' in keys: keys.remove('unknownCounter')
        for key in keys:
            text = "%s%s "%(text, prefix)
            text = self.keyAssist(text, types, items[key], newline=newline)
            text = "%s /\n"%(text)
        return text
    
    def makeMESH(self, meshes, meshTypes, meshOrder=False):
        text = ''
        meshList = list(meshes.keys())
        if 'unknownCounter' in meshList: meshList.remove('unknownCounter')
        if (meshOrder is not False): meshList = [meshList[x] for x in meshOrder]
        for key in meshList:
            text = "%s&MESH "%(text)
            text = self.keyAssist(text, meshTypes, meshes[key])
            text = "%s /\n"%(text)
        return text
    
    def makeRAMP(self, ramps):
        text = ''
        for key in list(ramps.keys()):
            ID = ramps[key]['ID']
            for F, T in zip(ramps[key]['F'], ramps[key]['T']):
                text = "%s&RAMP ID='%s', T = %0.4f, F = %0.4f, /\n"%(text, ID, T, F)
        return text
    
    def splitLineIntoKeys(self, line2):
        line = line2.replace('\n', ',').replace('\r', ',')
        while (',,' in line) or ('  ' in line):
            line = line.replace(',,', ',').replace('  ', ' ')    
            
        regex1 = r"\(\s*\d+\s*,\s*\d+\s*\)"
        regex2 = r"\(\s*\d+\s*:\s*\d+\s*\)"
        try:
            line = re.sub(regex1, regex2, line)
        except:
            pass
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
        while '' in updatedKeys:
            updatedKeys.remove('')
        for i, txt in enumerate(updatedKeys):
            while txt[0] == ' ':
                txt = txt[1:]
            updatedKeys[i] = txt
        for i, txt in enumerate(updatedKeys):
            while txt[-1] == ' ' or txt[-1] == ',':
                txt = txt[:-1]
            updatedKeys[i] = txt
        return updatedKeys
    