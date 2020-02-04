# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 08:17:40 2019

@author: JHodges
"""

import os
import scipy.spatial as scsp

#from fdsTools import fdsParser
#from fdsTools import smokeviewParser

def in_hull(p,hull):
    if not isinstance(hull,scsp.Delaunay):
        hull = scsp.Delaunay(hull)
    return hull.find_simplex(p)>=0

def pointsFromXB(XB,extend=[0,0,0]):
    pts = [[XB[0]-extend[0],XB[2]-extend[1],XB[4]-extend[2]],
           [XB[0]-extend[0],XB[2]-extend[1],XB[5]+extend[2]],
           [XB[0]-extend[0],XB[3]+extend[1],XB[4]-extend[2]],
           [XB[0]-extend[0],XB[3]+extend[1],XB[5]+extend[2]],
           [XB[1]+extend[0],XB[2]-extend[1],XB[4]-extend[2]],
           [XB[1]+extend[0],XB[2]-extend[1],XB[5]+extend[2]],
           [XB[1]+extend[0],XB[3]+extend[1],XB[4]-extend[2]],
           [XB[1]+extend[0],XB[3]+extend[1],XB[5]+extend[2]]]
    return pts

def checkDevices(fdsFile,smvFile=None):
    chid, devcs, obsts, vents, surfs, ramps, ctrls = fdsParser.parseFDSFile(fdsFile)
    if smvFile is None:
        smvFile = fdsFile.split(os.sep)
        smvFile[-1] = "%s.smv"%(chid)
        smvFile = os.sep.join(smvFile)
    grids, obsts, bndfs, surfs  = smokeviewParser.parseSMVFile(smvFile)
    
    obstHulls = []
    for i in range(0,obsts.shape[0]):
        obst = obsts[i,:]
        pts = pointsFromXB(obst[13:19])
        obstHull = scsp.Delaunay(pts)
        obstHulls.append(obstHull)
    
    for key in list(devcs.keys()):
        pt = devcs[key]['XYZ']
        if pt:
            for hull in obstHulls:
                if in_hull(pt,hull):
                    print("Device %s is within obstruction."%(key))
