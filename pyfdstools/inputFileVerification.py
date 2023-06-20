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
import scipy.spatial as scsp
from .smokeviewParser import parseSMVFile
from .utilities import pointsFromXB, in_hull

def checkDevices(fdsFile, smvFile=None):
    ''' This routine checks if any devices in an fds file are overlapping
    with an obstruction
    
    Inputs:
        fdsFile: fdsFileOperations with the different lines defined
        smvFile: filename of the smv file. If not included, determined by chid
    Outputs:
        devcNames: device names that are overlapping with an obstruction
    '''
    chid = fdsFile.head['ID']['CHID']
    devcs = fdsFile.devcs
    obsts = fdsFile.obsts
    surfs = fdsFile.surfs
    if smvFile is None:
        smvFile = fdsFile.split(os.sep)
        smvFile[-1] = "%s.smv"%(chid)
        smvFile = os.sep.join(smvFile)
    grids, obsts, bndfs, surfs, files, bndes = parseSMVFile(smvFile)
    
    obstHulls = []
    for i in range(0,obsts.shape[0]):
        obst = obsts[i,:]
        pts = pointsFromXB(obst[13:19])
        obstHull = scsp.Delaunay(pts)
        obstHulls.append(obstHull)
    
    devcNames = []
    for key in list(devcs.keys()):
        pt = devcs[key]['XYZ']
        if pt:
            for hull in obstHulls:
                if in_hull(pt,hull):
                    print("Device %s is within obstruction."%(key))
                    devcNames.append(key)
    devcNames = list(set(devcNames))
    return devcNames