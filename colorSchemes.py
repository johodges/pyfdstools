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
# EXAMPLES:
# See the examples subroutine for example operation.
#
#=========================================================================
# # IMPORTS
#=========================================================================
import numpy as np
from matplotlib.colors import ListedColormap

def getVTcolors():
    colors = np.array([[134, 31, 65],
                       [232, 119, 34],
                       [117, 120, 123],
                       [255, 255, 255],
                       [80, 133, 144],
                       [247, 234, 72],
                       [206, 0, 88],
                       [100, 38, 103],
                       [237, 139, 0],
                       [44, 213, 196],
                       [229, 225, 230],
                       [215, 210, 203],
                       [198, 70, 0]], dtype=np.float32)
    colors = colors/255
    return colors

def getJHcolors():
    colors = np.array([[0, 65, 101],
                       [229, 114, 0],
                       [136, 139, 141],
                       [170, 39, 44],
                       [119, 197, 213],
                       [161, 216, 132],
                       [255, 200, 69],
                       [101, 0, 65],
                       [0, 229, 114],
                       [141, 136, 139],
                       [44, 170, 39],
                       [213, 119, 197],
                       [132, 161, 216],
                       [69, 255, 200],
                       [65, 101, 0],
                       [114, 0, 229]], dtype=np.float32)
    colors = colors/255
    return colors

def buildSMVcolormap():
    newcmp = np.zeros((256,4))
    
    colors = np.array([[0,0,1,1],
              [0,1,1,1],
              [0,1,0,1],
              [1,1,0,1],
              [1,0,0,1],])
    colorInds = np.array([0, 64, 128, 192, 255])
    
    j = 0
    for i in range(0,len(newcmp)):
        if i == colorInds[j]:
            newcmp[i,:] = colors[j,:]
            j = j + 1
        else:
            m = (colors[j,:]-colors[j-1,:])/(colorInds[j]-colorInds[j-1])
            b = colors[j,:]-m*(colorInds[j]-colorInds[j-1])
            newcmp[i] = colors[j-1,:]+m*(i-colorInds[j-1])
    cmap = ListedColormap(newcmp)
    
    return cmap
