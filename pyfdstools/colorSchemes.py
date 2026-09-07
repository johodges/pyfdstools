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
# EXAMPLES:
# See the examples subroutine for example operation.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import numpy as np
from matplotlib.colors import ListedColormap

def getVTcolors():
    """Returns the Virginia Tech brand color sequence

    Returns
    -------
    array(13, 3)
        Array of RGB values in the range 0 to 1
    """

    colors = np.array([[134, 31, 65],
                       [232, 119, 34],
                       [117, 120, 123],
                       [200, 200, 200],
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
    """Returns an alternative categorical color sequence

    Returns
    -------
    array(16, 3)
        Array of RGB values in the range 0 to 1
    """

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

def buildSMVcolormap(percentile=None, width=None):
    """Builds the smokeview blue-cyan-green-yellow-red colormap

    Optionally a black band can be drawn across the colormap so that one
    value stands out, which plotSlice uses to highlight a threshold such
    as a tenability criterion.

    Parameters
    ----------
    percentile : float, optional
        Position of the highlight band as a fraction of the color scale,
        in the range 0 to 1. No band is drawn when omitted
    width : int, optional
        Half-width of the highlight band in colormap entries. Defaults
        to 2 when a percentile is given without one

    Returns
    -------
    matplotlib.colors.ListedColormap
        Colormap with 256 entries
    """

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
    if percentile is not None:
        # width is optional at the call site (plotSlice passes its
        # highlightWidth straight through, and that defaults to None),
        # so supply a default rather than raising a TypeError on the
        # subtraction below. Clamp the band to the colormap so that a
        # highlight near either end does not wrap around.
        if width is None:
            width = 2
        width = int(width)
        pind = int(percentile*newcmp.shape[0])
        lo = max(0, pind - width)
        hi = min(newcmp.shape[0], pind + width)
        if hi > lo:
            newcmp[lo:hi, :3] = 0.0
    cmap = ListedColormap(newcmp)
    
    return cmap
