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
import numpy as np
import matplotlib.pyplot as plt
import glob
import warnings
import zipfile
import os
import struct
import scipy.interpolate as scpi
import pandas as pd
from collections import defaultdict
from .fdsFileOperations import fdsFileOperations
from .utilities import getDatatypeByEndianness, getEndianness
from .utilities import getFileListFromZip, getFileList, getSmvFile
from .utilities import zopen, zreadlines
from .utilities import getFileListFromResultDir
from .utilities import getGridsFromXyzFiles, getAbsoluteGrid, rearrangeGrid
from .utilities import readXYZfile
from .colorSchemes import buildSMVcolormap
from .smokeviewParser import parseSMVFile
from collections.abc import Iterable

def _resolveResultDir(workingDir, chid):
    """Resolves the directory FDS actually wrote its output files to

    When the input file sets RESULTS_DIR on the &DUMP namelist, FDS
    writes its output into that subdirectory rather than alongside the
    input file. Reading it back therefore requires resolving the input
    file first. Zip archives are returned unchanged, since the archive
    path is used directly for lookups inside it.

    Parameters
    ----------
    workingDir : str
        Directory containing the FDS input file, or a zip archive
    chid : str
        FDS CHID of the case

    Returns
    -------
    str
        Directory the output files were written to
    """

    if '.zip' in workingDir:
        return workingDir
    fdsFiles = getFileList(workingDir, chid, 'fds')
    if len(fdsFiles) == 0:
        return workingDir
    fdsFile = fdsFileOperations()
    fdsFile.importFile(fdsFiles[0])
    if fdsFile.dump['ID'] is not False:
        if fdsFile.dump['ID']['RESULTS_DIR'] is not False:
            return (workingDir + os.sep
                    + fdsFile.dump['ID']['RESULTS_DIR'] + os.sep)
    return workingDir


def time2str(time, decimals=2):
    """Converts a timestamp to a string
    
    This subroutine converts a float timestamp to a string with the
    period replaced by an underscore.
    
    Parameters
    ----------
    time : float
        Timestamp
    decimals : int, optional
        Number of decimal places to include (default 2)
    
    Returns
    -------
    str
        Timestamp with decimal replaced by underscore
        number of patches
    """
    decStr = '{:.{prec}f}'.format(time - np.floor(time), prec=decimals)
    decStr = decStr.replace('0.','_')
    timeStr = '%08.0f%s'%(np.floor(time), decStr)
    return timeStr

def mesh2str(mesh):
    """Converts a mesh number to string
    
    This subroutine converts an integer mesh number to a string.
    FDS specifies mesh numbers as zero-padded integers with 4 digits.
    
    Parameters
    ----------
    mesh : int
        Integer mesh number
    
    Returns
    -------
    str
        Mesh number as a string zero-padded to 4 digits
    """
    
    meshStr = "%04.0f"%(mesh)
    return meshStr


def readP3Dfile(file):
    """Reads data from plot3D file
    
    This subroutine reads data from a plot3D file.
    TODO: Update this subroutine to use zopen
    
    Parameters
    ----------
    file : str
        String containing the path to a plot3D file
    
    Returns
    -------
    array(NX, NY, NZ, NT)
        Array containing float data in local coordinates for each time
    array()
        Array containing header information from plot3D file
    """
    
    f = zopen(file)
    data1 = f.read()
    f.close()
    header = np.frombuffer(data1, dtype=np.int32, count=5)
    _ = np.frombuffer(data1, dtype=np.float32, count=7)
    (nx, ny, nz) = (header[1], header[2], header[3])
    data = np.frombuffer(data1, dtype=np.float32, count=nx*ny*nz*5)
    data = np.reshape(data, (int(data.shape[0]/5),5), order='F')
    return data, header[1:-1]

def writeP3Dfile(file, data):
    """Writes data to plot3D file
    
    This subroutine writes data to a plot3D file.
    
    Parameters
    ----------
    file : str
        String containing the path to a plot3D file
    
    Returns
    -------
    array(NX, NY, NZ, NT)
        Array containing float data in local coordinates for each time
    array()
        Array containing header information from plot3D file
    """
    
    with open(file,'wb') as f:
        f.write(b'\x0c\x00\x00\x00')
        nx, ny, nz, v = data.shape
        (nx, ny, nz, v) = (int(nx), int(ny), int(nz), int(v))
        f.write(nx.to_bytes(4, 'little'))
        f.write(ny.to_bytes(4, 'little'))
        f.write(nz.to_bytes(4, 'little'))
        f.write(b'\x0c\x00\x00\x00')
        empty = np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.float32)
        empty.tofile(f)
        d1 = data.flatten(order='F')
        d = d1.tobytes()
        f.write(d)


def buildDataFile(meshStr, time):
    """Builds a plot3D datafile name
    
    Parameters
    ----------
    meshStr : mesh number string
    time : float timestamp
    
    Returns
    -------
    str
        String containing the name of the plot3D data file
    """
    
    dataStr = meshStr.replace('.xyz','_%s.q'%(time2str(time)))
    return dataStr
    
def printExtents(grid, data):
    """Prints the extents of plot3D grid and data
    
    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Array containing float global coordinates
    data : array(NX, NY, NZ, 5, NT)
        Array containing float data in local coordinates for each time
    """
    
    print("%0.4f < x < %0.4f"%(np.min(grid[:,0]),np.max(grid[:,0])))
    print("%0.4f < y < %0.4f"%(np.min(grid[:,1]),np.max(grid[:,1])))
    print("%0.4f < z < %0.4f"%(np.min(grid[:,2]),np.max(grid[:,2])))
    print("%0.4f < IBLK < %0.4f"%(np.min(grid[:,3]),np.max(grid[:,3])))
    
    print("")
    
    print("%0.4f < T < %0.4f"%(np.min(data[:,0]),np.max(data[:,0])))
    print("%0.4f < U < %0.4f"%(np.min(data[:,1]),np.max(data[:,1])))
    print("%0.4f < V < %0.4f"%(np.min(data[:,2]),np.max(data[:,2])))
    print("%0.4f < W < %0.4f"%(np.min(data[:,3]),np.max(data[:,3])))
    print("%0.4f < HRRPUV < %0.4f"%(
            np.min(data[:,4]),np.max(data[:,4])))
        

def findSliceLocation(grid, data, axis, value, plot3d=False):
    """Extracts slice location from absolute grid
    
    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Array containing grid coordinates
    data : array(NX, NY, NZ)
        Array containing float data for each grid coordinate
    axis : int
        Axis number to query
    value : float
        Axis value to query
    plot3d : bool
        Flag specifying whether data is plot3d or not
    
    Returns
    -------
    array(NX, NZ)
        Array containing the 2-D slice x coordinates
    array(NX, NZ)
        Array containing the 2-D slice z coordinates
    array(NX, NZ)
        Array containing the 2-D slice data
    """
    
    xGrid = grid[:, :, :, 0]
    yGrid = grid[:, :, :, 1]
    zGrid = grid[:, :, :, 2]
    if axis == 1: 
        nGrid = xGrid
        nValues = nGrid[:, 0, 0]
    if axis == 2: 
        nGrid = yGrid
        nValues = nGrid[0, :, 0]
    if axis == 3: 
        nGrid = zGrid
        nValues = nGrid[0, 0, :]
    ind = np.argmin(abs(nValues - value))
    if axis == 1:
        x = np.squeeze(yGrid[ind, :, :])
        z = np.squeeze(zGrid[ind, :, :])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[ind, :, :, :])
        else:
            data2 = np.squeeze(data[ind, :, :])
    elif axis == 2:
        x = np.squeeze(xGrid[:, ind, :])
        z = np.squeeze(zGrid[:, ind, :])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:, ind, :, :])
        else:
            data2 = np.squeeze(data[:, ind, :])
    elif axis == 3:
        x = np.squeeze(xGrid[:, :, ind])
        z = np.squeeze(yGrid[:, :, ind])
        if len(data.shape) == 5:
            data2 = np.squeeze(data[:, :, ind, :])
        else:
            data2 = np.squeeze(data[:, :, ind])
    if plot3d:
        T = np.squeeze(data2[:, :, 0])
        U = np.squeeze(data2[:, :, 1])
        V = np.squeeze(data2[:, :, 2])
        W = np.squeeze(data2[:, :, 3])
        HRR = np.squeeze(data2[:, :, 4])
        return x, z, T, U, V, W, HRR
    else:
        return x, z, data2

def plotSlice(x, z, data_slc, axis, fig=None, ax=None,
              cmap=None, figsize=None, fs=16, figsizeMult=4,
              qnty_mn=None, qnty_mx=None, cbarnumticks=10,
              tickDecimals=None,
              levels=None, cbarticks=None, clabel=None,
              highlightValue=None, highlightWidth=None,
              reverseXY=False, percentile=None,
              xmn=None, xmx=None, zmn=None, zmx=None,
              xlabel=None, zlabel=None,
              addCbar=True, fixXLims=True, fixZLims=True,
              title=None, contour=True, extend='both',
              linecontour=False,
              returnIm=False):
    """Plots a 2-D slice as a filled contour or image

    Parameters
    ----------
    x : array(N, M)
        First in-plane coordinate of each point
    z : array(N, M)
        Second in-plane coordinate of each point
    data_slc : array(N, M)
        Values to plot
    axis : int
        Axis normal to the slice (1 = x, 2 = y, 3 = z), used to pick the
        default axis labels
    fig : matplotlib.figure.Figure, optional
        Figure to draw into. A new figure is created when omitted
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. New axes are created when omitted
    cmap : str or matplotlib.colors.Colormap, optional
        Colormap to use. The smokeview colormap is used when omitted,
        as it is when the string 'SMV' is given
    figsize : tuple, optional
        Figure size. Derived from the aspect ratio of the slice when
        omitted
    fs : int, optional
        Font size (default 16)
    figsizeMult : int, optional
        Length in inches of the shorter figure axis (default 4)
    qnty_mn : float, optional
        Lower limit of the color scale. Data minimum when omitted
    qnty_mx : float, optional
        Upper limit of the color scale. Data maximum when omitted
    cbarnumticks : int, optional
        Number of colorbar ticks (default 10)
    tickDecimals : int, optional
        Number of decimals shown on the colorbar ticks
    levels : int or list, optional
        Number of contour levels, or explicit level values
    cbarticks : list, optional
        Explicit colorbar tick locations
    clabel : str, optional
        Label for the colorbar
    highlightValue : float, optional
        Value at which the smokeview colormap places its highlight band
    highlightWidth : float, optional
        Width of the highlight band
    reverseXY : bool, optional
        Swap the horizontal and vertical axes (default False)
    percentile : float, optional
        Position of the colormap highlight, in the range 0 to 1
    xmn : float, optional
        Lower horizontal axis limit. Data minimum when omitted
    xmx : float, optional
        Upper horizontal axis limit. Data maximum when omitted
    zmn : float, optional
        Lower vertical axis limit. Data minimum when omitted
    zmx : float, optional
        Upper vertical axis limit. Data maximum when omitted
    xlabel : str, optional
        Horizontal axis label. Derived from axis when omitted
    zlabel : str, optional
        Vertical axis label. Derived from axis when omitted
    addCbar : bool, optional
        Draw a colorbar (default True)
    fixXLims : bool, optional
        Apply xmn and xmx to the axes (default True)
    fixZLims : bool, optional
        Apply zmn and zmx to the axes (default True)
    title : str, optional
        Axes title
    contour : bool, optional
        Draw filled contours rather than an image (default True)
    extend : str, optional
        Which end of the color scale is extended: 'both', 'below' (also
        spelled 'min'), 'above' (also spelled 'max') or 'neither'
        (default 'both')
    linecontour : bool, optional
        Draw line contours rather than filled contours (default False)
    returnIm : bool, optional
        Also return the mappable produced by the plotting call

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    matplotlib.axes.Axes
        Axes containing the plot
    matplotlib.cm.ScalarMappable
        The mappable, returned only when returnIm is True
    """

    # 'below' and 'above' are pyfdstools spellings which describe which
    # end of the scale is extended. matplotlib names the same two
    # options 'min' and 'max'; passing the pyfdstools spelling straight
    # through raised a ValueError from contourf and colorbar.
    extendAliases = {'below': 'min', 'above': 'max'}
    mplExtend = extendAliases.get(extend, extend)
    if mplExtend not in ('neither', 'both', 'min', 'max'):
        raise ValueError(
            "extend must be one of 'neither', 'both', 'below'/'min' or "
            "'above'/'max'; received %r" % (extend,))

    if qnty_mn is None:
        qnty_mn = np.nanmin(data_slc)
    if qnty_mx is None:
        qnty_mx = np.nanmax(data_slc)
    if qnty_mn == qnty_mx:
        qnty_mn = 0
        qnty_mx = 1
    if highlightValue is not None:
        percentile = (highlightValue - qnty_mn) / (qnty_mx - qnty_mn)
    if xmn is None: xmn = x.min()
    if xmx is None: xmx = x.max()
    if zmn is None: zmn = z.min()
    if zmx is None: zmx = z.max()
    if (cmap is None) or (isinstance(cmap, str) and cmap == 'SMV'):
        cmap = buildSMVcolormap(
                percentile=percentile, width=highlightWidth)
    # An integer here would tell contourf to choose that many levels
    # spanning the *data*, which ignores qnty_mn and qnty_mx: the
    # colorbar would then span the data rather than the requested range,
    # and any tick outside the data would be clamped onto the extension
    # arrow, printing on top of its neighbours. Build the levels
    # explicitly so that the default matches what an explicit
    # levels=100 already did.
    if levels is None:
        levels = np.linspace(qnty_mn, qnty_mx, 100)
    elif isinstance(levels, Iterable):
        levels = np.array(levels)
    else:
        levels = np.linspace(qnty_mn, qnty_mx, levels)
    if cbarticks is None:
        cbarticks = np.linspace(qnty_mn, qnty_mx, cbarnumticks)
        if highlightValue is not None:
            cx = cbarticks[1] - cbarticks[0]
            cbarticks = list(cbarticks)
            cbarticks = [x for x in cbarticks if abs(highlightValue-x) > 0.25*cx]
            cbarticks.append(highlightValue)
            cbarticks = sorted(cbarticks)

    # Drop ticks outside the color scale. matplotlib does not discard
    # them, it clamps them to the end of the bar, so a caller passing a
    # tick list wider than the scale gets a stack of overlapping labels
    # on the extension arrow rather than a missing tick.
    tol = 1e-9 * max(1.0, abs(qnty_mx - qnty_mn))
    inRange = [t for t in np.asarray(cbarticks, dtype=float)
               if (t >= qnty_mn - tol) and (t <= qnty_mx + tol)]
    if len(inRange) < len(np.asarray(cbarticks).ravel()):
        dropped = len(np.asarray(cbarticks).ravel()) - len(inRange)
        print("Warning, %d colorbar tick(s) fall outside the color scale "
              "%0.4f to %0.4f and were dropped." % (dropped, qnty_mn, qnty_mx))
    cbarticks = inRange


    if reverseXY:
        (x1 , z1, d1) = (z, x, data_slc[:, :])
        zrange = xmx-xmn #x.max()-x.min()
        xrange = zmx-zmn #z.max()-z.min()
    else:
        (x1, z1, d1) = (x, z, data_slc[:, :])
        xrange = xmx-xmn #x.max()-x.min()
        zrange = zmx-zmn #z.max()-z.min()
    if figsize is None:
        if zrange > xrange:
            figsize = (figsizeMult, figsizeMult * zrange / xrange)
        else:
            figsize = (figsizeMult * xrange / zrange, figsizeMult)
    if (fig is None) or (ax is None):
        fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)
    if contour:
        if linecontour:
            im = ax.contour(
                x1, z1, d1, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx,
                levels=levels, extend=mplExtend)
        else:
            im = ax.contourf(
                x1, z1, d1, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx,
                levels=levels, extend=mplExtend)
    else:
        # Copy before masking: d1 is a view of the caller's array, and
        # writing NaN into it would corrupt the data they passed in.
        dout = d1[::-1, :].copy() # default behavior is extend above and below
        if mplExtend == 'min':
            dout[dout > qnty_mx] = np.nan
        elif mplExtend == 'max':
            dout[dout < qnty_mn] = np.nan
        im = ax.imshow(dout, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx, extent=[x1.min(), x1.max(),z1.min(), z1.max()])
    if addCbar:
        if tickDecimals is not None:
            fmt = "%"+".%d"%(tickDecimals)+"f"
            #fmt = FormatScalarFormatter("%"+fmtStr+"f")
        else:
            fmt = None
        cbar = fig.colorbar(
            im, cmap=cmap, extend=mplExtend, ticks=cbarticks, format=fmt)
        cbar.ax.tick_params(labelsize=fs)
        if clabel is not None:
            cbar.set_label(clabel, fontsize=fs)
    if xlabel is None:
        if abs(axis) == 1: xlabel = 'y (m)'
        if abs(axis) == 2: xlabel = 'x (m)'
        if abs(axis) == 3: xlabel = 'x (m)'
    if zlabel is None:
        if abs(axis) == 1: zlabel = 'z (m)'
        if abs(axis) == 2: zlabel = 'z (m)'
        if abs(axis) == 3: zlabel = 'y (m)'
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(zlabel, fontsize=fs)
    if fixXLims:
        ax.set_xlim(xmn, xmx)
    if fixZLims:
        ax.set_ylim(zmn, zmx)
    if title is not None:
        ax.set_title(title, fontsize=fs)
    ax.tick_params(labelsize=fs)
    #fig.tight_layout()
    if returnIm:
        return fig, ax, im
    else:
        return fig, ax

def _plot3dTimeFromName(name):
    """Recovers the timestamp encoded in a plot3D file name

    FDS names plot3D output files '<chid>_<mesh>_<time>.q', where the
    decimal point of the time is written as the letter p, for example
    'case001_1_30p02.q' for t = 30.02 s.

    Parameters
    ----------
    name : str
        Plot3D file name or path

    Returns
    -------
    float
        Timestamp encoded in the name
    """

    stem = os.path.basename(name).split('.q')[0]
    timeStr = stem.split('_')[-1]
    if 'p' in timeStr:
        whole, fraction = timeStr.split('p')
        return round(float(whole) + float(fraction)/(10.0**len(fraction)), 4)
    return float(timeStr)


def readPlot3Ddata(chid, resultDir, time, verbose=False):
    """Reads plot3D data for every mesh and maps it onto one grid

    Parameters
    ----------
    chid : str
        FDS CHID of the case
    resultDir : str
        Directory containing the FDS results, or a zip archive
    time : float
        Query time. The nearest plot3D output time is used
    verbose : bool, optional
        Print the extents of each mesh as it is read (default False)

    Returns
    -------
    array(NX, NY, NZ, 3)
        Absolute grid coordinates spanning every mesh
    array(NX, NY, NZ, 5)
        Plot3D data on the absolute grid. The last axis holds
        temperature, u, v, w and HRRPUV

    Raises
    ------
    FileNotFoundError
        If the case contains no xyz or plot3D files, or if a mesh has no
        plot3D file at any time
    """

    xyzFiles = getFileList(resultDir, chid, 'xyz')
    tFiles = getFileList(resultDir, chid, 'q')
    if len(xyzFiles) == 0:
        raise FileNotFoundError(
            "No xyz file for chid %s was found in %s. Add "
            "WRITE_XYZ=.TRUE. to the &DUMP namelist to have FDS write "
            "one." % (chid, resultDir))
    if len(tFiles) == 0:
        raise FileNotFoundError(
            "No plot3D (.q) file for chid %s was found in %s."
            % (chid, resultDir))

    grids = dict()
    datas = dict()

    for xyzFile in xyzFiles:
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s_"%(chid, mesh)

        # Select the plot3D file belonging to this mesh, then the time
        # nearest the one requested. Reading a single file for every
        # mesh, as earlier releases did, silently returned the wrong
        # mesh's data for multi-mesh cases.
        meshFiles = [x for x in tFiles
                     if os.path.basename(x).startswith(
                         '%s_%s_'%(chid, mesh) if mesh != '' else chid)]
        if len(meshFiles) == 0:
            raise FileNotFoundError(
                "No plot3D (.q) file was found for mesh %s of chid %s "
                "in %s." % (mesh, chid, resultDir))
        meshTimes = np.array([_plot3dTimeFromName(x) for x in meshFiles])
        dataFile = meshFiles[int(np.argmin(abs(meshTimes - time)))]

        grid, gridHeader = readXYZfile(xyzFile)
        data, dataHeader = readP3Dfile(dataFile)
        (nx, ny, nz) = (dataHeader[0], dataHeader[1], dataHeader[2])

        if verbose:
            printExtents(grid, data)

        xGrid, yGrid, zGrid = rearrangeGrid(grid)

        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid

        datas[meshStr] = np.reshape(data, (nx, ny, nz, 5), order='F')

    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    data_abs = np.full((xGrid_abs.shape[0],
                        xGrid_abs.shape[1],
                        xGrid_abs.shape[2],
                        5), np.nan)

    for key in list(grids.keys()):
        (xGrid, yGrid, zGrid) = (grids[key]['xGrid'],
                                 grids[key]['yGrid'],
                                 grids[key]['zGrid'])

        xloc = np.where(abs(xGrid_abs-xGrid[0,0,0]) == 0)[0][0]
        yloc = np.where(abs(yGrid_abs-yGrid[0,0,0]) == 0)[1][0]
        zloc = np.where(abs(zGrid_abs-zGrid[0,0,0]) == 0)[2][0]

        (NX, NY, NZ) = np.shape(xGrid)

        data_abs[xloc:xloc+NX, yloc:yloc+NY, zloc:zloc+NZ, :] = datas[key]

    return grid_abs, data_abs


def extractResultDirAndChidFromSlcfName(slcfFile):
    """Recovers the result directory and CHID from a slice file path

    Parameters
    ----------
    slcfFile : str
        Path to a slice file. FDS names these
        '<chid>_<mesh>_<quantity index>.sf'

    Returns
    -------
    str
        Directory containing the slice file
    str
        FDS CHID of the case
    """

    resultDir = os.sep.join(os.path.abspath(slcfFile).split(os.sep)[:-1])
    chid = '_'.join(os.path.abspath(slcfFile).split(os.sep)[-1].split('_')[:-2])
    return resultDir, chid

def _readSlcfFrames(f, timesSLCF, NX, NY, NZ, datatype, time, dt,
                    headerSize=142):
    """Reads the requested frames from an open slice file

    The three supported queries are handled here so that every slice
    reader treats them identically:

    * time is None                -> every frame in the file
    * time given, dt is None      -> the single frame nearest time
    * time and dt given           -> the mean of the frames whose
                                     timestamps fall within
                                     time +/- dt/2

    Parameters
    ----------
    f : file
        Binary file positioned immediately after the slice header
    timesSLCF : array(NT)
        Timestamps of every frame in the file
    NX : int
        Number of cells along the x-axis of the slice
    NY : int
        Number of cells along the y-axis of the slice
    NZ : int
        Number of cells along the z-axis of the slice
    datatype : numpy.dtype
        Datatype including byte order used to decode the frames
    time : float or None
        Query time
    dt : float or None
        Averaging window centred on time
    headerSize : int, optional
        Size of the slice file header in bytes (default 142)

    Returns
    -------
    array(NX+1, NY+1, NZ+1, NT)
        Array containing the requested frames
    list or array
        Timestamps corresponding to the returned frames
    """

    shape = (NX+1, NY+1, NZ+1)
    frameBytes = 4 * (5 + (NX+1) * (NY+1) * (NZ+1))

    if time is None:
        NT = len(timesSLCF)
        datas2 = np.zeros((NX+1, NY+1, NZ+1, NT), dtype=np.float32)
        for i in range(0, NT):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            if data is False:
                # The file is shorter than its time record implies.
                datas2 = datas2[:, :, :, :i]
                return datas2, timesSLCF[:i]
            datas2[:, :, :, i] = np.reshape(data, shape, order='F')
        return datas2, timesSLCF

    if dt is None:
        datas2 = np.zeros((NX+1, NY+1, NZ+1, 1), dtype=np.float32)
        i = int(np.argmin(abs(timesSLCF - time)))
        f.seek(i * frameBytes + headerSize, 0)
        t, data = readNextTime(f, NX, NY, NZ, datatype)
        datas2[:, :, :, 0] = np.reshape(data, shape, order='F')
        return datas2, [timesSLCF[i]]

    # Averaging window. Select by timestamp rather than by stepping a
    # fixed number of frames so that a non-uniform slice output
    # interval is handled correctly, and divide by the number of frames
    # actually accumulated.
    t1 = time - dt/2
    t2 = time + dt/2
    inds = np.where(np.logical_and(timesSLCF >= t1, timesSLCF <= t2))[0]
    if inds.size == 0:
        inds = np.array([int(np.argmin(abs(timesSLCF - time)))])

    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1), dtype=np.float64)
    ts = []
    for ind in inds:
        f.seek(int(ind) * frameBytes + headerSize, 0)
        t, data = readNextTime(f, NX, NY, NZ, datatype)
        if data is False:
            break
        datas2[:, :, :, 0] += np.reshape(data, shape, order='F')
        ts.append(float(np.squeeze(t)))
    if len(ts) > 0:
        datas2[:, :, :, 0] = datas2[:, :, :, 0] / float(len(ts))
        times = [(np.min(ts) + np.max(ts)) / 2]
    else:
        times = [time]
    return np.array(datas2, dtype=np.float32), times


def readSingleSlcfFile(slcfFile,
                       time=None, dt=None, saveTimesFile=False,
                       endianness="<"):
    """Reads the data from a single slice file

    Handles both 2-D and 3-D slices; the returned array always has the
    shape (NX+1, NY+1, NZ+1, NT), where the singleton axis of a 2-D
    slice is retained.

    Parameters
    ----------
    slcfFile : str
        Path to a slice file, or to a slice file inside a zip archive
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window centred on time. Ignored when time is omitted
    saveTimesFile : bool, optional
        If True the timestamps are cached alongside the slice file so
        that subsequent reads do not have to rescan it (default False)
    endianness : str, optional
        Byte order of the file, '<' or '>' (default '<')

    Returns
    -------
    list
        Six component list [iX, eX, iY, eY, iZ, eZ] of the slice extents
        in cell indices
    array(NX+1, NY+1, NZ+1, NT)
        Array containing the slice data
    list or array
        Timestamps corresponding to the returned data
    """

    if saveTimesFile:
        timesFile = slcfFile.replace('.sf','_times.csv')
    else:
        timesFile = None
    datatype = getDatatypeByEndianness(np.float32, endianness)
    timesSLCF = readSLCFtimes(slcfFile, timesFile, endianness=endianness)
    f = zopen(slcfFile)

    qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    datas2, times = _readSlcfFrames(
        f, timesSLCF, NX, NY, NZ, datatype, time, dt)
    f.close()
    lims = [iX, eX, iY, eY, iZ, eZ]
    return lims, datas2, times


def readSLCF3Ddata(chid, workingDir, quantityToExport,
                   time=None, dt=None, saveTimesFile=False, verbose=False,
                   axis=None, value=None):
    """Reads 3-D slice data for a quantity and maps it onto one grid

    Every mesh which wrote a 3-D slice of quantityToExport is read and
    interpolated onto the absolute grid spanning all meshes in the case.

    Parameters
    ----------
    chid : str
        FDS CHID of the case
    workingDir : str
        Directory containing the FDS results, or a zip archive
    quantityToExport : str
        FDS quantity to read, for example 'TEMPERATURE'
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window centred on time
    saveTimesFile : bool, optional
        Cache the slice timestamps next to each slice file
    verbose : bool, optional
        Print progress while reading each mesh
    axis : int, optional
        If given together with value, extract only the 2-D plane normal
        to this axis (1 = x, 2 = y, 3 = z), which avoids assembling the
        full 3-D array
    value : float, optional
        Coordinate of the extracted plane along axis

    Returns
    -------
    array
        Absolute grid coordinates. Shape (NX, NY, NZ, 3) for a full 3-D
        read, or (N, M, 2) when axis and value select a plane
    array
        Slice data on the absolute grid, with time as the last axis
    array
        Timestamps of the returned data
    str
        Units of the quantity as recorded in the slice file

    Notes
    -----
    Returns ``(False, False, False, False)`` when the case contains no
    3-D slice of the requested quantity.
    """

    endianness = getEndianness(workingDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    smvOutputs = parseSMVFile(getSmvFile(workingDir, chid))
    
    smv_grids = smvOutputs['grids']

    resultDir = _resolveResultDir(workingDir, chid)

    grids = defaultdict(bool)
    foundSomething = False
    tinds = []
    outputUnits = None

    for i in range(0, len(smv_grids)):
        if verbose: print("Starting grid %d"%(i+1))
        xs = smv_grids[i][0][:, 1]
        ys = smv_grids[i][1][:, 1]
        zs = smv_grids[i][2][:, 1]
        xGrid, yGrid, zGrid = np.meshgrid(xs, ys, zs)
        xGrid = np.swapaxes(xGrid, 0, 1)
        yGrid = np.swapaxes(yGrid, 0, 1)
        zGrid = np.swapaxes(zGrid, 0, 1)
        
        meshStr = "%s"%(chid) if len(smv_grids) == 1 else "%s_%d_"%(chid, i+1)
        #print(i, meshStr)
        if '.zip' in resultDir:
            slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s%s*.sf"%(resultDir, os.sep, meshStr))
        #print(i, "%s%s%s*.sf"%(resultDir, os.sep, meshStr), slcfFiles)
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        if (axis is not None) and (value is not None):
            skip = False
            if (axis == 1) and ((value < xs[0]) or (value > xs[-1])): skip = True
            if (axis == 2) and ((value < ys[0]) or (value > ys[-1])): skip = True
            if (axis == 3) and ((value < zs[0]) or (value > zs[-1])): skip = True
            if skip and verbose: print("Skipping mesh %d as selected plane (%d, %0.4f) not contained in mesh %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f"%(i+1, axis, value, xs[0], xs[-1], ys[0], ys[-1], zs[0], zs[-1]))
            if skip: continue
        grids[meshStr]['include'] = True
        times3D = []
        datas3D = []
        lims3D = []
        
        for slcfFile in slcfFiles:
            if saveTimesFile:
                timesFile = slcfFile.replace('.sf','_times.csv')
            else:
                timesFile = None
            timesSLCF = readSLCFtimes(
                slcfFile, timesFile, endianness=endianness)
            f = zopen(slcfFile)

            qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(
                f, endianness)
            correctQuantity = (qty == quantityToExport)
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            if correctQuantity and threeDimSlice:
                if verbose: print("Reading slice %s"%(slcfFile))
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                datas2, times = _readSlcfFrames(
                    f, timesSLCF, NX, NY, NZ, datatype, time, dt)
                lims3D.append([iX, eX, iY, eY, iZ, eZ])
                datas3D.append(datas2)
                times3D.append(times)
                outputUnits = uts
            f.close()

        if len(datas3D) == 0:
            # This mesh holds no 3-D slice of the requested quantity.
            # Mark it excluded so that the assembly loop below skips it
            # instead of indexing an empty list.
            if verbose:
                print("No 3-D slice of %s found for mesh %s"
                      % (quantityToExport, meshStr))
            grids[meshStr]['include'] = False
            continue

        foundSomething = True
        grids[meshStr]['datas3D'] = datas3D
        grids[meshStr]['lims3D'] = lims3D
        tinds.append(datas3D[0].shape[3])
    
    if not foundSomething:
        print("No 3D slice data found for qty %s"%(quantityToExport))
        print("Quantities available in this case:")
        readSLCFquantities(chid, workingDir, printInfo=True)
        return False, False, False, False

    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    tInd = int(np.nanmin(tinds))
    xshp, yshp, zshp = xGrid_abs.shape
    if (axis is not None) and (value is not None):
        if axis == 1: xshp = 1
        if axis == 2: yshp = 1
        if axis == 3: zshp = 1
    data_abs = np.zeros((xshp, yshp, zshp, tInd))
    data_abs[:, :, :, :] = np.nan
    
    abs_xs = xGrid_abs[:, 0, 0]
    abs_ys = yGrid_abs[0, :, 0]
    abs_zs = zGrid_abs[0, 0, :]
    
    for key in list(grids.keys()):
        x1 = grids[key]['xGrid'].flatten()
        y1 = grids[key]['yGrid'].flatten()
        z1 = grids[key]['zGrid'].flatten()
        
        if not grids[key]['include']: continue
        
        xi = np.where(np.logical_or(np.isclose(abs_xs, x1.min()), np.isclose(abs_xs,x1.max())))[0]
        yi = np.where(np.logical_or(np.isclose(abs_ys, y1.min()), np.isclose(abs_ys,y1.max())))[0]
        zi = np.where(np.logical_or(np.isclose(abs_zs, z1.min()), np.isclose(abs_zs,z1.max())))[0]
        
        xg, yg, zg = np.meshgrid(abs_xs[xi[0]:xi[1]+1], abs_ys[yi[0]:yi[1]+1], abs_zs[zi[0]:zi[1]+1])
        xg = np.swapaxes(xg, 0, 1)
        yg = np.swapaxes(yg, 0, 1)
        zg = np.swapaxes(zg, 0, 1)
        
        p = np.zeros((xg.flatten().shape[0],3))
        p[:, 0] = xg.flatten()
        p[:, 1] = yg.flatten()
        p[:, 2] = zg.flatten()
        
        xs = grids[key]['xGrid'][:, 0, 0]
        ys = grids[key]['yGrid'][0, :, 0]
        zs = grids[key]['zGrid'][0, 0, :]        
        d = np.zeros_like(data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, zi[0]:zi[1]+1, 0])
        for t in range(0, tInd):
            d[:, :, :] = np.nan
            if t < grids[key]['datas3D'][0].shape[3]:
                lims3D = grids[key]['lims3D']
                iX, eX, iY, eY, iZ, eZ = lims3D[0]
                #d[iX:eX+1,iY:eY+1,iZ:eZ+1] = grids[key]['datas3D'][0][:, :, :, t]
                d = grids[key]['datas3D'][0][:, :, :, t] #hybrid meshing
                #print(xs.shape, ys.shape, zs.shape, grids[key]['datas3D'][0][:, :, :, t].shape)
                #print(grids[key]['lims3D'])
                #d2 = scpi.interpn((xs, ys, zs), grids[key]['datas3D'][0][:, :, :, t], p, method='linear', fill_value=None, bounds_error=False)
                #print(d2r.shape)
                #print(data_abs.shape, data_abs[xi[0]:xi[1]+1, 0, zi[0]:zi[1]+1, t].shape)
                d2 = scpi.interpn((xs, ys, zs), d, p, method='linear', fill_value=None, bounds_error=False)
                d2r = np.reshape(d2, xg.shape)
                if (axis is not None) and (value is not None):
                    if axis == 1:
                        xind = np.argmin(abs(xs-value))
                        #d2 = scpi.interpn((ys, zs), d[xind, :, :], p[:,1:], method='linear', fill_value=None, bounds_error=False)
                        #d2r = np.reshape(d2, xg[0, :, :].shape)
                        data_abs[0, yi[0]:yi[1]+1, zi[0]:zi[1]+1, t] = d2r[xind, :, :]
                    if axis == 2: 
                        yind = np.argmin(abs(ys-value))
                        #d2 = scpi.interpn((xs, zs), d[:, yind, :], p[:,[0, 2]], method='linear', fill_value=None, bounds_error=False)
                        #print(d2.shape, xg.shape)
                        #d2r = np.reshape(d2, xg[:, 0, :].shape)
                        #print(d2r.shape, data_abs[0, yi[0]:yi[1]+1, zi[0]:zi[1]+1, t].shape)
                        data_abs[xi[0]:xi[1]+1, 0, zi[0]:zi[1]+1, t] = d2r[:, yind, :]
                    if axis == 3: 
                        zind = np.argmin(abs(zs-value))
                        #d2 = scpi.interpn((xs, ys), d[:, :, zind], p[:,:-1], method='linear', fill_value=None, bounds_error=False)
                        #d2r = np.reshape(d2, xg[:, :, 0].shape)
                        data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, 0, t] = d2r[:, :, zind]
                else:
                    data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, zi[0]:zi[1]+1, t] = d2r
    if (axis is not None) and (value is not None):
        if axis == 1: 
            grid_out = grid_abs[0, :, :, :][:, :, [1, 2]]
            data_out = data_abs[0, :, :, :]
        if axis == 2: 
            grid_out = grid_abs[:, 0, :, :][:, :, [0, 2]]
            data_out = data_abs[:, 0, :, :]
        if axis == 3: 
            grid_out = grid_abs[:, :, 0, :][:, :, [0, 1]]
            data_out = data_abs[:, :, 0, :]
    else:
        grid_out = grid_abs
        data_out = data_abs
    return grid_out, data_out, times3D[0], outputUnits


def readSLCF3DdataXYZ(chid, resultDir, quantityToExport,
                      time=None, dt=None, saveTimesFile=False):
    """Reads 3-D slice data using xyz files to build the mesh grids

    .. deprecated::
        Use :func:`readSLCF3Ddata`, which reads the mesh grids from the
        smokeview file and therefore does not require the case to have
        been run with WRITE_XYZ=.TRUE. on the &DUMP namelist. This
        routine is retained for backwards compatibility.

    Parameters
    ----------
    chid : str
        FDS CHID of the case
    resultDir : str
        Directory containing the FDS results, or a zip archive
    quantityToExport : str
        FDS quantity to read, for example 'TEMPERATURE'
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window centred on time
    saveTimesFile : bool, optional
        Cache the slice timestamps next to each slice file

    Returns
    -------
    array(NX, NY, NZ, 3)
        Absolute grid coordinates
    array(NX, NY, NZ, NT)
        Slice data on the absolute grid
    array(NT)
        Timestamps of the returned data
    """

    warnings.warn(
        "readSLCF3DdataXYZ is deprecated and will be removed in a future "
        "release; use readSLCF3Ddata instead.",
        DeprecationWarning, stacklevel=2)
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    if '.zip' in resultDir:
        xyzFiles = getFileListFromZip(resultDir, chid, 'xyz')
    else:
        xyzFiles = glob.glob("%s%s%s*.xyz"%(resultDir,os.sep, chid))
    #print(xyzFiles)
    grids = defaultdict(bool)
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid)
        
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s_"%(chid, mesh)
        if '.zip' in resultDir:
            slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s%s*.sf"%(resultDir, os.sep, meshStr))
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        times3D = []
        datas3D = []
        lims3D = []
        #print("%s%s_*.sf"%(resultDir, meshStr))
        for slcfFile in slcfFiles:
            if saveTimesFile:
                timesFile = slcfFile.replace('.sf','_times.csv')
            else:
                timesFile = None
            timesSLCF = readSLCFtimes(slcfFile, timesFile, endianness=endianness)
            times = []
            f = zopen(slcfFile)

            qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
            # Check if slice is correct quantity
            correctQuantity = (qty == quantityToExport)
            # Check if slice is 2-dimensional
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            #print(qty, quantityToExport)
            if correctQuantity and threeDimSlice:
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                # Check if slice is 3-D
                #print(slcfFile, qty, sName, uts, iX, eX, iY, eY, iZ, eZ)
                shape = (NX+1, NY+1, NZ+1)
                if time == None:
                    NT = len(timesSLCF)
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, NT))
                    for i in range(0, NT):
                        t, data = readNextTime(f, NX, NY, NZ, datatype)
                        data = np.reshape(data, shape, order='F')
                        datas2[:, :, :, i] = data
                    times = timesSLCF
                elif (time != None) and (dt == None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF-time))
                    f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
                    t, data = readNextTime(f, NX, NY, NZ, datatype)
                    data = np.reshape(data, shape, order='F')
                    datas2[:, :, :, 0] = data
                    times = [timesSLCF[i]]
                elif (time != None) and (dt != None):
                    datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
                    i = np.argmin(abs(timesSLCF - (time - dt/2)))
                    j = np.argmin(abs(timesSLCF - (time + dt/2)))
                    f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
                    data = True
                    for ii in range(i, j+1):
                        if data is not False:
                            t, data = readNextTime(f, NX, NY, NZ, datatype)
                            data = np.reshape(data, shape, order='F')
                            datas2[:, :, :, 0] += data
                    if j - i > 0:
                        datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
                    times = [timesSLCF[i]]
                lims3D.append([iX, eX, iY, eY, iZ, eZ])
                datas3D.append(datas2)
                times3D.append(times)
                outputUnits = uts
            f.close()
        grids[meshStr]['datas3D'] = datas3D
        grids[meshStr]['lims3D'] = lims3D
    
    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    if len(grids[list(grids.keys())[0]]['datas3D']) == 0:
        print("No 3D slice data found for qty %s"%(quantityToExport))
        return False, False, False, False
    tInd = grids[list(grids.keys())[0]]['datas3D'][0].shape[3]
    data_abs = np.zeros((xGrid_abs.shape[0],
                         xGrid_abs.shape[1],
                         xGrid_abs.shape[2],
                         tInd))
    data_abs[:, :, :, :] = np.nan
    
    abs_xs = xGrid_abs[:, 0, 0]
    abs_ys = yGrid_abs[0, :, 0]
    abs_zs = zGrid_abs[0, 0, :]
    
    for key in list(grids.keys()):
        x1 = grids[key]['xGrid'].flatten()
        y1 = grids[key]['yGrid'].flatten()
        z1 = grids[key]['zGrid'].flatten()
        
        xi = np.where(np.logical_or(np.isclose(abs_xs, x1.min()), np.isclose(abs_xs,x1.max())))[0]
        yi = np.where(np.logical_or(np.isclose(abs_ys, y1.min()), np.isclose(abs_ys,y1.max())))[0]
        zi = np.where(np.logical_or(np.isclose(abs_zs, z1.min()), np.isclose(abs_zs,z1.max())))[0]
        
        xg, yg, zg = np.meshgrid(abs_xs[xi[0]:xi[1]+1], abs_ys[yi[0]:yi[1]+1], abs_zs[zi[0]:zi[1]+1])
        xg = np.swapaxes(xg, 0, 1)
        yg = np.swapaxes(yg, 0, 1)
        zg = np.swapaxes(zg, 0, 1)
        
        p = np.zeros((xg.flatten().shape[0],3))
        p[:, 0] = xg.flatten()
        p[:, 1] = yg.flatten()
        p[:, 2] = zg.flatten()
        
        xs = grids[key]['xGrid'][:, 0, 0]
        ys = grids[key]['yGrid'][0, :, 0]
        zs = grids[key]['zGrid'][0, 0, :]        
        d = np.zeros_like(data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, zi[0]:zi[1]+1, 0])
        for t in range(0, tInd):
            d[:, :, :] = np.nan
            if t < grids[key]['datas3D'][0].shape[3]:
                lims3D = grids[key]['lims3D']
                iX, eX, iY, eY, iZ, eZ = lims3D[0]
                #d[iX:eX+1,iY:eY+1,iZ:eZ+1] = grids[key]['datas3D'][0][:, :, :, t]
                d = grids[key]['datas3D'][0][:, :, :, t] #hybrid meshing
                #print(xs.shape, ys.shape, zs.shape, grids[key]['datas3D'][0][:, :, :, t].shape)
                #print(grids[key]['lims3D'])
                #d2 = scpi.interpn((xs, ys, zs), grids[key]['datas3D'][0][:, :, :, t], p, method='linear', fill_value=None, bounds_error=False)
                d2 = scpi.interpn((xs, ys, zs), d, p, method='linear', fill_value=None, bounds_error=False)
                d2r = np.reshape(d2, xg.shape)
                data_abs[xi[0]:xi[1]+1, yi[0]:yi[1]+1, zi[0]:zi[1]+1, t] = d2r
    
    
    
    
    '''
    for key in list(grids.keys()):
        xGrid = grids[key]['xGrid']
        yGrid = grids[key]['yGrid']
        zGrid = grids[key]['zGrid']
        datas3D = grids[key]['datas3D']
        lims3D = grids[key]['lims3D']
        print("Starting Grid %s"%(key))
        for data, lim, times in zip(datas3D, lims3D, times3D):
            xloc_mn = np.where(np.isclose(
                    abs(xGrid_abs - xGrid[lim[0], 0, 0]),
                    0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(
                    abs(xGrid_abs - xGrid[lim[1], 0, 0]),
                    0, atol=1e-04))[0][0]
            yloc_mn = np.where(np.isclose(
                    abs(yGrid_abs - yGrid[0, lim[2], 0]),
                    0, atol=1e-04))[1][0]
            yloc_mx = np.where(np.isclose(
                    abs(yGrid_abs - yGrid[0, lim[3], 0]),
                    0, atol=1e-04))[1][0]
            zloc_mn = np.where(np.isclose(
                    abs(zGrid_abs - zGrid[0, 0, lim[4]]),
                    0, atol=1e-04))[2][0]
            zloc_mx = np.where(np.isclose(
                    abs(zGrid_abs - zGrid[0, 0, lim[5]]),
                    0, atol=1e-04))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            ANX = xloc_mx-xloc_mn + 1
            ANY = yloc_mx-yloc_mn + 1
            ANZ = zloc_mx-zloc_mn + 1
            if (NX != ANX) or (NY != ANY) or (NZ != ANZ):
                x = xGrid[lim[0]:lim[1]+1, 0, 0]
                y = yGrid[0, lim[2]:lim[3]+1, 0]
                z = zGrid[0, 0, lim[4]:lim[5]+1]
                
                xi = grid_abs[xloc_mn:xloc_mx+1, 
                              yloc_mn:yloc_mx+1, 
                              zloc_mn:zloc_mx+1,
                              0].flatten()
                yi = grid_abs[xloc_mn:xloc_mx+1,
                              yloc_mn:yloc_mx+1,
                              zloc_mn:zloc_mx+1,
                              1].flatten()
                zi = grid_abs[xloc_mn:xloc_mx+1,
                              yloc_mn:yloc_mx+1,
                              zloc_mn:zloc_mx+1,
                              2].flatten()
                
                x = np.round(x, decimals=4)
                y = np.round(y, decimals=4)
                z = np.round(z, decimals=4)
                
                xi = np.round(xi, decimals=4)
                yi = np.round(yi, decimals=4)
                zi = np.round(zi, decimals=4)
                
                xi[xi < np.min(x)] = np.min(x)
                xi[xi > np.max(x)] = np.max(x)
                yi[yi < np.min(y)] = np.min(y)
                yi[yi > np.max(y)] = np.max(y)
                zi[zi < np.min(z)] = np.min(z)
                zi[zi > np.max(z)] = np.max(z)
                
                tmpGrid = np.array([xi, yi, zi]).T
                for i in range(0, NT):
                    interpolator = scpi.RegularGridInterpolator(
                            (x, y, z), data[:, :, :, i])
                    data2 = interpolator(tmpGrid)
                    data2 = np.reshape(
                            data2, (ANX, ANY, ANZ), order='C')
                    try:
                        data_abs[xloc_mn:xloc_mx+1,
                                 yloc_mn:yloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading mesh %s at time %0.0f"%(key, i))
            else:
                try:
                    data_abs[xloc_mn:xloc_mx+1,
                             yloc_mn:yloc_mx+1,
                             zloc_mn:zloc_mx+1,
                             :] = data
                except:
                    try:
                        NTT = min([data_abs.shape[3], data.shape[3]])
                        data_abs[xloc_mn:xloc_mx+1,
                                 yloc_mn:yloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 :NTT] = data[:, :, :, :NTT]
                        print("Error loading mesh %s at time %0.0f"%(key, i))
                    except:
                        print("Error loading mesh %s at all times"%(key))
                        print(data.shape, data_abs[xloc_mn:xloc_mx+1, yloc_mn:yloc_mx+1, zloc_mn:zloc_mx+1, :].shape, NTT)
    '''
    return grid_abs, data_abs, times3D[0], outputUnits





def readSLCF2Ddata(chid, resultDir, quantityToExport,
                   time=None, dt=None):
    """Reads every 2-D slice of a quantity onto the absolute grid

    .. deprecated::
        Use :func:`query2dAxisValue`, which selects a single plane by
        axis and coordinate, reads the grid from the smokeview file
        rather than requiring xyz files, and handles cell-centered
        slices. This routine is retained for backwards compatibility.

    Parameters
    ----------
    chid : str
        FDS CHID of the case
    resultDir : str
        Directory containing the FDS results, or a zip archive
    quantityToExport : str
        FDS quantity to read, for example 'TEMPERATURE'
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window centred on time

    Returns
    -------
    array(NX, NY, NZ, 3)
        Absolute grid coordinates
    array(NX, NY, NZ, NT)
        Slice data mapped onto the absolute grid
    array(NT)
        Timestamps of the returned data
    """

    warnings.warn(
        "readSLCF2Ddata is deprecated and will be removed in a future "
        "release; use query2dAxisValue instead.",
        DeprecationWarning, stacklevel=2)
    slcfDir = _resolveResultDir(resultDir, chid)
    if '.zip' in resultDir:
        xyzFiles = getFileListFromZip(slcfDir, chid, 'xyz')
    else:
        xyzFiles = glob.glob("%s%s%s*.xyz"%(slcfDir, os.sep, chid))
    grids = defaultdict(bool)
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid)
        
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
        if '.zip' in slcfDir:
            slcfFiles = getFileListFromZip(slcfDir, chid, 'sf')
            slcfFiles = [x for x in slcfFiles if '%s'%(meshStr) in x]
        else:
            slcfFiles = glob.glob("%s%s%s_*.sf"%(slcfDir, os.sep, meshStr))
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
        
        datas2D = []
        lims2D = []
        coords2D = []
        times2D = []
        timesOut = np.zeros((0,))
        for slcfFile in slcfFiles:
            timesSLCF = readSLCFtimes(slcfFile, None, endianness)
            times = []
            f = zopen(slcfFile)
            
            qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
            # Check if slice is correct quantity
            correctQuantity = (qty == quantityToExport)
            # Check if slice is 2-dimensional
            threeDimSlice = (eX-iX > 0) and (eY-iY > 0) and (eZ-iZ > 0)
            if correctQuantity and not threeDimSlice:
                (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
                datas2, times = _readSlcfFrames(
                    f, timesSLCF, NX, NY, NZ, datatype, time, dt)
                lims2D.append([iX, eX, iY, eY, iZ, eZ])
                datas2D.append(datas2)
                coords2D.append([xGrid[iX, iY, iZ],
                                 yGrid[iX, iY, iZ],
                                 zGrid[iX, iY, iZ]])
                times2D.append(np.array(times))
                timesOut = np.array(times)
            f.close()
        grids[meshStr]['datas2D'] = datas2D
        grids[meshStr]['lims2D'] = lims2D
        grids[meshStr]['coords2D'] = coords2D
        
    grid_abs = getAbsoluteGrid(grids)
    xGrid_abs = grid_abs[:, :, :, 0]
    yGrid_abs = grid_abs[:, :, :, 1]
    zGrid_abs = grid_abs[:, :, :, 2]
    tInds = []
    for key in list(grids.keys()):
        datas2D = grids[key]['datas2D']
        if (datas2D is not False) and (len(datas2D) > 0):
            tInds.append(datas2D[0].shape[3])
    if len(tInds) == 0:
        raise ValueError(
            "No 2-D slice of %s was found in %s."
            % (quantityToExport, resultDir))
    tInd = int(np.min(tInds))
    data_abs = np.zeros((xGrid_abs.shape[0],
                         xGrid_abs.shape[1],
                         xGrid_abs.shape[2],
                         tInd))
    data_abs[:, :, :, :] = np.nan
    for key in list(grids.keys()):
        xGrid = grids[key]['xGrid']
        yGrid = grids[key]['yGrid']
        zGrid = grids[key]['zGrid']
        datas2D = grids[key]['datas2D']
        lims2D = grids[key]['lims2D']
        coords2D = grids[key]['coords2D']
        
        for data, coord in zip(datas2D, coords2D):
            xloc = np.where(np.isclose(
                    abs(xGrid_abs - coord[0]), 0, atol=1e-06))[0][0]
            yloc = np.where(np.isclose(
                    abs(yGrid_abs - coord[1]), 0, atol=1e-06))[1][0]
            zloc = np.where(np.isclose(
                    abs(zGrid_abs - coord[2]), 0, atol=1e-06))[2][0]
            (NX, NY, NZ, NT) = np.shape(data)
            NT = min([NT, data_abs.shape[3]])
            data_abs[xloc:xloc+NX,
                     yloc:yloc+NY,
                     zloc:zloc+NZ,
                     :NT] = data[:, :, :, :NT]
    return grid_abs, data_abs, timesOut







def extractPoint(point, grid, data):
    """Returns the data at the grid point nearest a coordinate

    Parameters
    ----------
    point : array-like
        Three component [x, y, z] coordinate to query
    grid : array(NX, NY, NZ, 3)
        Absolute grid coordinates
    data : array(NX, NY, NZ, N)
        Data on the absolute grid

    Returns
    -------
    array(N)
        Data at the nearest grid point. A warning is printed when the
        nearest point is more than 0.25 m away in total from the query
    """

    ind = np.argmin(np.sum(abs(grid-point),axis=3).flatten())
    ind = np.unravel_index(ind, grid[:,:,:,0].shape)
    x = grid[ind[0],ind[1],ind[2],0]
    y = grid[ind[0],ind[1],ind[2],1]
    z = grid[ind[0],ind[1],ind[2],2]
    d = data[ind[0],ind[1],ind[2],:]
    err = abs(point[0]-x)+abs(point[1]-y)+abs(point[2]-z)
    if err > 0.25:
        errPt = "%0.4f, %0.4f, %0.4f, is %0.4f"%(x, y, z, err)
        print("Warning, error for point %s."%(errPt))
    return d

def readNextTime(f, NX, NY, NZ, datatype):
    """Reads one timestep record from an open slice file

    Parameters
    ----------
    f : file
        Binary file positioned at the start of a timestep record
    NX : int
        Number of cells along the x-axis of the slice
    NY : int
        Number of cells along the y-axis of the slice
    NZ : int
        Number of cells along the z-axis of the slice
    datatype : numpy.dtype
        Float datatype including byte order

    Returns
    -------
    array(1)
        Timestamp of the frame
    array or bool
        Flat array of the frame values, or False if the record could not
        be read because the file ended early
    """

    _ = np.frombuffer(f.read(8), dtype=datatype)
    time = np.frombuffer(f.read(4), dtype=datatype)
    _ = np.frombuffer(f.read(8), dtype=datatype)
    try:
        data = np.frombuffer(f.read((NX+1)*(NY+1)*(NZ+1)*4), 
                             dtype=datatype)
    except:
        data = False
    return time, data

def readSLCFheader(f, endianness, byteSize=False):
    """Reads the 142 byte header of a slice file

    Parameters
    ----------
    f : file
        Binary file positioned at the start of the file
    endianness : str
        Byte order of the file, '<' or '>'
    byteSize : bool, optional
        Return the extents as a single six component tuple rather than
        as six separate values (default False)

    Returns
    -------
    str
        FDS quantity recorded in the file
    str
        Short name of the quantity
    str
        Units of the quantity
    tuple or six ints
        Slice extents [iX, eX, iY, eY, iZ, eZ] in cell indices, as one
        tuple when byteSize is True and as six values otherwise
    """

    data = f.read(142)
    header = data[:110]
    size = struct.unpack('%siiiiii'%(endianness), data[118:142])
    tmp = header.split(b'\x1e')
    quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
    shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
    units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')
    
    if byteSize:
        return quantity, shortName, units, size
    else:
        iX, eX, iY, eY, iZ, eZ = size
        return quantity, shortName, units, iX, eX, iY, eY, iZ, eZ

def readSLCFtimes(file, timesFile=None, endianness=None):
    """Reads the timestamps of every frame in a slice file

    Parameters
    ----------
    file : str
        Path to a slice file, or to a slice file inside a zip archive
    timesFile : str, optional
        Path of a csv cache of the timestamps. It is read when it exists
        and written after a scan when it does not
    endianness : str, optional
        Byte order of the file, '<' or '>'. Determined from the case's
        .end file when omitted

    Returns
    -------
    array(NT)
        Array of timestamps
    """

    if timesFile is not None:
        if os.path.exists(timesFile):
            times = np.loadtxt(timesFile, delimiter=',')
            return times
    if endianness is None:
        resultDir, chid = extractResultDirAndChidFromSlcfName(file)
        endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)

    f = zopen(file)
    qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    headerSize = 142
    frameFloats = (NX+1)*(NY+1)*(NZ+1) + 5
    frameBytes = 4 * frameFloats

    fileSize = None
    if '.zip' not in file:
        try:
            fileSize = os.path.getsize(file)
        except OSError:
            fileSize = None

    if fileSize is not None:
        # Seek to each frame's time record instead of reading the whole
        # file. A slice file is frequently several gigabytes, and only
        # four bytes per frame are needed here.
        numberOfFrames = int((fileSize - headerSize) // frameBytes)
        times = np.zeros(numberOfFrames, dtype=np.float32)
        for i in range(0, numberOfFrames):
            f.seek(headerSize + i*frameBytes + 8, 0)
            times[i] = np.frombuffer(f.read(4), dtype=datatype)[0]
    else:
        # Members of a zip archive are decompressed as a stream, so
        # seeking past data is no cheaper than reading it.
        data = f.read()
        remainder = len(data) % 4
        if remainder != 0:
            data = data[:-remainder]
        fullFile = np.frombuffer(data, dtype=datatype)
        times = np.array(fullFile[2::frameFloats])
    f.close()

    if timesFile is not None:
        np.savetxt(timesFile, times)
    return times

def readSLCFquantities(chid, workingDir, printInfo=False):
    """Lists the slice files written by a case and what they contain

    Parameters
    ----------
    chid : str
        FDS CHID of the case
    workingDir : str
        Directory containing the FDS results, or a zip archive
    printInfo : bool, optional
        Print each slice file as it is inspected (default False)

    Returns
    -------
    list
        Quantity recorded in each slice file
    list
        Path of each slice file
    list
        Six component extent [iX, eX, iY, eY, iZ, eZ] of each slice
    list
        Mesh string parsed from each slice file name
    list
        Whether each slice holds cell-centered data
    list
        Units of the quantity in each slice file
    """

    resultDir = _resolveResultDir(workingDir, chid)

    smvData = parseSMVFile(getSmvFile(workingDir, chid))
    
    '''
    try:
        smvFile = os.path.join(resultDir, '%s.smv'%(chid))
        smvData = parseSMVFile(smvFile)
    except:
        smvDir = os.sep.join(os.path.abspath(resultDir).split(os.sep)[:-1])
        smvFile = os.path.join(smvDir, '%s.smv'%(chid))
        smvData = parseSMVFile(smvFile)
    '''
    files = smvData['files']
    files2 = defaultdict(bool)
    files2['SLICES'] = defaultdict(bool)
    for key in list(files['SLICES'].keys()):
        n = os.path.basename(key)
        files2['SLICES'][n] = files['SLICES'][key]
    smvData['files'] = files2
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    
    if '.zip' in resultDir:
        slcfFiles = getFileListFromZip(resultDir, chid, 'sf')
        zip = zipfile.ZipFile(resultDir, 'r')
    else:
        slcfFiles = glob.glob("%s%s_*.sf"%(resultDir, chid))
    endianness = getEndianness(resultDir, chid)
    quantities = []
    dimensions = []
    meshes = []
    centers = []
    units = []
    all_files_from_smv = list(files['SLICES'].keys())
    #print(all_files_from_smv)
    #assert False, "Stopped"
    for file in slcfFiles:
        if '.zip' in resultDir:
            f = zip.open(file.split("%s%s"%('.zip',os.sep))[1])
        else:
            f = open(file, 'rb')
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
        quantities.append(qty)
        dimensions.append([iX, eX, iY, eY, iZ, eZ])
        n = file.split(chid)[-1].split('_')
        meshStr = n[-2]
        meshes.append(meshStr)
        #print(files['SLICES'][file.split(os.sep)[-1]])
        ftocheck = file.split(os.sep)[-1]
        if ftocheck not in all_files_from_smv:
            print("WARNING file %s not in smokeview slices"%(ftocheck))
            continue
        if printInfo:
            print(file)
            print(files['SLICES'])
            print(files['SLICES'][ftocheck])
        centers.append(files['SLICES'][file.split(os.sep)[-1]]['CELL_CENTERED'])
        units.append(uts)
        f.close()
    if '.zip' in resultDir:
        zip.close()
    return quantities, slcfFiles, dimensions, meshes, centers, units

def buildQTYstring(chid, resultDir, qty):
    """Finds the file name suffix FDS used for a slice quantity

    Parameters
    ----------
    chid : str
        FDS CHID of the case
    resultDir : str
        Directory containing the FDS results, or a zip archive
    qty : str
        FDS quantity, for example 'TEMPERATURE'

    Returns
    -------
    str
        Suffix used in the slice file names for this quantity

    Raises
    ------
    ValueError
        If the case contains no slice of the requested quantity
    """

    quantities, slcfFiles, dimensions, meshes, centers, units = \
        readSLCFquantities(chid, resultDir)
    inds = np.where([qty == x for x in quantities])[0]
    if inds.size == 0:
        raise ValueError(
            "Quantity %s not found in %s. Known quantities: %s"
            % (qty, resultDir, ', '.join(sorted(set(quantities)))))
    ind = int(inds[0])
    quantityStr = slcfFiles[ind].split('.sf')[0].split('_')[-1]
    return quantityStr

def getLimsFromGrid(grid):
    """Returns the bounding box of an absolute grid

    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Absolute grid coordinates

    Returns
    -------
    list
        Six component list [xmin, xmax, ymin, ymax, zmin, zmax]

    See Also
    --------
    pyfdstools.extractBoundaryData.getPatchLimsFromGrid : the unrelated
        routine which converts boundary patch cell indices to
        coordinates, and which carried this same name before v0.0.24
    """

    xGrid = grid[:, :, :, 0]
    yGrid = grid[:, :, :, 1]
    zGrid = grid[:, :, :, 2]
    
    (xmn, xmx) = (xGrid[:, 0, 0].min(), xGrid[:, 0, 0].max())
    (ymn, ymx) = (yGrid[0, :, 0].min(), yGrid[0, :, 0].max())
    (zmn, zmx) = (zGrid[0, 0, :].min(), zGrid[0, 0, :].max())
    
    return [xmn, xmx, ymn, ymx, zmn, zmx]

def visualizePlot3D(x, z, T, U, V, W, HRR,
                    qnty_mn=None, qnty_mx=None):
    """Plots the temperature field of a plot3D slice

    Parameters
    ----------
    x : array(N, M)
        First in-plane coordinate of each point
    z : array(N, M)
        Second in-plane coordinate of each point
    T : array(N, M)
        Temperature values, which are the values plotted
    U : array(N, M)
        Velocity component along x, accepted for symmetry with
        findSliceLocation and not currently plotted
    V : array(N, M)
        Velocity component along y
    W : array(N, M)
        Velocity component along z
    HRR : array(N, M)
        Heat release rate per unit volume
    qnty_mn : float, optional
        Lower limit of the color scale. Data minimum when omitted
    qnty_mx : float, optional
        Upper limit of the color scale. Data maximum when omitted

    See Also
    --------
    plotSlice : the general purpose slice plotting routine, which offers
        control over labels, limits, colormaps and output
    """

    cmap = buildSMVcolormap()
    
    xrange = x.max()-x.min()
    zrange = z.max()-z.min()
    
    if zrange > xrange:
        plt.figure(figsize=(12*xrange/zrange,12))
    else:
        plt.figure(figsize=(12,12*zrange/xrange))
    if qnty_mn is None:
        qnty_mn = np.nanmin(T)
    if qnty_mx is None:
        qnty_mx = np.nanmax(T)
    levels = np.linspace(qnty_mn, qnty_mx, 100)
    plt.contourf(x, z, T, cmap=cmap, vmin=qnty_mn, vmax=qnty_mx,
                 levels=levels, extend='both')
    plt.colorbar()


def read2dSliceFile(slcfFile, chid, time=None, dt=None, cen=False, grid=None):
    """Reads a single 2-D slice file and returns it with its coordinates

    Parameters
    ----------
    slcfFile : str
        Path to a slice file, or to one inside a zip archive
    chid : str
        FDS CHID of the case
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window. Combined with time it returns a single
        averaged frame; on its own it applies a running average of this
        width to every frame
    cen : bool, optional
        Shift the coordinates to cell centers, for a slice written with
        CELL_CENTERED=.TRUE. (default False)
    grid : dict, optional
        Mesh grid with the keys 'xGrid', 'yGrid' and 'zGrid'. Read from
        the mesh's xyz file when omitted

    Returns
    -------
    array(N, M)
        First in-plane coordinate of each point
    array(N, M)
        Second in-plane coordinate of each point
    array(N, M, NT)
        Slice values
    list or array
        Timestamps of the returned frames
    list
        Six component list of the slice bounds in coordinates

    Notes
    -----
    Returns None if the file holds a 3-D rather than a 2-D slice.
    """

    resultDir = os.sep.join(os.path.abspath(slcfFile).split(os.sep)[:-1])
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    timesSLCF = readSLCFtimes(slcfFile, None, endianness=endianness)
    times = []
    
    if grid == None:
        xyzFile = '%s%s'%('_'.join(slcfFile.split('_')[:-1]), '.xyz')
        grids = getGridsFromXyzFiles([xyzFile], chid)
        grid = grids[list(grids.keys())[0]]
    xGrid = grid['xGrid']
    yGrid = grid['yGrid']
    zGrid = grid['zGrid']
    
    if cen:
        dx = np.round((grid['xGrid'][-1, 0, 0] - grid['xGrid'][0, 0, 0]) / (grid['xGrid'][:, 0, 0].shape[0]-1), decimals=4)
        dy = np.round((grid['yGrid'][0, -1, 0] - grid['yGrid'][0, 0, 0]) / (grid['yGrid'][0, :, 0].shape[0]-1), decimals=4)
        dz = np.round((grid['zGrid'][0, 0, -1] - grid['zGrid'][0, 0, 0]) / (grid['zGrid'][0, 0, :].shape[0]-1), decimals=4)
    else:
        dx = 0
        dy = 0
        dz = 0
    f = zopen(slcfFile)
    
    qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
    headerSize = 142
    
    (NX, NY, NZ) = (eX - iX, eY - iY, eZ - iZ)
    if (NX == 0):
        slcf_axis = 1
    elif (NY == 0):
        slcf_axis = 2
    elif (NZ == 0):
        slcf_axis = 3
    else:
        slcf_axis = -1
    
    if slcf_axis < 0:
        print("Not a 2-D slice.")
        return None
    shape = (NX+1, NY+1, NZ+1)
    if (time == None) and (dt == None):
        NT = len(timesSLCF)
        datas2 = np.zeros((NX+1, NY+1, NZ+1, NT), dtype=np.float32)
        for i in range(0, NT):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, i] = data
        times = timesSLCF
    elif (time == None) and (dt != None):
        NT = len(timesSLCF)
        datas3 = np.zeros((NX+1, NY+1, NZ+1, NT))
        for i in range(0, NT):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas3[:, :, :, i] = data
        # Running mean over the window [t-dt/2, t+dt/2] for every
        # frame. timesSLCF is monotonic, so both window bounds advance
        # monotonically and the means can be formed from a prefix sum:
        # one pass over the data instead of re-averaging a slab of
        # frames for every timestep.
        firstInds = np.searchsorted(timesSLCF, timesSLCF - dt/2, side='left')
        lastInds = np.searchsorted(timesSLCF, timesSLCF + dt/2, side='right')
        counts = (lastInds - firstInds).astype(np.float64)
        counts[counts < 1] = 1.0
        if np.isnan(datas3).any():
            # A prefix sum would smear a single NaN across every
            # subsequent frame, so fall back to the direct nanmean.
            datas2 = np.zeros_like(datas3)
            for i in range(0, NT):
                datas2[:, :, :, i] = np.nanmean(
                    datas3[:, :, :, firstInds[i]:lastInds[i]], axis=3)
        else:
            cumulative = np.zeros(
                (NX+1, NY+1, NZ+1, NT+1), dtype=np.float64)
            np.cumsum(datas3, axis=3, out=cumulative[:, :, :, 1:])
            datas2 = (cumulative[:, :, :, lastInds]
                      - cumulative[:, :, :, firstInds]) / counts
        times = timesSLCF
    elif (time != None) and (dt == None):
        datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
        i = np.argmin(abs(timesSLCF-time))
        f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
        t, data = readNextTime(f, NX, NY, NZ, datatype)
        data = np.reshape(data, shape, order='F')
        datas2[:, :, :, 0] = data
        times = [timesSLCF[i]]
        NT = 1
    elif (time != None) and (dt != None):
        datas2 = np.zeros((NX+1, NY+1, NZ+1, 1))
        
        t1 = max([0, time-dt/2])
        t2 = min([timesSLCF[-1], time+dt/2])
        inds = np.where(np.logical_and(timesSLCF >= t1, timesSLCF <= t2))[0]
        
        #print(timesSLCF, t1, t2, timesSLCF[inds])
        ts = []
        for ind in inds:
            f.seek(ind * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)) + headerSize, 0)
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, 0] += data
            ts.append(t)
        datas2[:, :, :, 0] = datas2[:, :, :, 0] / float(len(inds))
        times = [(np.min(ts) + np.max(ts))/2]
        '''
        i = np.argmin(abs(timesSLCF - (time - dt/2)))
        j = np.argmin(abs(timesSLCF - (time + dt/2)))
        f.seek(i * 4 * (5 + (NX+1) * (NY+1) * (NZ+1)), 1)
        for ii in range(i, j+1):
            t, data = readNextTime(f, NX, NY, NZ, datatype)
            data = np.reshape(data, shape, order='F')
            datas2[:, :, :, 0] += data
        if j - i > 0:
            datas2[:, :, :, 0] = datas2[:, :, :, 0] / (j-i)
        times = [(timesSLCF[i] + timesSLCF[j])/2]
        
        print(i, timesSLCF[i], j, timesSLCF[j])
        '''
        NT = 1
    coords = [xGrid[iX, iY, iZ], xGrid[eX, eY, eZ],
              yGrid[iX, iY, iZ], yGrid[eX, eY, eZ],
              zGrid[iX, iY, iZ], zGrid[eX, eY, eZ]]
    f.close()
    if slcf_axis == 1:
        x = yGrid[iX, iY:eY+1, iZ:eZ+1]-dy/2
        z = zGrid[iX, iY:eY+1, iZ:eZ+1]-dz/2
        d = datas2[0, :, :, :]
    elif slcf_axis == 2:
        x = xGrid[iX:eX+1, iY, iZ:eZ+1]-dx/2
        z = zGrid[iX:eX+1, iY, iZ:eZ+1]-dz/2
        d = datas2[:, 0, :, :]
    elif slcf_axis == 3:
        x = xGrid[iX:eX+1, iY:eY+1, iZ]-dx/2
        z = yGrid[iX:eX+1, iY:eY+1, iZ]-dy/2
        d = datas2[:, :, 0, :]
    
    return x, z, d, times, coords

def getAxisAndValueFromXB(XB, grid, cen):
    """Determines which plane a slice lies in from its cell extents

    Parameters
    ----------
    XB : list
        Six component slice extent [iX, eX, iY, eY, iZ, eZ] in cell
        indices
    grid : dict
        Mesh grid with the keys 'xGrid', 'yGrid' and 'zGrid'
    cen : bool
        Whether the slice holds cell-centered data, in which case the
        coordinate is shifted by half a cell

    Returns
    -------
    int
        Axis normal to the slice (1 = x, 2 = y, 3 = z), or -1 when the
        extents describe a 3-D slice
    float
        Coordinate of the slice along that axis, or -1 for a 3-D slice
    """

    NX = XB[1] - XB[0]
    NY = XB[3] - XB[2]
    NZ = XB[5] - XB[4]
    
    dx = np.round((grid['xGrid'][-1, 0, 0] - grid['xGrid'][0, 0, 0]) / (grid['xGrid'][:, 0, 0].shape[0]-1), decimals=4)
    dy = np.round((grid['yGrid'][0, -1, 0] - grid['yGrid'][0, 0, 0]) / (grid['yGrid'][0, :, 0].shape[0]-1), decimals=4)
    dz = np.round((grid['zGrid'][0, 0, -1] - grid['zGrid'][0, 0, 0]) / (grid['zGrid'][0, 0, :].shape[0]-1), decimals=4)
    if (NX == 0):
        axis = 1
        if cen:
            value = grid['xGrid'][XB[0], 0, 0] - dx/2
        else:
            value = grid['xGrid'][XB[0], 0, 0]
    elif (NY == 0):
        axis = 2
        if cen:
            value = grid['yGrid'][0, XB[2], 0] - dy/2
        else:
            value = grid['yGrid'][0, XB[2], 0]
    elif (NZ == 0):
        axis = 3
        if cen:
            value = grid['zGrid'][0, 0, XB[4]] - dz/2
        else:
            value = grid['zGrid'][0, 0, XB[4]]
    else:
        axis = -1
        value = -1
    return axis, value

def query2dAxisValue(workingDir, chid, quantity, axis, value, time=None,
                     dt=None, atol=1e-8, printInfo=False, verbose=False):
    """Reads a 2-D slice at a given axis and coordinate

    Every mesh which wrote a 2-D slice of the requested quantity in the
    requested plane is read and mapped onto the absolute grid spanning
    all meshes in the case. Mesh grids are read from the smokeview file,
    so the case does not need to have been run with WRITE_XYZ=.TRUE.

    Parameters
    ----------
    workingDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case
    quantity : str
        FDS quantity to read, for example 'TEMPERATURE'
    axis : int
        Axis normal to the queried slice (1 = x, 2 = y, 3 = z)
    value : float
        Coordinate of the queried slice along axis
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window centred on time. When time is omitted, a
        running average of this width is applied to every frame
    atol : float, optional
        Tolerance used when matching value against the available slice
        coordinates (default 1e-8)
    printInfo : bool, optional
        Print each slice file as it is read (default False)
    verbose : bool, optional
        Print progress while reading each mesh (default False)

    Returns
    -------
    defaultdict
        Dictionary with the keys:
        'x'     - array(N, M) of the first in-plane coordinate
        'z'     - array(N, M) of the second in-plane coordinate
        'datas' - array(N, M, NT) of slice values
        'times' - timestamps of each frame
    str
        Units of the quantity as recorded in the slice file

    Notes
    -----
    Returns ``(None, None)`` and prints the slices which are available
    when the requested plane does not exist in the case.

    Examples
    --------
    >>> data, units = query2dAxisValue(
    ...     workingDir, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    >>> data['datas'].shape[-1]
    1
    """

    endianness = getEndianness(workingDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    smvOutputs = parseSMVFile(getSmvFile(workingDir, chid))
    
    smv_grids = smvOutputs['grids']
    smv_slcf = smvOutputs['files']['SLICES']
    
    resultDir = _resolveResultDir(workingDir, chid)

    grids = defaultdict(bool)
    for i in range(0, len(smv_grids)):
        if verbose: print("Starting grid %d"%(i+1))
        xs = smv_grids[i][0][:, 1]
        ys = smv_grids[i][1][:, 1]
        zs = smv_grids[i][2][:, 1]
        xGrid, yGrid, zGrid = np.meshgrid(xs, ys, zs)
        xGrid = np.swapaxes(xGrid, 0, 1)
        yGrid = np.swapaxes(yGrid, 0, 1)
        zGrid = np.swapaxes(zGrid, 0, 1)
        
        meshStr = "%s"%(chid) if len(smv_grids) == 1 else "%s_%d_"%(chid, i+1)
        #print(i, "%s%s%s*.sf"%(resultDir, os.sep, meshStr), slcfFiles)
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
    meshStrIsChid = False
    if len(smv_grids) == 1: meshStrIsChid = True
    
    grids_abs = getAbsoluteGrid(grids)
    if abs(axis) == 1:
        xAbs = grids_abs[0, :, :, 1]
        zAbs = grids_abs[0, :, :, 2]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    elif abs(axis) == 2:
        xAbs = grids_abs[:, 0, :, 0]
        zAbs = grids_abs[:, 0, :, 2]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    elif abs(axis) == 3:
        xAbs = grids_abs[:, :, 0, 0]
        zAbs = grids_abs[:, :, 0, 1]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    quantities, slcfFiles, dimensions, meshes, centers, units = readSLCFquantities(chid, workingDir, printInfo=printInfo)
    if printInfo:
        print(quantities)
        print(slcfFiles)
    
    datas = defaultdict(bool)
    foundSlice = False
    available_slices = []
    for qty, slcfFile, dim, cen, uts in zip(quantities, slcfFiles, dimensions, centers, units):
        if qty == quantity:
            n = slcfFile.split(chid)[-1].split('_')
            meshStr = n[-2]
            if meshStr == '': meshStr = chid
            tmp = [len(x) for x in list(grids.keys())]
            if np.max(tmp) < 4:
                meshStr = str(int(meshStr))
            if meshStrIsChid:
                meshStr = chid
            else:
                meshStr = chid + "_%s_"%(meshStr)
            #mesh = slcfFile.split(chid)[-1].split('.sf')[0].split('_')[-2]
            #meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
            #print(slcfFile, n, meshStr, grids.keys())
            slcf_axis, slcf_value = getAxisAndValueFromXB(dim, grids[meshStr], cen)
            
            available_slices.append([slcf_axis, slcf_value])
            if np.isclose(axis, abs(slcf_axis)) and np.isclose(slcf_value, value, atol=atol):
                if printInfo:
                    print("Reading %s"%(slcfFile), time, dt)
                x, z, d, times, coords = read2dSliceFile(slcfFile, chid, time=time, dt=dt, grid=grids[meshStr])
                if cen:
                    x = (x[1:, :]+x[:-1, :])/2
                    z = (z[:,1:]+z[:,:-1])/2
                    x = x[:,:-1]
                    z = z[:-1,:]
                    d = d[1:,1:]
                #print(x.shape, z.shape, d.shape)
                slcfName = '%s_%0.4f_%0.4f_%0.4f_%0.4f_%0.4f_%0.4f'%(qty, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5])
                datas[slcfName] = defaultdict(bool)
                datas[slcfName]['limits'] = coords
                datas[slcfName]['times'] = times
                datas[slcfName]['datas'] = d.copy()
                datas[slcfName]['x'] = x.copy()
                datas[slcfName]['z'] = z.copy()
                datas[slcfName]['center'] = cen
                foundSlice = True
                outUnits = uts
    if not foundSlice:
        print("Warning did not find a 2-D slice of %s on axis %d at value %0.4f"%(quantity, axis, value))
        print("Available slices of qty %s:"%(quantity))
        print('\tAxis\tValue')
        for slcf_axis, slcf_value in available_slices:
            print("\t%d\t%0.4f"%(slcf_axis, slcf_value))
        print("Note an axis of -1 indicates a 3-D slice which is not currently supported in this function.")
        return None, None        
    data_abs = np.zeros((xAbs.shape[0], xAbs.shape[1], len(times)), dtype=np.float32)
    data_abs_c = np.zeros((xAbs_c.shape[0], xAbs_c.shape[1], len(times)), dtype=np.float32)
    for slcfName in list(datas.keys()):
        x = datas[slcfName]['x']
        z = datas[slcfName]['z']
        d = datas[slcfName]['datas']
        c = datas[slcfName]['center']
        
        if c:
            xloc_mn = np.where(np.isclose(abs(xAbs_c - np.nanmin(x)), 0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xAbs_c - np.nanmax(x)), 0, atol=1e-04))[0][0]
            zloc_mn = np.where(np.isclose(abs(zAbs_c - np.nanmin(z)), 0, atol=1e-04))[1][0]
            zloc_mx = np.where(np.isclose(abs(zAbs_c - np.nanmax(z)), 0, atol=1e-04))[1][0]
        else:
            xloc_mn = np.where(np.isclose(abs(xAbs - np.nanmin(x)), 0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xAbs - np.nanmax(x)), 0, atol=1e-04))[0][0]
            zloc_mn = np.where(np.isclose(abs(zAbs - np.nanmin(z)), 0, atol=1e-04))[1][0]
            zloc_mx = np.where(np.isclose(abs(zAbs - np.nanmax(z)), 0, atol=1e-04))[1][0]
            
        (NX, NZ, NT) = np.shape(d)
        ANX = xloc_mx-xloc_mn + 1
        ANZ = zloc_mx-zloc_mn + 1
        
        if (NX != ANX) or (NZ != ANZ):
            if c:
                xi = xAbs_c[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
                zi = zAbs_c[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
            else:
                xi = xAbs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
                zi = zAbs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
            
            x = np.round(x, decimals=4)
            z = np.round(z, decimals=4)
            
            xi = np.round(xi, decimals=4)
            zi = np.round(zi, decimals=4)
            
            xi[xi < np.min(x)] = np.min(x)
            xi[xi > np.max(x)] = np.max(x)
            zi[zi < np.min(z)] = np.min(z)
            zi[zi > np.max(z)] = np.max(z)
            
            tmpGrid = np.array([xi, zi]).T
            for i in range(0, NT):
                interpolator = scpi.RegularGridInterpolator((x[:, 0], z[0, :]), d[:, :, i])
                data2 = interpolator(tmpGrid)
                data2 = np.reshape(data2, (ANX, ANZ), order='C')
                if c:
                    try:
                        data_abs_c[xloc_mn:xloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading %s at time %0.0f"%(slcfName, i))
                else:
                    try:
                        data_abs[xloc_mn:xloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading %s at time %0.0f"%(slcfName, i))
        else:
            try:
                data_abs[xloc_mn:xloc_mx+1,
                         zloc_mn:zloc_mx+1,
                         :] = d
            except:
                try:
                    NTT = min([data_abs.shape[2], d.shape[2]])
                    data_abs[xloc_mn:xloc_mx+1,
                             zloc_mn:zloc_mx+1,
                             :NTT] = d[:, :, :NTT]
                    print("Error loading %s at time %0.0f"%(slcfName, i))
                except:
                    print("Error loading %s at all times"%(slcfName))
                    print(d.shape, data_abs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1, :].shape, NTT)
                    
    data_abs_out = defaultdict(bool)
    data_abs_out['x'] = xAbs
    data_abs_out['z'] = zAbs
    data_abs_out['datas'] = data_abs
    data_abs_out['times'] = times
    return data_abs_out, outUnits

def query2dAxisValueXYZ(resultDir, chid, quantity, axis, value, time=None,
                        dt=None, atol=1e-8, printInfo=False):
    """Reads a 2-D slice using xyz files to build the mesh grids

    .. deprecated::
        Use :func:`query2dAxisValue`, which reads the mesh grids from
        the smokeview file and therefore does not require the case to
        have been run with WRITE_XYZ=.TRUE. on the &DUMP namelist. This
        routine is retained for backwards compatibility.

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case
    quantity : str
        FDS quantity to read, for example 'TEMPERATURE'
    axis : int
        Axis normal to the queried slice (1 = x, 2 = y, 3 = z)
    value : float
        Coordinate of the queried slice along axis
    time : float, optional
        Query time. Every frame is returned when omitted
    dt : float, optional
        Averaging window centred on time
    atol : float, optional
        Tolerance used when matching value against the available slice
        coordinates (default 1e-8)
    printInfo : bool, optional
        Print each slice file as it is inspected (default False)

    Returns
    -------
    defaultdict
        Dictionary with the keys 'x', 'z', 'datas' and 'times'
    str
        Units of the quantity
    """

    warnings.warn(
        "query2dAxisValueXYZ is deprecated and will be removed in a "
        "future release; use query2dAxisValue instead.",
        DeprecationWarning, stacklevel=2)
    xyzFiles = getFileListFromResultDir(resultDir, chid, 'xyz')
    if len(xyzFiles) == 0:
        xyzFiles = getFileListFromResultDir(resultDir+'*'+os.sep, chid, 'xyz')
        resultDir2 = os.path.dirname(xyzFiles[0]) + os.sep 
    else:
        resultDir2 = resultDir
    grids = getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = getAbsoluteGrid(grids)
    if abs(axis) == 1:
        xAbs = grids_abs[0, :, :, 1]
        zAbs = grids_abs[0, :, :, 2]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    elif abs(axis) == 2:
        xAbs = grids_abs[:, 0, :, 0]
        zAbs = grids_abs[:, 0, :, 2]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    elif abs(axis) == 3:
        xAbs = grids_abs[:, :, 0, 0]
        zAbs = grids_abs[:, :, 0, 1]
        xAbs_c = (xAbs[1:, :] + xAbs[:-1, :])/2
        xAbs_c = xAbs_c[:, :-1]
        zAbs_c = (zAbs[:,1:] + zAbs[:,:-1])/2
        zAbs_c = zAbs_c[:-1,:]
    quantities, slcfFiles, dimensions, meshes, centers, units = readSLCFquantities(chid, resultDir2, printInfo=printInfo)
    if printInfo:
        print(quantities)
        print(slcfFiles)
    
    datas = defaultdict(bool)
    foundSlice = False
    available_slices = []
    for qty, slcfFile, dim, cen, uts in zip(quantities, slcfFiles, dimensions, centers, units):
        if qty == quantity:
            n = slcfFile.split(chid)[-1].split('_')
            meshStr = n[-2]
            if meshStr == '': meshStr = chid
            tmp = [len(x) for x in list(grids.keys())]
            if np.max(tmp) < 4:
                meshStr = str(int(meshStr))
            #mesh = slcfFile.split(chid)[-1].split('.sf')[0].split('_')[-2]
            #meshStr = "%s"%(chid) if mesh == '' else "%s_%s"%(chid, mesh)
            
            slcf_axis, slcf_value = getAxisAndValueFromXB(dim, grids[meshStr], cen)
            
            available_slices.append([slcf_axis, slcf_value])
            if np.isclose(axis, abs(slcf_axis)) and np.isclose(slcf_value, value, atol=atol):
                if printInfo:
                    print("Reading %s"%(slcfFile), time, dt)
                x, z, d, times, coords = read2dSliceFile(slcfFile, chid, time=time, dt=dt)
                if cen:
                    x = (x[1:, :]+x[:-1, :])/2
                    z = (z[:,1:]+z[:,:-1])/2
                    x = x[:,:-1]
                    z = z[:-1,:]
                    d = d[1:,1:]
                #print(x.shape, z.shape, d.shape)
                slcfName = '%s_%0.4f_%0.4f_%0.4f_%0.4f_%0.4f_%0.4f'%(qty, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5])
                datas[slcfName] = defaultdict(bool)
                datas[slcfName]['limits'] = coords
                datas[slcfName]['times'] = times
                datas[slcfName]['datas'] = d.copy()
                datas[slcfName]['x'] = x.copy()
                datas[slcfName]['z'] = z.copy()
                datas[slcfName]['center'] = cen
                foundSlice = True
                outUnits = uts
    if not foundSlice:
        print("Warning did not find a 2-D slice of %s on axis %d at value %0.4f"%(quantity, axis, value))
        print("Available slices of qty %s:"%(quantity))
        print('\tAxis\tValue')
        for slcf_axis, slcf_value in available_slices:
            print("\t%d\t%0.4f"%(slcf_axis, slcf_value))
        print("Note an axis of -1 indicates a 3-D slice which is not currently supported in this function.")
        return None, None        
    data_abs = np.zeros((xAbs.shape[0], xAbs.shape[1], len(times)), dtype=np.float32)
    data_abs_c = np.zeros((xAbs_c.shape[0], xAbs_c.shape[1], len(times)), dtype=np.float32)
    for slcfName in list(datas.keys()):
        x = datas[slcfName]['x']
        z = datas[slcfName]['z']
        d = datas[slcfName]['datas']
        c = datas[slcfName]['center']
        
        if c:
            xloc_mn = np.where(np.isclose(abs(xAbs_c - np.nanmin(x)), 0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xAbs_c - np.nanmax(x)), 0, atol=1e-04))[0][0]
            zloc_mn = np.where(np.isclose(abs(zAbs_c - np.nanmin(z)), 0, atol=1e-04))[1][0]
            zloc_mx = np.where(np.isclose(abs(zAbs_c - np.nanmax(z)), 0, atol=1e-04))[1][0]
        else:
            xloc_mn = np.where(np.isclose(abs(xAbs - np.nanmin(x)), 0, atol=1e-04))[0][0]
            xloc_mx = np.where(np.isclose(abs(xAbs - np.nanmax(x)), 0, atol=1e-04))[0][0]
            zloc_mn = np.where(np.isclose(abs(zAbs - np.nanmin(z)), 0, atol=1e-04))[1][0]
            zloc_mx = np.where(np.isclose(abs(zAbs - np.nanmax(z)), 0, atol=1e-04))[1][0]
            
        (NX, NZ, NT) = np.shape(d)
        ANX = xloc_mx-xloc_mn + 1
        ANZ = zloc_mx-zloc_mn + 1
        
        if (NX != ANX) or (NZ != ANZ):
            if c:
                xi = xAbs_c[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
                zi = zAbs_c[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
            else:
                xi = xAbs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
                zi = zAbs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1].flatten()
            
            x = np.round(x, decimals=4)
            z = np.round(z, decimals=4)
            
            xi = np.round(xi, decimals=4)
            zi = np.round(zi, decimals=4)
            
            xi[xi < np.min(x)] = np.min(x)
            xi[xi > np.max(x)] = np.max(x)
            zi[zi < np.min(z)] = np.min(z)
            zi[zi > np.max(z)] = np.max(z)
            
            tmpGrid = np.array([xi, zi]).T
            for i in range(0, NT):
                interpolator = scpi.RegularGridInterpolator((x[:, 0], z[0, :]), d[:, :, i])
                data2 = interpolator(tmpGrid)
                data2 = np.reshape(data2, (ANX, ANZ), order='C')
                if c:
                    try:
                        data_abs_c[xloc_mn:xloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading %s at time %0.0f"%(slcfName, i))
                else:
                    try:
                        data_abs[xloc_mn:xloc_mx+1,
                                 zloc_mn:zloc_mx+1,
                                 i] = data2
                    except:
                        print("Error loading %s at time %0.0f"%(slcfName, i))
        else:
            try:
                data_abs[xloc_mn:xloc_mx+1,
                         zloc_mn:zloc_mx+1,
                         :] = d
            except:
                try:
                    NTT = min([data_abs.shape[2], d.shape[2]])
                    data_abs[xloc_mn:xloc_mx+1,
                             zloc_mn:zloc_mx+1,
                             :NTT] = d[:, :, :NTT]
                    print("Error loading %s at time %0.0f"%(slcfName, i))
                except:
                    print("Error loading %s at all times"%(slcfName))
                    print(d.shape, data_abs[xloc_mn:xloc_mx+1, zloc_mn:zloc_mx+1, :].shape, NTT)
                    
    data_abs_out = defaultdict(bool)
    data_abs_out['x'] = xAbs
    data_abs_out['z'] = zAbs
    data_abs_out['datas'] = data_abs
    data_abs_out['times'] = times
    return data_abs_out, outUnits

def renderSliceCsvs(data, chid, outdir):
    """Writes each frame of a 2-D slice to its own csv file

    Parameters
    ----------
    data : dict
        Slice dictionary as returned by query2dAxisValue, with the keys
        'x', 'z', 'datas' and 'times'
    chid : str
        FDS CHID of the case, used as the file name prefix
    outdir : str
        Directory the csv files are written to

    Notes
    -----
    Each file is named '<chid>_<time>.csv' and is written with the
    second in-plane coordinate as the row index and the first as the
    column headers.
    """

    times = data['times']
    xs = data['x'][:, 0]
    zs = data['z'][0, :]
    for i in range(0, len(times)):
        time = times[i]
        outFile = os.path.join(outdir, '%s_%0.4f.csv'%(chid, time))
        print("\t%d/%d: Rendering time %0.4f to file %s"%(i, len(times), time, outFile))
        d = pd.DataFrame(data['datas'][:, :, i].T, index=zs, columns=xs)
        d.to_csv(outFile)


def writeSLCFheader(f, quantity, shortName, units, size, endianness):
    """Writes the 142 byte header of a slice file

    Parameters
    ----------
    f : file
        Binary file open for writing
    quantity : str
        FDS quantity name
    shortName : str
        Short name of the quantity
    units : str
        Units of the quantity
    size : array(6)
        Slice extents [iX, eX, iY, eY, iZ, eZ] in cell indices
    endianness : str
        Byte order to write, '<' or '>'
    """

    sz = struct.pack('%s%0.0fi'%(endianness, len(size)), *size)
    qty = str.encode("{:<30}".format(quantity))
    sn = str.encode("{:<30}".format(shortName))
    un = str.encode("{:<30}".format(units))
    f.write(b'\x1e\x00\x00\x00')
    f.write(qty)
    f.write(b'\x1e\x00\x00\x00\x1e\x00\x00\x00')
    f.write(sn)
    f.write(b'\x1e\x00\x00\x00\x1e\x00\x00\x00')
    f.write(un)
    f.write(b'\x1e\x00\x00\x00\x18\x00\x00\x00')
    f.write(sz)
    f.write(b'\x18\x00\x00\x00')

def writeSLCFTime(f, time, data, endianness):
    """Writes one timestep record to a slice file

    Parameters
    ----------
    f : file
        Binary file positioned at the end of the previous record
    time : float
        Timestamp of the frame
    data : array
        Flat array of the frame values, in Fortran order
    endianness : str
        Byte order to write, '<' or '>'
    """

    f.write(b'\x04\x00\x00\x00')
    t = time.tobytes()
    f.write(t)
    f.write(b'\x04\x00\x00\x00')
    f.write(struct.pack('%si'%(endianness), data.shape[0]*4))
    d = data.tobytes()
    f.write(d)
    if data.shape[0]*4 <= 65535:
        f.write(struct.pack('%sHH'%(endianness), data.shape[0]*4,0))
    else:
        f.write(struct.pack('%sI'%(endianness), data.shape[0]*4))

def writeSlice(outFile, resultDir, chid, data, times, axis, val,
                       outQty, sName, uts, meshnum, smvFile=None, endianness="<", suffix=None):
    """Writes a derived 2-D slice as a slice file smokeview can read

    Used to add a computed quantity to an existing case, for example a
    radiative heat flux slice built from device output. Passing smvFile
    also registers the new slice in that smokeview file, without which
    smokeview will not display it.

    Parameters
    ----------
    outFile : str
        Name of the slice file to write, relative to resultDir
    resultDir : str
        Directory the slice file is written to
    chid : str
        FDS CHID of the case
    data : list
        List of NT frames, each an array(N, M) of values
    times : array(NT)
        Timestamps of each frame
    axis : int
        Axis normal to the slice (1 = x, 2 = y, 3 = z)
    val : int
        Cell index of the slice plane along axis
    outQty : str
        FDS quantity name to record
    sName : str
        Short name of the quantity
    uts : str
        Units of the quantity
    meshnum : int
        Mesh number the slice belongs to
    smvFile : str, optional
        Name of a smokeview file to append the slice record to
    endianness : str, optional
        Byte order to write, '<' or '>' (default '<')
    suffix : str, optional
        Suffix distinguishing this slice in the smokeview record
    """

    outPath = os.path.join(resultDir, outFile)
    smvPath = os.path.join(resultDir, smvFile)
    
    if axis == 1:
        (NX, NY, NZ) = (0, data[0].shape[0]-1, data[0].shape[1]-1)
        size = np.array([val, val, 0, NY, 0, NZ], dtype=np.int32)
        X = [val, val, 0, NY, 0, NZ]
    elif axis == 2:
        (NX, NY, NZ) = (data[0].shape[0]-1, 0, data[0].shape[1]-1)
        size = np.array([0, NX, val, val, 0, NZ], dtype=np.int32)
        X = [0, NX, val, val, 0, NZ]
    elif axis == 3:
        (NX, NY, NZ) = (data[0].shape[0]-1, data[0].shape[1]-1, 0)
        size = np.array([0, NX, 0, NY, val, val], dtype=np.int32)
        X = [0, NX, 0, NY, val, val]
        
    NT = len(times)
    f = zopen(outPath, 'wb')
    writeSLCFheader(f, outQty, sName, uts, size, endianness)
    shape2 = ((NX+1)*(NY+1)*(NZ+1),)
    for i in range(0, NT):
        data_out = np.reshape(data[i], shape2, order='F')
        writeSLCFTime(f, times[i], data_out, endianness)
    f.close()
    
    if smvFile is not None:
        writeSliceToSmv(smvPath, meshnum, X, outQty, outPath, sName, uts, suffix=suffix)

def writeSliceToSmv(file, meshNum, X, outQty, outFile, sName, uts, suffix):
    """Appends a slice record to a smokeview file

    Parameters
    ----------
    file : str
        Path to the smokeview file to append to
    meshNum : int
        Mesh number the slice belongs to
    X : list
        Six component slice extent in cell indices
    outQty : str
        FDS quantity name
    outFile : str
        Path of the slice file being registered
    sName : str
        Short name of the quantity
    uts : str
        Units of the quantity
    suffix : str
        Suffix distinguishing this slice from others of the same
        quantity
    """

    if suffix is None: suffix = "1 \n"
    with open(file, 'a') as f:
        f.write('SLCF     %0.0f # STRUCTURED &     %0.0f    %0.0f     %0.0f    %0.0f     %0.0f    %0.0f !      %s'%(meshNum, X[0], X[1], X[2], X[3], X[4], X[5], suffix))
        f.write(' %s\n'%(outFile))
        f.write(' %s\n'%(outQty))
        f.write(' %s\n'%(sName))
        f.write(' %s\n\n'%(uts))

def getAxisFromLims(lims):
    """Determines which plane a slice lies in from its cell extents

    Parameters
    ----------
    lims : list
        Six component slice extent [iX, eX, iY, eY, iZ, eZ] in cell
        indices

    Returns
    -------
    int
        Axis normal to the slice (1 = x, 2 = y, 3 = z), or -1 when the
        extents describe a 3-D slice
    int or None
        Cell index of the slice plane along that axis, or None for a
        3-D slice
    """

    iX, eX, iY, eY, iZ, eZ = lims
    (NX, NY, NZ) = (eX - iX+1, eY - iY+1, eZ - iZ+1)
    if (NX == 1):
        slcf_axis = 1
        val = eX
    elif (NY == 1):
        slcf_axis = 2
        val = eY
    elif (NZ == 1):
        slcf_axis = 3
        val = eZ
    else:
        slcf_axis = -1
        val = None
    return slcf_axis, val

def slcfTimeAverage(slcfFile, dt, outFile=None, outQty=None, outdt=None):
    """Time-averages one slice file and writes the result as a slice file

    Parameters
    ----------
    slcfFile : str
        Path to the slice file to average
    dt : float
        Averaging window in seconds
    outFile : str, optional
        Path the averaged slice file is written to. Derived from the
        input name when omitted
    outQty : str, optional
        Quantity name to record in the output. Derived from the input
        quantity and the window when omitted
    outdt : float, optional
        Output timestep. The input timestep is kept when omitted

    Returns
    -------
    str
        Path of the slice file which was written
    str
        Quantity name recorded in it
    """

    
    # Read the data
    resultDir, chid = extractResultDirAndChidFromSlcfName(slcfFile)
    endianness = getEndianness(resultDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    timesSLCF = readSLCFtimes(slcfFile, None, endianness=endianness)
    f = zopen(slcfFile)
    qty, sName, uts, size = readSLCFheader(f, endianness, byteSize=True)
    iX, eX, iY, eY, iZ, eZ = size
    (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
    NT = len(timesSLCF)
    datas2 = np.zeros((NX+1, NY+1, NZ+1, NT), dtype=np.float32)
    shape = (NX+1, NY+1, NZ+1)
    for i in range(0, NT):
        t, data = readNextTime(f, NX, NY, NZ, datatype)
        if data is not False:
            data2 = np.reshape(data, shape, order='F')
            datas2[:, :, :, i] = data2
    f.close()
    
    # Average the data
    if outdt is not None:
        outTimes = np.linspace(0, timesSLCF[-1], int((timesSLCF[-1]/outdt)+1), dtype=np.float32)
        
    else:
        outTimes = timesSLCF
    
    NT2 = outTimes.shape[0]
    data_avg = np.zeros((NX+1, NY+1, NZ+1, NT2), dtype=np.float32)
    
    for i in range(0, NT2):
        t = outTimes[i]
        j1 = np.argmin(abs(timesSLCF - t - dt/2))
        j0 = np.argmin(abs(timesSLCF - t + dt/2))
        #i0 = max([int(j - dt/dt2), 0])
        #i1 = min([int(j + dt/dt2), NT])
        #print(i, j0, j1, timesSLCF[j0], timesSLCF[j1])
        data_avg[:, :, :, i] = np.nanmean(datas2[:, :, :, j0:j1], axis=3)
    
    # Write the data
    if outFile is None:
        outFile = slcfFile.replace('.sf', '_avg.sf')
    if outQty is None:
        outQty = sName + ' (%0.0ds Avg)'%(dt)
    f = zopen(outFile, 'wb')
    writeSLCFheader(f, outQty, sName, uts, size)
    shape2 = ((NX+1)*(NY+1)*(NZ+1),)
    for i in range(0, NT2):
        data = np.reshape(data_avg[:,:,:,i], shape2, order='F')
        writeSLCFTime(f, outTimes[i], data)
    f.close()
    if outdt is None:
        print("Wrote %s with %0.1f average window"%(outFile, dt))
    else:
        print("Wrote %s with %0.1f average window at interval %0.1f"%(outFile, dt, outdt))
    header = qty, sName, uts, iX, eX, iY, eY, iZ, eZ
    return outFile, header
    

def slcfsTimeAverage(resultDir, chid, fdsQuantity, dt, outDir=None, outQty=None, outdt=None):
    """Time-averages every slice of a quantity across all meshes

    Writes one averaged slice file per input slice, plus a smokeview
    file which registers them, so that the averaged field can be opened
    in smokeview alongside the original.

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case
    fdsQuantity : str
        FDS quantity to average, for example 'TEMPERATURE'
    dt : float
        Averaging window in seconds
    outDir : str, optional
        Directory the averaged files are written to. Defaults to
        resultDir
    outQty : str, optional
        Quantity name to record in the output
    outdt : float, optional
        Output timestep. The input timestep is kept when omitted

    Returns
    -------
    list
        Paths of the averaged slice files
    str
        Quantity name recorded in them
    list
        Paths of the slice files which were averaged
    str
        Path of the smokeview file which registers the new slices
    """

    slcfFiles = getFileList(resultDir, chid, 'sf')
    filesWithQueriedQuantity = []
    endianness = getEndianness(resultDir, chid)
    for sliceFile in slcfFiles:
        f = zopen(sliceFile)
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = readSLCFheader(f, endianness)
        f.close()
        # Check if slice is correct quantity
        correctQuantity = (qty == fdsQuantity)
        if correctQuantity:
            filesWithQueriedQuantity.append(sliceFile)
    
    outFiles = []
    meshInfo = []
    for sliceFile in slcfFiles:
        print("Starting to avg file %s"%(sliceFile))
        outFile, header = slcfTimeAverage(sliceFile, dt, outFile=None, outQty=outQty, outdt=outdt)
        outFiles.append(outFile)
        print("wrote %s"%(outFile))
        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = header
        (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
        meshNum = int(sliceFile.split('_')[-2])
        meshInfo.append([meshNum, NX, NY, NZ])
    
    # Add new slices to smokeview file
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    linesSMV = zreadlines(smvFile)
    smvData = parseSMVFile(smvFile)
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    slices = files['SLICES']
    if outQty is None:
        outQty = sName + ' (%0.0ds Avg)'%(dt)
    
    for outFile, info in zip(outFiles, meshInfo):
        fname = outFile.replace('_avg.sf','.sf').split(os.sep)[-1]
        lineText = slices[fname]['LINETEXT']
        lineText = lineText.replace('.sf', '_avg.sf')
        lineText = lineText[:-4] + ('%s'%(int(lineText.replace('\n','').split()[-1])+1000)).ljust(4) + '\n'
        linesSMV.append(lineText.replace('.sf', '_avg.sf'))
        meshNum, NX, NY, NZ = info
        #linesSMV.append('SLCF     %0.0f # STRUCTURED &     0    %0.0f     0    %0.0f     0    %0.0f !      %0.0f\n'%(meshNum, NX, NY, NZ, meshNum))
        linesSMV.append(' %s\n'%(outFile.split(os.sep)[-1]))
        linesSMV.append(' %s\n'%(outQty))
        linesSMV.append(' %s\n'%(sName))
        linesSMV.append(' %s\n'%(uts))
    
    smvText = '\n'.join(linesSMV)
    smvText = smvText + '\n'
    smvText = smvText.replace('\n\n','\n')
    
    if outDir is None:
        newSmvFile = outFile.split(os.sep)[:-1]
        newSmvFile.append(smvFile.split(os.sep)[-1].replace('.smv','_avg.smv'))
        newSmvFile = os.sep.join(newSmvFile)
    else:
        newSmvFile = outDir.split(os.sep)
        newSmvFile.append(smvFile.split(os.sep)[-1].replace('.smv','_avg.smv'))
        newSmvFile = os.sep.join(newSmvFile)
        
    with open(newSmvFile, 'w') as f:
        f.write(smvText)
    print("wrote to %s"%(newSmvFile))
    
    return outFiles, outQty, filesWithQueriedQuantity, newSmvFile


def writeXYZfile(file, grid):
    """Writes grid to an xyz file
    
    This subroutine writes a grid to an xyz file.
    
    Parameters
    ----------
    file : str
        String containing the path to an xyz file
    
    grid : array(NX, NY, NZ, 3)
        Array containing float global coordinates
    """
    nx = np.unique(grid[:, 0]).shape[0]
    ny = np.unique(grid[:, 1]).shape[0]
    nz = np.unique(grid[:, 2]).shape[0]
    v = 4
    with open(file,'wb') as f:
        f.write(b'\x0c\x00\x00\x00')
        (nx, ny, nz, v) = (int(nx), int(ny), int(nz), int(v))
        f.write(nx.to_bytes(4, 'little'))
        f.write(ny.to_bytes(4, 'little'))
        f.write(nz.to_bytes(4, 'little'))
        f.write(b'\x0c\x00\x00\x00')
        #empty = np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.float32)
        empty = np.array([0], dtype=np.float32)
        empty.tofile(f)
        d1 = grid.flatten(order='F')
        d = d1.tobytes()
        f.write(d)
        byteStr = struct.pack('I', nx*ny*nz*v*4)
        f.write(byteStr)