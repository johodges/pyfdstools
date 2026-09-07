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
# This script has utilities used to throughout the package for tasks
# associated with visualization and filtering.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import numpy as np
import matplotlib.colors as pltc
import scipy.spatial as scsp
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
import os
import zipfile
import glob
import struct
from collections import defaultdict
from .colorSchemes import getVTcolors
import warnings

# numpy removed the np.trapz alias in numpy 2.0 in favour of
# np.trapezoid. Bind whichever name the installed numpy provides so that
# the package works on both numpy 1.x and numpy 2.x.
_trapezoid = getattr(np, 'trapezoid', None)
if _trapezoid is None:
    _trapezoid = np.trapz

def astFromGhf(ghf, h, e, Tgauge=20):
    """Calculates adiabatic surface temperature from gauge heat flux

    Solves the steady state energy balance at a heat flux gauge for the
    adiabatic surface temperature (AST):

        e*sigma*(AST^4 - Tgauge^4) + h*(AST - Tgauge) = ghf

    The quartic is solved analytically, so ghf, h and e may be scalars
    or arrays of matching shape.

    Parameters
    ----------
    ghf : float or array
        Gauge heat flux in kW/m2
    h : float or array
        Convective heat transfer coefficient in kW/m2-K
    e : float or array
        Emissivity of the gauge surface
    Tgauge : float, optional
        Gauge temperature in degrees Celsius (default 20)

    Returns
    -------
    float or array
        Adiabatic surface temperature in degrees Celsius

    Notes
    -----
    Releases before v0.0.24 mixed Celsius and Kelvin in the energy
    balance and returned the root without converting it back to
    Celsius, so a zero gauge heat flux did not return the gauge
    temperature. Values produced by this routine differ from those
    releases.
    """

    # sigma is the Stefan-Boltzmann constant expressed in kW/m2-K4 so
    # that it is consistent with h and ghf.
    sigma = 5.67e-11
    if np.any(np.asarray(h) > 1):
        print("Warning h > 1 kW/m2-K")

    # The energy balance is solved entirely in Kelvin. Releases before
    # v0.0.24 built the constant term from the gauge temperature in
    # Celsius while raising it to the fourth power in Kelvin, and
    # returned the Kelvin root without converting it back. The result
    # was that a zero gauge heat flux did not return the gauge
    # temperature: with h = 0.01 kW/m2-K and Tgauge = 20 C it returned
    # 57.6 rather than 20.
    Tgauge_K = Tgauge + 273.15

    a = e*sigma
    b = h
    c = -(a*Tgauge_K**4 + b*Tgauge_K + ghf)

    # Real positive root of a*T^4 + b*T + c = 0 in Kelvin, via the
    # resolvent of the depressed quartic.
    alpha = ((3**0.5)*(27*(a**2)*(b**4) - 256*(a**3)*(c**3))**0.5
             + 9*a*(b**2))**(1/3)
    beta = 4*((2/3)**(1/3))*c
    gamma = (18**(1/3))*a
    M = ((beta/alpha) + (alpha/gamma))**0.5

    Tast_K = (1/2)*(-M + ((2*b)/(a*M) - M**2)**0.5)
    return Tast_K - 273.15

def timeAverage2(data, times, window):
    """Applies an exponentially weighted running average to a series

    .. deprecated::
        Use :func:`timeAverage`, which resamples onto a uniform time
        base and applies a true boxcar average. This routine is retained
        for backwards compatibility with existing scripts.

    Parameters
    ----------
    data : array
        Array whose last axis is time
    times : array(NT)
        Array of timestamps corresponding to the last axis of data
    window : float
        Averaging window in seconds

    Returns
    -------
    array
        Array of the same shape as data containing the averaged values
    """

    warnings.warn(
        "timeAverage2 is deprecated and will be removed in a future "
        "release; use timeAverage instead.",
        DeprecationWarning, stacklevel=2)
    sz = data.shape
    dt = window/2
    tmax = len(times)
    if len(sz) == 1:
        data2 = np.zeros((1,sz[0]))
        data2[0, :] = data
    else:
        # Reshape data into array for time averaging
        pts = np.prod(sz[:-1])
        data2 = np.reshape(data, (pts,sz[-1]))
    
    # Time average the array
    data3 = np.zeros_like(data2) #data2.copy()
    
    data3[:, 0] = data2[:, 0]
    tmax = min([tmax, data2.shape[-1], data3.shape[-1]])
    for i in range(1, tmax):
        data_dt = times[i]-times[i-1]
        if data_dt > dt:
            data3[:, i] = data2[:, i]
        elif (times[i]-times[0]) < dt:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                data3[:,i] = (np.nanmean(data2[:,:i], axis=1)*(times[i-1]-times[0]) + data2[:,i]*data_dt)/(times[i]-times[0])
        else:
            data3[:,i] = (data3[:, i-1]*(dt-data_dt)+ data2[:, i]*data_dt)/dt
    data4 = np.reshape(data3, sz)
    return data4

def timeAverage(data, times, window, outdt=-1, smoothEnds=False,
                queryTime=-1):
    """Boxcar time-averages a data array over a moving window

    The data are first interpolated onto a uniform time base built from
    the smallest positive timestep found in times, then convolved with a
    rectangular window.

    Parameters
    ----------
    data : array(NX, NY, NT)
        Array whose last axis is time
    times : array(NT)
        Array of timestamps corresponding to the last axis of data
    window : float
        Averaging window in seconds
    outdt : float, optional
        Timestep of the returned series. A value <= 0 returns the series
        on the internal uniform time base (default -1)
    smoothEnds : bool, optional
        If True the first and last half-windows are averaged over the
        partial window available; if False they are copied from the
        interpolated series unchanged (default False)
    queryTime : float, optional
        If > 0, average over a single window centred on this time and
        return a single frame (default -1)

    Returns
    -------
    array
        Array containing the time-averaged data
    array or float
        Timestamps of the returned data, or queryTime when a single
        frame was requested
    """

    tmax = np.nanmax(times)
    tmin = np.nanmin(times)
    if window > (tmax-tmin):
        print("Warning, windowSize > time interval. Not time averaging.")
        return data, times
    
    dts = times[1:] - times[:-1]
    dt = np.nanmin(dts[dts > 0])
    
    if queryTime > 0:
        tmin = np.floor((queryTime-window/2)/dt)*dt
        tmax = np.ceil((queryTime+window/2)/dt)*dt
        data2 = np.zeros((data.shape[0], data.shape[1], 1))
        t1 = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
        for i in range(0, data.shape[0]):
            for j in range(0, data.shape[1]):
                v1 = np.interp(t1, times, data[i, j, :])
                data2[i, j, 0] = np.nanmean(v1)
        return data2, queryTime
    
    tmin = np.floor(tmin/dt)*dt
    tmax = np.ceil(tmax/dt)*dt
    
    t1 = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
    N = int(np.round(window/dt))
    N2 = int(N/2)
    f = np.zeros((N)) + 1
    f = f/f.sum()
    
    if data.shape[-1] < times.shape[-1]:
        print("Warning, data shape is less than timesteps, truncating.")
        times = times[:data.shape[-1]]
    
    data2 = np.zeros((data.shape[0], data.shape[1], int(len(t1)-N+1)))
    for i in range(0, data.shape[0]):
        for j in range(0, data.shape[1]):
            v1 = np.interp(t1, times, data[i, j, :])
            data2[i, j, :] = np.convolve(v1, f, mode='valid')
            if smoothEnds:
                for k in range(0, N2):
                    data2[i,j,k] = np.nanmean(v1[k:k+N2])
                for k in range(data2.shape[2]-N, data2.shape[2]):
                    data2[i,j,k] = np.nanmean(v1[k:k+N2])
            else:
                data2[i, j, :N2] = v1[:N2]
                data2[i ,j,-N2:] = v1[-N2:]
    
    if outdt > 0:
        tmin = np.floor(tmin/outdt)*outdt
        tmax = np.ceil(tmax/outdt)*outdt
        t2 = np.linspace(tmin, tmax, int((tmax-tmin)/outdt + 1))
        data3 = np.zeros((data.shape[0], data.shape[1], t2.shape[0]))
        for i in range(0, data.shape[0]):
            for j in range(0, data.shape[1]):
                data3[i, j, :] = np.interp(t2, t1[:data2.shape[2]], data2[i, j, :])
        return data3, t2
    else:
        return data2, t1[:data2.shape[2]]

def kalmanFilter(z, Q=1e-5, R=0.5**2):
    """Applies a scalar Kalman filter to a series of measurements

    Useful for smoothing noisy device output, such as a thermocouple
    trace, without the phase lag a boxcar average introduces.

    Parameters
    ----------
    z : array(N)
        Series of measurements to filter
    Q : float, optional
        Process variance. Larger values let the estimate follow the
        measurements more closely (default 1e-5)
    R : float, optional
        Measurement variance. Larger values smooth more heavily
        (default 0.25)

    Returns
    -------
    array(N)
        Filtered series
    """

    # intial parameters
    sz = z.shape[0] # size of array
    
    # allocate space for arrays
    xhat=np.zeros(sz)      # a posteri estimate of x
    P=np.zeros(sz)         # a posteri error estimate
    xhatminus=np.zeros(sz) # a priori estimate of x
    Pminus=np.zeros(sz)    # a priori error estimate
    K=np.zeros(sz)         # gain or blending factor
    
    # intial guesses
    xhat[0] = z[0]
    P[0] = 1.0
    
    for k in range(1,sz):
        # time update
        xhatminus[k] = xhat[k-1]
        Pminus[k] = P[k-1]+Q
    
        # measurement update
        K[k] = Pminus[k]/( Pminus[k]+R )
        xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
        P[k] = (1-K[k])*Pminus[k]
    
    return xhat

def smvVisual(obstructions, surfaces, namespace, fs=16, fig=None, ax=None,
              limits=[0, 15, 0, 8, 0, 5]):
    """Renders smokeview obstructions as a 3-D figure

    Parameters
    ----------
    obstructions : list
        List of obstructions parsed from a smokeview file
    surfaces : list
        List of surfaces parsed from a smokeview file, used for colors
    namespace : str
        Prefix used to build the saved figure name
    fs : int, optional
        Font size for the axis labels and ticks (default 16)
    fig : matplotlib.figure.Figure, optional
        Figure to draw into. A new figure is created when omitted
    ax : matplotlib.axes.Axes, optional
        3-D axes to draw into. New axes are created when omitted
    limits : list, optional
        Six component list of axis limits
        [xmin, xmax, ymin, ymax, zmin, zmax]

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the rendered obstructions
    matplotlib.axes.Axes
        Axes containing the rendered obstructions
    """

    if fig is None:
        fig = plt.figure(figsize=(12, 12))
    if ax is None:
        # Axes3D(fig) no longer attaches itself to the figure in
        # matplotlib >= 3.4; add_subplot is the supported spelling.
        ax = fig.add_subplot(projection='3d')

    for obst in obstructions:
        pts, colors = getPtsFromObst(obst, surfaces)
        for pt, color in zip(pts, colors):
            f = a3.art3d.Poly3DCollection(pt)
            f.set_color(color)
            f.set_edgecolor('k')
            #f.set_alpha(1)
            ax.add_collection3d(f)
    plt.xlim(limits[0],limits[1])
    plt.ylim(limits[2],limits[3])
    ax.set_zlim(limits[4],limits[5])
    plt.xlabel('x (m)',fontsize=fs)
    plt.ylabel('y (m)',fontsize=fs)
    ax.set_zlabel('z (m)',fontsize=fs)
    plt.tick_params(labelsize=fs)
    plt.savefig('%s_smvvisual.png'%(namespace),dpi=300)
    
    return fig, ax

def buildSMVgeometry(file):
    """Parses surfaces and obstructions from a smokeview file

    Parameters
    ----------
    file : str
        String containing the path to a smokeview file

    Returns
    -------
    list
        List of surfaces, each containing
        [name, Tign, emissivity, type, texture width, texture height,
         r, g, b, a]
    list
        List of obstructions, each a list of floats read from the
        smokeview file
    """

    with open(file,'r') as f:
        lines = f.readlines()
    inds = []
    for i in range(0,len(lines)):
        if lines[i][0] != ' ' and lines[i][0] != '-':
            inds.append(i)
    surfaces = []
    obstructions = []
    for ind in inds:
        if 'SURFACE' in lines[ind]:
            sname = ' '.join(lines[ind+1].split())
            (Tign,eps) = (lines[ind+2].split()[0],lines[ind+2].split()[1])
            (stype,t_width,t_height) = (lines[ind+3].split()[0],lines[ind+3].split()[1],lines[ind+3].split()[2])
            (c1,c2,c3,c4) = (lines[ind+3].split()[3],lines[ind+3].split()[4],lines[ind+3].split()[5],lines[ind+3].split()[6])
            surfaces.append([sname,Tign,eps,stype,t_width,t_height,c1,c2,c3,c4])
        if 'OBST' in lines[ind] and '_OBST' not in lines[ind]:
            nObst = int(lines[ind+1].split()[0])
            for i in range(0,nObst):
                obst = [float(x) for x in lines[ind+i+2].split()]
                obstructions.append(obst)
    return surfaces, obstructions

def getPtsFromObst(obst, surfaces):
    """Builds corner points and face colors for a smokeview obstruction

    Parameters
    ----------
    obst : list
        Obstruction record parsed from a smokeview file. Entries 0-5 are
        the bounding box and entries 7-12 are the surface indices of the
        y-, y+, x-, x+, z- and z+ faces
    surfaces : list
        List of surfaces parsed from a smokeview file

    Returns
    -------
    array(8, 3)
        Array containing the corner coordinates of the obstruction
    list
        List of six (r, g, b, a) tuples, one per face
    """

    pts = []
    colors = []
    pts = np.array([[obst[0],obst[2],obst[4]],
                   [obst[0],obst[2],obst[5]],
                   [obst[1],obst[2],obst[5]],
                   [obst[1],obst[2],obst[4]],
                   [obst[0],obst[3],obst[5]],
                   [obst[1],obst[3],obst[5]],
                   [obst[1],obst[3],obst[4]],
                   [obst[0],obst[3],obst[4]]])
    
    # y-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[1],obst[2],obst[4]),
    #            (obst[1],obst[2],obst[5]),(obst[0],obst[2],obst[5])])
    surf = surfaces[int(obst[7])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # y-positive surface
    #pts.append([(obst[0],obst[3],obst[4]),(obst[1],obst[3],obst[4]),
    #            (obst[1],obst[3],obst[5]),(obst[0],obst[3],obst[5])])
    surf = surfaces[int(obst[8])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # x-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[0],obst[2],obst[5]),
    #            (obst[0],obst[3],obst[5]),(obst[0],obst[3],obst[4])])
    surf = surfaces[int(obst[9])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # x-positive surface
    #pts.append([(obst[1],obst[2],obst[4]),(obst[1],obst[2],obst[5]),
    #            (obst[1],obst[3],obst[5]),(obst[1],obst[3],obst[4])])
    surf = surfaces[int(obst[10])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # z-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[1],obst[2],obst[4]),
    #            (obst[1],obst[3],obst[4]),(obst[0],obst[3],obst[4])])
    surf = surfaces[int(obst[11])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # z-positive surface
    #pts.append([(obst[0],obst[2],obst[5]),(obst[1],obst[2],obst[5]),
    #            (obst[1],obst[3],obst[5]),(obst[0],obst[3],obst[5])])
    surf = surfaces[int(obst[12])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    
    return pts, colors

def maxValueCSV(times, mPts, names, namespace):
    """Writes a time series of per-group maximum values to a csv file

    Parameters
    ----------
    times : array(NT)
        Array of timestamps
    mPts : array(NT, NG)
        Array whose rows correspond to times and whose columns
        correspond to different groups
    names : list
        List of NG group names used as column headers
    namespace : str
        Prefix used to build the output file name

    Returns
    -------
    str
        Name of the csv file which was written
    """

    numberOfGroups = mPts.shape[1]
    header = 'Time,'
    for i in range(0,numberOfGroups):
        name = names[i].replace(',','_')
        header = header+name+','
    header = header[:-1]+'\n'
    data = np.append(np.reshape(times,(times.shape[0],1)),mPts,axis=1)
    csvName = '%s.csv'%(namespace)
    print("Saving max value csv to %s"%(csvName))
    np.savetxt(csvName, data, delimiter=',', header=header)
    return '%s.csv'%(namespace)   

def maxValuePlot(times, mPts, names, figName, fs=16, lw=3, pcs=None, vName='',
                 yticks=None, xticks=None):
    """Plots a time series of per-group maximum values

    Parameters
    ----------
    times : array(NT)
        Array of timestamps
    mPts : array(NT, NG)
        Array whose rows correspond to times and whose columns
        correspond to different groups
    names : list
        List of NG group names used in the legend
    figName : str
        Path the figure is saved to
    fs : int, optional
        Font size (default 16)
    lw : int, optional
        Line width (default 3)
    pcs : list, optional
        List of colors, one per group. Generated when omitted
    vName : str, optional
        Label for the y-axis
    yticks : list, optional
        Explicit y-axis tick locations
    xticks : list, optional
        Explicit x-axis tick locations

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """

    numberOfGroups = mPts.shape[1]
    if pcs is None:
        pcs = getVTcolors()
        if len(pcs) < numberOfGroups: pcs = getPlotColors(numberOfGroups)
    fig = plt.figure(figsize=(12,8))
    for i in range(0,numberOfGroups):
        plt.plot(times,mPts[:,i],color=pcs[i],label=names[i],linewidth=lw)
    if yticks is not None: plt.yticks(yticks)
    if xticks is not None: plt.xticks(xticks)
    plt.legend(fontsize=fs)
    plt.xlabel('time (s)',fontsize=fs)
    plt.ylabel('%s'%(vName),fontsize=fs)
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    plt.savefig(figName, dpi=300)
    plt.show()
    return fig

def getPlotColors(numberOfGroups):
    """Builds a list of visually distinct random plot colors

    Colors whose channel sum falls outside [0.3, 2.7] are rejected so
    that neither very dark nor very light colors are produced.

    Parameters
    ----------
    numberOfGroups : int
        Number of colors to generate

    Returns
    -------
    list
        List of hex color strings
    """

    pcs = []
    for i in range(0, numberOfGroups):
        v = np.random.rand(3)
        su = np.sum(v)
        while (su > 2.7) or (su < 0.3):
            v = np.random.rand(3)
            su = np.sum(v)
        tmp = pltc.rgb2hex(v)
        pcs.append(tmp)
    return pcs

def pointsFromXB(XB,extend=[0,0,0]):
    ''' This routine builds a list of XYZ points from an obstruction XB
    
    Inputs:
        XB: Septuplet containing [xmin, xmax, ymin, ymax, zmin, zmax]
        extend: Float array containing amount to extend polygon along each axis
    Outputs:
        pts: List of corner points
    '''
    pts = [[XB[0]-extend[0],XB[2]-extend[1],XB[4]-extend[2]],
           [XB[0]-extend[0],XB[2]-extend[1],XB[5]+extend[2]],
           [XB[0]-extend[0],XB[3]+extend[1],XB[4]-extend[2]],
           [XB[0]-extend[0],XB[3]+extend[1],XB[5]+extend[2]],
           [XB[1]+extend[0],XB[2]-extend[1],XB[4]-extend[2]],
           [XB[1]+extend[0],XB[2]-extend[1],XB[5]+extend[2]],
           [XB[1]+extend[0],XB[3]+extend[1],XB[4]-extend[2]],
           [XB[1]+extend[0],XB[3]+extend[1],XB[5]+extend[2]]]
    return pts

def in_hull(p, hull):
    """Tests whether points fall inside a convex hull

    Parameters
    ----------
    p : array(N, D)
        Array of points to test
    hull : array(M, D) or scipy.spatial.Delaunay
        Points defining the hull, or a pre-computed triangulation

    Returns
    -------
    array(N)
        Boolean array which is True for points inside the hull
    """

    if not isinstance(hull,scsp.Delaunay):
        hull = scsp.Delaunay(hull)
    return hull.find_simplex(p)>=0

def pts2polygons(groups):
    """Builds convex hull polygons from groups of point sets

    Point sets which cannot be triangulated (for example a degenerate
    planar set) are reported and skipped.

    Parameters
    ----------
    groups : list
        List of groups, each of which is a list of point arrays

    Returns
    -------
    list
        List of groups, each containing a list of
        scipy.spatial.ConvexHull objects
    int
        Number of groups
    """

    polygons = []
    
    for group in groups:
        linkedPolygons = []
        for pts in group:
            try:
                linkedPolygons.append(scsp.ConvexHull(pts))
            except:
                print("Failed points:")
                print(pts)
        polygons.append(linkedPolygons)
        
    return polygons, len(polygons)

def getFileList(resultDir, chid, extension):
    """Lists the result files for a case with a given extension

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or the path to a zip
        archive containing them
    chid : str
        FDS CHID of the case
    extension : str
        File extension to search for, without the leading period

    Returns
    -------
    list
        List of matching file paths. Files inside an archive are
        returned as '<archive>.zip<os.sep><name inside the archive>'
    """

    if '.zip' in resultDir:
        files = getFileListFromZip(resultDir, chid, extension)
    else:
        path = os.path.join(resultDir, '%s*.%s'%(chid, extension))
        files = glob.glob(path)
    return files

def getFileListFromZip(filename, chid, extension):
    """Lists the files in a zip archive matching a chid and extension

    Parameters
    ----------
    filename : str
        Path to the zip archive
    chid : str
        FDS CHID of the case
    extension : str
        File extension to search for, without the leading period

    Returns
    -------
    list
        List of paths of the form
        '<archive>.zip<os.sep><name inside the archive>'
    """

    filelist = []
    with zipfile.ZipFile(filename, 'r') as zip:
        for info in zip.infolist():
            if info.filename.split('.')[-1] == extension:
                if chid in info.filename:
                    filelist.append("%s%s%s"%(filename, os.sep, info.filename))
    return filelist

def zreadlines(file):
    """Reads the lines of a text file which may live inside an archive

    Parameters
    ----------
    file : str
        Path to a text file, or to a file inside a zip archive using the
        '<archive>.zip<os.sep><name inside the archive>' convention

    Returns
    -------
    list
        List of strings, one per line, with line endings removed
    """

    f = zopen(file, readtype='r')
    lines = f.readlines()
    f.close()
    # Members of an archive come back as bytes. Strip the line endings
    # for both sources so that a file read from a directory and the same
    # file read from an archive produce identical lines.
    cleaned = []
    for line in lines:
        if isinstance(line, bytes):
            line = line.decode('utf-8', errors='replace')
        cleaned.append(line.replace('\r', '').replace('\n', ''))
    return cleaned

def zopen(file, readtype='rb'):
    """Opens a file which may live inside a zip archive

    Parameters
    ----------
    file : str
        Path to a file, or to a file inside a zip archive using the
        '<archive>.zip<os.sep><name inside the archive>' convention
    readtype : str, optional
        Mode used when opening a file which is not inside an archive
        (default 'rb'). Files inside an archive are always opened
        in binary mode

    Returns
    -------
    file
        Open file object
    """

    if '.zip' in file:
        zname = '%s.zip'%(file.split('.zip')[0])
        fname = file.split('.zip%s'%(os.sep))[1]
        zip = zipfile.ZipFile(zname, 'r')
        f = zip.open(fname)
    else:
        f = open(file, readtype)
    return f


def getSmvFile(resultDir, chid):
    """Returns the smokeview file for a case

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case

    Returns
    -------
    str
        Path to the smokeview file

    Raises
    ------
    FileNotFoundError
        If no smokeview file for the case is present. Indexing the
        result of getFileList directly raised a bare IndexError here,
        which gave the caller nothing to act on.
    """

    smvFiles = getFileList(resultDir, chid, 'smv')
    if len(smvFiles) == 0:
        raise FileNotFoundError(
            "No smokeview (.smv) file for chid %s was found in %s."
            % (chid, resultDir))
    return smvFiles[0]


def getEndianness(resultDir, chid):
    """Determines the byte order used to write a case's output files

    FDS writes a .end file whose second record is the integer 1. The
    byte order which decodes that record as 1 is the byte order of every
    other binary output file for the case. Little-endian is assumed when
    no .end file is present, which is correct for all mainstream
    platforms FDS is built for.

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case

    Returns
    -------
    str
        '<' for little-endian or '>' for big-endian, in the notation
        used by the struct and numpy modules
    """

    endFiles = getFileList(resultDir, chid, 'end')
    if len(endFiles) == 0:
        #print('Unable to find endianness file, %s.'%(os.path.join(resultDir, chid+'.end')))
        #print("Assuming little-endian.")
        return "<"
        
    f = zopen(endFiles[0])
    data = f.read()
    f.close()
    if struct.unpack('<i',data[4:8])[0] == 1:
        return '<'
    elif struct.unpack('>i',data[4:8])[0] == 1:
        return ">"
    else:
        print("Unable to determine endianness from file, %s."%(os.path.join(resultDir, chid+'.end')))
        print("Assuming little-endian.")
        return "<"

def getDatatypeByEndianness(datatype1, endianness):
    """Applies a byte order to a numpy datatype

    Parameters
    ----------
    datatype1 : type or numpy.dtype
        Base datatype, for example numpy.float32
    endianness : str
        '<' for little-endian or '>' for big-endian

    Returns
    -------
    numpy.dtype
        Datatype with the requested byte order applied

    Raises
    ------
    ValueError
        If endianness is neither '<' nor '>'
    """

    if endianness == '>':
        datatype2 = np.dtype(datatype1).newbyteorder('>')
    elif endianness == '<':
        datatype2 = np.dtype(datatype1).newbyteorder('<')
    else:
        raise ValueError(
            "Endianness must be '<' or '>'; received %s" % (endianness))
    return datatype2


def getFileListFromResultDir(resultDir, chid, ext):
    """Lists the result files for a case with a given extension

    This is a backwards-compatible alias for :func:`getFileList`, which
    it now delegates to. New code should call getFileList directly.

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case
    ext : str
        File extension to search for, without the leading period

    Returns
    -------
    list
        List of matching file paths
    """

    return getFileList(resultDir, chid, ext)
    

def getGridsFromXyzFiles(xyzFiles, chid):
    """Builds a dictionary of mesh grids from a list of xyz files

    Parameters
    ----------
    xyzFiles : list
        List of paths to xyz files written by FDS
    chid : str
        FDS CHID of the case, used to recover the mesh number from each
        file name

    Returns
    -------
    defaultdict
        Dictionary keyed by mesh string, each value a dictionary with
        'xGrid', 'yGrid' and 'zGrid' arrays of shape (NX, NY, NZ)
    """

    grids = defaultdict(bool)
    for xyzFile in xyzFiles:
        grid, gridHeader = readXYZfile(xyzFile)
        xGrid, yGrid, zGrid = rearrangeGrid(grid)
        mesh = xyzFile.split(chid)[-1].split('.xyz')[0].replace('_','')
        meshStr = "%s"%(chid) if mesh == '' else mesh
        
        grids[meshStr] = defaultdict(bool)
        grids[meshStr]['xGrid'] = xGrid
        grids[meshStr]['yGrid'] = yGrid
        grids[meshStr]['zGrid'] = zGrid
    return grids


def rearrangeGrid(grid, iblock=False):
    """Builds a meshgrid based on grid array
    
    Parameters
    ----------
    grid : array(NX, NY, NZ, 3)
        Array containing float global coordinates
    
    Returns
    -------
    array(NX, NY, NZ)
        Array containing float global coordinates for x-axis
    array(NX, NY, NZ)
        Array containing float global coordinates for y-axis
    array(NX, NY, NZ)
        Array containing float global coordinates for z-axis
    """
    
    xs = np.unique(grid[:,0])
    ys = np.unique(grid[:,1])
    zs = np.unique(grid[:,2])
    xGrid, yGrid, zGrid = np.meshgrid(xs, ys, zs)
    xGrid = np.swapaxes(xGrid, 0, 1)
    yGrid = np.swapaxes(yGrid, 0, 1)
    zGrid = np.swapaxes(zGrid, 0, 1)
    
    if iblock:
        iGrid = np.zeros_like(xGrid)
        for i in range(0, xGrid.shape[0]):
            for j in range(0, xGrid.shape[1]):
                for k in range(0, xGrid.shape[2]):
                    dist = abs(xs-xGrid[i, j, k]) + abs(ys-yGrid[i, j, k]) + abs(zs-zGrid[i, j, k])
                    ind = np.argmin(dist)
                    iGrid[i, j, k] = grid[ind, 3]
        return xGrid, yGrid, zGrid, iGrid
    else:
        return xGrid, yGrid, zGrid


def readXYZfile(file):
    """Reads points from an xyz file
    
    This subroutine reads grid coordinates from an xyz file. Note,
    xyz file can be generated in FDS by adding WRITE_XYZ=.TRUE. in the
    &DUMP namelist.
    
    Parameters
    ----------
    file : str
        String containing the path to an xyz file or xyz file in an
        archive
    
    Returns
    -------
    array(NX, NY, NZ, 3)
        Array containing float global coordinates
    array()
        Array containing header information from xyz file
    """
    
    try:
        f = zopen(file)
    except FileNotFoundError:
        tmp = file[:-4].split('_')
        meshStr = str(int(tmp[-1]))
        file2 = '%s_%s.xyz'%('_'.join(tmp[:-1]), meshStr)
        f = zopen(file2)
    header = struct.unpack('<iiiiif', f.read(24))
    (nx, ny, nz) = (header[1], header[2], header[3])
    data = np.frombuffer(f.read(nx*ny*nz*4*4), dtype=np.float32)
    grid = np.reshape(data, (int(data.shape[0]/4), 4),order='F')
    f.close()
    return grid, header[1:-1]
    

def getAbsoluteGrid(grids, makeUniform=False, decimals=4):
    """Builds absolute grid from defaultdict of local grids
    
    Parameters
    ----------
    grids : defaultdict
        Dictionary containing arrays(NX, NY, NZ) for each grid
    makeUniform : bool
        Boolean flag indicating whether the absolute grid should be
        interpolated to be rectangular at the most resolved mesh.
    Returns
    -------
    array(NX, NY, NZ, 3)
        Array containing the absolute grid coordinates
    """
    
    if makeUniform:
        mins = []
        maxs = []
        deltas = []
        for key in list(grids.keys()):
            xGrid = grids[key]['xGrid']
            yGrid = grids[key]['yGrid']
            zGrid = grids[key]['zGrid']
            mins.append([xGrid.min(), yGrid.min(), zGrid.min()])
            maxs.append([xGrid.max(), yGrid.max(), zGrid.max()])
            dx = np.round(xGrid[1, 0, 0] - xGrid[0, 0, 0], decimals=decimals)
            dy = np.round(yGrid[0, 1, 0] - yGrid[0, 0, 0], decimals=decimals)
            dz = np.round(zGrid[0, 0, 1] - zGrid[0, 0, 0], decimals=decimals)
            deltas.append([dx, dy, dz])
        absMins = np.min(mins, axis=0)
        absMaxs = np.max(maxs, axis=0)
        absDeltas = np.min(deltas, axis=0)
        
        Nx = int(np.round((absMaxs[0] - absMins[0]) / absDeltas[0]) + 1)
        Ny = int(np.round((absMaxs[1] - absMins[1]) / absDeltas[1]) + 1)
        Nz = int(np.round((absMaxs[2] - absMins[2]) / absDeltas[2]) + 1)
        
        xs = np.linspace(absMins[0], absMaxs[0], Nx)
        ys = np.linspace(absMins[1], absMaxs[1], Ny)
        zs = np.linspace(absMins[2], absMaxs[2], Nz)
        
        xGrid_abs, yGrid_abs, zGrid_abs = np.meshgrid(xs, ys, zs)
        xGrid_abs = np.swapaxes(xGrid_abs, 0, 1)
        yGrid_abs = np.swapaxes(yGrid_abs, 0, 1)
        zGrid_abs = np.swapaxes(zGrid_abs, 0, 1)
        
        grid_abs = np.zeros((xGrid_abs.shape[0],
                             xGrid_abs.shape[1],
                             xGrid_abs.shape[2],
                             3))
        grid_abs[:, :, :, 0] = xGrid_abs
        grid_abs[:, :, :, 1] = yGrid_abs
        grid_abs[:, :, :, 2] = zGrid_abs
    else:
        all_xs = []
        all_ys = []
        all_zs = []
        for key in list(grids.keys()):
            xs = grids[key]['xGrid'][:, 0, 0]
            ys = grids[key]['yGrid'][0, :, 0]
            zs = grids[key]['zGrid'][0, 0, :]
            all_xs.extend(np.round(xs, decimals=decimals))
            all_ys.extend(np.round(ys, decimals=decimals))
            all_zs.extend(np.round(zs, decimals=decimals))
        
        abs_xs = np.unique(all_xs)
        abs_ys = np.unique(all_ys)
        abs_zs = np.unique(all_zs)
        xGrid_abs, yGrid_abs, zGrid_abs = np.meshgrid(abs_xs, abs_ys, abs_zs)
        
        xGrid_abs = np.swapaxes(xGrid_abs, 0, 1)
        yGrid_abs = np.swapaxes(yGrid_abs, 0, 1)
        zGrid_abs = np.swapaxes(zGrid_abs, 0, 1)
        
        grid_abs = np.zeros((xGrid_abs.shape[0],
                             xGrid_abs.shape[1],
                             xGrid_abs.shape[2],
                             3))
        grid_abs[:, :, :, 0] = xGrid_abs
        grid_abs[:, :, :, 1] = yGrid_abs
        grid_abs[:, :, :, 2] = zGrid_abs
    return grid_abs


def getTwoZone(z, val, lowtohigh=True):
    """Reduces a vertical profile to an equivalent two-zone model

    Applies the two-zone reduction of a continuous vertical profile
    described in the FDS verification guide: the interface height is
    chosen so that the integrals of the profile and of its reciprocal
    are preserved, then the layer values are the averages above and
    below that height.

    Parameters
    ----------
    z : array(N)
        Array of elevations
    val : array(N)
        Array of profile values at each elevation
    lowtohigh : bool, optional
        Retained for backwards compatibility and ignored. The ordering
        of z is now detected from the data itself, so a profile given in
        either direction produces the same result (default True)

    Returns
    -------
    float
        Average value of the lower layer
    float
        Average value of the upper layer
    float
        Elevation of the interface between the two layers

    Notes
    -----
    Releases before v0.0.24 reversed an ascending profile when
    lowtohigh was True, which broke the interpolation and returned an
    upper layer value cooler than the lower layer. Values produced by
    this routine differ from those releases when the input was
    ascending.
    """

    z = np.asarray(z, dtype=float)
    val = np.asarray(val, dtype=float)

    # np.interp requires the sample points to increase, so normalise the
    # profile to ascending elevation whichever way the caller supplied
    # it. Reversing an already ascending profile, as releases before
    # v0.0.24 did when lowtohigh was True, fed np.interp a descending
    # x-array and produced an upper layer cooler than the lower one.
    if z[0] > z[-1]:
        z = z[::-1]
        val = val[::-1]

    val_low = val[0]
    H = z.max()
    H0 = z.min()
    tmpZ = np.linspace(0, H, 101)
    tmpV = np.interp(tmpZ, z, val)
    if np.isclose(np.nanmin(tmpV), np.nanmax(tmpV), atol=1e-1):
        zInt = H
        return val_low, val_low, zInt
    
    I1 = _trapezoid(tmpV, tmpZ)
    I2 = _trapezoid(1/tmpV, tmpZ)
    zInt = val_low*(I1*I2-H**2)/(I1+I2*val_low**2-2*val_low*H)
    
    zU = np.linspace(zInt, H, num=50)
    
    val_high_tmp = np.interp(zU, z, val)
    val_high = _trapezoid(val_high_tmp, zU)/(H-zInt)
    
    zL = np.linspace(0, zInt, num=50)
    
    val_low_tmp = np.interp(zL, z, val)
    val_low = _trapezoid(val_low_tmp, zL)/(zInt-H0)
    
    return val_low, val_high, zInt
