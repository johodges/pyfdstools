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
from .colorSchemes import getVTcolors

def kalmanFilter(z, Q=1e-5, R=0.5**2):
    # This subroutine applies a kalman filter to an input set of data.
    #
    #Inputs:
    #    z: series of data to be filtered
    #    Q: process variance
    #    R: measurement variance
    #Outputs:
    #    xhat: filtered series of data
    
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

def smvVisual(obstructions,surfaces,namespace,fs=16,fig=None,ax=None,
              limits=[0,15,0,8,0,5]):
    if fig is None: fig = plt.figure(figsize=(12,12))
    if ax is None: ax = a3.Axes3D(fig)
    
    for obst in obstructions:
        pts, colors = getPtsFromObst(obst,surfaces)
        print(pts)
        print(colors)
        for pt, color in zip(pts,colors):
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

def getPtsFromObst(obst,surfaces):
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
    '''  mPts rows correlated to times, columns correlated to different groups. '''
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

def maxValuePlot(times, mPts, names, namespace, fs=16, lw=3, pcs=None, vName='',
                 yticks=None, xticks=None):
    '''  mPts rows correlated to times, columns correlated to different groups. '''
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
    figName = '%s_maxTPlot.png'%(namespace)
    print("Saving max value figure to %s"%(figName))
    plt.savefig(figName, dpi=300)
    plt.show()
    return fig

def getPlotColors(numberOfGroups):
    pcs = []
    for i in range(0,numberOfGroups):
        pcs.append(pltc.rgb2hex(np.random.rand(3)))
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
    if not isinstance(hull,scsp.Delaunay):
        hull = scsp.Delaunay(hull)
    return hull.find_simplex(p)>=0

def pts2polygons(groups):
    '''
    Build polygons from series of points.
    '''
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
    if '.zip' in resultDir:
        files = getFileListFromZip(resultDir, chid, extension)
    else:
        path = os.path.join(resultDir, '%s*.%s'%(chid, extension))
        files = glob.glob(path)
    return files

def getFileListFromZip(filename, chid, extension):
    filelist = []
    with zipfile.ZipFile(filename, 'r') as zip:
        for info in zip.infolist():
            if info.filename.split('.')[-1] == extension:
                if chid in info.filename:
                    filelist.append("%s%s%s"%(filename, os.sep, info.filename))
    return filelist

def zreadlines(file):
    f = zopen(file, readtype='r')
    lines = f.readlines()
    if '.zip' in file:
        lines = [line.decode("utf-8").replace('\r','').replace('\n','') for line in lines]
    f.close()
    return lines

def zopen(file, readtype='rb'):
    if '.zip' in file:
        zname = '%s.zip'%(file.split('.zip')[0])
        fname = file.split('.zip%s'%(os.sep))[1]
        zip = zipfile.ZipFile(zname, 'r')
        f = zip.open(fname)
    else:
        f = open(file, readtype)
    return f

def getTwoZone(z, val, lowtohigh=True):
    if lowtohigh:
        z = z[::-1]
        val = val[::-1]
        val_low = val[-1]
    else:
        val_low = val[0]
    H = z.max()
    H0 = z.min()
    tmpZ = np.linspace(0, H, 101)
    tmpV = np.interp(tmpZ, z, val)
    
    I1 = np.trapz(tmpV, tmpZ)
    I2 = np.trapz(1/tmpV, tmpZ)
    zInt = val_low*(I1*I2-H**2)/(I1+I2*val_low**2-2*val_low*H)
    
    zU = np.linspace(zInt, H, num=50)
    
    val_high_tmp = np.interp(zU, z, val)
    val_high = np.trapz(val_high_tmp, zU)/(H-zInt)
    
    zL = np.linspace(0, zInt, num=50)
    
    val_low_tmp = np.interp(zL, z, val)
    val_low = np.trapz(val_low_tmp, zL)/(zInt-H0)
    
    return val_low, val_high, zInt