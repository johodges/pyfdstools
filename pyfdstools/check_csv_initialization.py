# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 12:05:09 2025

@author: jhodges
"""

import pyfdstools as fds
from preprocessing_initialization_greedy import greedy_uvwts_read, extract_contour, make_hull, in_hull
import numpy as np

if __name__ == "__main__":
    
    resultdir = "E:\\1JLH-NIST2023\\marcos_plume\\Thermo_Only_0p5cm_Jon\\fine_grid2\\"
    chid = "test"
    time = 'avg'
    T_axis=2
    T_value=0.15
    T_plot_mn=None
    T_plot_mx=None
    
    xs, ys, zs, xGrid1, yGrid1, zGrid1, Tf1, Uf1, Sf1 = greedy_uvwts_read(resultdir, chid, time)
    
    Tf = np.copy(Tf1)
    Uf = np.copy(Uf1)
    Sf = np.copy(Sf1)
    
    # Read smokeview file and fill each obst with nan
    smvFilePath = fds.getFileList(resultdir, chid, 'smv')[0]
    smvData = fds.parseSMVFile(smvFilePath)
    (smvGrids, smvObsts) = (smvData['grids'], smvData['obsts'])
    (smvBndfs, smvSurfs) = (smvData['bndfs'], smvData['surfs'])
    (smvFiles, bndes) = (smvData['files'], smvData['bndes'])
    for i in range(0, smvObsts.shape[0]):
        xs1 = smvObsts[i, 13:15]
        ys1 = smvObsts[i, 15:17]
        zs1 = smvObsts[i, 17:19]
        '''
        xs1 = smvObsts[i, :2]
        ys1 = smvObsts[i, 2:4]
        zs1 = smvObsts[i, 4:6]
        '''
        try:
            hull = make_hull(xs1, ys1, zs1)
        except:
            continue
        
        z1_ind = max([np.argmin(abs(zs-zs1[0]))-1,0])
        z2_ind = min([np.argmin(abs(zs-zs1[1]))+1,zs.shape[0]])
        
        x1_ind = max([np.argmin(abs(xs-xs1[0]))-1,0])
        x2_ind = min([np.argmin(abs(xs-xs1[1]))+1,xs.shape[0]])
        
        y1_ind = max([np.argmin(abs(ys-ys1[0]))-1,0])
        y2_ind = min([np.argmin(abs(ys-ys1[1]))+1,ys.shape[0]])
        
        for k in range(z1_ind, z2_ind):
            for j in range(y1_ind, y2_ind):
                for i in range(x1_ind, x2_ind):
                    p = [xs[i], ys[j], zs[k]]
                    if in_hull(p, hull):
                        Tf[i, j, k] = np.nan
                        Uf[i, j, k] = np.nan
                        Sf[i, j, k,:] = np.nan
    
    xx, yy, dd2 = extract_contour(T_axis, T_value, Tf, xs, ys, zs, xGrid1, yGrid1, zGrid1, True, None)
    fig2, ax2 = fds.plotSlice(xx, yy, dd2, T_axis, qnty_mn=T_plot_mn, qnty_mx=T_plot_mx)