#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example reads a device file of heat flux predictions and 
# generates a slice for visualization in smokeview
#
#=======================================================================
# # IMPORTS
#=======================================================================
import pyfdstools as fds
import os, argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

if __name__ == "__main__":
    print("Starting %s example"%(__file__))
    systemPath = os.path.dirname(os.path.abspath(__file__))
    chid = "hfg_slice"
    path = os.path.join(systemPath, "data", chid+".fds")
    working_dir = os.path.join(systemPath, "data", "%s.zip"%(chid))
    outdir = os.path.join(systemPath, "generated")
    
    smvFile = fds.getFileList(working_dir, chid, 'smv')[0]
    f = fds.zopen(smvFile, 'r')
    smvTxt = f.read()
    f.close()
    if type(smvTxt) is bytes:
        smvTxt = smvTxt.decode('utf-8')
    if '.zip' in smvFile:
        smvFileCustom = smvFile.replace('.smv','_custom.smv').split('.zip')[-1]
        smvFileCustom = os.path.basename(smvFileCustom)
    else:
        smvFileCustom = smvFile.replace('.smv','_custom.smv')
        smvFileCustom = os.path.basename(smvFileCustom)
    smvPathCustom = os.path.join(outdir, smvFileCustom)
    with open(smvPathCustom, 'w') as f:
        f.write(smvTxt)
    
    ref_slcf_qty = 'TEMPERATURE'
    ref_slcf_axis = 2
    devc_namespace = 'rhf01'
    devc_x0 = -0.95
    devc_dx = 0.1
    devc_y0 = -0.75
    devc_dy = 0.0
    
    quantities, slcfFiles, dimensions, meshes, centers, units = fds.readSLCFquantities(chid, working_dir)
    
    slcfFiles = [x for x,y in zip(slcfFiles, quantities) if ref_slcf_qty == y]
    
    d = fds.load_csv(working_dir, chid, '_devc')
    c = [x for x in d.columns if devc_namespace in x]
    d_x = [devc_x0+devc_dx*(float(x.split('-')[-1].replace('"',''))-1) for x in c]
    d_y = [devc_y0+devc_dy*(float(x.split('-')[-1].replace('"',''))-1) for x in c]
    d_z = [(float(x.split('-')[-2].replace('p','.'))) for x in c]
    
    d_x = np.round(d_x, decimals=2)
    d_y = np.round(d_y, decimals=2)
    d_z = np.round(d_z, decimals=2)
    
    x_unique = np.unique(d_x)
    y_unique = np.unique(d_y)
    z_unique = np.unique(d_z)
    
    x_grid, y_grid, z_grid = np.meshgrid(x_unique, y_unique, z_unique)
    d_grid = np.zeros((y_grid.shape[0], y_grid.shape[1], y_grid.shape[2], d.values.shape[0]))
    
    for i in range(0, len(d_x)):
        x_ind = np.argwhere(d_x[i] == x_unique)
        y_ind = np.argwhere(d_y[i] == y_unique)
        z_ind = np.argwhere(d_z[i] == z_unique)
        d_grid[y_ind, x_ind, z_ind, :] = d[c[i]].values
    
    d_grid = np.swapaxes(d_grid, 0, 1)
    x_grid = np.swapaxes(x_grid, 0, 1)
    z_grid = np.swapaxes(z_grid, 0, 1)
    
    if ref_slcf_axis == 1:
        xx = y_grid[0,:,:]
        zz = z_grid[0,:,:]
        dd = d_grid[0, :, :, :]
        xxx = y_unique
        zzz = z_unique
    elif ref_slcf_axis == 2:
        xx = x_grid[:,0,:]
        zz = z_grid[:,0,:]
        dd = d_grid[:, 0, :, :]
        xxx = x_unique
        zzz = z_unique
    elif ref_slcf_axis == 3:
        xx = x_grid[:,:,0]
        zz = y_grid[:,:,0]
        dd = d_grid[:, :, 0, :]
        xxx = x_unique
        zzz = y_unique
    #fds.plotSlice(xx, zz, dd[:, :, -1], ref_slcf_axis)
    
    d_time = d['Time'].values
    for slcfFile in slcfFiles:
        lims, _, _ = fds.readSingleSlcfFile(slcfFile)
        x, z, d, times, coords = fds.read2dSliceFile(slcfFile, chid)
        data = np.zeros_like(d)
        for i in range(0, len(d_time)):
            r = scipy.interpolate.RegularGridInterpolator((xxx, zzz), dd[:, :, i], method='linear', bounds_error=False, fill_value=0.0)
            data[:, :, i] = r((x, z))
        
        #fds.plotSlice(x, z, data[:, :, -1], ref_slcf_axis)
        data = np.squeeze(data)
        data2 = [np.array(data[:,:,ii], dtype=np.float32) for ii in range(0, data.shape[2])]
        axis, val = fds.getAxisFromLims(lims)
        meshnum = int(slcfFile.split('_')[-2])
        if '.zip' in slcfFile:
            outFile = slcfFile.replace('.sf','_custom.sf').split('.zip')[-1]
        else:
            outFile = slcfFile.replace('.sf','_custom.sf')
        outFile = os.path.basename(outFile)
        fds.writeSlice(outFile, outdir, chid, data2, times, axis, val,
                              'RADIATIVE HEAT FLUX GAS', 'RHF', 'kW/m2', meshnum, smvFile=smvFileCustom)