# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 16:08:10 2025

@author: jhodges
"""

import pyfdstools as fds
import numpy as np
import glob, os
import matplotlib.pyplot as plt
import scipy.spatial as scsp
import scipy.interpolate as scip

def in_hull(p, hull):
    return hull.find_simplex(p)>=0

def make_hull(xs1, ys1, zs1):
    pts = [[xs1[0], ys1[0], zs1[0]],
           [xs1[-1], ys1[0], zs1[0]],
           [xs1[0], ys1[-1], zs1[0]],
           [xs1[-1], ys1[-1], zs1[0]],
           [xs1[0], ys1[0], zs1[-1]],
           [xs1[-1], ys1[0], zs1[-1]],
           [xs1[0], ys1[-1], zs1[-1]],
           [xs1[-1], ys1[-1], zs1[-1]]
           ]
    hull = scsp.Delaunay(pts)
    return hull

def get_center_grid(grid1):
    grid2 = (grid1[:-1,:-1,:-1] + grid1[1:,1:,1:])/2
    return grid2

def get_flat_grid(grid):
    grid_flat = grid.flatten(order='F')
    return grid_flat

def grid_to_points(xGrid, yGrid, zGrid):
    xGrid_flat = get_flat_grid(xGrid)
    yGrid_flat = get_flat_grid(yGrid)
    zGrid_flat = get_flat_grid(zGrid)
    pts = np.zeros((xGrid_flat.shape[0], 3))
    pts[:, 0] = xGrid_flat
    pts[:, 1] = yGrid_flat
    pts[:, 2] = zGrid_flat
    return pts

def flatten_data(data):
    if len(data.shape) > 3:
        data_flat = np.zeros((data.shape[0]*data.shape[1]*data.shape[2], data.shape[3]))
        for i in range(0, data.shape[3]):
            data_flat[:, i] = data[:,:,:,i].flatten(order='F')
    else:
        data_flat = data.flatten(order='F')
    return data_flat

def interpolate(interp, pts, mn, mx):
    fine_flat = np.zeros((pts.shape[0], len(interp)))
    for i in range(0, len(interp)):
        tmp = interp[i](pts)
        tmp[tmp > mx[i]] = mx[i]
        tmp[tmp < mn[i]] = mn[i]
        fine_flat[:,i] = tmp
    return fine_flat

def write_csv_data(data_flat1, namespace, varspace, grid, centered):
    out_file = namespace.replace('#VARIABLE#',varspace)
    if len(data_flat1) == 1:
        data_flat = np.reshape(data_flat1, (data_flat1.shape[0],1))
    else:
        data_flat = data_flat1
    with open(out_file, 'w') as f:
        if centered:
            f.write('1,%d,1,%d,1,%d\n'%(grid.shape))
        else:
            f.write('0,%d,0,%d,0,%d\n'%(grid.shape))
        np.savetxt(f, data_flat, fmt='%0.13f', delimiter=',')

def find_uvw_times(resultdir, chid):
    uvw_files = glob.glob(os.path.join(resultdir + os.sep + chid + '_uvw*.csv'))
    uvw_times = [float(x.split('_')[-2].replace('t','')) for x in uvw_files]
    
    unique_times = list(set(uvw_times))
    unique_times.sort()
    return unique_times

def find_uvw_meshes(resultdir, chid):
    uvw_files = glob.glob(os.path.join(resultdir + os.sep + chid + '_uvw*.csv'))
    uvw_meshes = [int(x.split('_')[-1].replace('m','').replace('.csv','')) for x in uvw_files]
    unique_meshes = list(set(uvw_meshes))
    unique_meshes.sort()
    return unique_meshes

def get_smv_grids(resultdir, chid):
    smvFile = os.path.join(resultdir, chid + '.smv')
    smvOutputs = fds.parseSMVFile(smvFile)
    grids = smvOutputs['grids']
    return grids

def greedy_uvwts_read(resultdir, chid, time):
    unique_meshes = find_uvw_meshes(resultdir, chid)
    grids = get_smv_grids(resultdir, chid)
    
    # Load Coarse Results
    coarse_results = {}
    xs_all = []
    ys_all = []
    zs_all = []
    for i in range(0, len(unique_meshes)):
        xs1 = grids[i][0][:,1]
        ys1 = grids[i][1][:,1]
        zs1 = grids[i][2][:,1]
        
        xs2 = (xs1[:-1] + xs1[1:])/2
        ys2 = (ys1[:-1] + ys1[1:])/2
        zs2 = (zs1[:-1] + zs1[1:])/2
        
        xGrid1, yGrid1, zGrid1 = np.meshgrid(xs1, ys1, zs1, indexing='ij')
        
        # Load UVW
        last_file = resultdir + os.sep + chid + '_uvw_t%d_m%d.csv'%(time, unique_meshes[i])
        d = np.loadtxt(last_file, delimiter=',', skiprows=1)
        U = np.zeros((xGrid1.shape[0], xGrid1.shape[1], xGrid1.shape[2], d.shape[1]))
        for j in range(0, d.shape[1]):
            U[:, :, :, j] = np.reshape(d[:,j], shape=xGrid1.shape, order='F')
        
        # Load T
        last_file = resultdir + os.sep + chid + '_tmp_t%d_m%d.csv'%(time, unique_meshes[i])
        d = np.loadtxt(last_file, delimiter=',', skiprows=1)
        T = np.reshape(d, shape=xGrid1[1:,1:,1:].shape, order='F')
        
        # Load Spec
        last_file = resultdir + os.sep + chid + '_spec_t%d_m%d.csv'%(time, unique_meshes[i])
        d = np.loadtxt(last_file, delimiter=',', skiprows=1)
        S = np.zeros((T.shape[0], T.shape[1], T.shape[2], d.shape[1]))
        for j in range(0, d.shape[1]):
            S_tmp = np.reshape(d[:,j], shape=T.shape, order='F')
            S[:, :, :, j] = S_tmp
        
        # Make convex hull for checking if a fine mesh is in a coarse mesh
        hull = make_hull(xs1, ys1, zs1)
        
        # Store coarse results
        coarse_results[i] = {'U': U, 'T': T, 'S': S, 'hull': hull, 'xs1': xs1, 'xs2': xs2, 'ys1': ys1, 'ys2': ys2, 'zs1': zs1, 'zs2': zs2, 'xGrid1': xGrid1, 'yGrid1': yGrid1, 'zGrid1': zGrid1}
        
        xs_all.append(xs1)
        ys_all.append(ys1)
        zs_all.append(zs1)
    
    # Make abs mesh
    unique_xs = np.unique(xs_all)
    unique_ys = np.unique(ys_all)
    unique_zs = np.unique(zs_all)
    
    xGrid1f, yGrid1f, zGrid1f = np.meshgrid(unique_xs, unique_ys, unique_zs, indexing='ij')
    Tf = np.zeros_like(xGrid1f) + np.nan
    Sf = np.zeros((xGrid1f.shape[0], xGrid1f.shape[1], xGrid1f.shape[2], S.shape[3])) + np.nan
    Uf = np.zeros((xGrid1f.shape[0], xGrid1f.shape[1], xGrid1f.shape[2], U.shape[3])) + np.nan
    
    for i in range(0, len(unique_meshes)):
        xGrid1 = coarse_results[i]['xGrid1']
        yGrid1 = coarse_results[i]['yGrid1']
        zGrid1 = coarse_results[i]['zGrid1']
        xs1 = coarse_results[i]['xs1']
        xs2 = coarse_results[i]['xs2']
        ys1 = coarse_results[i]['ys1']
        ys2 = coarse_results[i]['ys2']
        zs1 = coarse_results[i]['zs1']
        zs2 = coarse_results[i]['zs2']
        T = coarse_results[i]['T']
        U = coarse_results[i]['U']
        S = coarse_results[i]['S']
        
        xind1 = np.argmin(abs(unique_xs-xs1[0]))
        xind2 = np.argmin(abs(unique_xs-xs1[-1]))
        yind1 = np.argmin(abs(unique_ys-ys1[0]))
        yind2 = np.argmin(abs(unique_ys-ys1[-1]))
        zind1 = np.argmin(abs(unique_zs-zs1[0]))
        zind2 = np.argmin(abs(unique_zs-zs1[-1]))
        
        Tf[xind1:xind2, yind1:yind2, zind1:zind2] = T
        Uf[xind1:xind2+1, yind1:yind2+1, zind1:zind2+1, :] = U
        Sf[xind1:xind2, yind1:yind2, zind1:zind2, :] = S
    
    return unique_xs, unique_ys, unique_zs, xGrid1f, yGrid1f, zGrid1f, Tf, Uf, Sf
    

def extract_contour(axis, value, data, xs, ys, zs, xGrid1, yGrid1, zGrid1, centered, dataInd):
    if axis == 1:
        ind = np.argmin(abs(xs - value))
        if dataInd == None:
            dd = data[ind, :, :]
        else:
            dd = data[ind, :, :, dataInd]
        if centered:
            xx = (yGrid1[ind, 1:, 1:] + yGrid1[ind, :-1, :-1])/2
            yy = (zGrid1[ind, 1:, 1:] + zGrid1[ind, :-1, :-1])/2
            dd = dd[:-1,:-1]
        else:            
            xx = yGrid1[ind, :, :]
            yy = zGrid1[ind, :, :]
    elif axis == 2:
        ind = np.argmin(abs(ys - value))
        if dataInd == None:
            dd = data[:, ind, :]
        else:
            dd = data[:, ind, :, dataInd]
        if centered:
            xx = (xGrid1[1:, ind, 1:] + xGrid1[:-1, ind, :-1])/2
            yy = (zGrid1[1:, ind, 1:] + zGrid1[:-1, ind, :-1])/2
            dd = dd[:-1,:-1]
        else:            
            xx = xGrid1[:, ind, :]
            yy = zGrid1[:, ind, :]
    elif axis == 3:
        ind = np.argmin(abs(zs - value))
        if dataInd == None:
            dd = data[:, :, ind]
        else:
            dd = data[:, :, ind, dataInd]
        if centered:
            xx = (yGrid1[1:, 1:, ind] + yGrid1[:-1, :-1, ind])/2
            yy = (zGrid1[1:, 1:, ind] + zGrid1[:-1, :-1, ind])/2
            dd = dd[:-1,:-1]
        else:            
            xx = xGrid1[:, :, ind]
            yy = zGrid1[:, :, ind]
    return xx, yy, dd

def fill_centered_ghost_cell(data):
    if len(data.shape) > 3:
        data[-1, :, :, :] = data[-2, :, :, :]
        data[:, -1, :, :] = data[:, -2, :, :]
        data[:, :, -1, :] = data[:, :, -2, :]
    else:
        data[-1, :, :] = data[-2, :, :]
        data[:, -1, :] = data[:, -2, :]
        data[:, :, -1] = data[:, :, -2]
    return data

def get_interp(data, xs, ys, zs):
    if len(data.shape) == 3:
        data = np.reshape(data, (data.shape[0], data.shape[1], data.shape[2], 1))
        
    interps, mns, mxs = [], [], []
    for i in range(0, data.shape[3]):
        interp = scip.RegularGridInterpolator((xs, ys, zs), data[:, :, :, i], bounds_error=False, fill_value=None)
        interps.append(interp)
        mns.append(np.nanmin(data[:, :, :, i]))
        mxs.append(np.nanmax(data[:, :, :, i]))
    return interps, mns, mxs

def get_pts_from_grid(xGrid1, yGrid1, zGrid1):
    xGrid2 = get_center_grid(xGrid1)
    yGrid2 = get_center_grid(yGrid1)
    zGrid2 = get_center_grid(zGrid1)
    pts1 = grid_to_points(xGrid1, yGrid1, zGrid1)
    pts2 = grid_to_points(xGrid2, yGrid2, zGrid2)
    return pts1, pts2

def greedy_upsample(coarse_resultdir, coarse_chid, fine_resultdir, fine_chid, time=-1,
                    show_uvw=False, uvw_axis=1, uvw_value=None, uvw_plot_mn=None, uvw_plot_mx=None, uplot_ind=0,
                    show_T=False, T_axis=1, T_value=None, T_plot_mn=None, T_plot_mx=None,
                    show_s=False, s_axis=1, s_value=None, s_plot_mn=None, s_plot_mx=None, splot_ind=0,
                    xmn=None, xmx=None, ymn=None, ymx=None, zmn=None, zmx=None):
    
    unique_times = find_uvw_times(coarse_resultdir, coarse_chid)
    unique_meshes = find_uvw_meshes(coarse_resultdir, coarse_chid)
    
    if time == -1:
        time1 = unique_times[-1]
    
    if time == None:
        time1 = unique_times[-1]
    print("Reading one uvwts file")
    xs, ys, zs, xGrid1, yGrid1, zGrid1, Tf, Uf, Sf = greedy_uvwts_read(coarse_resultdir, coarse_chid, time1)
    
    if time == None:
        numTimes = len(unique_times[:-1])
        for t in range(0, numTimes):
            time2 = unique_times[:-1][t]
            print("Reading time %0.4f uvwts file %0.0f / %0.0f"%(time2, t, numTimes))
            xs, ys, zs, xGrid1, yGrid1, zGrid1, Tf1, Uf1, Sf1 = greedy_uvwts_read(coarse_resultdir, coarse_chid, time2)
            Tf = Tf + Tf1
            Uf = Uf + Uf1
            Sf = Sf + Sf1
        Tf = Tf/len(unique_times)
        Uf = Uf/len(unique_times)
        Sf = Sf/len(unique_times)
        tname = "avg"
    else:
        tname = "%d"%(unique_times[-1])
    print("Finished reading uvwts files. Cleaning data.")
    Sf = fill_centered_ghost_cell(Sf)
    Tf = fill_centered_ghost_cell(Tf)
    
    # Read smokeview file and fill each obst with nan
    smvFilePath = fds.getFileList(coarse_resultdir, coarse_chid, 'smv')[0]
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
    
    if (show_uvw):
        print("Plotting visualization of UVW")
        xx, yy, dd1 = extract_contour(uvw_axis, uvw_value, Uf, xs, ys, zs, xGrid1, yGrid1, zGrid1, False, uplot_ind)
        fig1, ax1 = fds.plotSlice(xx, yy, dd1, uvw_axis, qnty_mn=uvw_plot_mn, qnty_mx=uvw_plot_mx)
    if (show_T):
        print("Plotting visualization of T")
        xx, yy, dd2 = extract_contour(T_axis, T_value, Tf, xs, ys, zs, xGrid1, yGrid1, zGrid1, True, None)
        fig2, ax2 = fds.plotSlice(xx, yy, dd2, T_axis, qnty_mn=T_plot_mn, qnty_mx=T_plot_mx)
    if (show_s):
        print("Plotting visualization of S")
        xx, yy, dd3 = extract_contour(s_axis, s_value, Sf, xs, ys, zs, xGrid1, yGrid1, zGrid1, True, splot_ind)
        fig3, ax3 = fds.plotSlice(xx, yy, dd3, s_axis, qnty_mn=s_plot_mn, qnty_mx=s_plot_mx)
    
    #assert False, "Stopped"
    
    #Tf[Tf < 20.01] = np.nan
    
    print("Finished cleaning data. Building interpolator functions.")
    
    T_interps, T_mns, T_mxs = get_interp(Tf, xs, ys, zs)
    U_interps, U_mns, U_mxs = get_interp(Uf, xs, ys, zs)
    S_interps, S_mns, S_mxs = get_interp(Sf, xs, ys, zs)
    print("Finished building interpolator functions.")
    
    # Loop through fine meshes
    fine_fdsfile = os.path.join(fine_resultdir + os.sep + fine_chid + '.fds')
    fdsFile = fds.fdsFileOperations()
    fdsFile.importFile(fine_fdsfile)
    
    mesh_keys = fdsFile.meshes
    mesh_keys.pop('unknownCounter')
    
    # count meshes for printing
    total_meshes = 0
    for key in mesh_keys:
        MULT_ID = fdsFile.meshes[key]["MULT_ID"]
        if MULT_ID is not False:
            mult = fdsFile.mult[MULT_ID]
            I_LOWER, J_LOWER, K_LOWER = 0, 0, 0
            if mult['I_LOWER'] is not False: I_LOWER = mult['I_LOWER']
            if mult['J_LOWER'] is not False: J_LOWER = mult['J_LOWER']
            if mult['K_LOWER'] is not False: K_LOWER = mult['K_LOWER']
            I_UPPER, J_UPPER, K_UPPER = 0, 0, 0
            if mult['I_UPPER'] is not False: I_UPPER = mult['I_UPPER']
            if mult['J_UPPER'] is not False: J_UPPER = mult['J_UPPER']
            if mult['K_UPPER'] is not False: K_UPPER = mult['K_UPPER']
            DX, DY, DZ = 0, 0, 0
            if mult['DX'] is not False: DX = mult['DX']
            if mult['DY'] is not False: DY = mult['DY']
            if mult['DZ'] is not False: DZ = mult['DZ']
            
            for K in range(K_LOWER, K_UPPER+1):
                for J in range(J_LOWER, J_UPPER+1):
                    for I in range(I_LOWER, I_UPPER+1):
                        total_meshes += 1
        else:
            total_meshes += 1
    
    mesh_counter = 1
    for key in mesh_keys:
        IJK = fdsFile.meshes[key]["IJK"]
        XB = fdsFile.meshes[key]["XB"]
        MULT_ID = fdsFile.meshes[key]["MULT_ID"]
        xs = np.linspace(XB[0], XB[1], IJK[0]+1)
        ys = np.linspace(XB[2], XB[3], IJK[1]+1)
        zs = np.linspace(XB[4], XB[5], IJK[2]+1)
        
        xGrid_fine, yGrid_fine, zGrid_fine = np.meshgrid(xs, ys, zs, indexing='ij')
        
        if MULT_ID is not False:
            mult = fdsFile.mult[MULT_ID]
            I_LOWER, J_LOWER, K_LOWER = 0, 0, 0
            if mult['I_LOWER'] is not False: I_LOWER = mult['I_LOWER']
            if mult['J_LOWER'] is not False: J_LOWER = mult['J_LOWER']
            if mult['K_LOWER'] is not False: K_LOWER = mult['K_LOWER']
            I_UPPER, J_UPPER, K_UPPER = 0, 0, 0
            if mult['I_UPPER'] is not False: I_UPPER = mult['I_UPPER']
            if mult['J_UPPER'] is not False: J_UPPER = mult['J_UPPER']
            if mult['K_UPPER'] is not False: K_UPPER = mult['K_UPPER']
            DX, DY, DZ = 0, 0, 0
            if mult['DX'] is not False: DX = mult['DX']
            if mult['DY'] is not False: DY = mult['DY']
            if mult['DZ'] is not False: DZ = mult['DZ']
            
            for K in range(K_LOWER, K_UPPER+1):
                for J in range(J_LOWER, J_UPPER+1):
                    for I in range(I_LOWER, I_UPPER+1):
                        print("Writing mesh number %0.0f/%0.0f."%(mesh_counter, total_meshes))
                        # Build fine grid
                        xGrid_fine1 = np.copy(xGrid_fine) + DX*I
                        yGrid_fine1 = np.copy(yGrid_fine) + DY*J
                        zGrid_fine1 = np.copy(zGrid_fine) + DZ*K
                        fine_pts1, fine_pts2 = get_pts_from_grid(xGrid_fine1, yGrid_fine1, zGrid_fine1)
                        xGrid_fine2 = get_center_grid(xGrid_fine1)
                        
                        # Write fine initialization files
                        namespace = fine_resultdir + os.sep + fine_chid + '_#VARIABLE#_t%s_m%d.csv'%(tname, mesh_counter)
                        
                        # Write temperature
                        fine_T_flat = interpolate(T_interps, fine_pts2, T_mns, T_mxs)
                        write_csv_data(fine_T_flat, namespace, 'tmp', xGrid_fine2, True)
                        
                        # Write velocity
                        fine_U_flat = interpolate(U_interps, fine_pts1, U_mns, U_mxs)
                        write_csv_data(fine_U_flat, namespace, 'uvw', xGrid_fine2, False)
                        
                        # Write species
                        fine_S_flat = interpolate(S_interps, fine_pts2, S_mns, S_mxs)
                        write_csv_data(fine_S_flat, namespace, 'spec', xGrid_fine2, True)
                        su = np.nansum(fine_S_flat,axis=1)
                        #print(su.max(), su.min())
                        mesh_counter += 1
        else:
            fine_pts1, fine_pts2 = get_pts_from_grid(xGrid_fine, yGrid_fine, zGrid_fine)
            xGrid_fine2 = get_center_grid(xGrid_fine)
            
            # Write fine initialization files
            namespace = fine_resultdir + os.sep + fine_chid + '_#VARIABLE#_t%d_m%d.csv'%(unique_times[-1], mesh_counter)
            
            # Write temperature
            fine_T_flat = interpolate(T_interps, fine_pts2, T_mns, T_mxs)
            write_csv_data(fine_T_flat, namespace, 'tmp', xGrid_fine2, True)
            
            # Write velocity
            fine_U_flat = interpolate(U_interps, fine_pts1, U_mns, U_mxs)
            write_csv_data(fine_U_flat, namespace, 'uvw', xGrid_fine2, False)
            
            # Write species
            fine_S_flat = interpolate(S_interps, fine_pts2, S_mns, S_mxs)
            write_csv_data(fine_S_flat, namespace, 'spec', xGrid_fine2, True)
            
            mesh_counter += 1


if __name__ == "__main__":
    
    coarse_resultdir = "50kw_coarse"
    coarse_chid = "50kw_les_methane_1cm_cold"
    
    fine_resultdir = "50kw_fine"
    fine_chid = "50kw_dns_methane_2mm_1024core"

    coarse_resultdir = "E:\\1JLH-NIST2023\\marcos_plume\\Thermo_Only_0p5cm_Jon\\"
    coarse_chid = "test"

    
    fine_resultdir = "E:\\1JLH-NIST2023\\marcos_plume\\Thermo_Only_0p5cm_Jon\\fine_grid2\\"
    fine_chid = "test"
    
    
    greedy_upsample(coarse_resultdir, coarse_chid, fine_resultdir, fine_chid, time=None,
                    show_uvw=True, uvw_axis=2, uvw_value=0.15, uvw_plot_mn=None, uvw_plot_mx=None, uplot_ind=2,
                    show_T=True, T_axis=2, T_value=0.15, T_plot_mn=None, T_plot_mx=None,
                    show_s=True, s_axis=2, s_value=0.15, s_plot_mn=None, s_plot_mx=None, splot_ind=0)
    
    