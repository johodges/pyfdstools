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
import os, shutil
from collections import defaultdict
from .utilities import zopen, getFileList
from .colorSchemes import buildSMVcolormap
from .smokeviewParser import parseSMVFile

def getBingeoms(resultDir, chid):
    """Reads the binary geometry (.bingeom) files written for a case

    Parameters
    ----------
    resultDir : str
        Directory containing the FDS results, or a zip archive
    chid : str
        FDS CHID of the case
    """

    bingeomFiles = getFileList(resultDir, chid, 'bingeom')
    for bingeomFile in bingeomFiles:
        f = zopen(bingeomFile, 'rb')
        data = f.read()
        f.close()
        

def getBndeQuantities(smvFile):
    """Lists the boundary element files a case wrote and what is in them

    Boundary element (.be) files hold boundary data on unstructured
    &GEOM surfaces, as opposed to the .bf files used for rectangular
    obstructions.

    Parameters
    ----------
    smvFile : str
        Path to the case's smokeview file, or to one inside a zip
        archive

    Returns
    -------
    defaultdict
        Dictionary keyed by boundary element file name, each holding
        'quantity', 'gridfile' (the geometry file the data map onto),
        'datafile', 'mesh' and 'var'
    """

    smvData = parseSMVFile(smvFile)
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    quantities = defaultdict(bool)
    file_quantities = defaultdict(bool)
    for bnde in bndes:
        if quantities[bnde[3]] is False:
            quantities[bnde[3]] = [bnde[1]]
        else:
            quantities[bnde[3]].append(bnde[1])
        file_quantities[bnde[1]] = defaultdict(bool)
        file_quantities[bnde[1]]['quantity'] = bnde[3]
        file_quantities[bnde[1]]['gridfile'] = bnde[2]
        file_quantities[bnde[1]]['datafile'] = bnde[1]
        file_quantities[bnde[1]]['mesh'] = bnde[0]
        file_quantities[bnde[1]]['var'] = bnde[4]
    return file_quantities

def readGcfFile(gcffile):
    """Reads a geometry connectivity (.gcf) file

    The gcf file holds the vertices and faces of the &GEOM surfaces a
    case defines, which the boundary element data are mapped onto.

    Parameters
    ----------
    gcffile : str
        Path to a gcf file, or to one inside a zip archive

    Returns
    -------
    array(NV, 3)
        Vertex coordinates
    array(NF, 3)
        One-based vertex indices of each triangular face
    array(NF)
        Surface index of each face
    """

    f = zopen(gcffile, 'rb')
    data = f.read()
    f.close()
    header = np.frombuffer(data[:76], dtype=np.int32)
    geom_type = header[1]
    formatVersion = header[4]
    firstFrameStatic = header[9]
    stime = header[12]
    nverts = header[15]
    nfaces = header[16]
    nvols = header[17]
    
    base = 80
    verts = np.frombuffer(data[base:base+nverts*3*4], dtype=np.float32)
    base = base + 8 + nverts*3*4
    faces = np.frombuffer(data[base:base+nfaces*3*4], dtype=np.int32)
    base = base + 8 + nfaces*3*4
    locations = np.frombuffer(data[base:base+nfaces*3*4], dtype=np.int32)
    base = base + 8 + nfaces*3*4
    float_placeholder1 = np.frombuffer(data[base:base+nfaces*3*4], dtype=np.float32)
    #base = base + 8 + nfaces*3*4
    #float_placeholder2 = np.frombuffer(data[base:base+nfaces*3*4], dtype=np.float32)
    
    (x, y, z) = (verts[::3], verts[1::3], verts[2::3])
    vertices = np.array([x, y, z]).T
    faces2 = faces.reshape((int(faces.shape[0]/3),3))
    return vertices, faces2, header

def readBeFile(befile):
    """Reads a boundary element (.be) file

    Parameters
    ----------
    befile : str
        Path to a be file, or to one inside a zip archive

    Returns
    -------
    array(NT)
        Timestamps of each frame
    array(NV, NT)
        Value at each geometry vertex for each frame
    array
        Raw integer header, which writeBeFile needs to write the file
        back out
    """

    f = zopen(befile, 'rb')
    data = f.read()
    f.close()
    header = np.frombuffer(data[:64], dtype=np.int32)
    one = header[1]
    version = header[4]
    stime = header[7]
    n_vert_s_vals = header[13]
    #n_vert_d_vals = header[11]
    n_face_s_vals = header[15]
    #n_face_d_vals = header[13]
    d = np.frombuffer(data[64:], dtype=np.float32)
    times =d[n_vert_s_vals+2::n_vert_s_vals+11]
    times = np.insert(times, stime, 0)
    
    vals = np.zeros((n_vert_s_vals, times.shape[0]), dtype=np.float32)
    loc = 0
    for i in range(0, times.shape[0]):
        vals[:, i] = d[loc:loc+n_vert_s_vals]
        loc = loc + n_vert_s_vals + 11
    
    return times, vals, header

def writeBeFile_old(out_file, header, out_data):
    """Writes a boundary element file in the pre-time-series layout

    .. deprecated::
        Use :func:`writeBeFile`, which writes the timestamps alongside
        the data. Kept for reading back files written by earlier
        versions of this package.

    Parameters
    ----------
    out_file : str
        Path the be file is written to
    header : array
        Raw integer header, as returned by readBeFile
    out_data : array(NV, NT)
        Value at each geometry vertex for each frame
    """

    outarray = header.tobytes()
    #outarray = outarray + np.array([np.min(out_data[:, 0]), np.max(out_data[:, 0])], dtype=np.float32).tobytes()
    #outarray = outarray + np.array([27, 27], dtype=np.float32).tobytes()
    outarray = outarray + out_data[:, 0].tobytes()
    #outarray = outarray + np.array(header[-1], dtype=np.int32).tobytes()
    outarray = outarray + b'\xc02\x00\x00\x04\x00\x00\x00' #np.array([0., 0.], dtype=np.float32).tobytes()
    for i in range(1, out_data.shape[1]):
        outarray = outarray + times[i].tobytes()
        outarray = outarray + b'\x04\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xb0\x0c\x00\x00\x10\x00\x00\x00\xc02\x00\x00'
        #outarray = outarray + np.array([0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float32).tobytes()
        outarray = outarray + out_data[:, i].tobytes()
        outarray = outarray + b'\xc02\x00\x00\x04\x00\x00\x00'
    outarray = outarray[:-4]
    
    with open(out_file, 'wb') as f:
        f.write(outarray)


def writeBeFile(out_file, header, out_data, times):
    """Writes a boundary element (.be) file smokeview can read

    Round-trips the header returned by readBeFile, so the usual way to
    build one is to read an existing be file, replace its values, and
    write it back out under a new name.

    Parameters
    ----------
    out_file : str
        Path the be file is written to
    header : array
        Raw integer header, as returned by readBeFile
    out_data : array(NV, NT)
        Value at each geometry vertex for each frame
    times : array(NT)
        Timestamps of each frame, as float32 so that they can be written
        directly into the record
    """

    outarray = header.tobytes()
    #outarray = outarray + np.array([np.min(out_data[:, 0]), np.max(out_data[:, 0])], dtype=np.float32).tobytes()
    #outarray = outarray + np.array([27, 27], dtype=np.float32).tobytes()
    outarray = outarray + out_data[:, 0].tobytes()
    outarray = outarray + np.array([header[-1],4], dtype=np.int32).tobytes()
    #outarray = outarray + b'\xc02\x00\x00\x04\x00\x00\x00' #np.array([0., 0.], dtype=np.float32).tobytes()
    for i in range(1, out_data.shape[1]):
        outarray = outarray + times[i].tobytes()
        outarray = outarray + b'\x04\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'
        outarray = outarray + np.array([header[-3]], dtype=np.int32).tobytes()
        outarray = outarray + b'\x10\x00\x00\x00'
        outarray = outarray + np.array([header[-1]], dtype=np.int32).tobytes()
        #outarray = outarray + np.array([0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float32).tobytes()
        outarray = outarray + out_data[:, i].tobytes()
        outarray = outarray + np.array([header[-1]], dtype=np.int32).tobytes() + b'\x04\x00\x00\x00'
        #b'\xc02\x00\x00\x04\x00\x00\x00'
    outarray = outarray[:-4]
    
    with open(out_file, 'wb') as f:
        f.write(outarray)

def appendNewBeFileToSMV(bnde, smvFile):
    """Registers a new boundary element file in a smokeview file

    .. warning::
        Not implemented. Writing a derived .be file is only half the
        job: smokeview will not display it until a matching BNDE record
        names it in the smokeview file. Until this is written, add that
        record by hand, or follow the pattern in
        extractBoundaryData.bndfsTimeAverage, which does the equivalent
        for rectangular boundary files.

    Parameters
    ----------
    bnde : list
        Boundary element record to add, in the layout parseBNDE returns:
        [mesh, data file, geometry file, quantity, variable number]
    smvFile : str
        Path to the smokeview file to add the record to
    """

    pass

if __name__ == "__main__":
    import trimesh
    from vedo import trimesh2vedo, show, Plotter
    resultDir = "E:\\projects\\1JLH-NIST2022\\geom_testing\\radiant_panel\\radiant_panel_plywood_geom_rot2\\"
    chid = "radiant_panel_plywood_geom_rot_cc"
    qty = 'WALL TEMPERATURE'
    qty = 'INCIDENT HEAT FLUX'
    
    resultDir = "E:\\projects\\1JLH-NIST2022\\geom_testing\\lattimer_tilted_wall_v1\\Lattimer_50_kW_00_degree\\"
    chid = "Lattimer_50_kW_00_degree"
    qty = 'GAUGE HEAT FLUX'
    
    #resultDir = "E:\\projects\\1JLH-NIST2022\\geom_testing\\radiant_panel\\radiant_panel_plywood_geom\\"
    #chid = "radiant_panel_plywood_geom"
    
    smvFiles = getFileList(resultDir, chid, 'smv')
    smvFile = smvFiles[0]
    outSmvFile = smvFile.replace('.smv','_mod.smv')
    shutil.copyfile(smvFile, outSmvFile)
    file_quantities = getBndeQuantities(smvFiles[0])
    
    gbfFiles = getFileList(resultDir, chid, 'gbf')
    beFiles = getFileList(resultDir, chid, 'be')
    allMeshes = []
    allValues = []
    for befile in beFiles:
        fname = os.path.abspath(befile).split(os.sep)[-1]
        if file_quantities[fname] is False:
            print("Warning %s not found in SMV file, may be expected if custom files generated already."%(fname))
            continue
        file_quantity = file_quantities[fname]['quantity']
        if file_quantity != qty:
            continue
        gname = file_quantities[fname]['gridfile']
        gcffile = befile.replace(fname, gname)
        
        vertices, faces2, header1 = readGcfFile(gcffile)
        times, vals, header2 = readBeFile(befile)
        
        '''
        arrival_times = np.zeros((vals.shape[0],times.shape[0]), dtype=np.float32)
        for i in range(0, arrival_times.shape[0]):
            inds = np.where(vals[i, :] > 380)[0]
            if len(inds) > 0:
                arrival_times[i, inds[0]:] = times[inds[0]]
                #print(times[inds[0]])
            else:
                pass
        '''
        #assert False, 'Stopped'
        
        #out_data = arrival_times
        out_data = vals
        out_file = befile.replace('.be','00.be')
        out_mesh = file_quantities[fname]['mesh']
        out_var = file_quantities[fname]['var']
        writeBeFile(out_file, header2, out_data, times)
        
        txt = 'BNDE  %s     %s\n'%(str(int(out_mesh)).rjust(4), str(int(out_var)).rjust(4))
        txt = txt + ' ' + fname.replace('.be','00.be') + '\n'
        txt = txt + ' ' + os.path.abspath(gcffile).split(os.sep)[-1] + '\n'
        txt = txt + ' ' + 'Ignition Time\n'
        txt = txt + ' ' + 't_ign\n'
        txt = txt + ' ' + 's\n\n'
        
        with open(outSmvFile, 'a') as f:
            f.write(txt)
        
        #print(arrival_times)
        print(befile, np.max(vals))
        
        #assert False, "Stopped"
        
        plt = Plotter(N=1, axes=3)
        mesh = trimesh.Trimesh(vertices=vertices,
                               faces=faces2-1,
                               process=False) 
        
        mesh2 = mesh.slice_plane(plane_origin=(0.0,0.1,0.0), plane_normal=(1.0,0.1,0.0))
        print(mesh2)
        
        #point_out = out_data[vertices[:, 0] > 0.0, -1]
        
        '''
        for facet in mesh.facets:
            mesh.visual.face_colors[facet] 
        mesh.show()
        '''
        allMeshes.append(mesh)
        allValues.append(out_data[:, -1])
    vmeshes = trimesh2vedo(allMeshes)
    
    vmin=4.5
    vmax=25.5
    n_colors=15
    
    percentile = None
    highlightWidth = None
    cmap = buildSMVcolormap(percentile=percentile, width=highlightWidth)
    cmaps = ('jet', 'PuOr', 'viridis')
    #scals = vmeshes[0].points()[:, 0]
    
    for i in range(0, len(allMeshes)):
        scals = vmeshes[i].cell_centers()[:, 0] + allValues[i] #+ vals[:, -2]
        vmeshes[i].cmap(cmap, scals, on='cells', n_colors=n_colors, vmin=vmin, vmax=vmax)
    pt = (0.475,0.1,1.1)
    vmeshes[0].add_scalarbar3d(s=[None,0.4], pos=pt, nlabels=n_colors+1)
    vmeshes[0].scalarbar.rotate_x(90, around=pt)
    #vmeshes[0].scalarbar.rotate_x(90).y(0.2)
        
    # add a 2D scalar bar to a mesh
    #vmeshes[0].addScalarBar(title="my scalarbar\nnumber #0", c='k')
    #vmeshes[0].add_scalarbar(title="my scalarbar\nnumber #0", c='k')
    
    # add 3D scalar bars
    #vmeshes[0].addScalarBar3D(pos=(1.2,2.1,-1.0), c='k')
    #vmeshes[0].add_scalarbar3d(pos=(1.2,2.1,-1.0), c='k')
    
    show(vmeshes, axes=1)
    plt.interactive().close()