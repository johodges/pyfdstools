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
import sys
print(sys.version_info, flush=True)

import pyfdstools as fds
import os
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import bpy


def exampleImportFile(fdsPath=None):
    if fdsPath is None:
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        fdsPath = os.path.join(systemPath, "examples", "case001.fds")
    fdsFile = fds.fdsFileOperations()
    print("\tStarting to import file from %s"%(fdsPath))
    fdsFile.importFile(fdsPath)
    print("\tFinished importing file.")
    return fdsFile


def exampleSaveFile(file=None, outdir=None, outname=None):
    if file is None:
        file = exampleImportFile()
    if outdir is None:
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        outdir = os.path.join(systemPath, "generated")
    if outname is None:
        location = os.path.join(outdir, '%s.fds'%(file.head['ID']['CHID']))
    print("\tStarting to save model.")
    file.saveModel(1, location, allowMeshSplitting=False)
    print("\tFinished saving model.")
    return location


def exampleErrorCalculation():
    data, keys = fds.readErrorTable()
    
    errorvalues = fds.calculatePercentile([100, 200, 300, 400, 500, 600], 'Plume Temperature', 0.95, fdsVersion='6.7.1')
    fig, ax1 = fds.plotPercentile(500, 'Plume Temperature', fdsVersion='6.7.1')
    return errorvalues


def exampleReadSlcf2dResults(resultDir=None, chid=None,
                             fdsQuantities = ['TEMPERATURE'],
                             tStart=0, tEnd=120):
    if (resultDir is None) and (chid is None):
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    datas = defaultdict(bool)
    
    for qty in fdsQuantities:
        grid, data, times = fds.readSLCF2Ddata(chid, resultDir, qty)
        tStartInd = np.argwhere(times >= tStart)[0][0]
        tEndInd = np.argwhere(times <= tEnd)[-1][0]
        data_tavg = np.nanmean(data[:, :, :, tStartInd:tEndInd], axis=3)
        datas[qty] = data_tavg.copy()
    datas['GRID'] = grid
    datas['TIMES'] = times
    
    return datas


def exampleExtract2dFromSlcf2d(datas,
                               fdsQuantities = ['TEMPERATURE'],
                               fdsUnits = ['C'],
                               axis=1, value=2.5,
                               qnty_mn=20, qnty_mx=150):
    datas2D = defaultdict(bool)    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        x, z, data_slc = fds.findSliceLocation(datas['GRID'], datas[qty], axis, value)
        datas2D[qty] = data_slc
        fig = fds.plotSlice(x, z, data_slc, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    datas2D['X'] = x
    datas2D['Z'] = z
    
    return datas2D, fig

def exampleReadSlcf3dResults(resultDir=None, chid=None,
                             fdsQuantities = ['TEMPERATURE'],
                             tStart=0, tEnd=120):
    if (resultDir is None) and (chid is None):
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    datas = defaultdict(bool)
    
    for qty in fdsQuantities:
        grid, data, times = fds.readSLCF3Ddata(chid, resultDir, qty)
        tStartInd = np.argwhere(times >= tStart)[0][0]
        tEndInd = np.argwhere(times <= tEnd)[-1][0]
        data_tavg = np.nanmean(data[:, :, :, tStartInd:tEndInd], axis=3)
        datas[qty] = data_tavg.copy()
    datas['GRID'] = grid
    datas['TIMES'] = times
    
    return datas

def exampleExtract2dFromSlcf3d(datas,
                               fdsQuantities = ['TEMPERATURE'],
                               fdsUnits = ['C'],
                               axis=2, value=4.4,
                               qnty_mn=20, qnty_mx=150):
    datas2D = defaultdict(bool)    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        x, z, data_slc = fds.findSliceLocation(datas['GRID'], datas[qty], axis, value)
        datas2D[qty] = data_slc
        fig = fds.plotSlice(x, z, data_slc, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    datas2D['X'] = x
    datas2D['Z'] = z
    
    return datas2D, fig

def exampleImportBndf(resultDir=None, chid=None, fdsFile=None,
                       fdsQuantities = ['WALL TEMPERATURE'],
                       fdsUnits = ['C'],
                       tStart=0, tEnd=120,
                       axis=-2, value=4.4,
                       qnty_mn=20, qnty_mx=100):
    if (resultDir is None) and (chid is None):
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    if fdsFile is None: fdsFilePath = fds.getFileList(resultDir, chid, 'fds')[0]
    
    datas, times = fds.queryBndf(resultDir, chid, fdsFilePath, fdsQuantities, fdsUnits, axis, value)
    tStartInd = np.argwhere(times >= tStart)[0][0]
    tEndInd = np.argwhere(times <= tEnd)[-1][0]
    
    for qty, unit in zip(fdsQuantities, fdsUnits):
        for mesh in list(datas[qty].keys()):
            data = datas[qty]['DATA']
            x = datas[qty]['X']
            z = datas[qty]['Z']
        meanData = np.mean(data[:, :, tStartInd:tEndInd], axis=2)
        fig = fds.plotSlice(x, z, meanData, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=[20, 40, 60, 80, 100, 120, 140])
    return datas, fig

def exampleExtractBndfMax(resultDir=None, chid=None, outDir=None,
                          fdsFile=None, smvFile=None, outputNamespace=None,
                          fdsQuantities = ['WALL TEMPERATURE'], fdsUnits = ['C'],
                          tStart=0, tEnd=120, tInt=1, tBand=3, orientations=[0],
                          axis=-2, value=4.4,
                          qnty_mn=20, qnty_mx=100,
                          yticks=[20, 50, 100, 150, 200, 250, 300, 350],
                          plot=True):
    if (resultDir is None) and (chid is None) and (outDir is None):
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
        outDir = os.path.join(systemPath, "generated")
    
    if fdsFile is None: fdsFilePath = fds.getFileList(resultDir, chid, 'fds')[0]
    if smvFile is None: smvFilePath = fds.getFileList(resultDir, chid, 'smv')[0]
    if outputNamespace is None: outputNamespace = "%s_max_"%(chid)
    
    outputName = os.path.join(outDir, outputNamespace)
    
    datas = fds.extractMaxBndfValues(fdsFilePath, smvFilePath, resultDir, chid, fdsQuantities,
                    tStart=tStart, tEnd=tEnd, tInt=tInt, tBand=tBand, orientations=orientations)
    
    figs = []
    for qty in fdsQuantities:
        times = datas[qty]['TIMES']
        mPts = datas[qty]['DATA']
        names = datas[qty]['NAMES']
        outName = "%s%s"%(outputName, qty)
        fds.maxValueCSV(times, mPts, names, outName)
        if plot:
            fig = fds.maxValuePlot(times, mPts, names, outName, vName=qty, yticks=yticks)
            figs.append(fig)
    return datas, figs


def example2dSliceToCsv(resultDir=None, outDir=None, chid=None,
                        quantity=None, unit=None, axis=None, value=None,
                        time=None, dt=None):
    if (resultDir is None) and (chid is None) and (outDir is None):
        try:
            import bpy
            systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
        except:
            systemPath = os.path.dirname(os.path.abspath(__file__))
        chid = "case001"        
        resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
        outDir = os.path.join(systemPath, "generated")
    try:
        os.mkdir(outDir)
    except:
        pass
    data = fds.query2dAxisValue(resultDir, chid, quantity, axis, value, time=time, dt=dt)
    fds.renderSliceCsvs(data, chid, outDir)
    fig, ax = fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], axis,
                        clabel="%s (%s)"%(quantity, unit))
    fig.savefig(os.path.join(outDir, '%s_%s_%0.0f_%0.4f_final_frame.png'%(chid, quantity, axis, value)))
    return data, fig


def runExamples():
    try:
        import bpy
        systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
    except:
        systemPath = os.path.dirname(os.path.abspath(__file__))
    exampleInputFdsFile = os.path.join(systemPath, "examples", "case001.fds")
    exampleOutputDir = os.path.join(systemPath, "generated")
    chid = "case001"
    resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    print("Importing model example", flush=True)
    file = exampleImportFile(exampleInputFdsFile)
    
    
    print("Saving model example", flush=True)
    exampleSaveFile(file=file, outdir=exampleOutputDir)
    
    
    print("Plotting error example", flush=True)
    exampleErrorCalculation()
    
    print("Importing SLCF2D results example", flush=True)
    datas = exampleReadSlcf2dResults(resultDir=resultDir, chid=chid)
    
    print("Extracting slice from SLCF2D results example", flush=True)
    datas2D, fig = exampleExtract2dFromSlcf2d(datas)
    
    print("Importing SLCF3D results example", flush=True)
    datas = exampleReadSlcf3dResults(resultDir=resultDir, chid=chid)
    
    print("Extracting 2-D slice from SLCF3D results example", flush=True)
    datas2D, fig = exampleExtract2dFromSlcf3d(datas)
    
    print("Importing BNDF results example", flush=True)
    datas, fig = exampleImportBndf(resultDir=resultDir, chid=chid)
    
    print("Extracting max value from BNDF results example", flush=True)
    datas, figs = exampleExtractBndfMax(resultDir=resultDir, chid=chid, outDir=exampleOutputDir)
    
    print("Rendering 2d slice to csv example.", flush=True)
    datas = example2dSliceToCsv(resultDir=resultDir, chid=chid, outDir=exampleOutputDir,
                                axis=1, value=2.45, time=30, dt=60, quantity='TEMPERATURE', unit='C')
    
    return systemPath, exampleInputFdsFile, exampleOutputDir, chid, resultDir

def create_custom_mesh(objname, px, py, pz):
 
    # Define arrays for holding data    
    myvertex = []
    myfaces = []

    # Create all Vertices

    # vertex 0
    mypoint = [(-1.0, -1.0, 0.0)]
    myvertex.extend(mypoint)

    # vertex 1
    mypoint = [(1.0, -1.0, 0.0)]
    myvertex.extend(mypoint)

    # vertex 2
    mypoint = [(-1.0, 1.0, 0.0)]
    myvertex.extend(mypoint)

    # vertex 3
    mypoint = [(1.0, 1.0, 0.0)]
    myvertex.extend(mypoint)

    # -------------------------------------
    # Create all Faces
    # -------------------------------------
    myface = [(0, 1, 3, 2)]
    myfaces.extend(myface)

    
    mymesh = bpy.data.meshes.new(objname)

    myobject = bpy.data.objects.new(objname, mymesh)

    bpy.context.scene.collection.objects.link(myobject)
    
    # Generate mesh data
    mymesh.from_pydata(myvertex, [], myfaces)
    # Calculate the edges
    mymesh.update(calc_edges=True)

    # Set Location
    myobject.location.x = px
    myobject.location.y = py
    myobject.location.z = pz

    return myobject



def create_custom_mesh2(objname, x0, x1, y0, y1, z0, z1, data):
 
    # Define arrays for holding data    
    myvertex = []
    myfaces = []
    
    # Create all Vertices

    # vertex 0
    mypoint = [(x0, y0, z0)]
    myvertex.extend(mypoint)

    # vertex 1
    mypoint = [(x1, y0, z0)]
    myvertex.extend(mypoint)

    # vertex 2
    mypoint = [(x0, y1, z1)]
    myvertex.extend(mypoint)

    # vertex 3
    mypoint = [(x1, y1, z1)]
    myvertex.extend(mypoint)

    # -------------------------------------
    # Create all Faces
    # -------------------------------------
    myface = [(0, 1, 3, 2)]
    myfaces.extend(myface)

    print("A", objname)
    mymesh = bpy.data.meshes.new(objname)
    print("B")
    myobject = bpy.data.objects.new(objname, mymesh)
    print("C")
    bpy.context.scene.collection.objects.link(myobject)
    print("D")
    # Generate mesh data
    mymesh.from_pydata(myvertex, [], myfaces)
    # Calculate the edges
    print("E")
    mymesh.update(calc_edges=True)

    # Set Location
    #myobject.location.x = px
    #myobject.location.y = py
    #myobject.location.z = pz
    print("F")
    myobject.location.x = 0
    myobject.location.y = 0
    myobject.location.z = 0
    
    tex = bpy.data.textures.new('%s-TEXTURE'%(objname), 'IMAGE')
    print(data.shape)
    mn = np.nanmin(data)
    mx = np.nanmax(data)
    d = (data - mn)/(mx - mn)
    d2 = np.zeros((d.shape[0], d.shape[1], 4))
    d2[:, :, 0] = d
    d2[:, :, 1] = d
    d2[:, :, 2] = d
    d2[:, :, 3] = 1
    print(d2.shape, d2.flatten().shape)
    
    image = bpy.data.images.new("%s-IMAGE"%(objname), width=d.shape[0], height=d.shape[1])
    image.pixels = d2.flatten()
    tex.image = image
    
    mat = bpy.data.materials.new('%s-DATA'%(objname))
    
    #bpy.ops.mesh.primitive_plane_add(size=2, enter_editmode=False, align='WORLD', location=(x0, y0, z0), scale=(1, 1, 1))
    #bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    print('%s-DATA'%(objname), flush=True)
    #bpy.ops.object.active_material.name = '%s-DATA'%(objname)
    #mat.texture_paint_images = [image]
    #mat.texture_slots.add()
    #ts = mat.texture_slots[0]
    #ts.texture = tex
    #ts.texture_coords = 'UV'
    #ts.uv_layer = 'default'
    #mtex.texture_coords = 'UV'
    
    myobject.data.materials = mat
    

    return myobject


if __name__ == '__main__':
    QTPLUGINPATH = os.getenv('QT_PLUGIN_PATH')
    print(QTPLUGINPATH)
    
    try:
        import bpy
        systemPath = os.path.dirname(bpy.context.space_data.text.filepath)
    except:
        systemPath = os.path.dirname(os.path.abspath(__file__))
    exampleInputFdsFile = os.path.join(systemPath, "examples", "case001.fds")
    exampleOutputDir = os.path.join(systemPath, "generated")
    chid = "case001"
    resultDir = os.path.join(systemPath, "examples", "%s.zip"%(chid))
    
    
    #datas, fig = exampleImportBndf(resultDir=resultDir, chid=chid)
    #print(datas['WALL TEMPERATURE'].keys(), flush=True)
    #print(datas.keys(), flush=True)
    
    '''
    datas, figs = exampleExtractBndfMax(resultDir=resultDir, chid=chid, outDir=exampleOutputDir, plot=False)
    times = datas['WALL TEMPERATURE']['TIMES']
    names = datas['WALL TEMPERATURE']['NAMES']
    ds = datas['WALL TEMPERATURE']['DATA']
    
    import pandas as pd
    print(ds.shape, flush=True)
    print(len(times), flush=True)
    print(len(names), flush=True)
    tmp = pd.DataFrame(ds, index=times, columns=names)
    print(tmp, flush=True)
    
    qty = 'WALL TEMPERATURE'
    for key in names:
        bpy.context.scene.collection.objects[key]['%s TIME'%(qty)] = tmp.index
        bpy.context.scene.collection.objects[key]['%s MAX'%(qty)] = tmp[key]
    '''
    
    
    '''
    systemPath, exampleInputFdsFile, exampleOutputDir, chid, resultDir = runExamples()
    '''
    
    '''
    print(bpy.context.scene.collection.objects['COUCH-BACK-XPOS']['WALL TEMPERATURE MAX'][:])
    print(bpy.context.scene.collection.objects['COUCH-BACK-XPOS']['WALL TEMPERATURE TIME'][:])
    '''
    
    
    qty = 'TEMPERATURE'
    axis = 2
    value = 4.4
    obst = 'Centerline'
    
    datas = exampleReadSlcf3dResults(resultDir=resultDir, chid=chid)
    datas2D, fig = exampleExtract2dFromSlcf3d(datas, axis=axis, value=value)
    
    TEMPERATURE = datas2D['TEMPERATURE']
    X = datas2D['X']
    Z = datas2D['Z']
    
    x0 = X[0, 0]
    x1 = X[-1, -1]
    y0 = value
    y1 = value
    z0 = Z[0, 0]
    z1 = Z[-1, -1]
    corners = [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)]
    #create_custom_mesh("Awesome_object", curloc[0], curloc[1], 0)
    objname = "%s_%s"%(obst, qty)
    print(objname)
    print(x0, x1, y0, y1, z0, z1, flush=True)
    print(TEMPERATURE.shape, flush=True)
    
    #from PIL import Image
    #import PIL
    from matplotlib import cm
    import matplotlib
    
    smv = fds.buildSMVcolormap()
    
    colors = smv.colors
    
    numColors = colors.shape[0]
    
    cmin = 20
    cmax = 200
    
    T2 = (TEMPERATURE - cmin)/(cmax - cmin)
    T2[T2 < 0] = 0
    T2[T2 > 1] = 1
    
    T3 = np.uint8(smv(T2)*255)
    print(T3)
    
    TEMPERATURE2 = np.zeros((TEMPERATURE.shape[0], TEMPERATURE.shape[1], 3))
    print(smv, flush=True)
    print(smv.colors, flush=True)
    matplotlib.image.imsave("E:\\projects\\customPythonModules\\pyfdstools\\temp.png", T3)
    #img = Image.fromarray(data)
    inputImages = [{"name":objname, "image":TEMPERATURE}]
    files = [{"name":"E:\\projects\\customPythonModules\\pyfdstools\\temp.png"}]
    
    bpy.ops.import_image.to_plane(files=files, location=(0,0,0), rotation=(0, 0, 0), force_reload=True, offset=False, shader='SHADELESS', relative=False)
#    (-0.537504, 0.853194, 0.784878)
    '''
    bpy.ops.transform.translate(value=(-0.537504, 0.853194, 0.784878), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False)
    '''

    
    #bpy.ops.mesh.primitive_plane_add(size=2, enter_editmode=False, align='WORLD', location=(x0, y0, z0), scale=(1, 1, 1))
    #bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    #myobj = create_custom_mesh2(objname, x0, x1, y0, y1, z0, z1, TEMPERATURE)
    #print(myobj.__dict__)
    #x, z, data_slc = fds.findSliceLocation(datas['GRID'], datas[qty], axis, value)
    
    
    
    
    '''
    print(x[:, 0], x[0, :])
    #corners = [(x[0, 0], value, z[0, 0])]
    #qty = 'TEMPERATURE'
    #grid, data, times = fds.readSLCF3Ddata(chid, resultDir, qty)
    
    #corners = 
    
    print(grid.shape, data.shape, times.shape)
    
    d = np.nanmean(datas['TEMPERATURE'], axis=2)
    import bpy
    context = bpy.context
    
    curloc = bpy.context.scene.cursor.location

    create_custom_mesh("Awesome_object", curloc[0], curloc[1], 0)
    '''