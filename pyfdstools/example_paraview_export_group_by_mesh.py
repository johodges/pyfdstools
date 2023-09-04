import pyfdstools as fds
import os
import numpy as np
import hashlib
import matplotlib.pyplot as plt
import stl
from collections import defaultdict

def obstToStl(resultDir, chid):
    smvFile = fds.getFileListFromResultDir(resultDir, chid, 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    (smvGrids, smvObsts) = (smvData['grids'], smvData['obsts'])
    (smvBndfs, smvSurfs) = (smvData['bndfs'], smvData['surfs'])
    (smvFiles, bndes) = (smvData['files'], smvData['bndes'])
    # Define the 12 triangles composing the cube
    faces = np.array([\
        [0,3,1],
        [1,3,2],
        [0,4,7],
        [0,7,3],
        [4,5,6],
        [4,6,7],
        [5,1,2],
        [5,2,6],
        [2,3,6],
        [3,7,6],
        [0,1,5],
        [0,5,4]])
    meshes = []
    for o in smvObsts:
        (xmn, xmx, ymn, ymx, zmn, zmx) = (o[0], o[1], o[2], o[3], o[4], o[5])
        # Define the 8 vertices of the cube
        vertices = np.array([[xmn, ymn, zmn],
                             [xmx, ymn, zmn],
                             [xmx, ymx, zmn],
                             [xmn, ymx, zmn],
                             [xmn, ymn, zmx],
                             [xmx, ymn, zmx],
                             [xmx, ymx, zmx],
                             [xmn, ymx, zmx]])
        # Create the mesh
        cube = stl.mesh.Mesh(np.zeros(faces.shape[0], dtype=stl.mesh.Mesh.dtype))
        for i, f in enumerate(faces):
            for j in range(3):
                cube.vectors[i][j] = vertices[f[j],:]
        meshes.append(cube)
    combined = stl.mesh.Mesh(np.concatenate([m.data for m in meshes]))
    combined.save(chid+'.stl') #,mode=stl.Mode.ASCII)

def writeVtkHeader(fname, t, ext):
    with open(fname+ext, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="%s" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n'%(t))

def writeVtkGrid(fname, x, y, z, gXB, lXB):
    with open(fname+'.vtr', 'a') as f:
        f.write('<RectilinearGrid WholeExtent="%d %d %d %d %d %d">\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
        f.write('<Piece Extent="%d %d %d %d %d %d">\n'%(lXB[0], lXB[1], lXB[2], lXB[3], lXB[4], lXB[5]))
        f.write('<Coordinates>\n')
        f.write('<DataArray type="Float64" Name="Array X-Axis" format="ascii" RangeMin="%0.1f" RangeMax="%0.1f">\n'%(x.min(), x.max()))
        for i in range(0, x.shape[0]): f.write('%0.1f '%(x[i]))
        f.write('\n')
        f.write('</DataArray>\n')
        f.write('<DataArray type="Float64" Name="Array Y-Axis" format="ascii" RangeMin="%0.1f" RangeMax="%0.1f">\n'%(y.min(), y.max()))
        for i in range(0, y.shape[0]): f.write('%0.1f '%(y[i]))
        f.write('\n')
        f.write('</DataArray>\n')
        f.write('<DataArray type="Float64" Name="Array Z-Axis" format="ascii" RangeMin="%0.1f" RangeMax="%0.1f">\n'%(z.min(), z.max()))
        for i in range(0, z.shape[0]): f.write('%0.1f '%(z[i]))
        f.write('\n')
        f.write('</DataArray>\n')
        f.write('</Coordinates>\n')

def writeVtkScalarsPoints(fname, qtys, datas):
    with open(fname+'.vtr', 'a') as f:
        f.write('<PointData>\n') # Scalars=%s>\n'%(",".join(['"%s"'%(q) for q in qtys])))
        for qty in qtys:
            f.write('<DataArray Name="%s" NumberOfComponents="1" type="Float64" format="ascii">\n'%(qty))
            for i in range(0, datas[qty].shape[0]): f.write('%0.1f '%(datas[qty][i]))
            f.write('\n</DataArray>\n')
        f.write('</PointData>\n')

def writeVtkScalarsCells(fname, qtys, datas):
    with open(fname+'.vtr', 'a') as f:
        f.write('<CellData>\n') # Scalars=%s>\n'%(",".join(['"%s"'%(q) for q in qtys])))
        for qty in qtys:
            f.write('<DataArray Name="%s" NumberOfComponents="1" type="Float64" format="ascii">\n'%(qty))
            for i in range(0, datas[qty].shape[0]): f.write('%0.1f '%(datas[qty][i]))
            f.write('\n</DataArray>\n')
        f.write('</CellData>\n')


def writeVtkTail(fname):
    with open(fname+'.vtr', 'a') as f:
        f.write('</Piece>\n')
        f.write('</RectilinearGrid>\n')
        f.write('</VTKFile>\n')

def writeVtkFile(fname, x, y, z, gXB, lXB, datas, cell=False):
    writeVtkHeader(fname, "RectilinearGrid", ".vtr")
    writeVtkGrid(fname, x, y, z, gXB, lXB)
    qtys = list(datas.keys())
    if cell:
        writeVtkScalarsCells(fname, qtys, datas)
    else:
        writeVtkScalarsPoints(fname, qtys, datas)
    writeVtkTail(fname)
    
def writeVtkTimeSeries(namespace, grid, series_data, times, gXB, lXB, cell=False):
    nx, ny, nz, _ = np.shape(grid)
    x = np.array(grid[:, 0, 0, 0], dtype='float64')
    y = np.array(grid[0, :, 0, 1], dtype='float64')
    z = np.array(grid[0, 0, :, 2], dtype='float64')
    
    for time in times:
        fname = namespace + '_%08d'%(time*1e3)
        datas = dict()
        for qty in list(series_data.keys()):
            tind = np.argmin(abs(series_data[qty]['times'] - time))
            datas[qty] = series_data[qty]['data'][:, :, :, tind].T.flatten()
        writeVtkFile(fname, x, y, z, gXB, lXB, datas, cell=cell)
        
def writeVtkPolyFile(fname, pieces):
    writeVtkHeader(fname, "PolyData", ".vtp")
    
    with open(fname + '.vtp', 'a') as f:
        f.write('<PolyData>\n')
    
    keys = list(pieces[0].keys())
    keys.remove('numberOfPolys')
    keys.remove('numberOfPoints')
    keys.remove('Points')
    keys.remove('Offset')
    keys.remove('Connectivity')
    pquantities = [key for key in keys if 'P-' in key]
    cquantities = [key for key in keys if 'C-' in key]
    for i in range(0, len(pieces)):
        piece = pieces[i]
        numberOfPolys = piece['numberOfPolys']
        numberOfPoints = piece['numberOfPoints']
        points = piece['Points']
        offsets = piece['Offset']
        connectivity = piece['Connectivity']
        
        with open(fname + '.vtp', 'a') as f:
            f.write('<Piece NumberOfPoints="%0d" NumberOfVerts="0" NumberOfLines="0" NumberOfPolys="%d">\n'%(numberOfPoints, numberOfPolys))
            f.write('    <Points>\n')
            f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n        ')
            for point in points:
                f.write('%0.1f '%(point))
            f.write('\n        </DataArray>\n    </Points>\n')
            if len(pquantities) > 0:
                f.write('    <PointData>\n')
                for qty in pquantities:
                    f.write('        <DataArray type="Float32" Name="%s" format="ascii">\n            '%(qty.replace('P-','')))
                    values = piece[qty]
                    for v in values:
                        f.write('%0.1f '%(v))
                    f.write('\n         </DataArray>\n')
                f.write('    </PointData>\n')
            if len(cquantities) > 0:
                f.write('    <CellData>\n')
                for qty in cquantities:
                    f.write('        <DataArray type="Float32" Name="%s" format="ascii">\n            '%(qty.replace('C-','')))
                    values = piece[qty]
                    for v in values:
                        f.write('%0.1f '%(v))
                    f.write('\n         </DataArray>\n')
                f.write('    </CellData>\n')
            f.write('    <Polys>\n')
            f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n            ')
            for c in connectivity:
                f.write('%d '%(c))
            f.write('\n        </DataArray>\n')
            f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n            ')
            for o in offsets:
                f.write('%d '%(o))
            f.write('\n        </DataArray>\n    </Polys>\n</Piece>\n')
    with open(fname + '.vtp', 'a') as f:
        f.write('</PolyData>\n</VTKFile>')
    
    return False

def writeVtkPointFile(fname, pieces):
    writeVtkHeader(fname, "PolyData", ".vtp")
    
    with open(fname + '.vtp', 'a') as f:
        f.write('<PolyData>\n')
    
    classes = list(pieces.keys())
    
    for i in range(0, len(classes)):
        piece = pieces[classes[i]]
        numberOfPoints = piece['numberOfPoints']
        points = piece['Points']
        quantities = list(piece['PointData'].keys())
        
        with open(fname + '.vtp', 'a') as f:
            f.write('<Piece NumberOfPoints="%0d" NumberOfVerts="0" NumberOfLines="0" NumberOfPolys="0">\n'%(numberOfPoints))
            #f.write('<Piece NumberOfPoints="%0d" NumberOfVerts="%0d" NumberOfLines="0" NumberOfPolys="0">\n'%(numberOfPoints, numberOfPoints))
            f.write('    <Points>\n')
            f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n        ')
            for point in points:
                f.write('%0.4f '%(point))
            f.write('\n        </DataArray>\n    </Points>\n')
            f.write('    <PointData>\n')
            for qty in quantities:
                f.write('        <DataArray type="Float32" Name="%s" format="ascii">\n            '%(qty))
                values = piece["PointData"][qty]
                for v in values:
                    f.write('%0.4f '%(v))
                f.write('\n         </DataArray>\n')
            f.write('    </PointData>\n')
            '''
            f.write('    <VertexData>\n')
            for qty in quantities:
                f.write('        <DataArray type="Float32" Name="%s" format="ascii">\n            '%(qty))
                values = piece["PointData"][qty]
                for v in values:
                    f.write('%0.4f '%(v))
                f.write('\n         </DataArray>\n')
            f.write('    </VertexData>\n')
            '''
            f.write('</Piece>\n')
    with open(fname + '.vtp', 'a') as f:
        f.write('</PolyData>\n</VTKFile>')

def writeVtkPolyTimeSeries(namespace, series_data, times):
    quantities = list(series_data.keys())
    qty = quantities[0]
    meshes = list(series_data[qty].keys())
    mesh = meshes[0]
    piece_times = series_data[qty][mesh]['times']
    piece_qty = len(series_data[qty][mesh]['pieces'])
    for time in times:
        fname = namespace + '_%08d'%(time*1e3)
        pieces = []
        tind = np.argmin(abs(piece_times - time))
        for mesh in meshes:
            qty = quantities[0]
            piece_qty = len(series_data[qty][mesh]['pieces'])
            for i in range(0, piece_qty):
                p = dict(series_data[qty][mesh]['pieces'][i])
                if p['PointData'] is not False:
                    p.pop('PointData')
                    for qty in list(series_data.keys()):
                        p['P-'+qty] = series_data[qty][mesh]['pieces'][i]['PointData'][:, tind]
                    pieces.append(p)
                if p['CellData'] is not False:
                    p.pop('CellData')
                    for qty in list(series_data.keys()):
                        p['C-'+qty] = series_data[qty][mesh]['pieces'][i]['CellData'][:, tind]
                    pieces.append(p)
        writeVtkPolyFile(fname, pieces)



def exportSl3dDataToVtk(chid, resultDir, outtimes=None):
    quantities, slcfFiles, dimensions, meshes, centers = fds.readSLCFquantities(chid, resultDir)
    
    uniqueMeshes = list(set(meshes))
    
    xyzFiles = fds.getFileListFromResultDir(resultDir, chid, 'xyz')
    grids = fds.getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = fds.getAbsoluteGrid(grids)
    
    gXB = [0, grids_abs.shape[0]-1, 0, grids_abs.shape[1]-1, 0, grids_abs.shape[2]-1]
    pieces=defaultdict(bool)
    for i in range(0, len(uniqueMeshes)):
        mesh = uniqueMeshes[i]
        xyzFile = False
        namespace = resultDir + os.sep + chid + '_sl3d_grid%04d'%(int(mesh))
        for j in range(0, len(slcfFiles)):
            if mesh != meshes[j]:
                continue
            slcfFile = slcfFiles[j]
            qty = quantities[j]
            
            lims, datas, times = fds.readSingleSlcfFile(slcfFile)
            
            if xyzFile is False:
                xyzFile = '_'.join(slcfFile.split('_')[:-1]) + '.xyz'
                grid, gridHeader = fds.readXYZfile(xyzFile)
                xGrid, yGrid, zGrid = fds.rearrangeGrid(grid)
                
                grid_combo = np.zeros((xGrid.shape[0], xGrid.shape[1], xGrid.shape[2],3))
                grid_combo[:, :, :, 0] = xGrid
                grid_combo[:, :, :, 1] = yGrid
                grid_combo[:, :, :, 2] = zGrid
                series_data = defaultdict(bool)
            if series_data[qty] is False:
                series_data[qty] = defaultdict(bool)
                series_data[qty]['data'] = datas
                series_data[qty]['times'] = times
                
                lims_xb = [xGrid[lims[0], 0, 0], xGrid[lims[1], 0, 0],
                           yGrid[0, lims[2], 0], yGrid[0, lims[3], 0],
                           zGrid[0, 0, lims[4]], zGrid[0, 0, lims[5]]]
                lXB     = [np.argmin(abs(grids_abs[:, 0, 0, 0]-lims_xb[0])),
                            np.argmin(abs(grids_abs[:, 0, 0, 0]-lims_xb[1])),
                            np.argmin(abs(grids_abs[0, :, 0, 1]-lims_xb[2])),
                            np.argmin(abs(grids_abs[0, :, 0, 1]-lims_xb[3])),
                            np.argmin(abs(grids_abs[0, 0, :, 2]-lims_xb[4])),
                            np.argmin(abs(grids_abs[0, 0, :, 2]-lims_xb[5]))]
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        writeVtkTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB)
        pieces[namespace] = defaultdict(bool)
        pieces[namespace]['source'] = namespace
        pieces[namespace]['times'] = times2
        pieces[namespace]['lims'] = lXB
    
    wrapper_namespace = resultDir + os.sep + chid + '_sl3d_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtr', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PRectilinearGrid">\n')
            f.write('    <PRectilinearGrid WholeExtent="%d %d %d %d %d %d" GhostLevel="0">\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
            f.write('        <PPointData>\n')
            for qty in list(set(quantities)):
                f.write('            <PDataArray type="Float64" Name="%s" NumberOfComponents="1"/>\n'%(qty))
            f.write('        </PPointData>\n')
            f.write('        <PCellData></PCellData>\n')
            f.write('        <PCoordinates>\n')
            f.write('            <PDataArray type="Float64" Name="Array X-Axis" NumberOfComponents="1"/>\n')
            f.write('            <PDataArray type="Float64" Name="Array Y-Axis" NumberOfComponents="1"/>\n')
            f.write('            <PDataArray type="Float64" Name="Array Z-Axis" NumberOfComponents="1"/>\n')
            f.write('        </PCoordinates>\n')
            for piece in list(pieces.keys()):
                lXB = pieces[piece]['lims']
                ext = "%d %d %d %d %d %d"%(lXB[0], lXB[1], lXB[2], lXB[3], lXB[4], lXB[5])
                f.write('         <Piece Extent="%s" Source="%s_%08d.vtr"/>\n'%(ext, piece.split(os.sep)[-1], time*1e3))
            f.write('    </PRectilinearGrid>\n')
            f.write('</VTKFile>\n')

def exportBndfDataToVtk(chid, resultDir, outtimes=None):
    # Boundary data
    quantities = fds.readBoundaryQuantities(resultDir, chid)
    uniqueQuantities = list(set(quantities))
    smvFile = fds.getFileList(resultDir, chid, 'smv')[0]
    bndfs = fds.getFileList(resultDir, chid, 'bf')
    fdsFile = fds.fdsFileOperations()
    fdsFile.importFile(fds.getFileList(resultDir, chid, 'fds')[0])
    meshes = list(fdsFile.meshes.keys())
    smvData = fds.parseSMVFile(smvFile)
    (smvGrids, smvObsts) = (smvData['grids'], smvData['obsts'])
    (smvBndfs, smvSurfs) = (smvData['bndfs'], smvData['surfs'])
    (smvFiles, bndes) = (smvData['files'], smvData['bndes'])
    bndf_dic = fds.linkBndfFileToMesh(meshes, bndfs, quantities)
    
    meshes = [x.split('_')[-2] for x in bndfs]
    uniqueMeshes = list(set(meshes))
    
    
    
    for i in range(0, len(uniqueMeshes)):
        series_data = dict()
        for qty in uniqueQuantities:
            series_data[qty] = dict()
        mesh = uniqueMeshes[i]
        for j in range(0, len(meshes)):
            if mesh != meshes[j]:
                continue
            file = bndfs[j]
            qty, shortName, units, npatch = fds.readBoundaryHeader(file)
            times, patches = fds.importBoundaryFile(file, gridNum=int(mesh)-1, grid=smvGrids)
            NT = len(times)
            pieces = []
            for patch in patches:
                lims = patch.lims
                data = patch.data
                ior = patch.orientation
                (N1, N2) = (data.shape[0], data.shape[1])
                if lims[0] == lims[1]: (NX, NY, NZ) = (1, N1, N2)
                if lims[2] == lims[3]: (NX, NY, NZ) = (N1, 1, N2)
                if lims[4] == lims[5]: (NX, NY, NZ) = (N1, N2, 1)
                
                numberOfPoints = (N1)*(N2)
                numberOfPolys = (N1-1) * (N2-1)
                
                Points = np.zeros((numberOfPoints*3,))
                PointData = np.zeros((numberOfPoints,NT))
                ConnectData = np.zeros((numberOfPoints*4,))
                
                if lims[0] == lims[1]:
                    y = np.linspace(lims[2], lims[3], NY)
                    z = np.linspace(lims[4], lims[5], NZ)
                    
                    counter = 0
                    for i in range(0, y.shape[0]):
                        for j in range(0, z.shape[0]):
                            Points[counter*3:(counter+1)*3] = [lims[0], y[i], z[j]]
                            counter += 1
                    counter = 0
                    for i in range(0, y.shape[0]):
                        for j in range(0, z.shape[0]):
                            PointData[counter, :] = data[i, j, :]
                            counter += 1
                    counter = 0
                    for i in range(0, z.shape[0]-1):
                        for j in range(0, y.shape[0]-1):
                            ConnectData[counter*4:(counter+1)*4] = [(z.shape[0])*j+i, (z.shape[0])*j+i+1, (z.shape[0])*(j+1)+i+1, (z.shape[0])*(j+1)+i]
                            counter += 1
                elif lims[2] == lims[3]:
                    x = np.linspace(lims[0], lims[1], NX)
                    z = np.linspace(lims[4], lims[5], NZ)
                    counter = 0
                    for i in range(0, x.shape[0]):
                        for j in range(0, z.shape[0]):
                            Points[counter*3:(counter+1)*3] = [x[i], lims[2], z[j]]
                            counter += 1
                    counter = 0
                    for i in range(0, x.shape[0]):
                        for j in range(0, z.shape[0]):
                            PointData[counter, :] = data[i, j, :]
                            counter += 1
                    counter = 0
                    for i in range(0, z.shape[0]-1):
                        for j in range(0, x.shape[0]-1):
                            ConnectData[counter*4:(counter+1)*4] = [(z.shape[0])*j+i, (z.shape[0])*j+i+1, (z.shape[0])*(j+1)+i+1, (z.shape[0])*(j+1)+i]
                            counter += 1
                    
                elif lims[4] == lims[5]:
                    x = np.linspace(lims[0], lims[1], NX)
                    y = np.linspace(lims[2], lims[3], NY)
                    counter = 0
                    for i in range(0, x.shape[0]):
                        for j in range(0, y.shape[0]):
                            Points[counter*3:(counter+1)*3] = [x[i], y[j], lims[4]]
                            counter += 1
                    counter = 0
                    for i in range(0, x.shape[0]):
                        for j in range(0, y.shape[0]):
                            PointData[counter, :] = data[i, j, :]
                            counter += 1
                    counter = 0
                    for i in range(0, y.shape[0]-1):
                        for j in range(0, x.shape[0]-1):
                            ConnectData[counter*4:(counter+1)*4] = [(y.shape[0])*j+i, (y.shape[0])*j+i+1, (y.shape[0])*(j+1)+i+1, (y.shape[0])*(j+1)+i]
                            counter += 1
                OffsetData = np.linspace(1, numberOfPolys, int(numberOfPolys))*4
                
                piece = dict()
                piece['numberOfPolys'] = numberOfPolys
                piece['numberOfPoints'] = numberOfPoints
                piece['Points'] = Points
                piece['PointData'] = PointData
                piece['CellData'] = False
                piece['Offset'] = OffsetData
                piece['Connectivity'] = ConnectData
                pieces.append(piece)
            series_data[qty][mesh] = dict()
            series_data[qty][mesh]['pieces'] = pieces
            series_data[qty][mesh]['times'] = times
        namespace = resultDir + os.sep + chid + '_bndf_%04d'%(int(mesh))
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        writeVtkPolyTimeSeries(namespace, series_data, times2)

    wrapper_namespace = resultDir + os.sep + chid + '_bndf_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtp', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PPolyData">\n')
            f.write('    <PPolyData GhostLevel="0">\n')
            f.write('        <PPointData>\n')
            for qty in uniqueQuantities:
                f.write('            <PDataArray type="Float32" Name="%s"/>\n'%(qty))
            f.write('        </PPointData>\n')
            f.write('        <PCellData></PCellData>\n')
            f.write('        <PPoints>\n')
            f.write('            <PDataArray type="Float32" NumberOfComponents="3"/>\n')
            f.write('        </PPoints>\n')
            for mesh in uniqueMeshes:
                source_name = chid + "_bndf_%04d_%08d.vtp"%(int(mesh), time*1e3)
                f.write('        <Piece Source="%s"/>\n'%(source_name))
            f.write('    </PPolyData>\n')
            f.write('</VTKFile>\n')
            
            
def exportPrt5DataToVtk(chid, resultDir, outtimes=None):
    # Particle data
    #particle_data, times_data = fds.importParticles(resultDir, chid)
    partFiles = fds.getFileList(resultDir, chid, 'prt5')
    meshes = [x.split('_')[-1].replace('.prt5','') for x in partFiles]
    uniqueMeshes = list(set(meshes))
    fname_written = defaultdict(bool)
    allQuantities = False
    for file in partFiles:
        particle_data, times = fds.importParticle(file)
        mesh = file.split('_')[-1].replace('.prt5','')
        if allQuantities is False:
            allQuantities = ['tag']
            pclasses = list(particle_data['classes'].keys())
            pclasses.remove('times')
            for pclass in pclasses:
                qtyNames = particle_data['classes'][pclass]['qtyNames']
                qtyUnits = particle_data['classes'][pclass]['qtyUnits']
                qids = [x + '(%s)'%(y) for x,y in zip(qtyNames, qtyUnits)]
                allQuantities.extend(qids)
            uniqueQuantities = list(set(allQuantities))
        
        times_data = defaultdict(bool)
        for time in times:
            times_data[time] = defaultdict(bool)
            times_data[time]['particles'] = []
        
        particles = list(particle_data['tags'].keys())
        for particle in particles:
            for t in particle_data['tags'][particle]['times']:
                times_data[t]['particles'].append(particle)
        
        particle_data['times'] = times
        
        namespace = resultDir + os.sep + chid + '_prt5_%04d'%(int(mesh))
        
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        
        for time in times: # enumerate([list(times_data.keys())[0]]):
            particles = times_data[time]['particles']
            pieces = defaultdict(bool)
            nParts = 0
            for particle in particles:
                pclass = particle_data['tags'][particle]['class']
                qtyNames = particle_data['classes'][pclass]['qtyNames']
                qtyUnits = particle_data['classes'][pclass]['qtyUnits']
                xyz = np.array(particle_data['tags'][particle]['xyz'])
                qids = [x + '(%s)'%(y) for x,y in zip(qtyNames, qtyUnits)]
                ptimes = np.array(particle_data['tags'][particle]['times'])
                
                if (ptimes.min() < time) and (ptimes.max() > time):
                    ptimeind = np.argmin(abs(ptimes-time))
                    nParts += 1
                    if pieces[pclass] is False:
                        pieces[pclass] = defaultdict(bool)
                        pieces[pclass]['numberOfPoints'] = 1
                        pieces[pclass]['Points'] = [xyz[ptimeind, :]]
                        pieces[pclass]['PointData'] = defaultdict(bool)
                        pieces[pclass]['PointData']['tag'] = [particle]
                        for qid in qids:
                            pieces[pclass]['PointData'][qid] = [particle_data['tags'][particle][qid][ptimeind]]
                    else:
                        pieces[pclass]['numberOfPoints'] += 1
                        pieces[pclass]['Points'].append(xyz[ptimeind, :])
                        pieces[pclass]['PointData']['tag'].append(particle)
                        for qid in qids:
                            pieces[pclass]['PointData'][qid].append(particle_data['tags'][particle][qid][ptimeind])
            classes = list(pieces.keys())
            for pclass in classes:
                pieces[pclass]['Points'] = np.array(pieces[pclass]['Points']).flatten()
                pieces[pclass]['xyz'] = np.array(pieces[pclass]['xyz']).flatten()
            
            for pclass in classes:
                if pieces[pclass]['numberOfPoints'] > 0:
                    pieces2 = dict()
                    pieces2[pclass] = pieces[pclass]
                    fname = namespace + '_' + pclass +'_%08d'%(time*1e3)
                    writeVtkPointFile(fname, pieces2)
                    fname_written[fname] = True
    for pclass in pclasses:
        wrapper_namespace = resultDir + os.sep + chid + '_prt5_out_' + pclass
        qtyNames = particle_data['classes'][pclass]['qtyNames']
        qtyUnits = particle_data['classes'][pclass]['qtyUnits']
        qids = [x + '(%s)'%(y) for x,y in zip(qtyNames, qtyUnits)]
        qids.extend(['tag'])
        for time in outtimes:
            fname = wrapper_namespace + '_%08d'%(time*1e3)
            with open(fname + '.pvtp', 'w') as f:
                f.write('<?xml version="1.0"?>\n')
                f.write('<VTKFile type="PPolyData">\n')
                f.write('    <PPolyData GhostLevel="0">\n')
                f.write('        <PPointData>\n')
                for qty in qids:
                    f.write('            <PDataArray type="Float32" Name="%s"/>\n'%(qty))
                f.write('        </PPointData>\n')
                f.write('        <PCellData></PCellData>\n')
                f.write('        <PPoints>\n')
                f.write('            <PDataArray type="Float32" NumberOfComponents="3"/>\n')
                f.write('        </PPoints>\n')
                for mesh in uniqueMeshes:
                    source_name = chid + "_prt5_%04d_%s_%08d.vtp"%(int(mesh), pclass, time*1e3)
                    if fname_written[resultDir + os.sep + source_name.replace('.vtp','')]:
                        f.write('        <Piece Source="%s"/>\n'%(source_name))
                f.write('    </PPolyData>\n')
                f.write('</VTKFile>\n')

def exportBndeDataToVtk(chid, resultDir, outtimes=None):
    smvFiles = fds.getFileList(resultDir, chid, 'smv')
    
    file_quantities = fds.getBndeQuantities(smvFiles[0])
    quantities = [file_quantities[key]['quantity'] for key in list(file_quantities.keys())]
    
    uniqueQuantities = list(set(quantities))
    
    gbfFiles = fds.getFileList(resultDir, chid, 'gbf')
    gcfFiles = fds.getFileList(resultDir, chid, 'gcf')
    beFiles = fds.getFileList(resultDir, chid, 'be')
    
    #meshes = [x.split('_')[-2] for x in beFiles]
    meshes = [x.split('_')[-1].replace('.gcf','') for x in gcfFiles]
    uniqueMeshes = list(set(meshes))
    outtimes = None
    
    for i in range(0, len(uniqueMeshes)):
        series_data = dict()
        for qty in uniqueQuantities:
            series_data[qty] = dict()
        mesh = uniqueMeshes[i]
        for j in range(0, len(meshes)):
            if mesh != meshes[j]:
                continue
            gcffile = gcfFiles[j]
            vertices, faces2, header1 = fds.readGcfFile(gcffile)
            numberOfPoints = vertices.shape[0]
            numberOfPolys = faces2.shape[0]
            Points = vertices.flatten() #np.zeros((numberOfPoints*3,))
            ConnectData = faces2.flatten()-1 #np.zeros((numberOfPolys*3,))
            OffsetData = np.linspace(1, numberOfPolys, int(numberOfPolys))*3
            
            beFiles = fds.getFileList(resultDir, chid+gcffile.split(chid)[-1].replace('.gcf','_'),'be')
            for befile in beFiles:
                times, vals, header2 = fds.readBeFile(befile)
                NT = times.shape[0]
                CellData = vals
                qty = file_quantities[chid + befile.split(chid)[-1]]['quantity']
                piece = dict()
                piece['numberOfPolys'] = numberOfPolys
                piece['numberOfPoints'] = numberOfPoints
                piece['Points'] = Points
                piece['PointData'] = False
                piece['CellData'] = CellData
                piece['Offset'] = OffsetData
                piece['Connectivity'] = ConnectData
                
                series_data[qty][mesh] = dict()
                series_data[qty][mesh]['pieces'] = [piece]
                series_data[qty][mesh]['times'] = times
        namespace = resultDir + os.sep + chid + '_bnde_%04d'%(int(mesh))
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        writeVtkPolyTimeSeries(namespace, series_data, times2)
    
    
    wrapper_namespace = resultDir + os.sep + chid + '_bnde_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtp', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PPolyData">\n')
            f.write('    <PPolyData GhostLevel="0">\n')
            f.write('        <PPointData></PPointData>\n')
            f.write('        <PCellData>\n')
            for qty in uniqueQuantities:
                f.write('            <PDataArray type="Float32" Name="%s"/>\n'%(qty))
            f.write('        </PCellData>\n')
            f.write('        <PPoints>\n')
            f.write('            <PDataArray type="Float32" NumberOfComponents="3"/>\n')
            f.write('        </PPoints>\n')
            for mesh in uniqueMeshes:
                source_name = chid + "_bnde_%04d_%08d.vtp"%(int(mesh), time*1e3)
                f.write('        <Piece Source="%s"/>\n'%(source_name))
            f.write('    </PPolyData>\n')
            f.write('</VTKFile>\n')

if __name__ == '__main__':
    
    chid = "case002"
    resultDir = "examples\\case002\\"
    
    '''
    (tmin, tmax) = (np.min(times), np.max(times))
    dt = 0.2
    nt = int(((tmax-tmin)/dt)+1)
    times_out = np.linspace(tmin, tmax, nt)
    '''
    
    exportSl3dDataToVtk(chid, resultDir, outtimes=None)
    exportBndfDataToVtk(chid, resultDir, outtimes=None)
    exportPrt5DataToVtk(chid, resultDir, outtimes=None)
    exportBndeDataToVtk(chid, resultDir, outtimes=None)
    
    #obstToStl(resultDir, chid)
    