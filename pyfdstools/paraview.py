from .extractS3D import extractS3dValues
from .extractBoundaryData import readBoundaryQuantities, linkBndfFileToMesh
from .extractBoundaryData import readBoundaryHeader, importBoundaryFile
from .extractGeomData import getBndeQuantities, readGcfFile, readBeFile
from .extractParticleData import importParticle
from .extractPlot3Ddata import readSLCFquantities, readSingleSlcfFile
from .fdsFileOperations import fdsFileOperations
from .utilities import getDatatypeByEndianness, getEndianness
from .utilities import getFileListFromZip, getFileList, zopen, zreadlines
from .utilities import getFileListFromResultDir
from .utilities import getGridsFromXyzFiles, getAbsoluteGrid, rearrangeGrid
from .utilities import readXYZfile
from .smokeviewParser import parseSMVFile
import os
import numpy as np
import stl
from collections import defaultdict
import evtk

def obstToStl(resultDir, chid, outDir=None):
    if outDir is None: outDir = resultDir
    smvFile = getFileListFromResultDir(resultDir, chid, 'smv')[0]
    smvData = parseSMVFile(smvFile)
    (grid, obst) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
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
    if len(obst) == 0:
        print("Warning no obst found in %s, not generating an STL"%(smvFile))
        return
    
    for o in obst:
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
    combined.save(os.path.join(outDir, chid+'.stl')) #,mode=stl.Mode.ASCII)

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

def writeVtkScalarsPoints(fname, qtys, datas, dtype):
    with open(fname, 'a') as f:
        f.write('<PointData>\n') # Scalars=%s>\n'%(",".join(['"%s"'%(q) for q in qtys])))
        for qty in qtys:
            f.write('<DataArray Name="%s" NumberOfComponents="1" type="%s" format="ascii">\n'%(qty, dtype))
            for i in range(0, datas[qty].shape[0]): 
                if 'Float' in dtype:
                    f.write('%0.1f '%(datas[qty][i]))
                elif 'UInt8' in dtype:
                    f.write('%0d '%(datas[qty][i]))
                else:
                    print("Warning dtype %s unknown"%(dtype))
            f.write('\n</DataArray>\n')
        f.write('</PointData>\n')

def writeVtkScalarsCells(fname, qtys, datas, dtype):
    with open(fname, 'a') as f:
        f.write('<CellData>\n') # Scalars=%s>\n'%(",".join(['"%s"'%(q) for q in qtys])))
        for qty in qtys:
            f.write('<DataArray Name="%s" NumberOfComponents="1" type="%s" format="ascii">\n'%(qty, dtype))
            for i in range(0, datas[qty].shape[0]):
                if 'Float' in dtype:
                    f.write('%0.1f '%(datas[qty][i]))
                elif 'UInt8' in dtype:
                    f.write('%d '%(datas[qty][i]))
                else:
                    print("Warning dtype %s unknown"%(dtype))
            f.write('\n</DataArray>\n')
        f.write('</CellData>\n')


def writeVtkTail(fname, ftype):
    with open(fname, 'a') as f:
        f.write('</Piece>\n')
        f.write('</%s>\n'%(ftype))
        f.write('</VTKFile>\n')

def writeVtkFile(fname, x, y, z, gXB, lXB, datas, cell=False, binary=True, dtype='Float64'):
    qtys = list(datas.keys())
    if binary:
        start = (gXB[0], gXB[2], gXB[4])
        if cell:
            evtk.hl.gridToVTK(fname, x, y, z, cellData = datas, start=start)
        else:
            evtk.hl.gridToVTK(fname, x, y, z, pointData = datas, start=start)
    else:
        writeVtkHeader(fname, "RectilinearGrid", ".vtr")
        writeVtkGrid(fname, x, y, z, gXB, lXB)
        if cell:
            writeVtkScalarsCells(fname+'.vtr', qtys, datas, dtype=dtype)
        else:
            writeVtkScalarsPoints(fname+'.vtr', qtys, datas, dtype=dtype)
        writeVtkTail(fname+'.vtr', 'RectilinearGrid')
    
def writeVtkTimeSeries(namespace, grid, series_data, times, gXB, lXB, cell=False, binary=True, dtype='Float32'):
    nx, ny, nz, _ = np.shape(grid)
    x = np.array(grid[:, 0, 0, 0], dtype='float64')
    y = np.array(grid[0, :, 0, 1], dtype='float64')
    z = np.array(grid[0, 0, :, 2], dtype='float64')
    
    for time in times:
        fname = namespace + '_%08d'%(time*1e3)
        datas = dict()
        for qty in list(series_data.keys()):
            tind = np.argmin(abs(series_data[qty]['times'] - time))
            tmp = series_data[qty]['data'][:, :, :, tind].T.flatten()
            if dtype == 'Float32':
                tmp = np.array(tmp, dtype='float32')
            elif dtype == 'Float64':
                tmp = np.array(tmp, dtype='float64')
            elif dtype == 'UInt8':
                tmp = np.array(tmp, dtype='uint8')
            else:
                print("Warning dtype %s unknown"%(dtype))
            datas[qty] = series_data[qty]['data'][:, :, :, tind].T.flatten()
        writeVtkFile(fname, x, y, z, gXB, lXB, datas, cell=cell, binary=binary, dtype=dtype)

def writeVtkImageFile(fname, x, y, z, gXB, lXB, datas, cell=False, binary=True, dtype='Float64'):
    qtys = list(datas.keys())
    dx = np.median(x[1:]-x[:-1])
    dy = np.median(y[1:]-y[:-1])
    dz = np.median(z[1:]-z[:-1])
    spacing = (dx, dy, dz)
    origin = (x.min(), y.min(), z.min())
    if binary:
        start = (gXB[0], gXB[2], gXB[4])
        if cell:
            evtk.hl.imageToVTK(fname, origin, spacing, cellData = datas, start=start)
        else:
            evtk.hl.imageToVTK(fname, origin, spacing, pointData = datas, start=start)
    else:
        writeVtkHeader(fname, "ImageData", ".vti")
        with open(fname+'.vti', 'a') as f:
            f.write('<ImageData WholeExtent="%d %d %d %d %d %d"\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
            f.write(' Origin="%0.4f %0.4f %0.4f" Spacing="%0.4f %0.4f %0.4f">\n'%(origin[0], origin[1], origin[2], dx, dy, dz))
            f.write('<Piece Extent="%d %d %d %d %d %d">\n'%(lXB[0], lXB[1], lXB[2], lXB[3], lXB[4], lXB[5]))
            
        #writeVtkGrid(fname, x, y, z, gXB, lXB)
        if cell:
            writeVtkScalarsCells(fname+'.vti', qtys, datas, dtype=dtype)
        else:
            writeVtkScalarsPoints(fname+'.vti', qtys, datas, dtype=dtype)
        writeVtkTail(fname+'.vti', 'ImageData')
        
def writeVtkImageTimeSeries(namespace, grid, series_data, times, gXB, lXB, cell=False, binary=True, dtype='Float32'):
    nx, ny, nz, _ = np.shape(grid)
    x = np.array(grid[:, 0, 0, 0], dtype='float64')
    y = np.array(grid[0, :, 0, 1], dtype='float64')
    z = np.array(grid[0, 0, :, 2], dtype='float64')
    
    for time in times:
        fname = namespace + '_%08d'%(time*1e3)
        datas = dict()
        for qty in list(series_data.keys()):
            tind = np.argmin(abs(series_data[qty]['times'] - time))
            tmp = series_data[qty]['data'][:, :, :, tind].T.flatten()
            if dtype == 'Float32':
                tmp = np.array(tmp, dtype='float32')
            elif dtype == 'Float64':
                tmp = np.array(tmp, dtype='float64')
            elif dtype == 'UInt8':
                tmp = np.array(tmp, dtype='uint8')
            else:
                print("Warning dtype %s unknown"%(dtype))
            datas[qty] = series_data[qty]['data'][:, :, :, tind].T.flatten()
        writeVtkImageFile(fname, x, y, z, gXB, lXB, datas, cell=cell, binary=binary, dtype=dtype)

def writeVtkPolyFile(fname, pieces, binary=True):
    if binary:
        keys = list(pieces[0].keys())
        keys.remove('numberOfPolys')
        keys.remove('numberOfPoints')
        keys.remove('Points')
        keys.remove('Offset')
        keys.remove('Connectivity')
        pquantities = [key for key in keys if 'P-' in key]
        cquantities = [key for key in keys if 'C-' in key]
        w = evtk.hl.VtkFile(fname, evtk.vtk.VtkPolyData)
        w.openGrid()
        data_to_append = []
        append_type = []
        for i in range(0, len(pieces)):
            piece = pieces[i]
            npolys = piece['numberOfPolys']
            npoints = piece['numberOfPoints']
            points = piece['Points']
            x = np.ascontiguousarray(points[::3])
            y = np.ascontiguousarray(points[1::3])
            z = np.ascontiguousarray(points[2::3])
            offsets = piece['Offset']
            connectivity = piece['Connectivity']
            
            w.openPiece(npolys=npolys, npoints=npoints)
            
            w.openElement("Points")
            w.addData("points", (x, y, z))
            w.closeElement("Points")
            w.openElement("Polys")
            w.addData("connectivity", connectivity)
            w.addData("offsets", offsets)
            w.closeElement("Polys")
            
            pointData = dict()
            for qty in pquantities:
                pointData[qty.replace('P-','')] = np.ascontiguousarray(piece[qty])
            cellData = dict()
            for qty in cquantities:
                cellData[qty.replace('C-','')] = np.ascontiguousarray(piece[qty])
            evtk.hl._addDataToFile(w, pointData=pointData, cellData=cellData)
            
            w.closePiece()
            data_to_append.append((x, y, z))
            data_to_append.append(connectivity)
            data_to_append.append(offsets)
            append_type.append(1)
            append_type.append(1)
            append_type.append(1)
            data_to_append.append((pointData, cellData))
            append_type.append(2)
        w.closeGrid()
        for i in range(0, len(data_to_append)):
            if append_type[i] == 1:
                w.appendData(data_to_append[i])
            elif append_type[i] == 2:
                evtk.hl._appendDataToFile(w, pointData=data_to_append[i][0], cellData=data_to_append[i][1])
            else:
                pass
        w.save()
    else:
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

def writeVtkPointFile(fname, pieces, binary=True):
    classes = list(pieces.keys())
    
    if binary:
        w = evtk.hl.VtkFile(fname, evtk.vtk.VtkPolyData)
        w.openGrid()
        data_to_append = []
        append_type = []
        for i in range(0, len(classes)):
            piece = pieces[classes[i]]
            
            npoints = piece['numberOfPoints']
            points = piece['Points']
            quantities = list(piece['PointData'].keys())
            x = np.ascontiguousarray(points[::3], dtype=np.float64)
            y = np.ascontiguousarray(points[1::3], dtype=np.float64)
            z = np.ascontiguousarray(points[2::3], dtype=np.float64)
            
            w.openPiece(npoints=npoints)
            
            w.openElement("Points")
            w.addData("points", (x, y, z))
            w.closeElement("Points")
            
            pointData = dict()
            for qty in quantities:
                pointData[qty] = np.ascontiguousarray(piece['PointData'][qty], dtype=np.float64)
            evtk.hl._addDataToFile(w, pointData=pointData,cellData={})
            
            w.closePiece()
            data_to_append.append((x, y, z))
            append_type.append(1)
            data_to_append.append(pointData)
            append_type.append(2)
        w.closeGrid()
        for i in range(0, len(data_to_append)):
            if append_type[i] == 1:
                w.appendData(data_to_append[i])
            elif append_type[i] == 2:
                evtk.hl._appendDataToFile(w, pointData=data_to_append[i], cellData={})
            else:
                pass
        w.save()
    else:
        writeVtkHeader(fname, "PolyData", ".vtp")
        
        with open(fname + '.vtp', 'a') as f:
            f.write('<PolyData>\n')
            
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

def writeVtkPolyTimeSeries(namespace, series_data, times, binary=True):
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
        writeVtkPolyFile(fname, pieces, binary=binary)



def exportSl3dDataToVtk(chid, resultDir, outtimes=None, outDir=None, binary=True, dtype=None, ftype='ImageData'):
    if outDir is None: outDir = resultDir
    quantities, slcfFiles, dimensions, meshes, centers = readSLCFquantities(chid, resultDir)
    twoDslice = [True if (dim[0] == dim[1]) or (dim[2] == dim[3]) or (dim[4] == dim[5]) else False for dim in dimensions]
    slcfs = getFileList(resultDir, chid, 'sf')
    if len(slcfs) == 0:
        print("No sf files found in %s with chid %s"%(resultDir, chid))
        return None
    if np.all(twoDslice):
        print("No 3d sf files found in %s with chid %s"%(resultDir, chid))
        return None
    if dtype is None:
        if binary:
            dtype='Float64'
        else:
            dtype='Float32'
            
    if ftype == 'ImageData':
        ftype = 'ImageData'
        pftype = 'PImageData'
        fext = '.vti'
        pext = '.pvti'
    elif ftype == 'RectilinearGrid':
        ftype = 'RectilinearGrid'
        pftype = 'PRectilinearGrid'
        fext = '.vtr'
        pext = '.pvtr'
    
    uniqueMeshes = sorted(list(set(meshes)))
    
    xyzFiles = getFileListFromResultDir(resultDir, chid, 'xyz')
    grids = getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = getAbsoluteGrid(grids)
    
    gXB = [0, grids_abs.shape[0]-1, 0, grids_abs.shape[1]-1, 0, grids_abs.shape[2]-1]
    pieces=defaultdict(bool)
    threeD_quantities = []
    
    # Need to fix 3-d slices that do not fill the mesh
    
    #for g in sorted(list(grids.keys())):
    #    print(grids[g]['zGrid'].shape)
    #print(gXB)
    #assert False, "Stopped"
    (xmin,xmax,ymin,ymax,zmin,zmax) = (1e12,-1e12,1e12,-1e12,1e12,-1e12)
    for i in range(0, len(uniqueMeshes)):
        mesh = uniqueMeshes[i]
        xyzFile = False
        namespace = outDir + os.sep + chid + '_sl3d_grid%04d'%(int(mesh))
        for j in range(0, len(slcfFiles)):
            if mesh != meshes[j]:
                continue
            if twoDslice[j]:
                continue
            slcfFile = slcfFiles[j]
            qty = quantities[j]
            threeD_quantities.append(qty)
            lims, datas, times = readSingleSlcfFile(slcfFile)
            if dtype == 'Float64':
                datas = np.array(datas, dtype='float64')
            elif dtype == 'Float32':
                datas = np.array(datas, dtype='float32')
            elif dtype == 'UInt8':
                datas = np.array(datas, dtype='uint8')
            if xyzFile is False:
                xyzFile = '_'.join(slcfFile.split('_')[:-1]) + '.xyz'
                grid, gridHeader = readXYZfile(xyzFile)
                xGrid, yGrid, zGrid = rearrangeGrid(grid)
                
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
            (xmin, xmax) = (min([xmin, xGrid.min()]), max([xmax, xGrid.max()]))
            (ymin, ymax) = (min([ymin, yGrid.min()]), max([ymax, yGrid.max()]))
            (zmin, zmax) = (min([zmin, zGrid.min()]), max([zmax, zGrid.max()]))
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        if ftype == 'RectilinearGrid':
            writeVtkTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB, binary=binary, dtype=dtype)
        elif ftype == 'ImageData':
            writeVtkImageTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB, binary=binary, dtype=dtype)
        pieces[namespace] = defaultdict(bool)
        pieces[namespace]['source'] = namespace
        pieces[namespace]['times'] = times2
        pieces[namespace]['lims'] = lXB
    
    x = np.array(grid_combo[:, 0, 0, 0], dtype='float64')
    y = np.array(grid_combo[0, :, 0, 1], dtype='float64')
    z = np.array(grid_combo[0, 0, :, 2], dtype='float64')
    dx = np.median(x[1:]-x[:-1])
    dy = np.median(y[1:]-y[:-1])
    dz = np.median(z[1:]-z[:-1])
    spacing = (dx, dy, dz)
    origin = (xmin, ymin, zmin)
    wrapper_namespace = outDir + os.sep + chid + '_sl3d_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + pext, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="%s">\n'%(pftype))
            if pftype == 'PRectilinearGrid':
                f.write('    <PRectilinearGrid WholeExtent="%d %d %d %d %d %d" GhostLevel="0">\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
            elif pftype == 'PImageData':
                f.write('    <PImageData WholeExtent="%d %d %d %d %d %d" GhostLevel="0" \n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
                f.write('                Origin="%0.4f %0.4f %0.4f" Spacing="%0.4f %0.4f %0.4f">\n'%(origin[0], origin[1], origin[2], dx, dy, dz))
            f.write('        <PPointData>\n')
            for qty in sorted(list(set(threeD_quantities))):
                f.write('            <PDataArray type="%s" Name="%s" NumberOfComponents="1"/>\n'%(dtype, qty))
            f.write('        </PPointData>\n')
            f.write('        <PCellData></PCellData>\n')
            f.write('        <PCoordinates>\n')
            f.write('            <PDataArray type="Float64" Name="Array X-Axis" NumberOfComponents="1"/>\n')
            f.write('            <PDataArray type="Float64" Name="Array Y-Axis" NumberOfComponents="1"/>\n')
            f.write('            <PDataArray type="Float64" Name="Array Z-Axis" NumberOfComponents="1"/>\n')
            f.write('        </PCoordinates>\n')
            for piece in sorted(list(pieces.keys())):
                lXB = pieces[piece]['lims']
                ext = "%d %d %d %d %d %d"%(lXB[0], lXB[1], lXB[2], lXB[3], lXB[4], lXB[5])
                f.write('         <Piece Extent="%s" Source="%s_%08d%s"/>\n'%(ext, piece.split(os.sep)[-1], time*1e3, fext))
            f.write('    </%s>\n'%(pftype))
            f.write('</VTKFile>\n')

def exportSl2dDataToVtk(chid, resultDir, outtimes=None, outDir=None, binary=True, dtype=None, ftype='ImageData'):
    if outDir is None: outDir = resultDir
    quantities, slcfFiles, dimensions, meshes, centers = readSLCFquantities(chid, resultDir)
    twoDslice = [True if (dim[0] == dim[1]) or (dim[2] == dim[3]) or (dim[4] == dim[5]) else False for dim in dimensions]
    slcfs = getFileList(resultDir, chid, 'sf')
    if len(slcfs) == 0:
        print("No sf files found in %s with chid %s"%(resultDir, chid))
        return None
    if ~np.any(twoDslice):
        print("No 2d sf files found in %s with chid %s"%(resultDir, chid))
        return None
    if dtype is None:
        if binary:
            dtype='Float64'
        else:
            dtype='Float32'
            
    if ftype == 'ImageData':
        ftype = 'ImageData'
        pftype = 'PImageData'
        fext = '.vti'
        pext = '.pvti'
    elif ftype == 'RectilinearGrid':
        ftype = 'RectilinearGrid'
        pftype = 'PRectilinearGrid'
        fext = '.vtr'
        pext = '.pvtr'
    
    xyzFiles = getFileListFromResultDir(resultDir, chid, 'xyz')
    grids = getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = getAbsoluteGrid(grids)
    
    outs = parseSMVFile(os.path.join(resultDir, chid + '.smv'))
    terrainHeight = 0
    outFiles = defaultdict(bool)
    for i in range(len(slcfFiles)):
        if twoDslice[i]:
            dim = dimensions[i]
            grid = grids[str(meshes[i])]
            if (dim[0] == dim[1]): axis = 1; value = grid['xGrid'][dim[0],0,0]
            if (dim[2] == dim[3]): axis = 2; value = grid['yGrid'][0,dim[2],0]
            if (dim[4] == dim[5]): axis = 3; value = grid['zGrid'][0,0,dim[4]]
            slice_type = outs['files']['SLICES'][os.path.basename(slcfFiles[i])]['LINETEXT'].split()[0]
            if slice_type == 'SLCT':
                if terrainHeight == 0:
                    terrainHeight = value
                else:
                    value = terrainHeight
            if outFiles[axis] is False:
                outFiles[axis] = defaultdict(bool)
            if outFiles[axis][value] is False:
                outFiles[axis][value] = defaultdict(bool)
            if outFiles[axis][value][meshes[i]] is False:
                outFiles[axis][value][meshes[i]] = defaultdict(bool)
            outFiles[axis][value][meshes[i]][quantities[i]] = slcfFiles[i]
            #print(os.path.basename(slcfFiles[i]), axis, value, slice_type)
            #assert False, "Stopped"
    
    gXB2 = [0, grids_abs.shape[0]-1, 0, grids_abs.shape[1]-1, 0, grids_abs.shape[2]-1]
    (xmin,xmax,ymin,ymax,zmin,zmax) = (1e12,-1e12,1e12,-1e12,1e12,-1e12)
    for axis in sorted(list(outFiles.keys())):
        for value in sorted(list(outFiles[axis].keys())):
            pieces=defaultdict(bool)
            quantities = []
            gXB = [x for x in gXB2]
            if axis == 1: gXB[0] = 1; gXB[1] = 1
            if axis == 2: gXB[2] = 1; gXB[3] = 1
            if axis == 3: gXB[4] = 1; gXB[5] = 1
            for mesh in sorted(list(outFiles[axis][value].keys())):
                xyzFile = False
                namespace = (outDir + os.sep + chid +'slcf_%d_%.6f_%06d'%(axis,value, int(mesh))).replace('.','_')
                for qty in sorted(list(outFiles[axis][value][mesh].keys())):
                    slcfFile = outFiles[axis][value][mesh][qty]
                    quantities.append(qty)
                    #print(axis, value, mesh, os.path.basename(file))
                    
                    lims, datas, times = readSingleSlcfFile(slcfFile)
                    if dtype == 'Float64':
                        datas = np.array(datas, dtype='float64')
                    elif dtype == 'Float32':
                        datas = np.array(datas, dtype='float32')
                    elif dtype == 'UInt8':
                        datas = np.array(datas, dtype='uint8')
                    if xyzFile is False:
                        xyzFile = '_'.join(slcfFile.split('_')[:-1]) + '.xyz'
                        grid, gridHeader = readXYZfile(xyzFile)
                        xGrid, yGrid, zGrid = rearrangeGrid(grid)
                        
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
                        (xmin, xmax) = (min([xmin, grids_abs[lXB[0],0,0,0]]), max([xmax, grids_abs[lXB[1],0,0,0]]))
                        (ymin, ymax) = (min([ymin, grids_abs[0,lXB[2],0,1]]), max([ymax, grids_abs[0,lXB[3],0,1]]))
                        (zmin, zmax) = (min([zmin, grids_abs[0,0, lXB[4],2]]), max([zmax, grids_abs[0,0, lXB[5],2]]))
                        if axis == 1: lXB[0] = 1; lXB[1] = 1
                        if axis == 2: lXB[2] = 1; lXB[3] = 1
                        if axis == 3: lXB[4] = 1; lXB[5] = 1
                
                if outtimes is None:
                    times2 = times
                    outtimes = times2
                else:
                    times2 = outtimes
                if ftype == 'RectilinearGrid':
                    writeVtkTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB, binary=binary, dtype=dtype)
                elif ftype == 'ImageData':
                    writeVtkImageTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB, binary=binary, dtype=dtype)
                pieces[namespace] = defaultdict(bool)
                pieces[namespace]['source'] = namespace
                pieces[namespace]['times'] = times2
                pieces[namespace]['lims'] = lXB
                
            x = np.array(grid_combo[:, 0, 0, 0], dtype='float64')
            y = np.array(grid_combo[0, :, 0, 1], dtype='float64')
            z = np.array(grid_combo[0, 0, :, 2], dtype='float64')
            dx = np.median(x[1:]-x[:-1])
            dy = np.median(y[1:]-y[:-1])
            dz = np.median(z[1:]-z[:-1])
            spacing = (dx, dy, dz)
            origin = (xmin, ymin, zmin)
            wrapper_namespace = (outDir + os.sep + chid + '_slcf_out_%d_%.6f'%(axis,value)).replace('.','_')
            for time in outtimes:
                fname = wrapper_namespace + '_%08d'%(time*1e3)
                with open(fname + pext, 'w') as f:
                    f.write('<?xml version="1.0"?>\n')
                    f.write('<VTKFile type="%s">\n'%(pftype))
                    if pftype == 'PRectilinearGrid':
                        f.write('    <PRectilinearGrid WholeExtent="%d %d %d %d %d %d" GhostLevel="0">\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
                    elif pftype == 'PImageData':
                        f.write('    <PImageData WholeExtent="%d %d %d %d %d %d" GhostLevel="0" \n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
                        f.write('                Origin="%0.4f %0.4f %0.4f" Spacing="%0.4f %0.4f %0.4f">\n'%(origin[0], origin[1], origin[2], dx, dy, dz))
                    f.write('        <PPointData>\n')
                    for qty in sorted(list(set(quantities))):
                        f.write('            <PDataArray type="%s" Name="%s" NumberOfComponents="1"/>\n'%(dtype, qty))
                    f.write('        </PPointData>\n')
                    f.write('        <PCellData></PCellData>\n')
                    f.write('        <PCoordinates>\n')
                    f.write('            <PDataArray type="Float64" Name="Array X-Axis" NumberOfComponents="1"/>\n')
                    f.write('            <PDataArray type="Float64" Name="Array Y-Axis" NumberOfComponents="1"/>\n')
                    f.write('            <PDataArray type="Float64" Name="Array Z-Axis" NumberOfComponents="1"/>\n')
                    f.write('        </PCoordinates>\n')
                    for piece in sorted(list(pieces.keys())):
                        lXB = pieces[piece]['lims']
                        ext = "%d %d %d %d %d %d"%(lXB[0], lXB[1], lXB[2], lXB[3], lXB[4], lXB[5])
                        f.write('         <Piece Extent="%s" Source="%s_%08d%s"/>\n'%(ext, piece.split(os.sep)[-1], time*1e3, fext))
                    f.write('    </%s>\n'%(pftype))
                    f.write('</VTKFile>\n')

def exportBndfDataToVtk(chid, resultDir, outtimes=None, outDir=None, binary=True):
    if outDir is None: outDir = resultDir
    # Boundary data
    quantities = readBoundaryQuantities(resultDir, chid)
    uniqueQuantities = list(set(quantities))
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    bndfs = getFileList(resultDir, chid, 'bf')
    if len(bndfs) == 0:
        print("No bf files found in %s with chid %s"%(resultDir, chid))
        return None
    fdsFile = fdsFileOperations()
    fdsFile.importFile(getFileList(resultDir, chid, 'fds')[0])
    meshes = list(fdsFile.meshes.keys())
    smvData = parseSMVFile(smvFile)
    (smvGrids, smvObsts) = (smvData['grids'], smvData['obsts'])
    (smvBndfs, smvSurfs) = (smvData['bndfs'], smvData['surfs'])
    (smvFiles, bndes) = (smvData['files'], smvData['bndes'])
    bndfs = [os.path.join(resultDir, x[1]) for x in smvBndfs]
    bndf_dic = linkBndfFileToMesh(meshes, bndfs, quantities)
    
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
            qty, shortName, units, npatch = readBoundaryHeader(file)
            times, patches = importBoundaryFile(file, gridNum=int(mesh)-1, grid=smvGrids)
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
        namespace = outDir + os.sep + chid + '_bndf_%04d'%(int(mesh))
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        writeVtkPolyTimeSeries(namespace, series_data, times2, binary=binary)

    wrapper_namespace = outDir + os.sep + chid + '_bndf_out'
    if binary:
        dataprecision = "Float64"
    else:
        dataprecision = "Float32"
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtp', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PPolyData">\n')
            f.write('    <PPolyData GhostLevel="0">\n')
            f.write('        <PPointData>\n')
            for qty in uniqueQuantities:
                f.write('            <PDataArray type="%s" Name="%s"/>\n'%(dataprecision, qty))
            f.write('        </PPointData>\n')
            f.write('        <PCellData></PCellData>\n')
            f.write('        <PPoints>\n')
            f.write('            <PDataArray type="%s" NumberOfComponents="3"/>\n'%(dataprecision))
            f.write('        </PPoints>\n')
            for mesh in uniqueMeshes:
                source_name = chid + "_bndf_%04d_%08d.vtp"%(int(mesh), time*1e3)
                f.write('        <Piece Source="%s"/>\n'%(source_name))
            f.write('    </PPolyData>\n')
            f.write('</VTKFile>\n')
            
            
def exportPrt5DataToVtk(chid, resultDir, outtimes=None, outDir=None, binary=True):
    if outDir is None: outDir = resultDir
    # Particle data
    partFiles = getFileList(resultDir, chid, 'prt5')
    if len(partFiles) == 0:
        print("No prt5 files found in %s with chid %s"%(resultDir, chid))
        return None
    meshes = [x.split('_')[-1].replace('.prt5','') for x in partFiles]
    uniqueMeshes = sorted(list(set(meshes)))
    fname_written = defaultdict(bool)
    allQuantities = False
    for file in partFiles:
        particle_data, times = importParticle(file)
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
            uniqueQuantities = sorted(list(set(allQuantities)))
        
        times_data = defaultdict(bool)
        for time in times:
            times_data[time] = defaultdict(bool)
            times_data[time]['particles'] = []
        
        particles = sorted(list(particle_data['tags'].keys()))
        for particle in particles:
            for t in particle_data['tags'][particle]['times']:
                times_data[t]['particles'].append(particle)
        
        particle_data['times'] = times
        
        namespace = outDir + os.sep + chid + '_prt5_%04d'%(int(mesh))
        
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        
        for time in outtimes: # enumerate([list(times_data.keys())[0]]):
            t_ind = np.argmin(abs(times-time))
            particles = times_data[times[t_ind]]['particles']
            pieces = defaultdict(bool)
            nParts = 0
            for particle in particles:
                pclass = particle_data['tags'][particle]['class']
                qtyNames = particle_data['classes'][pclass]['qtyNames']
                qtyUnits = particle_data['classes'][pclass]['qtyUnits']
                xyz = np.array(particle_data['tags'][particle]['xyz'])
                qids = [x + '(%s)'%(y) for x,y in zip(qtyNames, qtyUnits)]
                ptimes = np.array(particle_data['tags'][particle]['times'])
                
                if (ptimes.min() <= time) and (ptimes.max() >= time):
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
            classes = sorted(list(pieces.keys()))
            for pclass in classes:
                pieces[pclass]['Points'] = np.array(pieces[pclass]['Points']).flatten()
                pieces[pclass]['xyz'] = np.array(pieces[pclass]['xyz']).flatten()
            
            for pclass in classes:
                if pieces[pclass]['numberOfPoints'] > 0:
                    pieces2 = dict()
                    pieces2[pclass] = pieces[pclass]
                    fname = namespace + '_' + pclass +'_%08d'%(time*1e3)
                    writeVtkPointFile(fname, pieces2, binary=binary)
                    fname_written[fname] = True
    if binary:
        dataprecision = "Float64"
    else:
        dataprecision = "Float32"
    for pclass in pclasses:
        wrapper_namespace = outDir + os.sep + chid + '_prt5_out_' + pclass
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
                    f.write('            <PDataArray type="%s" Name="%s"/>\n'%(dataprecision, qty))
                f.write('        </PPointData>\n')
                f.write('        <PCellData></PCellData>\n')
                f.write('        <PPoints>\n')
                f.write('            <PDataArray type="%s" NumberOfComponents="3"/>\n'%(dataprecision))
                f.write('        </PPoints>\n')
                for mesh in uniqueMeshes:
                    source_name = chid + "_prt5_%04d_%s_%08d.vtp"%(int(mesh), pclass, time*1e3)
                    if fname_written[outDir + os.sep + source_name.replace('.vtp','')]:
                        f.write('        <Piece Source="%s"/>\n'%(source_name))
                f.write('    </PPolyData>\n')
                f.write('</VTKFile>\n')

def exportBndeDataToVtk(chid, resultDir, outtimes=None, outDir=None, binary=True):
    if outDir is None: outDir = resultDir
    smvFiles = getFileList(resultDir, chid, 'smv')
    
    file_quantities = getBndeQuantities(smvFiles[0])
    quantities = [file_quantities[key]['quantity'] for key in list(file_quantities.keys())]
    
    uniqueQuantities = list(set(quantities))
    
    gbfFiles = getFileList(resultDir, chid, 'gbf')
    gcfFiles = getFileList(resultDir, chid, 'gcf')
    beFiles = getFileList(resultDir, chid, 'be')
    
    if len(gcfFiles) == 0:
        print("No gcf files found in %s with chid %s"%(resultDir, chid))
        return None
    
    if len(beFiles) == 0:
        print("No be files found in %s with chid %s"%(resultDir, chid))
        return None
    
    #meshes = [x.split('_')[-2] for x in beFiles]
    meshes = [x.split('_')[-1].replace('.gcf','') for x in gcfFiles]
    uniqueMeshes = list(set(meshes))
    #outtimes = None
    if binary:
        dataprecision = "Float64"
    else:
        dataprecision = "Float32"
    
    for i in range(0, len(uniqueMeshes)):
        series_data = dict()
        for qty in uniqueQuantities:
            series_data[qty] = dict()
        mesh = uniqueMeshes[i]
        for j in range(0, len(meshes)):
            if mesh != meshes[j]:
                continue
            gcffile = gcfFiles[j]
            vertices, faces2, header1 = readGcfFile(gcffile)
            numberOfPoints = vertices.shape[0]
            numberOfPolys = faces2.shape[0]
            Points = vertices.flatten() #np.zeros((numberOfPoints*3,))
            ConnectData = faces2.flatten()-1 #np.zeros((numberOfPolys*3,))
            OffsetData = np.linspace(1, numberOfPolys, int(numberOfPolys))*3
            
            beFiles = getFileList(resultDir, chid+gcffile.split(chid)[-1].replace('.gcf','_'),'be')
            for befile in beFiles:
                times, vals, header2 = readBeFile(befile)
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
        namespace = outDir + os.sep + chid + '_bnde_%04d'%(int(mesh))
        if outtimes is None:
            times2 = times
            outtimes = times2
        else:
            times2 = outtimes
        writeVtkPolyTimeSeries(namespace, series_data, times2, binary=binary)
    
    
    wrapper_namespace = outDir + os.sep + chid + '_bnde_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtp', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PPolyData">\n')
            f.write('    <PPolyData GhostLevel="0">\n')
            f.write('        <PPointData></PPointData>\n')
            f.write('        <PCellData>\n')
            for qty in uniqueQuantities:
                f.write('            <PDataArray type="%s" Name="%s"/>\n'%(dataprecision, qty))
            f.write('        </PCellData>\n')
            f.write('        <PPoints>\n')
            f.write('            <PDataArray type="%s" NumberOfComponents="3"/>\n'%(dataprecision))
            f.write('        </PPoints>\n')
            for mesh in uniqueMeshes:
                source_name = chid + "_bnde_%04d_%08d.vtp"%(int(mesh), time*1e3)
                f.write('        <Piece Source="%s"/>\n'%(source_name))
            f.write('    </PPolyData>\n')
            f.write('</VTKFile>\n')

def exportS3dDataToVtk(chid, resultDir, outtimes=None, binary=False, decode=False, dtype='UInt8'):
    values, times = extractS3dValues(resultDir, chid, decode=decode)
    s3dfiles = getFileList(resultDir, chid, 's3d')
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    smvData = parseSMVFile(smvFile)
    (grids, obsts) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    
    meshes = [int(files['SMOKF3D'][x.split(os.sep)[-1]]['LINETEXT'].split()[1])-1 for x in s3dfiles]
    uniqueMeshes = list(set(meshes))
    
    quantities = [files['SMOKF3D'][x.split(os.sep)[-1]]['QUANTITY'] for x in s3dfiles]
    uniqueQuantities = list(set(quantities))
    
    xyzFiles = getFileListFromResultDir(resultDir, chid, 'xyz')
    grids = getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = getAbsoluteGrid(grids)
    
    gXB = [0, grids_abs.shape[0]-1, 0, grids_abs.shape[1]-1, 0, grids_abs.shape[2]-1]
    pieces=defaultdict(bool)
    for i in range(0, len(uniqueMeshes)):
        mesh = uniqueMeshes[i]
        xyzFile = False
        namespace = resultDir + os.sep + chid + '_s3d_grid%04d'%(int(mesh))
        for j in range(0, len(s3dfiles)):
            if mesh != meshes[j]:
                continue
            s3dfile = s3dfiles[j]
            qty = quantities[j]
            
            datas = values[mesh][qty]
            datas = np.swapaxes(np.swapaxes(np.swapaxes(datas, 0, 1),1,2),2,3)
            if dtype == 'Float64':
                datas = np.array(datas, dtype='float64')
            elif dtype == 'Float32':
                datas = np.array(datas, dtype='float32')
            elif dtype == 'UInt8':
                datas = np.array(datas, dtype='uint8')
            #lims, datas, times = fds.readSingleSlcfFile(slcfFile)
            
            if xyzFile is False:
                xyzFile = '_'.join(s3dfile.split('_')[:-1]) + '.xyz'
                grid, gridHeader = readXYZfile(xyzFile)
                xGrid, yGrid, zGrid = rearrangeGrid(grid)
                
                grid_combo = np.zeros((xGrid.shape[0], xGrid.shape[1], xGrid.shape[2],3))
                grid_combo[:, :, :, 0] = xGrid
                grid_combo[:, :, :, 1] = yGrid
                grid_combo[:, :, :, 2] = zGrid
                series_data = defaultdict(bool)
            if series_data[qty] is False:
                series_data[qty] = defaultdict(bool)
                series_data[qty]['data'] = datas
                series_data[qty]['times'] = times
                
                lims_xb = [xGrid[0, 0, 0], xGrid[-1, 0, 0],
                           yGrid[0, 0, 0], yGrid[0, -1, 0],
                           zGrid[0, 0, 0], zGrid[0, 0, -1]]
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
        if decode:
            dtype = 'Float32'
        else:
            dtype = 'UInt8'
        writeVtkTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB, binary=binary, dtype=dtype)
        pieces[namespace] = defaultdict(bool)
        pieces[namespace]['source'] = namespace
        pieces[namespace]['times'] = times2
        pieces[namespace]['lims'] = lXB
    if outtimes is None:
        return
    wrapper_namespace = resultDir + os.sep + chid + '_s3d_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtr', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PRectilinearGrid">\n')
            f.write('    <PRectilinearGrid WholeExtent="%d %d %d %d %d %d" GhostLevel="0">\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
            f.write('        <PPointData>\n')
            for qty in uniqueQuantities:
                f.write('            <PDataArray type="%s" Name="%s" NumberOfComponents="1"/>\n'%(dtype, qty))
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
    #exportBndfDataToVtk(chid, resultDir, outtimes=None)
    #exportPrt5DataToVtk(chid, resultDir, outtimes=None)
    #exportBndeDataToVtk(chid, resultDir, outtimes=None)
    exportS3dDataToVtk(chid, resultDir, outtimes=None)
    
    #obstToStl(resultDir, chid)
    
