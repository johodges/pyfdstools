import pyfdstools as fds
import os
import numpy as np
import hashlib
import matplotlib.pyplot as plt
import stl
from collections import defaultdict

def obstToStl(resultDir, chid):
    smvFile = fds.getFileListFromResultDir(resultDir, chid, 'smv')[0]
    grid, obst, bndfs, surfs, files = fds.parseSMVFile(smvFile)
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
    
def writeVtkTimeSeries(namespace, grid, series_data, times, cell=False):
    nx, ny, nz, _ = np.shape(grid)
    gXB = [0, nx-1, 0, ny-1, 0, nz-1]
    lXB = [0, nx-1, 0, ny-1, 0, nz-1]
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
    quantities = keys
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
            f.write('    <PointData>\n')
            for qty in quantities:
                f.write('        <DataArray type="Float32" Name="%s" format="ascii">\n            '%(qty))
                values = piece[qty]
                for v in values:
                    f.write('%0.1f '%(v))
                f.write('\n         </DataArray>\n')
            f.write('    </PointData>\n')
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
                p.pop('PointData')
                for qty in list(series_data.keys()):
                    p[qty] = series_data[qty][mesh]['pieces'][i]['PointData'][:, tind]
                pieces.append(p)
        writeVtkPolyFile(fname, pieces)

def exportSl3dDataToVtk(chid, resultDr):
    # Slice data
    quantities, slcfFiles, dimensions, meshes, centers = fds.readSLCFquantities(chid, resultDir)
    
    uniqueQuantities = list(set(quantities))
    series_data = dict()
    tmax = -1e6
    tmin = 1e6
    for qty in uniqueQuantities:
        grid, data, times = fds.readSLCF3Ddata(chid, resultDir, qty)
        if grid is not False:
            series_data[qty] = dict()
            series_data[qty]['times'] = times
            series_data[qty]['grid'] = grid
            series_data[qty]['data'] = data
            tmax = max([tmax, np.max(times)])
            tmin = min([tmin, np.min(times)])
    
    grids = dict()
    grid_hashes = []
    for qty in list(series_data.keys()):
        sha256 = hashlib.sha256()
        sha256.update(series_data[qty]['grid'])
        gridhash = sha256.hexdigest()
        if gridhash not in grid_hashes:
            grid_hashes.append(gridhash)
            grids[gridhash] = series_data[qty]['grid']
        series_data[qty]['gridhash'] = gridhash
    
    dt = 0.2
    nt = int(((tmax-tmin)/dt)+1)
    times = np.linspace(tmin, tmax, nt)
    
    for i, gridhash in enumerate(grid_hashes):
        grid = grids[gridhash]
        namespace = resultDir + os.sep + chid + '_sl3d_grid%04d'%(i)
        writeVtkTimeSeries(namespace, grid, series_data, times)


def exportBndfDataToVtk(chid, resultDir):
    # Boundary data
    quantities = fds.readBoundaryQuantities(resultDir, chid)
    uniqueQuantities = list(set(quantities))
    smvFile = fds.getFileList(resultDir, chid, 'smv')[0]
    bndfs = fds.getFileList(resultDir, chid, 'bf')
    fdsFile = fds.fdsFileOperations()
    fdsFile.importFile(fds.getFileList(resultDir, chid, 'fds')[0])
    meshes = list(fdsFile.meshes.keys())
    smvGrids, smvObsts, smvBndfs, smvSurfs, smvFiles = fds.parseSMVFile(smvFile)
    bndf_dic = fds.linkBndfFileToMesh(meshes, bndfs, quantities)
    
    series_data = dict()
    for qty in uniqueQuantities:
        bndfs_points = bndf_dic[qty]
        series_data[qty] = dict()
        for file, mesh in bndfs_points:
            mesh = int(mesh)
            times, patches = fds.importBoundaryFile(file, gridNum=mesh, grid=smvGrids)
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
                piece['Offset'] = OffsetData
                piece['Connectivity'] = ConnectData
                pieces.append(piece)
            series_data[qty][mesh] = dict()
            series_data[qty][mesh]['pieces'] = pieces
            series_data[qty][mesh]['times'] = times
    namespace = resultDir + os.sep + chid + '_bndf'
    writeVtkPolyTimeSeries(namespace, series_data, times)
    


def exportPrt5DataToVtk(chid, resultDir):
    # Particle data
    particle_data, times_data = fds.importParticles(resultDir, chid)
    
    times = list(times_data.keys())
    namespace = resultDir + os.sep + chid + '_prt5'
    
    for time in times: # enumerate([list(times_data.keys())[0]]):
        fname = namespace + '_%08d'%(time*1e3)
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
        
        if nParts > 0:
            writeVtkPointFile(fname, pieces)
        


if __name__ == '__main__':
    
    chid = "case002"
    resultDir = "case002\\"
    
    exportSl3dDataToVtk(chid, resultDir)
    exportBndfDataToVtk(chid, resultDir)
    exportPrt5DataToVtk(chid, resultDir)
    
    obstToStl(resultDir, chid)
    