import pyfdstools as fds
import os
import numpy as np
import hashlib
import matplotlib.pyplot as plt
import stl
from collections import defaultdict

if __name__ == '__main__':
    
    chid = "case002"
    resultDir = "examples\\case002\\"
    binary = False
    
    '''
    (tmin, tmax) = (np.min(times), np.max(times))
    dt = 0.2
    nt = int(((tmax-tmin)/dt)+1)
    times_out = np.linspace(tmin, tmax, nt)
    '''
    
    #fds.exportSl3dDataToVtk(chid, resultDir, outtimes=None, binary=True, dtype='Float64')
    #exportBndfDataToVtk(chid, resultDir, outtimes=None)
    #exportPrt5DataToVtk(chid, resultDir, outtimes=None)
    #exportBndeDataToVtk(chid, resultDir, outtimes=None)
    
    if binary:
        fds.exportS3dDataToVtk(chid, resultDir, outtimes=None, binary=True, dtype='Float64')
    else:
        fds.exportS3dDataToVtk(chid, resultDir, outtimes=None, binary=False, dtype='UInt8')
    
    #fds.exportS3dDataToVtk(chid, resultDir, outtimes=None, binary=False)
    
    '''
    outtimes = None
    values, times = fds.extractS3dValues(resultDir, chid)
    s3dfiles = fds.getFileList(resultDir, chid, 's3d')
    smvFile = fds.getFileList(resultDir, chid, 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    (grids, obsts) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    
    meshes = [int(files['SMOKF3D'][x.split(os.sep)[-1]]['LINETEXT'].split()[1])-1 for x in s3dfiles]
    uniqueMeshes = list(set(meshes))
    
    quantities = [files['SMOKF3D'][x.split(os.sep)[-1]]['QUANTITY'] for x in s3dfiles]
    uniqueQuantities = list(set(quantities))
    
    xyzFiles = fds.getFileListFromResultDir(resultDir, chid, 'xyz')
    grids = fds.getGridsFromXyzFiles(xyzFiles, chid)
    grids_abs = fds.getAbsoluteGrid(grids)
    
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
            #lims, datas, times = fds.readSingleSlcfFile(slcfFile)
            
            if xyzFile is False:
                xyzFile = '_'.join(s3dfile.split('_')[:-1]) + '.xyz'
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
        writeVtkTimeSeries(namespace, grid_combo, series_data, times2, lXB, lXB)
        pieces[namespace] = defaultdict(bool)
        pieces[namespace]['source'] = namespace
        pieces[namespace]['times'] = times2
        pieces[namespace]['lims'] = lXB
    
    wrapper_namespace = resultDir + os.sep + chid + '_s3d_out'
    for time in outtimes:
        fname = wrapper_namespace + '_%08d'%(time*1e3)
        with open(fname + '.pvtr', 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="PRectilinearGrid">\n')
            f.write('    <PRectilinearGrid WholeExtent="%d %d %d %d %d %d" GhostLevel="0">\n'%(gXB[0], gXB[1], gXB[2], gXB[3], gXB[4], gXB[5]))
            f.write('        <PPointData>\n')
            for qty in uniqueQuantities:
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
    '''
    
    
    fds.obstToStl(resultDir, chid)
    