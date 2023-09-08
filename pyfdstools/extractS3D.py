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
import matplotlib.pyplot as plt
import glob
import zipfile
import os
import struct
import scipy.interpolate as scpi
import pandas as pd
from collections import defaultdict
from .utilities import getDatatypeByEndianness, getEndianness
from .utilities import getFileListFromZip, getFileList, zopen, zreadlines
from .utilities import getFileListFromResultDir
from .utilities import getGridsFromXyzFiles, getAbsoluteGrid, rearrangeGrid
from .utilities import readXYZfile
from .colorSchemes import buildSMVcolormap
from .smokeviewParser import parseSMVFile
from itertools import groupby

def readS3dTime(data, time_ind, datashape, float32, int32):
    time = np.frombuffer(data[time_ind:time_ind+4], dtype=float32)[0]
    nchars_in = np.frombuffer(data[time_ind+12:time_ind+16], dtype=int32)[0]
    nchars_out = np.frombuffer(data[time_ind+16:time_ind+20], dtype=int32)[0]
    
    i = 0
    mark = np.uint8(255)
    out_pos = 0
    decoded_data = np.empty((datashape[0]*datashape[1]*datashape[2]), dtype=np.uint8)
    while i < nchars_out:
        # Credit to fdsreader for this parsing
        # https://github.com/FireDynamics/fdsreader/blob/master/fdsreader/smoke3d/smoke3d.py
        if data[time_ind+28+i] == mark:
            value = data[time_ind+28+i+1]
            repeats = data[time_ind+28+i+2]
            i += 3
        else:
            value = data[time_ind+28+i]
            repeats = 1
            i += 1
        decoded_data[out_pos:out_pos+repeats] = value
        out_pos += repeats
    return time, decoded_data.reshape(datashape, order='F'), nchars_out

def readS3dFile(file):
    f = zopen(file, 'rb')
    data = f.read()
    f.close()
    if struct.unpack('<i',data[4:8])[0] == 1:
        endianness = "<"
    elif struct.unpack('>i',data[4:8])[0] == 1:
        endianness = ">"
    else:
        print("Warning unable to determine endianness, assuming little endian")
        endianness = "<"
    float32 = np.dtype(np.float32).newbyteorder(endianness)
    int32 = np.dtype(np.int32).newbyteorder(endianness)
    header = np.frombuffer(data[:44], dtype=int32)
    (NX, NY, NZ) = (header[4]+1, header[6]+1, header[8]+1)
    data_shape = (NX, NY, NZ)
    ordered_data = []
    times = []
    
    time_ind = 44
    
    while time_ind < len(data):
        time, data_out, nchars_out = readS3dTime(data, time_ind, data_shape, float32, int32)
        ordered_data.append(data_out)
        times.append(time)
        time_ind = time_ind+36+nchars_out
    times = np.array(times, dtype=np.float32)
    ordered_data = np.array(ordered_data, dtype=np.uint8)
    return times, ordered_data

def extractS3dValues(resultDir, chid, decode=True):
    values = defaultdict(bool)
    s3dfiles = getFileList(resultDir, chid, 's3d')
    smvFile = getFileList(resultDir, chid, 'smv')[0]
    smvData = parseSMVFile(smvFile)
    (grids, obsts) = (smvData['grids'], smvData['obsts'])
    (bndfs, surfs) = (smvData['bndfs'], smvData['surfs'])
    (files, bndes) = (smvData['files'], smvData['bndes'])
    
    for s3dfile in s3dfiles:
        times, ordered_data = readS3dFile(s3dfile)
        file = s3dfile.split(os.sep)[-1]
        quantity = files['SMOKF3D'][file]['QUANTITY']
        f_ordered_data = np.array(ordered_data, dtype=np.float32)
        meshNum = int(files['SMOKF3D'][file]['LINETEXT'].split()[1])-1
        if decode:
            if quantity == 'SOOT DENSITY':
                grid = grids[meshNum]
                dx = np.median(grid[0][1:,1]-grid[0][:-1,1])
                dy = np.median(grid[1][1:,1]-grid[1][:-1,1])
                dz = np.median(grid[2][1:,1]-grid[2][:-1,1])
                dd = (dx*dy*dz)**(1/3)
                MAX_SMV = 8700
                MASS_EXTINCTION_COEFFICIENT = MAX_SMV
                factor = -1*(MASS_EXTINCTION_COEFFICIENT) * dd
                val_fds = np.log(-1*((f_ordered_data / 254) - 1.0))/factor
            if quantity == 'HRRPUV':
                MAX_SMV = 1200
                val_fds = (f_ordered_data/254)*MAX_SMV
            if quantity == 'TEMPERATURE':
                MAX_SMV = 2000
                TMPA = 20 + 273.15
                TMPM = 273.15
                val_fds = (f_ordered_data/254)*(MAX_SMV-(TMPA-TMPM)) + (TMPA-TMPM)
        else:
            val_fds = f_ordered_data
        if values[meshNum] is False:
            values[meshNum] = defaultdict(bool)
        values[meshNum][quantity] = val_fds
    return values, times

def encode_list(s_list):
    return [[len(list(group)), key] for key, group in groupby(s_list)]

def encodeS3dData(data, quantity, dx):
    if quantity == 'SOOT DENSITY':
        MAX_SMV = 8700
        MASS_EXTINCTION_COEFFICIENT = MAX_SMV
        factor = -1*(MASS_EXTINCTION_COEFFICIENT) * dx
        encoded_data = np.array(np.round((1-np.exp(data*factor))*254), dtype=np.uint8)
    if quantity == 'HRRPUV':
        MAX_SMV = 1200
        encoded_data = np.array(np.round((data/MAX_SMV)*254), dtype=np.uint8)
    if quantity == 'TEMPERATURE':
        MAX_SMV = 2000
        TMPA = 20 + 273.15
        TMPM = 273.15
        TMP_MIN = TMPA-TMPM
        encoded_data = np.array(np.round(254*(data - TMP_MIN)/(MAX_SMV-TMP_MIN)), dtype=np.uint8)
    return encoded_data

def writeS3dFile(file, times, data_out, quantity, dx):
    endianness = "<"
    float32 = np.dtype(np.float32).newbyteorder(endianness)
    int32 = np.dtype(np.int32).newbyteorder(endianness)
    
    f = zopen(file, 'wb')
    f.write(np.array([32, 1, 0], dtype=int32).tobytes())
    NT, NX, NY, NZ = data_out.shape
    f.write(np.array([0, NX-1, 0, NY-1, 0, NZ-1], dtype=int32).tobytes())
    f.write(np.array([32], dtype=int32).tobytes())
    
    nchars_in = NX*NY*NZ
    encoded_data = encodeS3dData(data_out, quantity, dx)
    for i, time in enumerate(times):
        f.write(np.array([4], dtype=int32).tobytes())
        f.write(times[i].tobytes())
        f.write(np.array([4, 8], dtype=int32).tobytes())
        
        #times[i].tobytes()+np.array([4], dtype=int32).tobytes()+
        array_data = encoded_data[i, :, :, :].reshape((nchars_in), order='F')
        encoded_list = encode_list(array_data)
        
        # Hack job to calculate number of chars before going through the looping logic
        count = [x[0] for x in encoded_list]
        counter = np.sum([x for x in count if x <4])
        count = [x for x in count if x > 3]
        counter = counter + np.sum([np.floor(x/254) for x in count])*3
        count = [x-np.floor(x/254)*254 for x in count]
        counter = counter + np.sum([x for x in count if x <4])
        count = [x for x in count if x > 3]
        nchars_out = counter + len(count)*3
        
        f.write(np.array([nchars_in, nchars_out, 8, nchars_out], dtype=int32).tobytes())
        tmp = b""
        for en in encoded_list:
            #if en[0] > 254:
            while en[0] > 254:
                #tmp = tmp + b'\xff' + en[1].tobytes() + b'\xfe'
                f.write(b'\xff' + en[1].tobytes() + b'\xfe')
                en[0] = en[0]-254
            if en[0] > 3:
                #tmp = tmp + b'\xff' + en[1].tobytes() + np.uint8(en[0]).tobytes()
                f.write(b'\xff' + en[1].tobytes() + np.uint8(en[0]).tobytes())
                en[0] = 0
            while en[0] > 0:
                #tmp = tmp + en[1].tobytes()
                f.write(en[1].tobytes())
                en[0] = en[0] - 1
        f.write(np.array([nchars_out], dtype=int32).tobytes())
        
    f.close()
    
def writeS3dFile_debug(file, debug_file, times, data_out, quantity, dx):
    endianness = "<"
    float32 = np.dtype(np.float32).newbyteorder(endianness)
    int32 = np.dtype(np.int32).newbyteorder(endianness)
    
    f = zopen(debug_file, 'rb')
    debug_data = f.read()
    f.close()
    
    f = zopen(file, 'wb')
    f.write(np.array([32, 1, 0], dtype=int32).tobytes())
    NT, NX, NY, NZ = data_out.shape
    f.write(np.array([0, NX-1, 0, NY-1, 0, NZ-1], dtype=int32).tobytes())
    f.write(np.array([32], dtype=int32).tobytes())
    
    nchars_in = NX*NY*NZ
    encoded_data = encodeS3dData(data_out, quantity, dx)
    f.close()
    for i, time in enumerate(times):
        f = zopen(file, 'ab')
        f.write(np.array([4], dtype=int32).tobytes())
        f.write(times[i].tobytes())
        f.write(np.array([4, 8], dtype=int32).tobytes())
        f.close()
        
        #times[i].tobytes()+np.array([4], dtype=int32).tobytes()+
        array_data = encoded_data[i, :, :, :].reshape((nchars_in), order='F')
        encoded_list = encode_list(array_data)
        
        # Hack job to calculate number of chars before going through the looping logic
        count = [x[0] for x in encoded_list]
        counter = np.sum([x for x in count if x <4])
        count = [x for x in count if x > 3]
        counter = counter + np.sum([np.floor(x/254) for x in count])*3
        count = [x-np.floor(x/254)*254 for x in count]
        counter = counter + np.sum([x for x in count if x <4])
        count = [x for x in count if x > 3]
        nchars_out = counter + len(count)*3
        
        f = zopen(file, 'ab')
        f.write(np.array([nchars_in, nchars_out, 8, nchars_out], dtype=int32).tobytes())
        f.close()
        tmp = b""
        for en in encoded_list:
            #if en[0] > 254:
            while en[0] > 254:
                #tmp = tmp + b'\xff' + en[1].tobytes() + b'\xfe'
                f = zopen(file, 'ab')
                f.write(b'\xff' + en[1].tobytes() + b'\xfe')
                f.close()
                f = zopen(file, 'rb')
                data = f.read()
                f.close()
                l = len(data)
                if debug_data[:l] != data:
                    assert False, "Stopped"
                en[0] = en[0]-254
            if en[0] > 3:
                #tmp = tmp + b'\xff' + en[1].tobytes() + np.uint8(en[0]).tobytes()
                f = zopen(file, 'ab')
                f.write(b'\xff' + en[1].tobytes() + np.uint8(en[0]).tobytes())
                f.close()
                f = zopen(file, 'rb')
                data = f.read()
                f.close()
                l = len(data)
                if debug_data[:l] != data:
                    assert False, "Stopped"
                en[0] = 0
            while en[0] > 0:
                #tmp = tmp + en[1].tobytes()
                f = zopen(file, 'ab')
                f.write(en[1].tobytes())
                f.close()
                f = zopen(file, 'rb')
                data = f.read()
                f.close()
                l = len(data)
                if debug_data[:l] != data:
                    assert False, "Stopped"
                en[0] = en[0] - 1
        f = zopen(file, 'ab')
        f.write(np.array([nchars_out], dtype=int32).tobytes())
        f.close()
    