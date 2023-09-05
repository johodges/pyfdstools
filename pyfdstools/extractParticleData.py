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
from collections import defaultdict
from .utilities import zopen, zreadlines, getFileList

def importParticles(resultDir, chid):
    partFiles = getFileList(resultDir, chid, 'prt5')
    particle_data = defaultdict(bool)
    particle_data['tags'] = defaultdict(bool)
    particle_data['classes'] = defaultdict(bool)
    allTimes = []
    for file in partFiles:
        p_data, times = importParticle(file)
        for key in list(p_data['tags'].keys()):
            particle_data['tags'][key] = p_data['tags'][key]
        allTimes.append(times)
        classes = list(p_data['classes'].keys())
        classes.remove('times')
        for key in classes:
            if particle_data['classes'][key] is False:
                particle_data['classes'][key] = defaultdict(bool)
                particle_data['classes'][key]['nPart'] = np.array(p_data['classes'][key]['nPart'])
                particle_data['classes'][key]['nQty'] = p_data['classes'][key]['nQty']
                particle_data['classes'][key]['qtyNames'] = p_data['classes'][key]['qtyNames']
                particle_data['classes'][key]['qtyUnits'] = p_data['classes'][key]['qtyUnits']
            else:
                #print(key)
                #print("dict", particle_data['classes'][key]['nPart'])
                #print("add", np.array(p_data['classes'][key]))
                particle_data['classes'][key]['nPart'] += np.array(p_data['classes'][key]['nPart'])
                
    allTimes = np.unique(np.squeeze(allTimes))
    
    times_data = defaultdict(bool)
    for time in allTimes:
        times_data[time] = defaultdict(bool)
        times_data[time]['particles'] = []
    
    particles = list(particle_data['tags'].keys())
    for particle in particles:
        for t in particle_data['tags'][particle]['times']:
            times_data[t]['particles'].append(particle)
    
    particle_data['times'] = allTimes
    
    return particle_data, times_data

def importParticle_meta(file):
    txt = zreadlines(file.replace('.prt5','.prt5.bnd'))
    times = [float(line.split()[0]) for line in txt if line[1] != ' ']
    f = zopen(file, 'rb')
    data = f.read()
    f.close()
    header = np.frombuffer(data[:32], dtype=np.int32)
    endianness = header[1]
    fdsVersion = header[4]
    nPart_class = header[7]
    
    counter = 40
    particle_data = dict()
    for i in range(0, nPart_class):
        pid = 'N-%04d'%(i)
        particle_data[pid] = dict()
        n_part_info = np.frombuffer(data[counter:counter+8], np.int32)
        nQty = n_part_info[0]
        particle_data[pid]['nQty'] = nQty
        counter = counter+8
        qtyNames = []
        qtyUnits = []
        for j in range(0, nQty):
            counter = counter+8
            qtyName = data[counter:counter+30].decode('utf-8')
            counter = counter+38
            qtyUnit = data[counter:counter+30].decode('utf-8')
            counter = counter+30
            qtyNames.append(qtyName.strip())
            qtyUnits.append(qtyUnit.strip())
        counter = counter+8
        particle_data[pid]['qtyNames'] = qtyNames
        particle_data[pid]['qtyUnits'] = qtyUnits
        particle_data[pid]['nPart'] = []
        particle_data[pid]['times'] = []
    
    particle_dict = defaultdict(bool)
    particle_dict['classes'] = particle_data
    particle_dict['tags'] = defaultdict(bool)
    particle_data['times'] = []
    for time in times:
        time2 = np.frombuffer(data[counter:counter+4], np.float32)[0]
        #tid = 'T=%0.8f'%(time)
        counter = counter + 12
        if abs(time-time2) > 0.01:
            print(np.frombuffer(data[counter-32:counter+32], np.float32))
            print("Error imported time %0.4f does not match meta data time %0.4f"%(time2, time))
            assert False, "Stopped"
        
        for i in range(0, nPart_class):
            pid = 'N-%04d'%(i)
            #particle_data[pid][tid] = dict()
            nPart = np.frombuffer(data[counter:counter+16], np.int32)[0]
            counter = counter+12
            nQty = particle_data[pid]['nQty']
            qtyNames = particle_data[pid]['qtyNames']
            qtyUnits = particle_data[pid]['qtyUnits']
            #print('nPart_class: ', i, 'npart', nPart)
            
            particle_data[pid]['nPart'].append(nPart)
            particle_data[pid]['times'].append(time)
            if nPart > 0:
                #particle_data[pid][tid]['nPart'] = nPart
                
                xyzs = np.frombuffer(data[counter:counter+12*nPart], np.float32)
                counter = counter + 12*nPart + 8
                
                tags = np.frombuffer(data[counter:counter+4*nPart], np.int32)
                counter = counter + 4*nPart + 0
                
                if nQty > 0:
                    counter = counter + 8
                    qtyValues = np.frombuffer(data[counter:counter+4*nPart*nQty], np.float32)
                    counter = counter + 4*nPart*nQty
                else:
                    qtyValues = []
                    counter = counter + 0
                counter = counter + 8 
                '''
                print('\t', time, i, ) 
                print('\t\txyz:', xyzs) 
                print('\t\ttags:', tags) 
                print('\t\tqtyNames:', qtyNames) 
                print('\t\tqtyValues:',qtyValues)
                '''
                for j in range(0, nPart):
                    ppid = pid + "-%04d"%(j)
                    #particle_data[pid][tid][ppid] = dict()
                    xyz = xyzs[j::nPart]
                    tag = tags[j]
                    qtyValue = qtyValues[j::nPart]
                    
                    if particle_dict['tags'][tag] is not False:
                        particle_dict['tags'][tag]['times'].append(time)
                        particle_dict['tags'][tag]['xyz'].append(xyz)
                        for k in range(0, nQty):
                            qid = qtyNames[k] + '(%s)'%(qtyUnits[k])
                            try:
                                particle_dict['tags'][tag][qid].append(qtyValue[k])
                            except IndexError:
                                print(particle_dict['tags'][tag], qid)
                                print(k, qtyValue)
                                assert False, "Stopped"
                    else:
                        particle_dict['tags'][tag] = defaultdict(bool)
                        particle_dict['tags'][tag]['times'] = [time]
                        particle_dict['tags'][tag]['xyz'] = [xyz]
                        for k in range(0, nQty):
                            qid = qtyNames[k] + '(%s)'%(qtyUnits[k])
                            particle_dict['tags'][tag][qid] = [qtyValue[k]]
                        particle_dict['tags'][tag]['class'] = pid
                    #particle_data[pid][tid][ppid]['xyz'] = xyz
                counter = counter + 0
            else:
                if nQty > 0:
                    counter = counter + 8
                counter = counter + 16
                pass
        counter = counter + 0
    
    return particle_dict, times

def importParticle(file):
    f = zopen(file, 'rb')
    data = f.read()
    f.close()
    header = np.frombuffer(data[:32], dtype=np.int32)
    endianness = header[1]
    fdsVersion = header[4]
    nPart_class = header[7]
    
    counter = 40
    particle_data = dict()
    for i in range(0, nPart_class):
        pid = 'N-%04d'%(i)
        particle_data[pid] = dict()
        n_part_info = np.frombuffer(data[counter:counter+8], np.int32)
        nQty = n_part_info[0]
        particle_data[pid]['nQty'] = nQty
        counter = counter+8
        qtyNames = []
        qtyUnits = []
        for j in range(0, nQty):
            counter = counter+8
            qtyName = data[counter:counter+30].decode('utf-8')
            counter = counter+38
            qtyUnit = data[counter:counter+30].decode('utf-8')
            counter = counter+30
            qtyNames.append(qtyName.strip())
            qtyUnits.append(qtyUnit.strip())
        counter = counter+8
        particle_data[pid]['qtyNames'] = qtyNames
        particle_data[pid]['qtyUnits'] = qtyUnits
        particle_data[pid]['nPart'] = []
        particle_data[pid]['times'] = []
    
    particle_dict = defaultdict(bool)
    particle_dict['classes'] = particle_data
    particle_dict['tags'] = defaultdict(bool)
    particle_data['times'] = []
    times = []
    while counter < len(data):
        time = np.frombuffer(data[counter:counter+4], np.float32)[0]
        times.append(time)
        counter = counter + 12
        
        for i in range(0, nPart_class):
            pid = 'N-%04d'%(i)
            #particle_data[pid][tid] = dict()
            nPart = np.frombuffer(data[counter:counter+16], np.int32)[0]
            counter = counter+12
            nQty = particle_data[pid]['nQty']
            qtyNames = particle_data[pid]['qtyNames']
            qtyUnits = particle_data[pid]['qtyUnits']
            #print('nPart_class: ', i, 'npart', nPart)
            
            particle_data[pid]['nPart'].append(nPart)
            particle_data[pid]['times'].append(time)
            if nPart > 0:
                #particle_data[pid][tid]['nPart'] = nPart
                
                xyzs = np.frombuffer(data[counter:counter+12*nPart], np.float32)
                counter = counter + 12*nPart + 8
                
                tags = np.frombuffer(data[counter:counter+4*nPart], np.int32)
                counter = counter + 4*nPart + 0
                
                if nQty > 0:
                    counter = counter + 8
                    qtyValues = np.frombuffer(data[counter:counter+4*nPart*nQty], np.float32)
                    counter = counter + 4*nPart*nQty
                else:
                    qtyValues = []
                    counter = counter + 0
                counter = counter + 8 
                '''
                print('\t', time, i, ) 
                print('\t\txyz:', xyzs) 
                print('\t\ttags:', tags) 
                print('\t\tqtyNames:', qtyNames) 
                print('\t\tqtyValues:',qtyValues)
                '''
                for j in range(0, nPart):
                    ppid = pid + "-%04d"%(j)
                    #particle_data[pid][tid][ppid] = dict()
                    xyz = xyzs[j::nPart]
                    tag = tags[j]
                    qtyValue = qtyValues[j::nPart]
                    
                    if particle_dict['tags'][tag] is not False:
                        particle_dict['tags'][tag]['times'].append(time)
                        particle_dict['tags'][tag]['xyz'].append(xyz)
                        for k in range(0, nQty):
                            qid = qtyNames[k] + '(%s)'%(qtyUnits[k])
                            try:
                                particle_dict['tags'][tag][qid].append(qtyValue[k])
                            except IndexError:
                                print(particle_dict['tags'][tag], qid)
                                print(k, qtyValue)
                                assert False, "Stopped"
                    else:
                        particle_dict['tags'][tag] = defaultdict(bool)
                        particle_dict['tags'][tag]['times'] = [time]
                        particle_dict['tags'][tag]['xyz'] = [xyz]
                        for k in range(0, nQty):
                            qid = qtyNames[k] + '(%s)'%(qtyUnits[k])
                            particle_dict['tags'][tag][qid] = [qtyValue[k]]
                        particle_dict['tags'][tag]['class'] = pid
                    #particle_data[pid][tid][ppid]['xyz'] = xyz
                counter = counter + 0
            else:
                if nQty > 0:
                    counter = counter + 8
                counter = counter + 16
                pass
        counter = counter + 0
    times = np.array(times, dtype=np.float32)
    return particle_dict, times