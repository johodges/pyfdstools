#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example extracts a boundary data from a queried orientation and
# location
#
#=======================================================================
# # IMPORTS
#=======================================================================
import pyfdstools as fds
import os, argparse
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    print("Starting %s example"%(__file__))
    systemPath = os.path.dirname(os.path.abspath(__file__))
    chid = "case001"
    path = os.path.join(systemPath, "data", chid+".fds")
    working_dir = os.path.join(systemPath, "data", "%s.zip"%(chid))
    outdir = os.path.join(systemPath, "generated")
    axis, value = -2, 4.4
    time, dt = -1, -1
    quantities = ['WALL TEMPERATURE']
    cbarticks = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    qnty_mn, qnty_mx = [20, 1000]
    
    parser = argparse.ArgumentParser(prog='Dump 2D slice to csv',
                                     description='Example routine to dump a 2D slice to a csv',
                                     epilog='Help goes here')
    parser.add_argument('--chid', help='CHID from FDS simulation', type=str, default=chid)
    parser.add_argument('--fds_file_name', help='Filename for FDS input file if different from <chid>.fds', type=str, default=chid)
    parser.add_argument('--working_dir', help='Directory containing FDS simulation results', type=str, default=working_dir)
    parser.add_argument('--quantities', help='quantity to query', type=str, nargs='*', default=quantities)
    parser.add_argument('--axis', help='axis to query (1=x, 2=y, 3=z)', type=int, default=axis)
    parser.add_argument('--value', help='axis value to query', type=float, default=value)
    parser.add_argument('--time', help='time to query, if -1 output all times', type=float, default=time)
    parser.add_argument('--dt', help='averaging window to apply, if -1 do not average', type=float, default=dt)
    parser.add_argument('--outdir', help='path to output directory', type=str, default=outdir)
    parser.add_argument('--qntymn', help='minimum value for plot', type=float, default=qnty_mn)
    parser.add_argument('--qntymx', help='minimum value for plot', type=float, default=qnty_mx)
    parser.add_argument('--cbarticks', help='ticks for plotting colorbar', type=float, nargs='*', default=cbarticks)
    
    args = parser.parse_args()
    
    chid = args.chid
    quantities = args.quantities
    axis = args.axis
    value = args.value
    time = args.time
    dt = args.dt
    qnty_mn = args.qntymn
    qnty_mx = args.qntymx
    fds_file_name = args.fds_file_name
    cbarticks = args.cbarticks
    
    if fds_file_name[:-4].lower() != '.fds': fds_file_name = fds_file_name + '.fds'
    fdsFilePath = os.path.join(working_dir, fds_file_name)
    
    print("\t1. Reading boundary data.")
    datas, times = fds.queryBndf(working_dir, chid, fdsFilePath, quantities, axis, value)
    
    if time == -1: time = times[-1]
    
    for i in range(0, len(quantities)):
        qty = quantities[i]
        print("\t2. (%d/%d) Plotting %s."%(i, len(quantities), qty))
        data = datas[qty]['DATA']
        x = datas[qty]['X']
        z = datas[qty]['Z']
        unit = datas[qty]['UNITS']
        
        if dt > 0:
            print("\t3. Averaging boundary data over time interval.")
            tStart = time-dt/2
            tEnd = time+dt/2
            tStartInd = np.argwhere(times >= tStart)[0][0]
            tEndInd = np.argwhere(times <= tEnd)[-1][0]
            meanData = np.mean(data[:, :, tStartInd:tEndInd], axis=2)
        else:
            print("\t3. No time averaging requested, interpolating nearest times to query time")
            tInd1 = np.argmin(abs(times-time))
            d1 = data[:, :, tInd1]
            if times[tInd1] > time:
                tInd2 = max([tInd1-1,0])
            else:
                tInd2 = min([tInd1+1,data.shape[2]-1])
            d2 = data[:, :, tInd2]
            trange = abs(times[tInd2]-times[tInd1])
            if trange > 0:
                weight1 = abs(times[tInd1]-time)/trange
                weight2 = abs(times[tInd2]-time)/trange
                meanData = d1*weight1 + d2*weight2
            else:
                meanData = d1
            
        print("\t4. Plotting requested frame")
        fig, ax = fds.plotSlice(x, z, meanData, axis,
                            qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                            clabel="%s (%s)"%(qty, unit),
                            cbarticks=cbarticks)
        outfile = os.path.join(outdir, '%s_%s_%0.0f_%0.4f_t%0.1f.png'%(chid, qty, axis, value, time))
        print("\t4. Saving frame of extracted boundary data to file %s."%(outfile))
        fig.savefig(outfile)
    plt.show()
    