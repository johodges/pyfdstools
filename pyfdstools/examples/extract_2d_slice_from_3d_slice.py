#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example extracts a queried 2-D slice from 3-D slice data.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import pyfdstools as fds
import os, argparse
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("Starting %s example"%(__file__))
    systemPath = os.path.dirname(os.path.abspath(__file__))
    chid = "case001"
    path = os.path.join(systemPath, "data", chid+".fds")
    working_dir = os.path.join(systemPath, "data", "%s.zip"%(chid))
    outdir = os.path.join(systemPath, "generated")
    axis, value = 2, 4.4
    time, dt = 30, 60
    quantity ='TEMPERATURE'
    
    parser = argparse.ArgumentParser(prog='Dump 2D slice to csv',
                                     description='Example routine to dump a 2D slice to a csv',
                                     epilog='Help goes here')
    parser.add_argument('--chid', help='CHID from FDS simulation', type=str, default=chid)
    parser.add_argument('--working_dir', help='Directory containing FDS simulation results', type=str, default=working_dir)
    parser.add_argument('--quantity', help='quantity to query', type=str, default=quantity)
    parser.add_argument('--axis', help='axis to query (1=x, 2=y, 3=z)', type=int, default=axis)
    parser.add_argument('--value', help='axis value to query', type=float, default=value)
    parser.add_argument('--time', help='time to query, if -1 output all times', type=float, default=time)
    parser.add_argument('--dt', help='averaging window to apply, if -1 do not average', type=float, default=dt)
    parser.add_argument('--outdir', help='path to output directory', type=str, default=outdir)
    parser.add_argument('--qntymn', help='minimum value for plot', type=float, default=None)
    parser.add_argument('--qntymx', help='minimum value for plot', type=float, default=None)
    
    args = parser.parse_args()
    
    chid = args.chid
    quantity = args.quantity
    axis = args.axis
    value = args.value
    time = args.time
    dt = args.dt
    qnty_mn = args.qntymn
    qnty_mx = args.qntymx
    
    if time == -1: time = None
    if dt == -1: dt = None
    
    print("\t1. Reading 3D slice data.")
    grid, data, times, unit = fds.readSLCF3Ddata(chid, working_dir, quantity, time=time, dt=dt)
    
    print("\t2. Finding queried 2D slice in 3D slice data.")
    x, z, data_slc = fds.findSliceLocation(grid, data[:,:,:,-1], axis, value)
    
    print("\t3. Plotting final frame of extracted 2D slice.")
    fig, ax = fds.plotSlice(x, z, data_slc, axis,
                        qnty_mn=qnty_mn, qnty_mx=qnty_mx,
                        clabel="%s (%s)"%(quantity, unit))
    outfile = os.path.join(outdir, '%s_%s_%0.0f_%0.4f_final_frame.png'%(chid, quantity, axis, value))
    print("\t4. Saving final frame of extracted 2D slice to file %s."%(outfile))
    fig.savefig(outfile)
    plt.show()
    
    
    