#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example dumps a queried 2d slice to a csv file.
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
    resultDir = os.path.join(systemPath, "data", "%s.zip"%(chid))
    outdir = os.path.join(systemPath, "generated")
    axis, value = 1, 2.55
    time, dt = 30, 60
    quantity ='TEMPERATURE'
    
    parser = argparse.ArgumentParser(prog='Dump 2D slice to csv',
                                     description='Example routine to dump a 2D slice to a csv',
                                     epilog='Help goes here')
    parser.add_argument('--chid', help='CHID from FDS simulation', type=str, default=chid)
    parser.add_argument('--quantity', help='quantity to query', type=str, default=quantity)
    parser.add_argument('--axis', help='axis to query (1=x, 2=y, 3=z)', type=int, default=axis)
    parser.add_argument('--value', help='axis value to query', type=float, default=value)
    parser.add_argument('--time', help='time to query, if -1 output all times', type=float, default=time)
    parser.add_argument('--dt', help='averaging window to apply, if -1 do not average', type=float, default=dt)
    parser.add_argument('--outdir', help='path to output directory', type=str, default=outdir)
    
    args = parser.parse_args()
    
    chid = args.chid
    quantity = args.quantity
    axis = args.axis
    value = args.value
    time = args.time
    dt = args.dt
    
    if time == -1: time = None
    if dt == -1: dt = None
    
    data, unit = fds.query2dAxisValue(resultDir, chid, quantity, axis, value, time=time, dt=dt)
    
    fds.renderSliceCsvs(data, chid, outdir)
    
    fig, ax = fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], axis,
                        clabel="%s (%s)"%(quantity, unit))
    fig.savefig(os.path.join(outdir, '%s_%s_%0.0f_%0.4f_final_frame.png'%(chid, quantity, axis, value)))
    plt.show()