#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example extracts all boundary data for named obstructions and
# outputs the maximum value for the queried quantity with time.
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
    qnty_mn, qnty_mx = [0, 1000]
    tStart, tEnd = 0, 120
    tInt, tBand = 1, 3
    yticks = np.linspace(0, 1000, 11)
    
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
    parser.add_argument('--orientations', help='axis orientations to include [0, -1, 1, -2, 2, -3, 3], 0 indicates all', type=int, nargs='*', default=[0])
    parser.add_argument('--yticks', help='ticks for plotting colorbar', type=float, nargs='*', default=yticks)
    
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
    orientations = args.orientations
    yticks = args.yticks
    
    if fds_file_name[:-4].lower() != '.fds': fds_file_name = fds_file_name + '.fds'
    fdsFilePath = os.path.join(working_dir, fds_file_name)
    
    smvFilePath = fds.getFileList(working_dir, chid, 'smv')[0]
    print("\t1. Reading boundary data.")
    datas = fds.extractMaxBndfValues(fdsFilePath, smvFilePath, working_dir, chid, quantities,
                    tStart=tStart, tEnd=tEnd, tInt=tInt, tBand=tBand, orientations=orientations)
    
    for i in range(0, len(quantities)):
        qty = quantities[i]
        
        times = datas[qty]['TIMES']
        mPts = datas[qty]['DATA']
        names = datas[qty]['NAMES']
        unit = datas[qty]['UNITS']
        
        outfile = os.path.join(outdir, '%s_%s_max_with_time'%(chid, qty))
        
        print("\t2. (%d/%d) Dumping %s to %s."%(i, len(quantities), qty, outfile+'.csv'))
        fds.maxValueCSV(times, mPts, names, outfile)
        print("\t2. (%d/%d) Plotting %s and saving to %s"%(i, len(quantities), qty, outfile+'.png'))
        fig = fds.maxValuePlot(times, mPts, names, outfile+'.png', vName="%s (%s)"%(qty, unit), yticks=yticks)
    
    plt.show()
    