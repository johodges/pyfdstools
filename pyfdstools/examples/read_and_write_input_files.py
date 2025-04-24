#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example demonstrates the functionality to read and write FDS
# input files using the pyfdstools library.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import pyfdstools as fds
import os, argparse

if __name__ == "__main__":
    print("Starting %s example"%(__file__))
    systemPath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(systemPath, "data", "case001.fds")
    outdir = os.path.join(systemPath, "generated")
    precision = 4
    
    parser = argparse.ArgumentParser(prog='ImportFile',
                                     description='Example routine to import an FDS input file into memory',
                                     epilog='Help goes here')
    parser.add_argument('--path', help='path to input file', type=str, default=path)
    parser.add_argument('--outdir', help='path to output directory', type=str, default=outdir)
    parser.add_argument('--precision', help='precision to use for numeric outputs', type=int, default=precision)
    args = parser.parse_args()
    
    if not os.path.exists(outdir): os.makedirs(outdir)
    print("\t1. Starting to import file from %s"%(path))
    file = fds.fdsFileOperations()
    file.importFile(args.path)
    print("\t2. Finished importing file.")
    
    print("\t3. Generating text output from read data.")
    text = file.generateFDStext(precision=args.precision)
    #print(text)
    
    location = os.path.join(outdir, '%s.fds'%(file.head['ID']['CHID']))
    
    file.saveModel(1, location)
    print("\t4. Input file written to: %s"%(location))
    