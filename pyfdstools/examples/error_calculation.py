#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This example plots the uncertainty statistics for a parameter based
# on the FDS validation statistics.
#
#=======================================================================
# # IMPORTS
#=======================================================================
import pyfdstools as fds
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":
    print("Starting %s example"%(__file__))
    data, keys = fds.readErrorTable()
    
    values = [100, 200, 300, 400, 500, 600]
    meanValue = 500
    quantity = 'Plume Temperature'
    percentile = 0.95
    fdsVersion = '6.7.1'
    
    parser = argparse.ArgumentParser(prog='Error Calculation',
                                     description='Example routine to apply FDS uncertainty to values',
                                     epilog='Help goes here')
    parser.add_argument('--values', help='list of values', type=float, nargs='*', default=values)
    parser.add_argument('--quantity', help='quantity for uncertainty statistics', type=str, default=quantity)
    parser.add_argument('--percentile', help='percentile to compute for each of values', type=float, default=percentile)
    parser.add_argument('--fdsVersion', help='fds version for uncertainty statistics', type=str, default=fdsVersion)
    parser.add_argument('--meanValue', help='single value to show distribution of uncertainty', type=float, default=meanValue)
    
    args = parser.parse_args()
    
    values = args.values
    quantity = args.quantity
    percentile = args.percentile
    fdsVersion = args.fdsVersion
    meanValue = args.meanValue
    
    errorvalues = fds.calculatePercentile(values, quantity, percentile, fdsVersion=fdsVersion)
    print("\t1. %0.4f percentiles"%(percentile))
    print("\t2. Values:      %s"%(','.join(['%0.1f'%(x) for x in values])))
    print("\t3. Percentiles: %s"%(','.join(['%0.1f'%(x) for x in errorvalues])))
    fig, ax1 = fds.plotPercentile(meanValue, quantity, fdsVersion=fdsVersion)
    plt.show()
    