# -*- coding: utf-8 -*-
"""
Created on Thu May 30 11:34:24 2019

@author: jhodges
"""

import matplotlib.pyplot as plt
import scipy.stats as scst
import numpy as np
import pandas as pd
import os

def readErrorTable(fdsVersion='6.7.1'):
    dataDir = os.sep.join(__file__.split(os.sep)[:-1])
    data = pd.read_csv("%s//fdsErrorTables//%s.csv"%(dataDir, fdsVersion), index_col=0)
    #data = pd.read_csv("fdsErrorTables//%s.csv"%(fdsVersion), index_col=0)
    keys = list(data.index.values)
    return data, keys

def calculatePercentile(values, quantity, percentile, fdsVersion='6.7.1'):
    data, quantities = readErrorTable(fdsVersion=fdsVersion)
    
    if (quantity in quantities):
        fdsBias = data.loc[quantity]['Bias']
        fdsSigma = data.loc[quantity]['sigmaM']
        values = np.array(values)
        errorValues = np.zeros_like(values)
        for i in range(0, values.shape[0]):
            mu = values[i]/fdsBias
            sigma = values[i]/fdsBias*fdsSigma
            x = np.linspace(mu-3*sigma,mu+3*sigma, 101)
            y_cdf = scst.norm.cdf(x, mu, sigma)
            ind = np.where(y_cdf > percentile)[0][0]
            errorValues[i] = x[ind]
        return errorValues
        
    else:
        print("Quantity '%s' not known."%(quantity))
        print("Known Quantities:")
        for qty in quantities:
            print("\t%s"%(qty))
    
def plotPercentile(value, quantity, fdsVersion='6.7.1'):
    data, quantities = readErrorTable(fdsVersion=fdsVersion)
    if (quantity in quantities):
        fdsBias = data.loc[quantity]['Bias']
        fdsSigma = data.loc[quantity]['sigmaM']
        mu = value/fdsBias
        sigma = value/fdsBias*fdsSigma
        x = np.linspace(mu-3*sigma,mu+3*sigma, 101)
        y_pdf = scst.norm.pdf(x, mu, sigma)
        y_cdf = scst.norm.cdf(x, mu, sigma)
        
        fs = 16
        lw = 3
        fig = plt.figure(figsize=(12,6))
        ax1 = fig.add_subplot(111)
        
        ax1.plot([value, value], [y_pdf.min(), y_pdf.max()],'--k', linewidth=lw, label='Predicted')
        ax1.plot(x, y_pdf, label='PDF', color=[134/255,31/255,65/255,1.0], linewidth=lw,)
        ax1.set_xlabel('%s'%(quantity),fontsize=fs)
        ax1.set_ylabel('PDF Probabilty Density',fontsize=fs, color=[134/255,31/255,65/255,1.0])
        ax1.tick_params('y', colors=[134/255,31/255,65/255,1.0], labelsize=fs)
        ax1.tick_params('x', labelsize=fs)
        
        ax2 = ax1.twinx()
        ax2.plot(x, y_cdf, label='CDF', color=[232/255,119/255,34/255,1.0], linewidth=lw,)
        ax2.tick_params('y', colors=[232/255,119/255,34/255,1.0], labelsize=fs)
        ax2.set_ylabel('CDF Probabilty',fontsize=fs, color=[232/255,119/255,34/255,1.0],rotation=270, labelpad=20)
        
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        
        ax2.legend(lines+lines2, labels+labels2, fontsize=fs)
        plt.tight_layout()
        return fig, ax1
        
    else:
        print("Quantity '%s' not known."%(quantity))
        print("Known Quantities:")
        for qty in quantities:
            print("\t%s"%(qty))
    
    
def getQuantities(fdsVersion='6.7.1'):
    data = pd.read_csv("fdsErrorTables//%s.csv"%(fdsVersion))
    keys = list(data.keys())
    return keys
    
if __name__ == "__main__":
    
    data, keys = readErrorTable()
    
    errorvalues = calculatePercentile([100, 200, 300, 400, 500, 600], 'Plume Temperature', 0.95, fdsVersion='6.7.1')
    fig, ax1 = plotPercentile(500, 'Plume Temperature', fdsVersion='6.7.1')