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
import matplotlib.pyplot as plt
import scipy.stats as scst
import numpy as np
import pandas as pd
import glob
import os

from .colorSchemes import getVTcolors

def readErrorTable(fdsVersion='6.7.1'):
    """Reads the model error table for an FDS version

    The tables are the model bias and relative standard deviation
    published for each validated quantity in the FDS Validation Guide.
    They are shipped with pyfdstools in the fdsErrorTables directory.

    Parameters
    ----------
    fdsVersion : str, optional
        FDS version whose error table is read (default '6.7.1'). The
        tables available are named after the versions they came from

    Returns
    -------
    pandas.DataFrame
        Table indexed by quantity, with the columns published in the
        validation guide including 'Bias' and 'sigmaM'
    list
        List of the quantity names in the table

    Raises
    ------
    FileNotFoundError
        If no error table is shipped for the requested FDS version
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))
    file = os.path.abspath("%s%sfdsErrorTables%s%s.csv"%(dir_path, os.sep, os.sep, fdsVersion))
    if not os.path.exists(file):
        available = sorted(
            os.path.basename(x).replace('.csv', '')
            for x in glob.glob(os.path.join(
                dir_path, 'fdsErrorTables', '*.csv')))
        raise FileNotFoundError(
            "No FDS error table for version %s. Available versions: %s"
            % (fdsVersion, ', '.join(available)))
    data = pd.read_csv(file, index_col=0)
    keys = list(data.index.values)
    return data, keys

def calculatePercentile(values, quantity, percentile, fdsVersion='6.7.1'):
    """Applies the FDS model uncertainty to predicted values

    Each predicted value is corrected for the model bias and treated as
    the mean of a normal distribution whose standard deviation is the
    published relative standard deviation. The value at the requested
    percentile of that distribution is returned, which is the basis of
    the probabilistic approach described in the FDS Validation Guide.

    Parameters
    ----------
    values : array-like
        Predicted values from the model
    quantity : str
        Quantity name as it appears in the error table
    percentile : float
        Percentile to evaluate, in the range 0 to 1
    fdsVersion : str, optional
        FDS version whose error table is used (default '6.7.1')

    Returns
    -------
    array
        Values at the requested percentile, one per input value. A zero
        prediction has zero spread under this multiplicative model and
        is returned bias-corrected as it is

    Raises
    ------
    ValueError
        If the quantity is not in the error table for this version

    Notes
    -----
    Releases before v0.0.24 evaluated the distribution on a 101 point
    grid rather than inverting it, so results were about 0.4% high, and
    raised IndexError for a zero value or for a percentile beyond about
    0.9987.
    """

    data, quantities = readErrorTable(fdsVersion=fdsVersion)

    if (quantity in quantities):
        fdsBias = data.loc[quantity]['Bias']
        fdsSigma = data.loc[quantity]['sigmaM']

        values = np.asarray(values, dtype=float)
        mu = values/fdsBias
        sigma = np.abs(mu)*fdsSigma

        # scipy's percent point function inverts the normal CDF
        # directly. Earlier releases evaluated the CDF on 101 points
        # spanning mu +/- 3 sigma and took the first above the
        # percentile, which was 0.4% high, raised IndexError for any
        # percentile beyond about 0.9987 (where the answer lies outside
        # three standard deviations), and raised IndexError again for a
        # zero value, where sigma is zero and the CDF is undefined.
        errorValues = np.array(mu, dtype=float)
        spread = sigma > 0
        errorValues[spread] = scst.norm.ppf(
            percentile, mu[spread], sigma[spread])
        # A zero prediction has zero spread under this multiplicative
        # uncertainty model, so it is returned bias-corrected as it is.
        return errorValues

    raise ValueError(
        "Quantity '%s' is not in the FDS %s error table. Known "
        "quantities: %s" % (quantity, fdsVersion, ', '.join(quantities)))


def plotPercentile(value, quantity, fdsVersion='6.7.1', colors=None):
    """Plots the model uncertainty distribution for a predicted value

    Parameters
    ----------
    value : float
        Predicted value from the model
    quantity : str
        Quantity name as it appears in the error table
    fdsVersion : str, optional
        FDS version whose error table is used (default '6.7.1')
    colors : array-like, optional
        Colors for the PDF and CDF curves. A default sequence is used
        when omitted

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    matplotlib.axes.Axes
        Axes carrying the probability density curve. The cumulative
        curve is drawn on a twin axis

    Raises
    ------
    ValueError
        If the quantity is not in the error table for this version
    """

    data, quantities = readErrorTable(fdsVersion=fdsVersion)
    if colors is None:
        colors = getVTcolors()
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
        ax1.plot(x, y_pdf, label='PDF', color=colors[0], linewidth=lw,)
        ax1.set_xlabel('%s'%(quantity),fontsize=fs)
        ax1.set_ylabel('PDF Probabilty Density',fontsize=fs, color=colors[0])
        ax1.tick_params('y', colors=colors[0], labelsize=fs)
        ax1.tick_params('x', labelsize=fs)
        
        ax2 = ax1.twinx()
        ax2.plot(x, y_cdf, label='CDF', color=colors[1], linewidth=lw,)
        ax2.tick_params('y', colors=colors[1], labelsize=fs)
        ax2.set_ylabel('CDF Probabilty',fontsize=fs, color=colors[1], rotation=270, labelpad=20)
        
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        
        ax2.legend(lines+lines2, labels+labels2, fontsize=fs)
        plt.tight_layout()
        return fig, ax1

    raise ValueError(
        "Quantity '%s' is not in the FDS %s error table. Known "
        "quantities: %s" % (quantity, fdsVersion, ', '.join(quantities)))


def getQuantities(fdsVersion='6.7.1'):
    """Lists the quantities in an FDS error table

    Parameters
    ----------
    fdsVersion : str, optional
        FDS version whose error table is read (default '6.7.1')

    Returns
    -------
    list
        List of the quantity names in the table

    Notes
    -----
    Releases before v0.0.24 read the table through a path relative to
    the working directory, so this only worked when called from inside
    the installed package directory, and returned the table's column
    names rather than its quantities.
    """

    data, keys = readErrorTable(fdsVersion=fdsVersion)
    return keys
    
