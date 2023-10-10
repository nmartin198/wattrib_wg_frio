# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_ProcCalib_Results
   :platform: Windows, Linux
   :synopsis: Module that provides processing to produce the calibration results

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

The functions in this module are designed to be executed after all of the
realizations for a simulation are completed. They output the results used
to calibrate the weather generator formulation for a particular basin.

"""
# Copyright and License
"""
Copyright 2023 Southwest Research Institute

Module Author: Nick Martin <nick.martin@alumni.stanford.edu>

This file is part of a custom weather generator framework with extreme events, 
hereafter WG Framework.

WG Framework is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WG Framework is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with WG Framework.  If not, see <https://www.gnu.org/licenses/>.

"""

# imports
import os
import pandas as pd
import numpy as np
import EAAWG_Inputs as WGI

# parameters
PROC_SPEI_START_TS = pd.Timestamp( 2030, 1, 1, 0 )
PROC_CN_START_TS = pd.Timestamp( 2031, 1, 1, 0 )
"""Climate Normals start for processing"""
PROC_CN_END_TS = WGI.END_DATE
"""Climate Normals end for processing"""
#CN2020_Ann_Ave_P = 953.9  # Blanco
#CN2020_Ann_Ave_P = 925.2  # Cibolo
CN2020_Ann_Ave_P = 782.0  # Frio
#CN2020_Ann_Ave_P = 862.9  # Guadalupe
#CN2020_Ann_Ave_P = 931.1  # Med-Cib
#CN2020_Ann_Ave_P = 886.7  # Medina
#CN2020_Ann_Ave_P = 669.9  # Nueces
#CN2020_Ann_Ave_P = 835.6  # Sab-Med
#CN2020_Ann_Ave_P = 808.7  # Sabinal
"""Annual average precipitation depth from 1991 through 2020 for this basin"""
DROUGHT_TARGETS = { 1 : [0.149, -218.774],
                    2 : [0.149, -220.737],
                    3 : [0.149, -282.874],
                    4 : [0.149, -360.204],
                    5 : [0.149, -464.334],
                    6 : [0.149, -508.962],
                    7 : [0.149, -539.548],
                    8 : [0.149, -554.424],
                    9 : [0.149, -527.857],
                   10 : [0.149, -418.466],
                   11 : [0.149, -327.735],
                   12 : [0.149, -270.022], }
"""Cumulative deficit with corresponding cumulative probability by month for
this basin"""


#-----------------------------------------------------------------------
# functions
def probDistLLogis( paramDict, npArray ):
    """Uses generalized logistic probability distribution to estimate cumulative
    probilities for each value in the Numpy array, npArray.
    
    Args:
        paramDict (dict): dictionary with best-fit parameter values for a 
                log-logisitic distribution. Must have keys: "k", "scale",
                "loc" which are the 3 required parameters
        npArray (np.ndarray): array from time series of monthly, rolling
                average values
    
    Returns:
        retArray (np.ndarray): cumulative probabilies for each npArray value
    """
    shape = paramDict["k"]
    location = paramDict["loc"]
    scale = paramDict["scale"]
    if shape == 0.0:
        # this is the special case of a logistic distribution with 2 params
        y = ( npArray - location ) / scale
    else:
        # this is the general case of the log-logistic distribution
        takeLogArray = 1.0 - ( shape * ( npArray - location ) / scale )
        useLogArray = np.where( takeLogArray <= 0.0, 1e-7, takeLogArray )
        y = ( -1.0 * ( 1.0 / shape ) ) * np.log( useLogArray )
    # end if
    retArray = 1.0 / ( 1.0 + np.exp( -1.0 * y  ) )
    # return
    return retArray


def ppfLLogis( paramDict, npArray ):
    """Uses generalized logistic probability distribution to estimate cumulative
    probilities for each value in the Numpy array, npArray.
    
    Args:
        paramDict (dict): dictionary with best-fit parameter values for a 
                log-logisitic distribution. Must have keys: "k", "scale",
                "loc" which are the 3 required parameters
        npArray (np.ndarray): array from time series of cumulative probability
                values. All values must be greater than zero and less than or equal to one.
    
    Returns:
        retArray (np.ndarray): cumulative deficit values corresponding to cumulative
                               probabilities
    """
    shape = paramDict["k"]
    location = paramDict["loc"]
    scale = paramDict["scale"]
    if shape == 0.0:
        # this is the special case of a logistic distribution with 2 params
        retArray = location - ( scale *  np.log( ( 1.0 - npArray ) / npArray ) )
    else:
        # this is the general case of the log-logistic distribution
        retArray = location + ( scale * ( ( 1.0 - np.power( ( ( 1.0 - npArray ) / npArray ), shape ) ) / shape ) )
    # end if
    # return
    return retArray


def estimatellogparams( npArray ):
    """Estimate the parameters of a log-logistic distribution from an
    array of annual values.
    
    Estimate is done using L-moments and the "Generalized logistic distribution".
    This distribtion is a reparameterized version of the log-logistic
    distribution of Ahmad et al. (1988). Estimation is done using 
    the equations and procedure in Appendix A.7 of "Regional Frequency
    Analysis", Hosking and Wallis (1997)
    
    To estimate the distribution parameters (shape, scale, and location),
    the L-moments l1, l2, and t3 need to be calculated. These three
    L-moments can be estimated from the first three, sample weighted
    probability moments (b0, b1, and b2).
    
    Args:
        npArray (np.ndarray): Numpy, 1D array
    
    Returns:
        log-logistic parameters in dictionary, D:
            D["k"]: k or shape
            D["scale"]: alpha or scale
            D["loc"]: Eta or location
    """
    # imports
    import math
    # don't do any checking for type and assume that will always
    #  be Numpy ndarray for single argument
    totLen = len( npArray )
    # need a sorted array in increasing order
    srtAr = np.sort( npArray )
    # calculate sample probability weighted moments: b0, b1, b2
    b0 = srtAr.mean()
    b1 = 0.0
    for iI in range(2, totLen + 1):
        b1 += ( ( iI - 1 ) / ( totLen - 1 ) ) * srtAr[iI-1]
    # end for
    b1 = b1 / totLen
    b2 = 0.0
    for iI in range( 3, totLen + 1 ):
        b2 += ( ( ( iI - 1 ) * ( iI - 2 ) ) / ( ( totLen - 1 ) * (totLen - 2 ) ) ) * srtAr[iI-1]
    # end for
    b2 = b2 / totLen
    # calculate sample L-moments: l1, l2, t3
    l1 = b0
    l2 = (2.0 * b1 ) - b0
    l3 = ( 6.0 * b2 ) -  ( 6.0 * b1 ) + b0
    t3 = l3 / l2
    # estimate the distribution parameters
    shape = -1.0 * t3
    scale = ( l2 * math.sin( shape * math.pi ) ) / ( shape * math.pi )
    location = l1 - ( scale * ( ( 1.0 / shape ) - ( math.pi / math.sin( shape * math.pi ) ) ) )
    retDict = { "k" : shape,
                "scale" : scale,
                "loc" : location, }
    # return
    return retDict


def processOutputsForPEST( NumReal ):
    """Collate and process outputs for PEST analyses
    
    There are two main results employed for calibration:
        1. average annual precipitation depth
        2. cumulative deficit corresponding to the 5X cumulative probability
    

    Parameters
    ----------
    NumReal : int
        Number of realizations in this simulation.

    Returns
    -------
    None.

    """
    # imports
    from copy import deepcopy
    # globals
    global PROC_SPEI_START_TS, PROC_CN_START_TS, PROC_CN_END_TS
    global DROUGHT_TARGETS
    # parameters
    # locals
    AnnAvePreList = list()
    AnnAvePETList = list()
    AnnAveDefList = list()
    RealsList = [ x for x in range( 1, NumReal+1 ) ]
    MonthsList = [ x for x in range(1, 13, 1) ]
    MonthTrackDict = dict()
    TargCDefDict = dict()
    # start
    copyList = list()
    for mon in MonthsList:
        MonthTrackDict[mon] = deepcopy( copyList )
    # end for
    # need to loop through all of the realizations to collate
    for RealNum in RealsList:
        # get the name and location for serialized DataFrame for this real
        H0FileName = "%s_R%d_DF.pickle" % (WGI.OUT_LABEL, RealNum)
        H0OutFP = os.path.normpath( os.path.join( WGI.OUT_DIR, WGI.OUT_SUB_DIR, 
                                                  H0FileName ) )
        # read in the DataFrame
        H0DF = pd.read_pickle( H0OutFP, compression='zip' )
        MonH0DF = H0DF[["Def_mm"]].loc[PROC_SPEI_START_TS:PROC_CN_END_TS].resample( 'MS', ).sum()
        H0DF = H0DF.loc[PROC_CN_START_TS:PROC_CN_END_TS].copy()
        AnnH0DF = H0DF[["Precip_mm", "ETo_mm", "Def_mm"]].resample( 'AS' ).sum()
        # process annual averages 
        AnnAvePreList.append( AnnH0DF["Precip_mm"].mean() )
        AnnAvePETList.append( AnnH0DF["ETo_mm"].mean() )
        AnnAveDefList.append( AnnH0DF["Def_mm"].mean() )
        # process to get collections of months to make SPEI distributions
        D3DF = MonH0DF.rolling(window=3,).sum()
        D3DF = D3DF.loc[PROC_CN_START_TS:].copy()
        D3DF["Month"] = D3DF.index.month
        for mon in MonthsList:
            m3Mon = D3DF[D3DF["Month"] == mon].copy()
            a3Mon = m3Mon["Def_mm"].to_numpy(dtype=np.float32)
            MonthTrackDict[mon].extend( a3Mon.tolist() )
        # end month for
    # end realization for
    # make the annual averages DataFrame
    DataDict = { "Precip_mm" : np.array( AnnAvePreList, dtype=np.float32 ),
                 "ETo_mm" : np.array( AnnAvePETList, dtype=np.float32 ),
                 "Def_mm" : np.array( AnnAveDefList, dtype=np.float32 ), }
    CoAnnAvesDF = pd.DataFrame( index=RealsList, data=DataDict )
    # process our collections of monthly, 3-month cumulative deficit values.
    for mon in MonthsList:
        a3Mon = np.array( MonthTrackDict[mon], dtype=np.float32 )
        # fit 'generalized logistic' distributions to this array
        lD3Mon = estimatellogparams( a3Mon )
        estCProbTarget = probDistLLogis( lD3Mon, np.array( [DROUGHT_TARGETS[mon][1]], 
                                                           dtype=np.float32 ) )
        #if np.isnan( estCProbTarget[0] ):
        #    meanMoDef = a3Mon.mean()
        #    meanDef = CoAnnAvesDF["Def_mm"].mean()
        #    if ( meanDef >= 0.0 ) or ( meanMoDef >= 0.0 ):
        #        cCDFEst = 1.0
        #    else:
        #        cCDFEst = 0.0
        #    # end if
        #else:
        #    cCDFEst = float( estCProbTarget[0] )
        ## end if
        cCDFEst = float( estCProbTarget[0] )
        TargCDefDict[mon] = cCDFEst
    # end for
    # now output
    OutFiler = os.path.normpath( os.path.join( WGI.CUR_DIR, "CollOuts.dat" ) )
    with open( OutFiler, 'w' ) as OF:
        OF.write( "Collated weather generator output for %d realizations \n\n" 
                  % NumReal)
        OF.write( "Annual average precipitation depth in mm: %6.1f \n" % 
                  CoAnnAvesDF["Precip_mm"].mean() )
        OF.write( "Annual average ETo depth in mm: %6.1f \n" % 
                  CoAnnAvesDF["ETo_mm"].mean() )
        OF.write( "Annual average deficit depth in mm: %6.1f \n\n" % 
                  CoAnnAvesDF["Def_mm"].mean() )
        for mon in MonthsList:
            OF.write( "Month %d - 3-month cumulative deficit (mm) %6.1f; "
                      "target cumulative probability %5.3f; "
                      "actual cumulative probability %5.3f \n" %
                      ( mon, DROUGHT_TARGETS[mon][1], DROUGHT_TARGETS[mon][0],
                        TargCDefDict[mon] ) )
        # end for
    # end with and close file
    # return
    return


#EOF