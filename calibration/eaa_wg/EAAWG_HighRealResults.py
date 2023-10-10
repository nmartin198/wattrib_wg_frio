# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_HighRealResults
   :platform: Windows, Linux
   :synopsis: Provides logic for handling one realization per process worker

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

Stores results for each realization. For precipitation, store daily precip 
depth. For other parameters, store a single value per day for the entire
domain. Functions are also provided for calculating ETo and Deficit.
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
import pandas as pd
import numpy as np
import EAAWG_Inputs as WGI


#------------------------------------------------------------------------
TMAX_IND = 0
TMIN_IND = 1
PRE_IND = 2
TOT_VALS = 3

# Program, global data structures
H0_REAL = None
"""Single realization output for H0 case"""

#--------------------------------------------------------------------------
# python functions
def createSimStructures(TotDays):
    """Convenience method to create and set the simulation structures
    """
    # imports
    # globals
    global H0_REAL, TOT_VALS
    # start of function
    H0_REAL = np.zeros( (TotDays, TOT_VALS ), dtype=np.float32 )
    # end


def assignDryDep( tIndex ):
    """Assign dry precipitation depth (i.e. 0) to the correct time interval
    
    Args:
        tIndex (int): current time interval
    
    """
    # imports
    # start
    global H0_REAL, PRE_IND
    # start function
    H0_REAL[tIndex, PRE_IND] = 0.0
    # end of function

def assignWetDep( tIndex, pVal):
    """Assign precipitation depth (i.e. 0) to the H0 data structure.
    
    Args:
        tIndex (int): current time interval
        pVal (float): current precipitation depth for the day
    
    """
    # globals
    global H0_REAL, PRE_IND
    # start
    H0_REAL[tIndex, PRE_IND] = pVal
    # end of function

def assignTemp( tIndex, MaxT, MinT ):
    """Convenience function to assign simulated temperature values.
    
    Args:
        tIndex (int): the time index
        MaxT (float): max temp
        MinT (float): min temp
    
    """
    # globals
    global H0_REAL, TMAX_IND, TMIN_IND
    # start
    H0_REAL[tIndex, TMAX_IND] = MaxT
    H0_REAL[tIndex, TMIN_IND] = MinT
    # end

def outputRealResults(RealNum, DT_INDEX):
    """Output the results for the current realization. Use Pandas DataFrames
    and pickles.
    
    Args:
        RealNum (int): current realization number.
        DT_INDEX (pd.DateTimeIndex): index for all outputs
    """
    # imports
    import pandas as pd
    from os import path
    # globals
    global H0_REAL, PRE_IND, TMAX_IND, TMIN_IND, TOT_VALS
    # start
    # file names
    H0FileName = "%s_R%d_DF.pickle" % (WGI.OUT_LABEL, RealNum)
    H0OutFP = path.normpath( path.join( WGI.OUT_DIR, WGI.OUT_SUB_DIR, 
                                        H0FileName ) )
    # make our DataFrames
    H0DDict = { "Tmax_C" : H0_REAL[:,TMAX_IND],
                "Tmin_C" : H0_REAL[:,TMIN_IND],
                "Precip_mm" : H0_REAL[:,PRE_IND], }
    # end of for
    H0DF = pd.DataFrame( index=DT_INDEX, data=H0DDict )
    H0DF.to_pickle( H0OutFP, compression='zip' )
    # end
    return

def outputWSResults(RealNum, DT_INDEX, TotDays):
    """Output the watershed results for the current realizationn. Use Pandas DataFrames
    and pickles. Uses area average of precipitation grid cells to calc the WS precip.
    
    Args:
        RealNum (int): current realization number.
        DT_INDEX (pd.DateTimeIndex): index for all outputs
        TotDays (int): total number of days in realization
    """
    # imports
    from os import path
    # globals
    global H0_REAL, H1_REAL, PRE_START_IND, TMAX_IND, TMIN_IND
    # start
    # file names
    #H0FileName = "WS_%s_R%d_DF.pickle" % (WGI.OUT_LABEL, RealNum)
    H0FileName = "%s_R%d_DF.pickle" % (WGI.OUT_LABEL, RealNum)
    H0OutFP = path.normpath( path.join( WGI.OUT_DIR, WGI.OUT_SUB_DIR, 
                                        H0FileName ) )
    # make our DataFrame
    H0DDict = { "Tmax_C" : H0_REAL[:,TMAX_IND],
                "Tmin_C" : H0_REAL[:,TMIN_IND],
                "Precip_mm" : H0_REAL[:,PRE_IND], }
    TAve = 0.5 * ( H0_REAL[:,TMAX_IND] + H0_REAL[:,TMIN_IND] )
    H0DDict["Tave_C"] = TAve
    H0DF = pd.DataFrame( index=DT_INDEX, data=H0DDict )
    H0DF["ETo_mm"] = calcPET_HS( DT_INDEX, H0DF )
    H0DF["Def_mm"] = H0DF["Precip_mm"] - H0DF["ETo_mm"]
    # write out all of our waterbalance related DataFrames
    H0DF.to_pickle( H0OutFP, compression='zip' )
    # end
    return

def calcPET_HS( DT_INDEX, H0DF ):
    """Calculate PET in mm using Hargreaves-Samani

    Args:
        DT_INDEX (pd.DateTimeIndex): index for all outputs
        H0DF (pd.DataFrame): dataframe with all the values
    
    Returns:
        ETo_mmd (np.array): PET depths per day in mm
    """
    # imports
    import math
    from EAAWG_Inputs import LAT_DEG
    # start of function
    # get the TAve array
    TAve = H0DF["Tave_C"].to_numpy(dtype=np.float32)
    MaxTDaily = H0DF[["Tmax_C"]].copy()
    MinTDaily = H0DF[["Tmin_C"]].copy()
    MonMaxTDaily = MaxTDaily.resample('MS').mean()
    MonMinTDaily = MinTDaily.resample('MS').mean()
    H0DF["MonDelta_T"] = 0.0
    for indx, row in H0DF.iterrows():
        curMonTS = pd.Timestamp( indx.year, indx.month, 1, 0, )
        H0DF.at[indx,"MonDelta_T"] = ( MonMaxTDaily.at[curMonTS,"Tmax_C"] - 
                                       MonMinTDaily.at[curMonTS,"Tmin_C"] )
    # end for
    Delta_T = H0DF["MonDelta_T"].to_numpy(dtype=np.float32)
    Delta_T = np.where( Delta_T < 1.0, 1.0, Delta_T )
    useDelta_T = np.power( Delta_T, 0.5 )
    # solar rad calcs
    DayOYr = DT_INDEX.dayofyear.to_numpy()
    SDec_rad = 0.4093 * np.sin( ( ( ( 2.0 * math.pi ) / 365.0 ) * DayOYr ) - 1.405 )
    SunS_rad = np.arccos( -1.0 * math.tan(math.radians(LAT_DEG)) * np.tan(SDec_rad) )
    RelDEtoS = 1.0 + 0.033 * np.cos( ( ( 2.0 * math.pi ) / 365.0 ) * DayOYr )
    #MaxDayHrs = (24.0/math.pi) * SunS_rad
    S_o_mmd = 15.392 * RelDEtoS * ( ( SunS_rad * math.sin( math.radians(LAT_DEG) ) * 
                                      np.sin( SDec_rad ) ) + 
                                   ( math.cos( math.radians(LAT_DEG) ) * 
                                     np.cos( SDec_rad ) * np.sin( SunS_rad ) ) ) 
    # now for the PET calc
    ETo_mmd = 0.0023 * S_o_mmd * useDelta_T * ( TAve + 17.8 )
    # return
    return ETo_mmd

def adjustETo( Ks, ETo, Precip, RRDay=1.0 ):
    """ Function to adjust ETo using a crop coefficient and if there was precipitation.
    Designed to be used in a apply( lambda row: )
    
    Args:
        Ks (float): crop coefficient
        ETo (float): current day's ETo in depth units
        Precip (float): current day's precip in depth units
        RRDay (float): reduction factor for ET on rainy days
    
    Returns:
        float: calculated potential evapotranspiration (PET)
        
    """
    PET = ETo * Ks
    if Precip > 0.0:
        PET = PET * RRDay
    # return
    return PET

def cleanAllEnd():
    """Convenience method to clean up at the end"""
    global H0_REAL
    # start
    H0_REAL = None

# EOF