# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_OtherWeather
   :platform: Windows, Linux
   :synopsis: Provides the logic for the other, non-precipitation, weather 
              parameters

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

Handles the calculation of daily values using helper functions from other 
modules. Also sends the daily values elsewhere for archive.

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
import numpy as np
import EAAWG_Inputs as WGI

# numpy set err
#np.seterr(all='raise')

# module level parameters
NUM_OTHER = 2
"""Number of other weather parameters in this case we are only dealing with
Min and max temperature"""
NUM_DAYS_YR = 366
"""Number of days in the year, maximum number"""
SIGMA_THRESH = 4.0
"""Sigma multiplier threshold for standard deviations"""
MIN_DAILY_DELTA = 4.0
"""The minimum allowable daily difference between maximum and minimum temps"""

# module level variables
A_DATA = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
"""A array for calculating or projecting the daily residual or error term."""
B_DATA = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
"""B array for calculating or projecting the daily residual or error term."""
M0 = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
"""M0 array for calculating A and B matrices. Not currently used"""
M1 = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
"""M1 array for calculating A and B matrices. Not currently used"""
WET_TMAX_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Average wet state daily Tmax, Fourier smoothed"""
WET_TMAX_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Daily wet state standard deviation of Tmax, Fourier smoothed"""
WET_TMIN_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Average daily wet state Tmin, Fourier smoothed"""
WET_TMIN_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Daily wet state standard deviation of Tmin, Fourier smoothed"""
DRY_TMAX_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Average dry state daily Tmax, Fourier smoothed"""
DRY_TMAX_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Daily dry state standard deviation of Tmax, Fourier smoothed."""
DRY_TMIN_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Average daily dry state Tmin, Fourier smoothed"""
DRY_TMIN_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
"""Daily dry state standard deviation of Tmin, Fourier smoothed."""
EPS_STD_NORMAL = list()
"""The list of standard normal error variates"""
EPS_NORM_SAMP = list()
"""The list of samplers for the standard normal error variates."""
EPSI_2 = np.ones( (1, NUM_OTHER), dtype=np.float64 )
"""Tracker for actual epsilon values"""
CHI_M0 = np.ones( (1, NUM_OTHER), dtype=np.float64 )
"""Current day Chi matrix"""
CHI_L1M0 = np.ones( (1, NUM_OTHER), dtype=np.float64 )
"""Previous day Chi matrix"""

#--------------------------------------------------------------------------
# functions
def populateMmats():
    """Convenience function to read in and populate the M0 and M1 matrices
    
    """
    # globals
    global M0, M1
    # start of function
    # do the M0
    # end of with block
    M0[0,0] = WGI.OW_M0_IN.at['rho_1X', 'rho_X1']
    M0[1,0] = WGI.OW_M0_IN.at['rho_2X', 'rho_X1']
    M0[0,1] = WGI.OW_M0_IN.at['rho_1X', 'rho_X2']
    M0[1,1] = WGI.OW_M0_IN.at['rho_2X', 'rho_X2']
    # end of with block
    M1[0,0] = WGI.OW_M1_IN.at['rho_1X', 'rho_X1_L1']
    M1[1,0] = WGI.OW_M1_IN.at['rho_2X', 'rho_X1_L1']
    M1[0,1] = WGI.OW_M1_IN.at['rho_1X', 'rho_X2_L1']
    M1[1,1] = WGI.OW_M1_IN.at['rho_2X', 'rho_X2_L1']
    # end of function
    return

def constructArrays():
    """Convenience method to construct calculation arrays. 

    """
    # imports
    # globals
    global WET_TMAX_AVE, WET_TMAX_STD, WET_TMIN_AVE, WET_TMIN_STD
    global DRY_TMAX_AVE, DRY_TMAX_STD, DRY_TMIN_AVE, DRY_TMIN_STD 
    global A_DATA, B_DATA
    # Now read in our values/files and construct
    WAveDF = WGI.OW_WET_AVE
    WET_TMAX_AVE[:] = np.array( WAveDF['Tmax_C'], dtype=np.float64 )
    WET_TMAX_AVE[:] = WET_TMAX_AVE[:] + WGI.AVE_WET_TMAX_ADD
    WET_TMIN_AVE[:] = np.array( WAveDF['Tmin_C'], dtype=np.float64 )
    WET_TMIN_AVE[:] = WET_TMIN_AVE[:] + WGI.AVE_WET_TMIN_ADD
    WStdDF = WGI.OW_WET_STD
    WET_TMAX_STD[:] = np.array( WStdDF['Tmax_C'], dtype=np.float64 )
    WET_TMIN_STD[:] = np.array( WStdDF['Tmin_C'], dtype=np.float64 )
    DAveDF = WGI.OW_DRY_AVE
    DRY_TMAX_AVE[:] = np.array( DAveDF['Tmax_C'], dtype=np.float64 )
    DRY_TMAX_AVE[:] = DRY_TMAX_AVE[:] + WGI.AVE_DRY_TMAX_ADD
    DRY_TMIN_AVE[:] = np.array( DAveDF['Tmin_C'], dtype=np.float64 )
    DRY_TMIN_AVE[:] = DRY_TMIN_AVE[:] + WGI.AVE_DRY_TMIN_ADD
    DStdDF = WGI.OW_DRY_STD
    DRY_TMAX_STD[:] = np.array( DStdDF['Tmax_C'], dtype=np.float64 )
    DRY_TMIN_STD[:] = np.array( DStdDF['Tmin_C'], dtype=np.float64 )
    # next do the A and B data
    A_DATA[:, :] = np.array( WGI.A_DATA_LIST, dtype=np.float64 )
    B_DATA[:, :] = np.array( WGI.B_DATA_LIST, dtype=np.float64 )
    # end of function
    return

def setupDistsSamples(seed_std_norm=None):
    """Setup the distributions and samplers for the white noise term.

    KWargs:
        seed_std_norm (int): standard normal sampler seed

    """
    # imports
    import EAAWG_StdNormal as WGSN
    # globals
    global EPS_STD_NORMAL, EPS_NORM_SAMP
    # start of function
    EPS_STD_NORMAL.append( WGSN.StdNormal() )
    EPS_NORM_SAMP.append( WGSN.ErrorTSampler(seed=seed_std_norm) )
    EPS_STD_NORMAL.append( WGSN.StdNormal() )
    EPS_NORM_SAMP.append( WGSN.ErrorTSampler(seed=seed_std_norm) )
    # end of for
    # end of function
    return

def updateTracker():
    """Update the epsilon value tracker array"""
    # globals
    global EPSI_2, EPS_STD_NORMAL, EPS_NORM_SAMP
    # start of function
    EPSI_2[0, 0] = EPS_STD_NORMAL[0].ranval1( EPS_NORM_SAMP[0].ranstate )
    EPSI_2[0, 1] = EPS_STD_NORMAL[1].ranval1( EPS_NORM_SAMP[1].ranstate )
    # now make some checks
    EPSI_2 = np.where( np.isnan( EPSI_2 ), 0.25, EPSI_2 )
    EPSI_2 = np.where( np.isinf( EPSI_2 ), 0.25, EPSI_2 )
    # end of function
    return

def calculateUpdate( curIndex, DayOYear, state ):
    """Calculate the current day and update
    
    Args:
        curIndex (int): current date index
        DayOYear (int): day of the year
        state (string): wet or dry state
    
    """
    # imports
    import EAAWG_HighRealResults as WGHRR
    # globals
    global WET_TMAX_AVE, WET_TMAX_STD, WET_TMIN_AVE, WET_TMIN_STD 
    global DRY_TMAX_AVE, DRY_TMAX_STD, DRY_TMIN_AVE, DRY_TMIN_STD
    global CHI_M0, MIN_DAILY_DELTA
    # start of function
    np.seterr(all='raise')
    if state == WGI.WET_STATE:
        try:
            cMaxT0 = ( (CHI_M0[0,0] * WET_TMAX_STD[(DayOYear - 1)]) + 
                        WET_TMAX_AVE[(DayOYear - 1)] )
        except:
            cMaxT0 = WET_TMAX_AVE[(DayOYear - 1)]
        try:
            cMinT0 = ( (CHI_M0[0,1] * WET_TMIN_STD[(DayOYear - 1)]) + 
                        WET_TMIN_AVE[(DayOYear - 1)] )
        except:
            cMinT0 = WET_TMIN_AVE[(DayOYear - 1)]
    else:
        try:
            cMaxT0 = ( (CHI_M0[0,0] * DRY_TMAX_STD[(DayOYear - 1)]) + 
                       DRY_TMAX_AVE[(DayOYear - 1)] )
        except:
            cMaxT0 = DRY_TMAX_AVE[(DayOYear - 1)]
        try:
            cMinT0 = ( (CHI_M0[0,1] * DRY_TMIN_STD[(DayOYear - 1)]) + 
                        DRY_TMIN_AVE[(DayOYear - 1)] )
        except:
            cMinT0 = DRY_TMIN_AVE[(DayOYear - 1)]
    # toggle back
    np.seterr(all='warn')
    acMaxT0 = np.array( cMaxT0, dtype=np.float64 )
    acMaxT0 = np.where( np.isnan( acMaxT0 ), 
                        DRY_TMAX_AVE[(DayOYear - 1)], acMaxT0 )
    acMaxT0 = np.where( np.isinf( acMaxT0 ), 
                        DRY_TMAX_AVE[(DayOYear - 1)], acMaxT0 )
    acMinT0 = np.array( cMinT0, dtype=np.float64 )
    acMinT0 = np.where( np.isnan( acMinT0 ), 
                        DRY_TMIN_AVE[(DayOYear - 1)], acMinT0 )
    acMinT0 = np.where( np.isinf( acMinT0 ), 
                        DRY_TMIN_AVE[(DayOYear - 1)], acMinT0 )
    passMaxT = float( acMaxT0 )
    passMinT = float( acMinT0 )
    if ( passMaxT - passMinT ) < MIN_DAILY_DELTA:
        passMaxT = passMinT + MIN_DAILY_DELTA
    # end if
    # now are ready to update
    WGHRR.assignTemp( curIndex, passMaxT, passMinT )
    # end of function
    return

def rollOverChis():
    """Convenience function to roll over the Chi squared values
    """
    # local imports
    from copy import deepcopy
    # globals
    global CHI_M0, CHI_L1M0, NUM_OTHER
    # start of function
    for iI in range(NUM_OTHER):
        CHI_L1M0[0, iI] = deepcopy( CHI_M0[0, iI] )
    # end for
    # end of function
    return

def calcCHI0( ):
    """Calculate the current Chi no lag
    
    """
    # globals
    global CHI_M0, CHI_L1M0, EPSI_2, A_DATA, B_DATA, SIGMA_THRESH
    # calculations
    # first set our A and B
    A_0 = A_DATA.copy()
    B_0 = B_DATA.copy()
    np.seterr(all='raise')
    try:
        CHI_M0[:,:] = (np.matmul(CHI_L1M0, A_0) + np.matmul( EPSI_2, B_0 ))[:,:]
    except:
        CHI_M0[:,:] = 1.0
    # do some checks
    np.seterr(all='warn')
    # first check to make sure that we actually got the calculation done
    CHI_M0 = np.where( np.isnan( CHI_M0 ), 1.0, CHI_M0 )
    CHI_M0 = np.where( np.isinf( CHI_M0 ), 1.0, CHI_M0 )
    # next we need to provide a sigma threshold
    CHI_M0 = np.where( CHI_M0 > SIGMA_THRESH, SIGMA_THRESH, CHI_M0 )
    CHI_M0 = np.where( CHI_M0 < (-1.0*SIGMA_THRESH), (-1.0*SIGMA_THRESH), 
                       CHI_M0 )
    # end of function
    return

def cleanAllEnd():
    """Convenience method to clean or delete all trackers at the end """
    global A_DATA, B_DATA, M0, M1, WET_TMAX_AVE, WET_TMAX_STD, WET_TMIN_AVE
    global WET_TMIN_STD, DRY_TMAX_AVE, DRY_TMAX_STD, DRY_TMIN_AVE 
    global DRY_TMIN_STD, EPS_STD_NORMAL, EPS_NORM_SAMP, EPSI_2
    global CHI_M0, CHI_L1M0
    # set to none
    A_DATA = None
    B_DATA = None
    M0 = None
    M1 = None
    WET_TMAX_AVE = None
    WET_TMAX_STD = None
    WET_TMIN_AVE = None
    WET_TMIN_STD = None
    DRY_TMAX_AVE = None
    DRY_TMAX_STD = None
    DRY_TMIN_AVE = None
    DRY_TMIN_STD = None
    EPS_STD_NORMAL = None
    EPS_NORM_SAMP = None
    EPSI_2 = None
    CHI_M0 = None
    CHI_L1M0 = None
    # end
    return

def setAllBegin():
    """Convenience method to reset all at the beginning of a realization """
    global NUM_OTHER, NUM_DAYS_YR, SIGMA_THRESH, A_DATA, B_DATA, M0, M1
    global WET_TMAX_AVE, WET_TMAX_STD, WET_TMIN_AVE, WET_TMIN_STD 
    global DRY_TMAX_AVE, DRY_TMAX_STD, DRY_TMIN_AVE, DRY_TMIN_STD
    global EPS_STD_NORMAL, EPS_NORM_SAMP, EPSI_2, CHI_M0, CHI_L1M0
    # now do the setting
    A_DATA = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
    B_DATA = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
    M0 = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
    M1 = np.ones( (NUM_OTHER, NUM_OTHER), dtype=np.float64 )
    WET_TMAX_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    WET_TMAX_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    WET_TMIN_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    WET_TMIN_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    DRY_TMAX_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    DRY_TMAX_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    DRY_TMIN_AVE = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    DRY_TMIN_STD = np.ones( NUM_DAYS_YR, dtype=np.float64 )
    EPS_STD_NORMAL = list()
    EPS_NORM_SAMP = list()
    EPSI_2 = np.ones( (1, NUM_OTHER), dtype=np.float64 )
    CHI_M0 = np.ones( (1, NUM_OTHER), dtype=np.float64 )
    CHI_L1M0 = np.ones( (1, NUM_OTHER), dtype=np.float64 )
    # end
    return

# EOF