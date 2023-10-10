"""
.. module:: EAAWGmp
   :platform: Windows, Linux
   :synopsis: Outermost wrapper and main of the stochastic weather generator 
              with multiprocessing

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

python EAAWGmp.py 10 --num_real 1000

10 workers and 1000 realizations. Chunk size sets the number that go to a 
worker at a time

Main will simulate from START_REAL to START_REAL + num_real of realizations. 
The random seed is set using the realization number so that can break the 
simulation into chunks of realizations and have reproducable results.

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

import argparse
from multiprocessing import Pool
import os
import sys
import site
site.addsitedir( os.getcwd() )
import EAAWG_ProcCalib_Results as RProc

START_REAL = 1
#START_REAL = 7001
"""Starting realization - only used when want to split realizations into 
multiple runs"""
DEF_NUM_REALIZATIONS = 100
"""Default number of realizations"""
# definition of our 'default seed values'
PDEPTH_DEF_SEED = int( 21342 )
"""Precipitation depth sampler, default seed"""
STD_NORM_DEF_SEED = int( 31344 )
"""Standard normal variate sampler, default seed"""
WET_STA_DEF_SEED = int( 41446 )
"""The wet state default seed"""
DRY_STA_DEF_SEED = int( 51548 )
"""The wet state default seed"""
EVENT_RECUR_DEF_SEED = int( 62379 )
"""Event recurrence interval default seed"""
EVENT_MAG_DEF_SEED = int( 73871 )
"""Event magnitude default seed"""
CHUNK_SIZE = 5
"""Chunk size to use with the mp module
This is effectively the number of realization sent to a worker process at one 
time."""


#-----------------------------------------------------------------------
# custom functions
def delPrevOutputs( ):
    """Delete the previous model outputs before starting a new run.

    Parameters
    ----------
    None.

    Returns
    -------
    None.

    """
    # imports
    import glob
    import EAAWG_Inputs as WGI
    # globals
    # parameters
    # locals
    # start
    # get the file path for outputs
    H0OutFP = os.path.normpath( os.path.join( WGI.OUT_DIR, WGI.OUT_SUB_DIR, "*") )
    allOldFiles = glob.glob( H0OutFP )
    for oldF in allOldFiles:
        try:
            os.remove(oldF)
        except OSError as e:
            print("Error: %s : %s" % (oldF, e.strerror))
        # end try
    # end for
    # return
    return


def WG_Worker_Main( RealNum, SNSeed, PDSeed, WSLSeed, DSLSeed, EVRecurSeed,
                    EVMagSeed ):
    """ Main functionality to run a single realization

    Args:
        RealNum (int): the current realization number
        SNSeed (int): base seed for the standard normal sampler
        PDSeed (int): the precipitation depth sampler seed
        WSLSeed (int): wet state spell length sampling seed
        DSLSeed (int): dry state spell length sampling seed
        EVRecurSeed(int): event recurrence interval sampling seed
        EVMagSeed(int): event magnitude sampling seed

    Returns:
        int. The return code::
            0 -- Success!
            1 -- Failure, generic

    """
    # imports
    import pandas as pd
    import EAAWG_Inputs as WGI
    import EAAWG_Dists_Samples as WGDS
    import EAAWG_PrecipDepth as WGPD
    import EAAWG_OtherWeather as WGOW
    import EAAWG_HighRealResults as WGHRR
    import EAAWG_Events as EXEV
    # 
    # set our local seeds
    pdSampSeed = PDSeed + RealNum
    sndSampSeed = SNSeed + RealNum
    wetSSampSeed = WSLSeed + RealNum
    drySSampSeed = DSLSeed + RealNum
    evrSSampSeed = EVRecurSeed + RealNum
    evmSSampSeed = EVMagSeed + RealNum
    # get our start and end
    start_date = WGI.START_DATE
    end_date = WGI.END_DATE
    # calculate our total days and set up our time index for output
    timeD = end_date - start_date
    TOTAL_DAYS = int( timeD.days + 1 )
    DT_INDEX = pd.date_range( start=start_date, end=end_date, freq='D' )
    # with the time index set, then we have everything that we need to
    # set-up our structures
    WGOW.constructArrays()
    # set-up the sampling and other things for standard weather generator
    WGDS.setDistributions( pdSampSeed, wetSSampSeed, drySSampSeed )
    WGDS.setTrackers( start_date.month )
    WGOW.setupDistsSamples(seed_std_norm=sndSampSeed)
    WGOW.updateTracker()
    # set-up sampling and custom events
    EXEV.setEvents(evrSSampSeed, evmSSampSeed)
    EXEV.setTrackers()
    # get our starting state sampler
    StarterSamp = WGPD.PrecipSampler(pd_sample_seed=(pdSampSeed - 1))
    # get our time index in a list for an iterator
    TimesList = DT_INDEX.to_pydatetime().tolist()
    # no loop at the realization level
    curMonth = start_date.month
    # now get the starting state
    TestVal = StarterSamp.getSingleVal()
    if TestVal > 0.5:
        h0State = WGI.WET_STATE
        h0remdur = WGDS.ST_WETSPELL[curMonth]
    else:
        h0State = WGI.DRY_STATE
        h0remdur = WGDS.ST_DRYSPELL[curMonth]
    # end if
    # get the event information
    NumCustomEvents = len( WGI.EVENT_KEYS )
    EvTrigTimes = [ EXEV.NEXT_EVENT_DICT[x][0] for x in WGI.EVENT_KEYS ]
    EvTrigMags = [ EXEV.NEXT_EVENT_DICT[x][1] for x in WGI.EVENT_KEYS ]
    # now create/set our realization tracking array
    WGHRR.createSimStructures(TOTAL_DAYS)
    # set the event processing vars
    eventTriggered = False
    eventTrigMag = 0.0
    # inner loop over times
    for jJ in range(TOTAL_DAYS):
        cTime = TimesList[jJ]
        # get the current month
        curMonth = cTime.month
        # sample all every time step
        WGDS.sampleAll( curMonth )
        WGOW.updateTracker()
        # get the current day of the year
        curDayoYr = cTime.timetuple().tm_yday                
        # now that everything is sampled check our state and if wet
        # then we get a precip depth
        if h0State == WGI.WET_STATE:
            if h0remdur <= 0:
                # then need to change state to dry
                h0State = WGI.DRY_STATE
                h0remdur = WGDS.ST_DRYSPELL[curMonth]
                WGHRR.assignDryDep( jJ )
            else:
                # check if any events should be triggered
                for eE in range(NumCustomEvents):
                    if EvTrigTimes[eE] < jJ:
                        eventTriggered = True
                        eventTrigMag += EvTrigMags[eE]
                        # now need to sample our new events
                        EXEV.updateTriggeredEvent( eE, cTime, jJ )
                    # end if
                # end for event
                # assign the sampled precipitation depth
                if eventTriggered:
                    pVal = eventTrigMag
                    # now update the event that tracked
                    EvTrigTimes = [ EXEV.NEXT_EVENT_DICT[x][0] for x in WGI.EVENT_KEYS ]
                    EvTrigMags = [ EXEV.NEXT_EVENT_DICT[x][1] for x in WGI.EVENT_KEYS ]
                    eventTriggered = False
                    eventTrigMag = 0.0
                else: 
                    pVal = WGDS.ST_PDEPTH[curMonth]
                # end if event
                WGHRR.assignWetDep( jJ, pVal )
            # end if a wet day
        else:
            # then the H0 branch is currently dry
            # check if time to toggle states
            if h0remdur <= 0:
                h0State = WGI.WET_STATE
                h0remdur = WGDS.ST_WETSPELL[curMonth]
                pVal = WGDS.ST_PDEPTH[curMonth]
                WGHRR.assignWetDep( jJ, pVal )
            else:
                WGHRR.assignDryDep( jJ )
        # do the other parameters
        WGOW.calcCHI0( )
        # now update our temp values
        WGOW.calculateUpdate( jJ, curDayoYr, h0State )
        # decrement counters before moving on
        h0remdur -= 1
        # copy our array
        WGOW.rollOverChis()
    # end of time for loop
    # now output the realization
    #WGHRR.outputRealResults( RealNum, DT_INDEX)
    WGHRR.outputWSResults( RealNum, DT_INDEX, TOTAL_DAYS )
    #EXEV.outputEventTracking( RealNum, DT_INDEX )
    EXEV.outputBigEvents( RealNum, DT_INDEX )
    # end of realizations loop
    # end
    return 0

if __name__ == "__main__":
    # use the command line processor so that can tell how many processes or cores to use
    parser = argparse.ArgumentParser(description='Project description')
    parser.add_argument(
        'nbr_workers', type=int, help='Number of workers e.g. 1, 2, 4, 8')
    parser.add_argument(
        '--num_real',
        type=int,
        default=DEF_NUM_REALIZATIONS,
        help='Number of realizations for simulation')
    # parse the command line arguments received
    args = parser.parse_args()
    # extract our arguments
    num_proc = args.nbr_workers
    num_real = args.num_real
    # remove any old files
    delPrevOutputs( )
    # output
    print("Using %d processes for %d realizations" % ( num_proc, num_real))
    print("Simulate realizations %d through %d" % (START_REAL, 
                                            ( START_REAL + num_real ) - 1) )
    # now check what our number of realizations are ...
    if num_real < 2:
        # this is the run once case
        tRes = WG_Worker_Main( 1, STD_NORM_DEF_SEED, PDEPTH_DEF_SEED, 
                               WET_STA_DEF_SEED, DRY_STA_DEF_SEED, 
                               EVENT_RECUR_DEF_SEED, EVENT_MAG_DEF_SEED )
        results = [ tRes ]
    else:
        # create our list of tuples to use for the mapping
        AllArgs = [ ( int(x), STD_NORM_DEF_SEED, PDEPTH_DEF_SEED, 
                     WET_STA_DEF_SEED, DRY_STA_DEF_SEED, EVENT_RECUR_DEF_SEED, 
                     EVENT_MAG_DEF_SEED) for x in 
                    range(START_REAL, START_REAL + num_real, 1) ]
        with Pool(processes=num_proc) as pool:
            results = pool.starmap( WG_Worker_Main, AllArgs, 
                                    chunksize=CHUNK_SIZE )
        # end of with block
    # check the results
    if num_real < 5:
        if results[0] == 0:
            SuccessMessage = "Finished %d realizations successfully" % num_real
            RProc.processOutputsForPEST(num_real)
            print("%s" % SuccessMessage)
            sys.exit(0)
        else:
            ErrorMessage = "Errors or issues in realizations!!!"
            print("%s" % ErrorMessage)
            sys.exit(-1)
        # end if
    else:
        TotFails = sum( results )
        if TotFails == 0:
            SuccessMessage = "Finished %d realizations successfully" % num_real
            RProc.processOutputsForPEST(num_real)
            print("%s" % SuccessMessage)
            sys.exit(0)
        else:
            NumSuccess = num_real - TotFails
            ErrorMessage = "Finished %d successful runs out of %d realizations" % (NumSuccess, num_real)
            print("%s" % ErrorMessage)
            sys.exit(-2)
        # end if
    # end if
    # end of main


# EOF