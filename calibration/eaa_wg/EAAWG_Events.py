# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_Events
   :platform: Windows, Linux
   :synopsis: Contains a daily "event" class. Events have an arrival time 
              (or really interarrival time) and an occurrence magnitude. 
              Consequently, they are composed Poisson distribution for 
              interarrival times and a uniform distribution for event 
              magnitude.

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>


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
from scipy import stats as scstats

# global parameters
ACT_EVENT_DICT = None
"""Dictionary for actual or instantiated event objects. Each event has two
random samplers, one for magnitude and one for recurrence, and has one 
ExtremeEvent object."""
TTRIG_EVENT_TRACK_DICT = None
"""Dictionary for tracking the events that are or could be triggered."""
NEXT_EVENT_DICT = None
"""Dictionary that stores the next event time and magnitude."""


# class definitions
class CreateEventError(Exception):
    """Custom exception error for when something happens with an event
    """
    def __init__(self, arg):
        self.args = arg


class ExtremeEvent(object):
    """Extreme event object
    An event is composed of both an arrival time (handled stochastically with
    a Poisson distribution) and an event magnitude (handled stochastically 
    with a Uniform distribution)."""
    
    def __init__(self, aveRecurYrs, lowMag, highMag, name="Custom Extreme Event" ):
        """Override of  initialization, creates a custom extreme event object.

        Parameters
        ----------
        aveRecurYrs : int
            Expected recurrence interval in years. Typically, 10-, 25-, 50-, or
            100-years.
        lowMag : float
            Minimum event magnitude for daily precipitation depth in mm.
        highMag : float
            Maximum event magnitude for daily precipitation depth in mm.
        name : str, optional
            Descriptive name for object. The default is "Custom Extreme Event".

        Raises
        ------
        CreateEventError
            Custom error handling.

        Returns
        -------
        None.

        """
        super().__init__()
        # do some checks and assign our properties
        if not isinstance( lowMag, float ):
            ErrorMsg = "lowMag must be a float!!!"
            raise CreateEventError( ErrorMsg )
        if not isinstance( lowMag, float ):
            ErrorMsg = "highMag must be a float!!!"
            raise CreateEventError( ErrorMsg )
        if not isinstance( aveRecurYrs, int ):
            ErrorMsg = "aveRecurYrs must be an int!!!"
            raise CreateEventError( ErrorMsg )
        if lowMag >= highMag:
            ErrorMsg = "highMag must be greater than lowMag !!!"
            raise CreateEventError( ErrorMsg )
        if aveRecurYrs < 2:
            ErrorMsg = "aveRecurYrs must be >= 2!!!"
            raise CreateEventError( ErrorMsg )
        self.name = name
        self.poisson = scstats.poisson( float(aveRecurYrs) )
        scaler = highMag - lowMag
        self.uniform = scstats.uniform( loc=lowMag, scale=scaler )
        return
    
    def nextEventDays( self, rstate ):
        """Select the arrival time for the next event.
        
        The Poisson distribution produces the arrival time in years. This is
        converted to days and should be considered days in future from 
        occurrence.

        Parameters
        ----------
        rstate : np.random_state
            Random state object that is the sampler.

        Returns
        -------
        eventOffsetDays : float
            Decimal days until the event occurs.

        """
        offsetYrs = self.poisson.rvs( size=1, random_state=rstate )
        eventOffsetDays = float( offsetYrs[0] ) * 365.25
        return eventOffsetDays
    
    def nextEventDailyMag( self, rstate ):
        """Sample the magnitude in mm for the next event.
        
        A scaled and shifted Uniform distribution is used to sample the event 
        magnitude.

        Parameters
        ----------
        rstate : np.random_state
            Random state object that is the sampler.

        Returns
        -------
        eventMag_mm : float
            Event magnitude in millimeters for the next event.

        """
        evMag = self.uniform.rvs( size=1, random_state=rstate )
        eventMag_mm = float( evMag[0] )
        return eventMag_mm


class EventRecurSampler(object):
    """An event recurrence probability sampler. 
    
    Use this to track the random state and produce reproduceable sampling 
    sequences. This uses the numpy random.uniform distribution to draw the 
    random numbers.
    """
    
    def __init__( self, evrecur_sample_seed=None ):
        """Default initialization method
        
        KWargs:
            evrecur_sample_seed (int): the seed to use for the sampler

        """
        super().__init__()
        self.ranstate = np.random.RandomState(seed=evrecur_sample_seed)
    
    def getSingleVal( self, ):
        """Produce a single value from the random sampler from a uniform
        distribution between [0.0 - 1.0]
        """
        return self.ranstate.uniform(low=0.0, high=1.0)
    
    def getSampleArray( self, N ):
        """Produce a sample of random numbers between 0.0 and 1.0 of size N.
        Will return a numpy array"""
        return self.ranstate.uniform(low=0.0, high=1.0, size=N)


class EventMagSampler(object):
    """An event magnitude probability sampler. 
    
    Use this to track the random state and produce reproduceable sampling 
    sequences. This uses the numpy random.uniform distribution to draw the 
    random numbers.
    """
    
    def __init__( self, evmag_sample_seed=None ):
        """Default initialization method
        
        KWargs:
            evmag_sample_seed (int): the seed to use for the sampler

        """
        super().__init__()
        self.ranstate = np.random.RandomState(seed=evmag_sample_seed)
    
    def getSingleVal( self, ):
        """Produce a single value from the random sampler from a uniform
        distribution between [0.0 - 1.0]
        """
        return self.ranstate.uniform(low=0.0, high=1.0)
    
    def getSampleArray( self, N ):
        """Produce a sample of random numbers between 0.0 and 1.0 of size N.
        Will return a numpy array"""
        return self.ranstate.uniform(low=0.0, high=1.0, size=N)


# custom functions for dealing with event objects
def setEvents( recurSeed, magSeed ):
    """Initial event setup

    Parameters
    ----------
    recurSeed : int
        Event recurrence interval sampling seed.
    magSeed : int
        Event magnitude sampling seed.

    Returns
    -------
    None.

    """
    # imports
    from EAAWG_Inputs import EVENT_DICT, EVENT_KEYS
    # global
    global ACT_EVENT_DICT
    # parameters
    # locals
    EvTrackDict = dict()
    evCnt = 1
    # start
    for evName in EVENT_KEYS:
        evVals = EVENT_DICT[evName]
        curEvent = ExtremeEvent(evVals[0], evVals[1][0], evVals[1][1], 
                                name=evName)
        curRecurSampler = EventRecurSampler(evrecur_sample_seed=(recurSeed+evCnt))
        curMagSampler = EventMagSampler(evmag_sample_seed=(magSeed+evCnt))
        EvTrackDict[evName] = [ curEvent, curRecurSampler, curMagSampler ]
        evCnt += 1
    # end of custom event for
    ACT_EVENT_DICT = EvTrackDict
    # return
    return


def setTrackers():
    """Setup and populate the initial events and tracking information

    Returns
    -------
    None.

    """
    # imports
    from EAAWG_Inputs import EVENT_KEYS
    # globals
    global ACT_EVENT_DICT, NEXT_EVENT_DICT, TTRIG_EVENT_TRACK_DICT
    # parameters
    # locals
    toBeTrigDict = dict()
    trackingDict = dict()
    # start
    for evName in EVENT_KEYS:
        evRecSamp = ACT_EVENT_DICT[evName][1]
        evMagSamp = ACT_EVENT_DICT[evName][2]
        curExEvent = ACT_EVENT_DICT[evName][0]
        curTrigTime = curExEvent.nextEventDays(evRecSamp.ranstate)
        curTrigMag = curExEvent.nextEventDailyMag(evMagSamp.ranstate)
        trackingDict[evName] = list()
        toBeTrigDict[evName] = [ curTrigTime, curTrigMag ]
    # end event for
    NEXT_EVENT_DICT = toBeTrigDict
    TTRIG_EVENT_TRACK_DICT = trackingDict
    # return
    return


def updateTriggeredEvent( evCnt, curTime, curElapsDays ):
    """Update an event that has been triggered

    Parameters
    ----------
    evCnt : int
        Index for event name in event dictionaries.
    curTime : pd.Timestamp
        Timestep for implementation of triggered event
    curElapsDays : int
        Current number of elapsed days into the simulation

    Returns
    -------
    None

    """
    # imports
    from EAAWG_Inputs import EVENT_KEYS
    # globals
    global ACT_EVENT_DICT, NEXT_EVENT_DICT, TTRIG_EVENT_TRACK_DICT
    # parameters
    # locals
    # start
    evRecSamp = ACT_EVENT_DICT[EVENT_KEYS[evCnt]][1]
    evMagSamp = ACT_EVENT_DICT[EVENT_KEYS[evCnt]][2]
    curExEvent = ACT_EVENT_DICT[EVENT_KEYS[evCnt]][0]
    curTrigTime = curExEvent.nextEventDays(evRecSamp.ranstate)
    curTrigMag = curExEvent.nextEventDailyMag(evMagSamp.ranstate)
    TTRIG_EVENT_TRACK_DICT[EVENT_KEYS[evCnt]].append( [curTime, 
                                NEXT_EVENT_DICT[EVENT_KEYS[evCnt]][1] ] )
    nextTrigTime = float( curElapsDays ) + curTrigTime
    NEXT_EVENT_DICT[EVENT_KEYS[evCnt]][0] = nextTrigTime
    NEXT_EVENT_DICT[EVENT_KEYS[evCnt]][1] = curTrigMag
    # return
    return


def outputEventTracking( RealNum, DT_INDEX ):
    """Output event information for this realization

    Parameters
    ----------
    RealNum : int
        Current realization index.
    DT_INDEX : pd.DateTimeIndex
        Index for all outputs

    Returns
    -------
    None.

    """
    # imports
    from os import path
    import pandas as pd
    from EAAWG_Inputs import EVENT_KEYS, OUT_LABEL, OUT_DIR, OUT_SUB_DIR
    # globals
    global NEXT_EVENT_DICT, TTRIG_EVENT_TRACK_DICT
    # parameters
    # locals
    # start
    EndTS = DT_INDEX[len(DT_INDEX)-1]
    StartTS = DT_INDEX[0]
    # prepare text file for output
    OFileName = "%s_R%d_EventsSummary.txt" % (OUT_LABEL, RealNum)
    OutFP = path.normpath( path.join( OUT_DIR, OUT_SUB_DIR, OFileName ) )
    with open( OutFP, 'w' ) as OF:
        OF.write("Event summary for realization %d \n" % RealNum )
        OF.write("    Simulation from %s through %s \n\n" % 
                  (StartTS.strftime("%Y-%m-%d"), EndTS.strftime("%Y-%m-%d")))
        for evName in EVENT_KEYS:
            OF.write("    %s events \n" % evName )
            trigList = TTRIG_EVENT_TRACK_DICT[evName]
            if len(trigList) > 0:
                for levList in trigList:
                    OF.write("        %s  %5.1f mm \n" % 
                             ( levList[0].strftime("%Y-%m-%d"), levList[1] ) )
                # end for event
            else:
                nextEventList = NEXT_EVENT_DICT[evName]
                futureTS = StartTS + pd.Timedelta(days=nextEventList[0])
                OF.write("        No events triggered, first event %s  %5.1f mm \n" 
                         % ( futureTS.strftime("%Y-%m-%d"), nextEventList[1] ) )
            # end if event triggered
            OF.write("\n")
        # end for custom event
    # end with and close file
    # return
    return

def outputBigEvents( RealNum, DT_INDEX ):
    """Output event information for this realization

    Parameters
    ----------
    RealNum : int
        Current realization index.
    DT_INDEX : pd.DateTimeIndex
        Index for all outputs

    Returns
    -------
    None.

    """
    # imports
    from os import path
    import pandas as pd
    from EAAWG_Inputs import EVENT_KEYS, OUT_LABEL, OUT_DIR, OUT_SUB_DIR
    # globals
    global NEXT_EVENT_DICT, TTRIG_EVENT_TRACK_DICT
    # parameters
    Key50 = "50-year"
    Key100 = "100-year"
    # locals
    # start
    EndTS = DT_INDEX[len(DT_INDEX)-1]
    StartTS = DT_INDEX[0]
    # prepare text file for output
    # examine our events to see if need to output
    ShouldOutput = False
    if Key50 in EVENT_KEYS:
        if len( TTRIG_EVENT_TRACK_DICT[Key50] ) > 0:
            ShouldOutput = True
        # end inner if
    # end outer if
    if Key100 in EVENT_KEYS:
        if len( TTRIG_EVENT_TRACK_DICT[Key100] ) > 0:
            ShouldOutput = True
        # end inner if
    # end outer if
    if ShouldOutput:
        OFileName = "%s_R%d_EventsSummary.txt" % (OUT_LABEL, RealNum)
        OutFP = path.normpath( path.join( OUT_DIR, OUT_SUB_DIR, OFileName ) )
        with open( OutFP, 'w' ) as OF:
            OF.write("Event summary for realization %d \n" % RealNum )
            OF.write("    Simulation from %s through %s \n\n" % 
                      (StartTS.strftime("%Y-%m-%d"), EndTS.strftime("%Y-%m-%d")))
            for evName in EVENT_KEYS:
                OF.write("    %s events \n" % evName )
                trigList = TTRIG_EVENT_TRACK_DICT[evName]
                if len(trigList) > 0:
                    for levList in trigList:
                        OF.write("        %s  %5.1f mm \n" % 
                                 ( levList[0].strftime("%Y-%m-%d"), levList[1] ) )
                    # end for event
                else:
                    nextEventList = NEXT_EVENT_DICT[evName]
                    futureTS = StartTS + pd.Timedelta(days=nextEventList[0])
                    OF.write("        No events triggered, first event %s  %5.1f mm \n" 
                             % ( futureTS.strftime("%Y-%m-%d"), nextEventList[1] ) )
                # end if event triggered
                OF.write("\n")
            # end for custom event
        # end with and close file
    # end if output
    # return
    return


def cleanAllEnd():
    """Convenience method to clean/delete all samplers at end

    """
    global ACT_EVENT_DICT, NEXT_EVENT_DICT, TTRIG_EVENT_TRACK_DICT
    # now set to None
    ACT_EVENT_DICT = None
    NEXT_EVENT_DICT = None
    TTRIG_EVENT_TRACK_DICT = None
    # end
    return


def setAllBegin():
    """Convenience method to seta all at beginning

    """
    global ACT_EVENT_DICT, NEXT_EVENT_DICT, TTRIG_EVENT_TRACK_DICT
    # now set to None
    ACT_EVENT_DICT = dict()
    NEXT_EVENT_DICT = dict()
    TTRIG_EVENT_TRACK_DICT = dict()
    # return
    return


#EOF