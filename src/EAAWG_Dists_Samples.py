# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_Dists_Samples
   :platform: Windows, Linux
   :synopsis: Contains collections of instantiated distributions and samplers.

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

Contains instantiated distributions and samplers

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

# project imports
import EAAWG_PrecipDepth as WGPD
import EAAWG_SpellLength as WGSL
import EAAWG_Inputs as WGI

# Distribution dictionaries
# spell distribution dictionaries
DDRY_SPELL_DISTS = dict()
DWET_SPELL_DISTS = dict()
# precip depth distribution dictionaries
DP_DEPTH_DISTS = dict()
"""Dictionary of precipitation depth distributions for
data periods, applies to both pathways.
"""

# Sampling dictionaries
# spell samplers
DDRY_SPELL_SAMP = dict()
DWET_SPELL_SAMP = dict()
# precip depth samplers
DP_DEPTH_SAMP = dict()
"""Dictionary of precipitation depth sampler for data periods,
applies to both pathways.
"""

# sample trackers
ST_DRYSPELL = dict()
ST_WETSPELL = dict()
ST_PDEPTH = dict()
"""Dictionary of tracked sample values"""

#-----------------------------------------------------------------------
# convenience set-up functions
def setTrackers( mon ):
    """Set-up our sample tracking dictionaries

    Parameters
    ----------
    mon : int
        current month.

    Returns
    -------
    None.

    """
    # globals
    global DDRY_SPELL_DISTS, DWET_SPELL_DISTS, DP_DEPTH_DISTS, DDRY_SPELL_SAMP 
    global DWET_SPELL_SAMP, DP_DEPTH_SAMP, ST_DRYSPELL, ST_WETSPELL, ST_PDEPTH
    # start
    MonthInts = list( range(1, 13, 1) )
    DrySYrDict = {}
    WetSYrDict = {}
    for jJ in MonthInts:
        dSamp = DDRY_SPELL_SAMP[jJ]
        DrySYrDict[jJ] = DDRY_SPELL_DISTS[jJ].ranval1(dSamp.ranstate)
        wSamp = DWET_SPELL_SAMP[jJ]
        WetSYrDict[jJ] = DWET_SPELL_DISTS[jJ].ranval1(wSamp.ranstate)
    # end of month for
    ST_DRYSPELL = DrySYrDict
    ST_WETSPELL = WetSYrDict
    # precip depth
    PDMonDict = {}
    for kK in MonthInts:
        pdSamp = DP_DEPTH_SAMP[kK]
        PDMonDict[kK] = DP_DEPTH_DISTS[kK].ranval1( pdSamp.ranstate, mon )
    # end of month for
    ST_PDEPTH = PDMonDict
    # end
    return

def setDistributions( pdSampSeed, wetSSampSeed, drySSampSeed ):
    """Go through and create the objects for all distributions.
    The parameters are specified in the input module. A parallel dictionary structure of
    distributions is created for the inputs.

    Args:
        pdSampSeed (int): precipitation depth sampling seed
        wetSSampSeed (int): wet state sampling seed
        drySSampSeed (int): dry state sampling seed

    """
    # imports
    # globals
    global DDRY_SPELL_DISTS, DWET_SPELL_DISTS, DP_DEPTH_DISTS, DDRY_SPELL_SAMP 
    global DWET_SPELL_SAMP, DP_DEPTH_SAMP
    # start
    MonthInts = list( range(1, 13, 1) )
    DrySYrDict = dict()
    WetSYrDict = dict()
    DrySmpDict = dict()
    WetSmpDict = dict()
    for jJ in MonthInts:
        drydist = WGSL.NegBinomial( WGI.DRY_SPELL_PARAMS[jJ][0],
                                    WGI.DRY_SPELL_PARAMS[jJ][1],
                                    location=WGI.DRY_SPELL_PARAMS[jJ][2],
                                    name="Dry spell, Month %d" % jJ )
        DrySYrDict[jJ] = drydist
        wetdist = WGSL.NegBinomial( WGI.WET_SPELL_PARAMS[jJ][0],
                                    WGI.WET_SPELL_PARAMS[jJ][1],
                                    location=WGI.WET_SPELL_PARAMS[jJ][2],
                                    name="Wet spell, Month %d" % jJ )
        WetSYrDict[jJ] = wetdist
        drysamp = WGSL.DryStateSampler(dry_state_seed=drySSampSeed)
        DrySmpDict[jJ] = drysamp
        wetsamp = WGSL.WetStateSampler(wet_state_seed=wetSSampSeed)
        WetSmpDict[jJ] = wetsamp
    # end of month for
    DDRY_SPELL_DISTS = DrySYrDict
    DWET_SPELL_DISTS = WetSYrDict
    DDRY_SPELL_SAMP = DrySmpDict
    DWET_SPELL_SAMP = WetSmpDict
    # end of data period for
    # next do the data period precipitation depth
    PDMonDict = dict()
    PSMonDict = dict()
    for kK in MonthInts:
        depdist = WGPD.Gamma2PCont( WGI.PRE_DEPTH_PARAMS[kK][0],
                                    WGI.PRE_DEPTH_PARAMS[kK][1], 
                                    loc=WGI.PRE_DEPTH_PARAMS[kK][2],
                                    scale=WGI.PRE_DEPTH_PARAMS[kK][3],
                                    name="Precip depth, Month %d" % kK )
        PDMonDict[kK] = depdist
        depsamp = WGPD.PrecipSampler(pd_sample_seed=pdSampSeed)
        PSMonDict[kK] = depsamp
    # end of month for
    DP_DEPTH_DISTS = PDMonDict
    DP_DEPTH_SAMP = PSMonDict
    # end
    return

def sampleAll( mon ):
    """Sample all distributions for every time step of every realization

    Parameters
    ----------
    mon : int
        Current month index.

    Returns
    -------
    None.

    """
    # global
    global DDRY_SPELL_DISTS, DWET_SPELL_DISTS, DP_DEPTH_DISTS, DDRY_SPELL_SAMP 
    global DWET_SPELL_SAMP, DP_DEPTH_SAMP, ST_DRYSPELL, ST_WETSPELL, ST_PDEPTH
    # start
    MonthInts = list( range(1, 13, 1) )
    for jJ in MonthInts:
        dSamp = DDRY_SPELL_SAMP[jJ]
        wSamp = DWET_SPELL_SAMP[jJ]
        ST_DRYSPELL[jJ] = DDRY_SPELL_DISTS[jJ].ranval1(dSamp.ranstate)
        ST_WETSPELL[jJ] = DWET_SPELL_DISTS[jJ].ranval1(wSamp.ranstate)
    # end of month for
    # precip depth
    for kK in MonthInts:
        pdSamp = DP_DEPTH_SAMP[kK]
        ST_PDEPTH[kK] = DP_DEPTH_DISTS[kK].ranval1( pdSamp.ranstate, mon )
    # end of month for
    # end
    return

def cleanAllEnd():
    """Convenience method to clean/delete all samplers at end

    """
    global DDRY_SPELL_DISTS, DWET_SPELL_DISTS, DP_DEPTH_DISTS, DDRY_SPELL_SAMP 
    global DWET_SPELL_SAMP, DP_DEPTH_SAMP, ST_DRYSPELL, ST_WETSPELL, ST_PDEPTH
    # now set to None
    DDRY_SPELL_DISTS = None
    DWET_SPELL_DISTS = None
    DDRY_SPELL_SAMP = None
    DWET_SPELL_SAMP = None
    DP_DEPTH_DISTS = None
    DP_DEPTH_SAMP = None
    ST_WETSPELL = None
    ST_PDEPTH = None
    ST_DRYSPELL = None
    # end
    return

def setAllBegin():
    """Convenience method to set all at beginning

    """
    global DDRY_SPELL_DISTS, DWET_SPELL_DISTS, DP_DEPTH_DISTS, DDRY_SPELL_SAMP 
    global DWET_SPELL_SAMP, DP_DEPTH_SAMP, ST_DRYSPELL, ST_WETSPELL, ST_PDEPTH
    # now set to None
    DDRY_SPELL_DISTS = dict()
    DWET_SPELL_DISTS = dict()
    DDRY_SPELL_SAMP = dict()
    DWET_SPELL_SAMP = dict()
    DP_DEPTH_DISTS = dict()
    DP_DEPTH_SAMP = dict()
    ST_WETSPELL = dict()
    ST_PDEPTH = dict()
    ST_DRYSPELL = dict()
    # end
    return

#EOF