# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_PrecipDepth
   :platform: Windows, Linux
   :synopsis: Contains logic for determining and selecting the precipitation 
              depth for a day

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
from EAAWG_Inputs import MON_MAX_PP

# global parameters
WD_THRESH = 0.255
"""Wet dry threshold as used in this study in millimeters. It provides the
minimum possible precipitation and so bounds our PDFs and CDFs"""



class CreateDistError(Exception):
    """Custom exception error for when something happens with a distribution
    """
    def __init__(self, arg):
        self.args = arg


class Gamma2PCont(object):
    """Gamma distribution object, 2 paramater gamma, for use in the Weather 
    Generator. Gamma distributions are used for wet day precipitation depths.
    In this form of the 2 paramater gamma, a and c are shape parameters. a 
    must be greater than zero, and c must not be equal to zero. This 
    generalized form of the gamma also allows for specification of location 
    (loc) and scale parameters."""
    
    def __init__(self, a, c, loc=0, scale=1.0, name="Custom 2 parameter gamma" ):
        """Override of default initialization method"""
        super().__init__()
        # do some checks and assign our properties
        if not isinstance( a, float ):
            ErrorMsg = "a must be a float!!!"
            raise CreateDistError( ErrorMsg )
        if not isinstance( c, float ):
            ErrorMsg = "c must be a float!!!"
            raise CreateDistError( ErrorMsg )
        if a <= 0.0:
            ErrorMsg = "a must be greater than zero!!!"
            raise CreateDistError( ErrorMsg )
        if c == 0.0:
            ErrorMsg = "c must not be equal to zero!!!"
            raise CreateDistError( ErrorMsg )
        self.name = name
        self.gamma = scstats.gengamma( a, c, loc=loc, scale=scale )
        return
    
    def ranval1( self, rstate, mon, maxes=MON_MAX_PP, thresh=WD_THRESH ):
        """With the specified val between 0.0 and 1.0, which is essentially
        a probability, return the corresponding value from the cdf.
        
        Args:
            rstate (np.random_state): the random state for the sampler
            mon (int): current month index
            maxes (dict): dictionary of maximum daily precipitation depth in
                          mm by month
            thresh (float): minimum precipitation depth in mm
            
        Returns:
            numd (int): the sampled wet day precipitation depth
            
        """
        numda = self.gamma.rvs( size=1, random_state=rstate )
        numd = numda[0]
        if numd > maxes[mon]:
            numd = maxes[mon]
        # end if
        if numd < thresh:
            numd = thresh
        # end if
        return numd
    
    def ranArray( self, asize, rstate, mon, maxes=MON_MAX_PP, thresh=WD_THRESH ):
        """Sample asize values and put them in an array
        
        Args:
            asize (int) = size of the return array
            rstate (np.random_state): the random state for the sampler
            mon (int) = current month index
            maxes (dict): dictionary of maximum daily precipitation depth in 
                          mm by month
            thresh (float): minimum precipitation depth in mm
            
        Returns:
            numda (np.array): the corresponding array of wet day precip depths
            
        """
        numda = self.gamma.rvs( size=asize, random_state=rstate )
        numda = np.where( numda > maxes[mon], maxes[mon], numda )
        numda = np.where( numda < thresh, thresh, numda )
        return numda


class PrecipSampler(object):
    """A precipitation depth probability sampler. 
    
    Use this to track the random
    state and produce reproduceable sampling sequences. This uses the numpy
    random.uniform distribution to draw the random numbers.
    """
    
    def __init__( self, pd_sample_seed=None ):
        """Default initialization method
        
        KWargs:
            pd_sample_seed (int): the seed to use for the sampler

        """
        super().__init__()
        self.ranstate = np.random.RandomState(seed=pd_sample_seed)
    
    def getSingleVal( self, ):
        """Produce a single value from the random sampler from a uniform
        distribution between [0.0 - 1.0]
        """
        return self.ranstate.uniform(low=0.0, high=1.0)
    
    def getSampleArray( self, N ):
        """Produce a sample of random numbers between 0.0 and 1.0 of size N.
        Will return a numpy array"""
        return self.ranstate.uniform(low=0.0, high=1.0, size=N)


#EOF