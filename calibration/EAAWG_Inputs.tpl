ptf $
# -*- coding: utf-8 -*-
"""
.. module:: EAAWG_Inputs
   :platform: Windows, Linux
   :synopsis: Provides input specification and shared module information

.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

Provides the inputs module where all input values and shared parameters can
be entered. This functionality can be replaced with a GUI in the future. Also 
holds the shared  data structures.

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

import pandas as pd
import pickle
import os

#------------------------------------------------------------------------
# Global simulation settings - these generally need to be changed
CUR_DIR = os.getcwd()
OUT_DIR = os.path.normpath( os.path.join( CUR_DIR, "Results" ) )
"""Location for model outputs"""
OUT_LABEL = r'Frio'
"""Label to use for outputting files to OUT_DIR"""
OUT_SUB_DIR = "Simulated"
"""Output subdirectory"""
START_DATE = pd.Timestamp(2024, 1, 1, 0, )
"""Starting time for production of the stochastic synthetic time series"""
END_DATE = pd.Timestamp( 2060, 12, 31, 23, 59, )
"""Ending time for production of the stochastic synthetic time series"""
K_c = 1.00
"""Crop coefficient to go from ETo to PET. Coefficient derived by comparing
calculated ETo to independently calculated/measured PET from "station" data.
In this implementation, the crop coefficient is not used."""
#LAT_DEG = 30.018 # Blanco
#LAT_DEG = 29.741 # Cibolo
LAT_DEG = 29.678 # Frio
#LAT_DEG = 29.985 # Guadalupe
#LAT_DEG = 29.627 # Med-Cib
#LAT_DEG = 29.729 # Medina
#LAT_DEG = 29.701 # Nueces
#LAT_DEG = 29.574 # Sab-Med
#LAT_DEG = 29.617 # Sabinal
"""Latitude in degrees for the basin centroid"""

# parameters - usually do not need to be changed
WET_STATE = "wet"
"""Keyword for wet state"""
DRY_STATE = "dry"
"""Keyword for dry state"""

#------------------------------------------------------------------------
# Other weather parameter model files
bas = "Frio"
InFiler = os.path.normpath( os.path.join( CUR_DIR, "Inputs", 
                        "OWeath_Smooth_WetAve_1981-2010_DictDF.pkl" ) )
with open( InFiler, 'rb' ) as IF:
    InDictWA = pickle.load( IF )
# end with          
OW_WET_AVE = InDictWA[bas]
"""Average wet day quantities by day of the year. Contains Tmax, Tmin, and Tave"""
InFiler = os.path.normpath( os.path.join( CUR_DIR, "Inputs", 
                        "OWeath_Smooth_DryAve_1981-2010_DictDF.pkl" ) )
with open( InFiler, 'rb' ) as IF:
    InDictDA = pickle.load( IF )
# end with
OW_DRY_AVE = InDictDA[bas]
"""Average dry day quantities by day of the year. Contains Tmax, Tmin, and Tave"""
InFiler = os.path.normpath( os.path.join( CUR_DIR, "Inputs", 
                        "OWeath_Smooth_WetStd_1981-2010_DictDF.pkl" ) )
with open( InFiler, 'rb' ) as IF:
    InDictWS = pickle.load( IF )
# end with
OW_WET_STD = InDictWS[bas]
"""Average wet day standard deviations for Tmax, Tmin, and Tave"""
InFiler = os.path.normpath( os.path.join( CUR_DIR, "Inputs", 
                        "OWeath_Smooth_DryStd_1981-2010_DictDF.pkl" ) )
with open( InFiler, 'rb' ) as IF:
    InDictDS = pickle.load( IF )
# end with
OW_DRY_STD = InDictDS[bas]
"""Average dry day standard deviations for Tmax, Tmin, and Tave"""

# other weather parameters correlation matrices
A_DATA_LIST = [ [ 0.56219295, 0.2050754 ],
                [ -0.03543818, 0.68296039 ], ]
"""A matrix for Tmax and Tmin"""
InFiler = os.path.normpath( os.path.join( CUR_DIR, "Inputs", 
                        "OWeath_Rho0_1991-2020_DFDict.pkl" ) )
with open( InFiler, 'rb' ) as IF:
    InDictM0 = pickle.load( IF )
# end with
OW_M0_IN = InDictM0[bas]
"""M0 matrix for calculating daily error or residual term"""
B_DATA_LIST = [ [ 0.74093011, 0.0 ],
                [ 0.17681369, 0.72149121 ], ]
""" B matrix for Tmax and Tmin"""
InFiler = os.path.normpath( os.path.join( CUR_DIR, "Inputs", 
                        "OWeath_Rho1_1991-2020_DFDict.pkl" ) )
with open( InFiler, 'rb' ) as IF:
    InDictM1 = pickle.load( IF )
# end with
OW_M1_IN = InDictM1[bas]
"""M1 matrix for calculating daily error or residual term"""

#------------------------------------------------------------------------
# distribution specifications
# For spell length distributions assume that one distribution applies to
#  the entire study area
DRY_SPELL_PARAMS = { 1 : [ $  m1_drysl_N   $, $ m1_drysl_P    $, 2.0 ], 
                     2 : [ $  m2_drysl_N   $, $ m2_drysl_P    $, 2.0 ], 
                     3 : [ $  m3_drysl_N   $, $ m3_drysl_P    $, 2.0 ], 
                     4 : [ $  m4_drysl_N   $, $ m4_drysl_P    $, 2.0 ], 
                     5 : [ $  m5_drysl_N   $, $ m5_drysl_P    $, 2.0 ], 
                     6 : [ $  m6_drysl_N   $, $ m6_drysl_P    $, 2.0 ], 
                     7 : [ $  m7_drysl_N   $, $ m7_drysl_P    $, 2.0 ], 
                     8 : [ $  m8_drysl_N   $, $ m8_drysl_P    $, 2.0 ], 
                     9 : [ $  m9_drysl_N   $, $ m9_drysl_P    $, 2.0 ], 
                    10 : [ $  m10_drysl_N  $, $ m10_drysl_P   $, 2.0 ], 
                    11 : [ $  m11_drysl_N  $, $ m11_drysl_P   $, 2.0 ], 
                    12 : [ $  m12_drysl_N  $, $ m12_drysl_P   $, 2.0 ], }
"""Dictionary of monthly parameters for negative binomial distributions to 
simulate dry spell durations. Parameters are N, P, and location. N must be 
greater than zero; P must be greater than zero and less than or equal to 1."""
WET_SPELL_PARAMS = { 1 : [ $  m1_wetsl_N   $, $ m1_wetsl_P    $, 1.0 ], 
                     2 : [ $  m2_wetsl_N   $, $ m2_wetsl_P    $, 1.0 ], 
                     3 : [ $  m3_wetsl_N   $, $ m3_wetsl_P    $, 1.0 ], 
                     4 : [ $  m4_wetsl_N   $, $ m4_wetsl_P    $, 1.0 ], 
                     5 : [ $  m5_wetsl_N   $, $ m5_wetsl_P    $, 1.0 ], 
                     6 : [ $  m6_wetsl_N   $, $ m6_wetsl_P    $, 1.0 ], 
                     7 : [ $  m7_wetsl_N   $, $ m7_wetsl_P    $, 1.0 ], 
                     8 : [ $  m8_wetsl_N   $, $ m8_wetsl_P    $, 1.0 ], 
                     9 : [ $  m9_wetsl_N   $, $ m9_wetsl_P    $, 1.0 ], 
                    10 : [ $  m10_wetsl_N  $, $ m10_wetsl_P   $, 1.0 ], 
                    11 : [ $  m11_wetsl_N  $, $ m11_wetsl_P   $, 1.0 ], 
                    12 : [ $  m12_wetsl_N  $, $ m12_wetsl_P   $, 1.0 ], }
"""Dictionary of monthly parameters for negative binomial distributions to 
simulate wet spell durations. Parameters are N, P, and location. N must be 
greater than zero; P must be greater than zero and less than or equal to 1."""

# Precipitation depth distribution parameters by month.
PRE_DEPTH_PARAMS = { 1 : [ $  m1_pdep_a    $, $  m1_pdep_c    $, 0.255, $  m1_pdep_sca  $ ],
                     2 : [ $  m2_pdep_a    $, $  m2_pdep_c    $, 0.255, $  m2_pdep_sca  $ ],
                     3 : [ $  m3_pdep_a    $, $  m3_pdep_c    $, 0.255, $  m3_pdep_sca  $ ],
                     4 : [ $  m4_pdep_a    $, $  m4_pdep_c    $, 0.255, $  m4_pdep_sca  $ ],
                     5 : [ $  m5_pdep_a    $, $  m5_pdep_c    $, 0.255, $  m5_pdep_sca  $ ],
                     6 : [ $  m6_pdep_a    $, $  m6_pdep_c    $, 0.255, $  m6_pdep_sca  $ ],
                     7 : [ $  m7_pdep_a    $, $  m7_pdep_c    $, 0.255, $  m7_pdep_sca  $ ],
                     8 : [ $  m8_pdep_a    $, $  m8_pdep_c    $, 0.255, $  m8_pdep_sca  $ ],
                     9 : [ $  m9_pdep_a    $, $  m9_pdep_c    $, 0.255, $  m9_pdep_sca  $ ],
                    10 : [ $  m10_pdep_a   $, $  m10_pdep_c   $, 0.255, $  m10_pdep_sca $ ],
                    11 : [ $  m11_pdep_a   $, $  m11_pdep_c   $, 0.255, $  m11_pdep_sca $ ],
                    12 : [ $  m12_pdep_a   $, $  m12_pdep_c   $, 0.255, $  m12_pdep_sca $ ], }
"""Dictionary of monthly parameters for the generalized, two parameter gamma 
distributions. The distribution parameters are 0 - a, 1 - c, 2 - loc, and
3 - scale. a must be greater than zero, and c must not be equal to zero."""

# maximum daily precipitation depth that can be used for each calendar month.
MON_MAX_PP = { 1 : 23.5,
               2 : 22.9,
               3 : 33.4,
               4 : 31.2,
               5 : 39.6,
               6 : 36.3,
               7 : 38.6,
               8 : 35.0,
               9 : 46.1,
               10 : 46.1,
               11 : 35.7,
               12 : 24.3, }
"""Maximum daily precipitation depth to allow for sampling from 2 parameter
gamma distributions. Us the 0.95 cumulative probability from the wet day
depth distributions for this."""

AVE_WET_TMAX_ADD = $ tmax_wet_add  $
"""Multiplier for average wet Tmax"""
AVE_DRY_TMAX_ADD = $ tmax_dry_add  $
"""Multiplier for average dry Tmax"""
AVE_WET_TMIN_ADD = $ tmin_wet_add  $
"""Multiplier for average wet Tmin"""
AVE_DRY_TMIN_ADD = $ tmin_dry_add  $
"""Multiplier for average dry Tmin"""

#-------------------------------------------------
# event dictionary
#    this one is for Frio
EVENT_DICT = { "2-year" : [ int(2), [96.0, $ ev_2yr_max    $],], 
               "5-year" : [ int(5), [141.0, $ ev_5yr_max    $],],
               "10-year" : [ int(10), [179.0, $ ev_10yr_max   $],],
               "25-year" : [ int(25), [236.0, $ ev_25yr_max   $],],
               "50-year" : [ int(50), [285.0, $ ev_50yr_max   $],],
               "100-year" : [ int(100), [343.0, $ ev_100yr_max  $],], }
EVENT_KEYS = list( EVENT_DICT.keys() )


#EOF