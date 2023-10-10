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
DRY_SPELL_PARAMS = { 1 : [ 3.915162388655806, .2977526333113043, 2.0 ],
                     2 : [ 5.312731046936296, .1988934314137431, 2.0 ],
                     3 : [ 3.603636540370255, .3428403644909424, 2.0 ],
                     4 : [ 2.704911745656023, .3295790193540502, 2.0 ],
                     5 : [ 6.953055983608310, .3176279916535475, 2.0 ],
                     6 : [ 5.058893398911076, .3010013234726136, 2.0 ],
                     7 : [ 6.095997869386332, .3669773538213991, 2.0 ],
                     8 : [ 3.177489210319657, .5094441058940438, 2.0 ],
                     9 : [ 3.735472701463134, .2370562126500331, 2.0 ],
                    10 : [ 4.606282010303014, .2057183456563516, 2.0 ],
                    11 : [ 3.000000000000000, .1761191459571701, 2.0 ],
                    12 : [ 3.756479818110243, .4118150622443279, 2.0 ], }
"""Dictionary of monthly parameters for negative binomial distributions to
simulate dry spell durations. Parameters are N, P, and location. N must be
greater than zero; P must be greater than zero and less than or equal to 1."""
WET_SPELL_PARAMS = { 1 : [ 4.510865915360880, .6451501587119356, 1.0 ],
                     2 : [ 2.475100135849106, .6696768747074860, 1.0 ],
                     3 : [ 1.675569413084325, .4654584992804364, 1.0 ],
                     4 : [ 2.556299700526773, .6365853084394536, 1.0 ],
                     5 : [ 3.469188391831542, .6069615384053114, 1.0 ],
                     6 : [ 3.115250111464483, .3393072067630571, 1.0 ],
                     7 : [ 1.830988852381821, .4184422242753315, 1.0 ],
                     8 : [ 1.714407941451173, .6526271751106262, 1.0 ],
                     9 : [ 2.355907294215807, .3353929549624517, 1.0 ],
                    10 : [ 2.381103082918876, .5459329374264738, 1.0 ],
                    11 : [ 2.138504952307077, .6512751306213146, 1.0 ],
                    12 : [ 2.440332865483639, .6921787755056502, 1.0 ], }
"""Dictionary of monthly parameters for negative binomial distributions to
simulate wet spell durations. Parameters are N, P, and location. N must be
greater than zero; P must be greater than zero and less than or equal to 1."""
 
# Precipitation depth distribution parameters by month.
PRE_DEPTH_PARAMS = { 1 : [ .7837115409758860, 1.347633111170500, 0.255, 5.610609853852632 ],
                     2 : [ 1.062803150829027, 1.341118445314213, 0.255, 6.391470021107692 ],
                     3 : [ 1.100735454047688, 1.272354996142382, 0.255, 8.348663906132668 ],
                     4 : [ 1.101843127926503, 1.210980653693744, 0.255, 6.991410943221600 ],
                     5 : [ 1.176806796894472, 1.534848341733679, 0.255, 7.508790937219044 ],
                     6 : [ 1.494046016720638, 1.372390726371365, 0.255, 10.52645546925398 ],
                     7 : [ 1.048404563404255, 1.491669130947536, 0.255, 9.731293267451120 ],
                     8 : [ .8286900863678252, 1.618721168571233, 0.255, 8.629731575543266 ],
                     9 : [ 1.101651079078155, 1.426794294882264, 0.255, 10.01386843532882 ],
                    10 : [ 1.536794107157059, 1.294664433648325, 0.255, 10.00000000000000 ],
                    11 : [ .7814800192944810, 1.396406230690366, 0.255, 6.756062486679516 ],
                    12 : [ .7951754020048494, 1.369403585110688, 0.255, 5.540832599808860 ], }
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
 
AVE_WET_TMAX_ADD = 5.098810943474146
"""Multiplier for average wet Tmax"""
AVE_DRY_TMAX_ADD = 6.941519835725408
"""Multiplier for average dry Tmax"""
AVE_WET_TMIN_ADD = .8523552116658896
"""Multiplier for average wet Tmin"""
AVE_DRY_TMIN_ADD = .8422783217775328
"""Multiplier for average dry Tmin"""
 
#-------------------------------------------------
# event dictionary
#    this one is for Frio
EVENT_DICT = { "2-year" : [ int(2), [96.0, 153.7479198867500],],
               "5-year" : [ int(5), [141.0, 247.6322773540412],],
               "10-year" : [ int(10), [179.0, 255.7052662952905],],
               "25-year" : [ int(25), [236.0, 290.9002966072379],],
               "50-year" : [ int(50), [285.0, 415.4495994011296],],
               "100-year" : [ int(100), [343.0, 498.1388101352216],], }
EVENT_KEYS = list( EVENT_DICT.keys() )
 
 
#EOF
