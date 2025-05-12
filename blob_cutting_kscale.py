###############################################
# cutting out cloud blobs in kscale model data
###############################################

import math as maths

import cartopy.crs as ccrs
import intake
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import easygems.healpix as egh

import warnings

# Suppress specific FutureWarnings matching the message pattern when using cat[...].to_dask()
warnings.filterwarnings(
    "ignore",
    message=".*The return type of `Dataset.dims` will be changed.*",
    category=FutureWarning,
)

cat = intake.open_catalog('https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml')['online']

# Show UM simulation keys:
#[key for key in cat if key.startswith('um_')]

# um_Africa_km4p4_RAL3P3_n1280_GAL9_nest - Africa LAM convection permitting nested GAL9 
# um_CTC_km4p4_RAL3P3_n1280_GAL9_nest - cyclic tropical channel convection permitting nested GAL9
# um_SAmer_km4p4_RAL3P3_n1280_GAL9_nest - South America LAM convection permitting nested GAL9
# um_SEA_km4p4_RAL3P3_n1280_GAL9_nest - South East Asia LAM convection permitting nested GAL9
# um_glm_n1280_CoMA9_TBv1p2 - global 10km parameter mass flux scheme
# um_glm_n1280_GAL9 - global 10km parameter 6A scheme
# um_glm_n2560_RAL3p3 - global 5km convection permitting 

sim = cat['um_glm_n2560_RAL3p3']

# tim = PT1H or PT3H
# zoom 0 - 10

ds = sim(zoom=8, time='PT1H').to_dask()

# list variable key
#list(ds)

# print information about specific variable
#print(ds.rlut)
