# Code for indetifying deep convective initations in kscale models
# outputs text file list of CIs
# -40C initations cooling at least 30 degrees per hour with no pre-existing cold cloud within 50km

######################
# Modules            #
######################

import math as maths
import cartopy.crs as ccrs
from scipy.ndimage import label
import intake
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import healpy
import sys
import matplotlib as mpl
import easygems.healpix as egh
import warnings

# Suppress specific FutureWarnings matching the message pattern when using cat[...].to_dask()
warnings.filterwarnings(
    "ignore",
    message=".*The return type of `Dataset.dims` will be changed.*",
    category=FutureWarning,
)

######################
# Functions          #
######################

def get_nn_lon_lat_index(nside, lons, lats):
    """for subsetting HEALPix ICON out onto regular lat/lon grid"""
    lons2, lats2 = np.meshgrid(lons, lats)
    return xr.DataArray(
        healpy.ang2pix(nside, lons2, lats2, nest=True, lonlat=True),
        coords=[("lat", lats), ("lon", lons)],
    )

def olr_to_bt(olr):
    """Application of Stefan-Boltzmann law - converts outgoing longwave to cloud top temperature"""
    sigma = 5.670373e-8
    tf = (olr/sigma)**0.25
    #Convert from bb to empirical BT (degC) - Yang and Slingo, 2001
    a = 1.228
    b = -1.106e-3
    Tb = (-a + np.sqrt(a**2 + 4*b*tf))/(2*b)
    return Tb - 273.15

########################
# Model Data           #
########################

cat = intake.open_catalog('https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml')['online']

zoom = 8 # ~9 km

idx = get_nn_lon_lat_index(
    2**zoom, np.linspace(27, 47, 200), np.linspace(-8, 7, 150))  # East Africa domain 10km

globsim = cat['um_glm_n2560_RAL3p3']
ds = globsim(zoom=zoom, time='PT1H').to_dask()

ref_lon_lat = ds.rlut.isel(time=1740,cell=idx)
ref = ref_lon_lat.values
datestr = str(ref_lon_lat.coords['time'].values)
lon1D = ref_lon_lat.coords['lon'].values
lat1D = ref_lon_lat.coords['lat'].values
lon2D, lat2D = np.meshgrid(lon1D, lat1D)
yr = datestr[0:4]
mn = datestr[5:7]
dy = datestr[8:10]
hr = datestr[11:13]
mt = datestr[14:16]

#########################
# blob cutting          #
#########################

fi = open('East-Africa_CIs_April_1200-1800LT.txt','w')

count = 0

for ddex in range(1737,2457,24):   # looping over days in April 2020 for East Africa (UTC+3)

  for tdex in range(0,6):     # looping over hours 1200 - 1800 LT

    sdex = ddex + tdex

    # current timestep
  
    rlut_lon_lat_c = ds.rlut.isel(time=sdex, cell=idx) # outgoing longwave radition OLR#
    rlut_c = rlut_lon_lat_c.values
    Tb_c = olr_to_bt(rlut_c) # convert OLR to brightness temperature

    datestr = str(rlut_lon_lat_c.coords['time'].values)
    cdate = datestr[0:10]
    ctim = datestr[11:16]

    #previous timestep
  
    rlut_lon_lat_p = ds.rlut.isel(time=sdex-1, cell=idx) # outgoing longwave radition OLR#
    rlut_p = rlut_lon_lat_p.values
    Tb_p = olr_to_bt(rlut_p) # convert OLR to brightness temperature

    Tb_c [Tb_c > -40] = 0
  
    labels, numL = label(Tb_c)
    u, inv = np.unique(labels, return_inverse=True)
    n=np.bincount(inv)
    goodinds= u[n > 0]

    for gi in goodinds:
        if gi == 0:
            continue
                        
        inds = np.where(labels == gi)
        
        latmax, latmin = lat2D[inds].max(), lat2D[inds].min()
        lonmax, lonmin = lon2D[inds].max(), lon2D[inds].min()

        marg = 0.5
        i, j = np.where( (lon2D>lonmin-marg) & (lon2D<lonmax+marg) & (lat2D>latmin-marg) & (lat2D<latmax+marg) )
        blat = lat2D[i.min():i.max()+1, j.min():j.max()+1]
        blon = lon2D[i.min():i.max()+1, j.min():j.max()+1]
        bcloud_c = Tb_c[i.min():i.max()+1, j.min():j.max()+1]
        bcloud_p = Tb_p[i.min():i.max()+1, j.min():j.max()+1]

        # filter pre-existing cold cloud within 50km
        
        if np.min(bcloud_p) < -40:
            continue

        Tdiff = np.min(bcloud_c) - np.min(bcloud_p)

        # filter for rapidly cooling pixels within 50km
        
        if Tdiff > -20:
            continue

        # indentify lat / lon of coldest pixel

        CTTmin = np.min(bcloud_c)
        mindex_x = np.where(bcloud_c == bcloud_c.min())[0][0]
        mindex_y = np.where(bcloud_c == bcloud_c.min())[1][0]

        clat = blat[mindex_x,mindex_y]
        clon = blon[mindex_x,mindex_y]

        print(stdate, cdate, ctim, clon, clat, sdex, file=fi)

print(count,"storm cases")

fi.close()
