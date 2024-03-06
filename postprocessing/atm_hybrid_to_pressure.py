#!/usr/bin/python
"""Vertical interpolation from atmospheric hybrid to pressure levels.

This script will read data on hybrid levels and interpolate them linearly onto 
pressure levels in log-pressure coordinates.

Contact: j.dejong3@uu.nl
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import sys
import glob
import numpy as np
from numba import float32, float64, guvectorize

lvls_p = sorted([200,250,300,500,700,900]) # new pressure levels in hPa
vars_interp = ['U'] # list of variables to interpolate
inputdir = "/home/jasperdj/nwo2021025/archive/lres_b.e10.B2000_CAM5.f09_g16.feedforward_2050.001/atm/hist/"
outputdir = "/home/jasperdj/esmlunch/"
inputfiles = sorted(glob.glob(inputdir+"lres_b.e10.B2000_CAM5.f09_g16.feedforward_2050.001.cam2.h0.2000-??.nc"))
outputfile = "lres.f09_g16.U.200001-200912.nc"


@guvectorize(
    "(float64[:], float64[:], float64[:], float32[:])",
    " (n), (n), (m) -> (m)",
    nopython=True,
)
def interp1d_gu(f, x, xi, out):
    """Interpolate one-dimensional field f(x) to xi in ln(x) coordinates."""
    i, imax, x0, f0 = 0, len(xi), x[0], f[0]
    while xi[i]<x0 and i < imax:
        out[i] = np.nan      
        i = i + 1 
    for x1,f1 in zip(x[1:], f[1:]):
        while xi[i] <= x1 and i < imax:
            out[i] = (f1-f0)/np.log(x1/x0)*np.log(xi[i]/x0)+f0
            i = i + 1
        x0, f0 = x1, f1
    while i < imax:
        out[i] = np.nan
        i = i + 1


def main():
	"""Open dataset, assign new pressure coordinate and interpolate.
	
	Data should be sorted such that pressure and lev (core dim) are increasing.
	"""
    ds = xr.open_mfdataset(inputfiles, parallel=True) # need a recent version of xarray 
    pres = (ds.hyam * ds.P0 + ds.hybm * ds.PS) # full 3D pressure  
    ds = ds.assign_coords({'plev':np.array(lvls_p, dtype='float64')*100})
    ds.plev.attrs.update({'units':'Pa'})
    dsi = xr.apply_ufunc(
		interp1d_gu, ds[['plev',*vars_interp]], pres, ds.plev, 
		input_core_dims=[['lev'],['lev'],['plev']], 
		output_core_dims=[['plev']], 
		exclude_dims=set(('lev',)), 
		dask='allowed', 
		on_missing_core_dim='drop'
    )
    dsi.to_netcdf(os.path.join(outputdir,outputfile))
    

if __name__ == '__main__':
	main()
