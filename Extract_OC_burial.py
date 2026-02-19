from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import date, datetime, timedelta
import os
path_spm = os.environ.get('PATHSPM')
path = os.environ.get('PATH_OUTPUT')
path_fig = os.environ.get('PATH_FIG')
import sys
sys.path.extend([path_spm])
from plotbcz import *
from matplotlib import cm
cmap2 = cm.jet
import cmocean.cm as cm
cmap2 = cm.solar_r

# script plots carbon sequestration over the NS domain
# it used to be bdtrc_02 bioirrigating below 10cm horizon
# or a settling flux for the semi-labile carbon
# but thanks to the implemented routine into ROMS code the extraction got easier

yrange = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2011,2012,2013,2014,2015,   2021,2022,2023]
burial_year = []
temp_year = []

for y in yrange:
	nc0 = Dataset(path+'/Hindcast_CE2COAST_HIS_%s_2c_atm3.nc' %(y), 'r', format='NETCDF4')
	nc1 = Dataset(path+'/Hindcast_CE2COAST_AVG_%s_2c_atm3.nc' %(y), 'r', format='NETCDF4')

	if y == yrange[0]:
		seq = np.zeros(np.shape(nc0.variables['h']))

	lats,lons,mask = nc0.variables['lat_rho'][:],nc0.variables['lon_rho'][:],nc0.variables['mask_rho'][:]

	b = np.nanmean(np.nansum(nc0.variables['boxcon_04'][:],axis=1),axis=0) * 86400 * 30
	burial_year.append(np.nansum(b))
	temp_year.append(np.nanmean(nc1.variables['temp'][:,0]))
	seq += b

seq = seq / len(yrange)

fig, ax = plt.subplots(figsize=(16, 10))

m1 = Basemap(projection='merc',llcrnrlat=51,urcrnrlat=59.0,llcrnrlon=-3.2,urcrnrlon=10.0,lat_0=3,resolution='i', ax=ax)
x2, y2 = m1(lons, lats)
cax = make_axes_locatable(ax).append_axes("right", size=0.4, pad=0.15)
bounds = [0,10,20,30,40,50,60]

m1.drawcoastlines()
seq[seq>bounds[-1]] = bounds[-1]-0.1
norm = mpl.colors.BoundaryNorm(bounds, cmap2.N)

CS2 = m1.contourf(x2,y2,seq,clevs=bounds,cmap=cmap2,norm=norm)
cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap2, boundaries = bounds, orientation='vertical',norm=norm)
fig.savefig(path_fig + "/Carbon_sequestration_clima.png", dpi=300, bbox_inches='tight')

