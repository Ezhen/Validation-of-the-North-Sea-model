from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import date, datetime, timedelta
import os
path_spm = os.environ.get('PATHSPM')
path = os.environ.get('PATH_OUTPUT')
path_fig = os.environ.get('PATH_FIG')
import sys
sys.path.extend([path_spm])
from plotbcz import *
from matplotlib import cm
import matplotlib as mpl
import cmocean.cm as cm
cmap2 = cm.algae
bounds=[0,1,2,3,4,5,7,10]
bounds=[0,0.5,1,2,3,4,5]
norm = mpl.colors.BoundaryNorm(bounds, cmap2.N)
import matplotlib
matplotlib.use('Agg')

# the script plot surface chlorophyll of ROMS simualtions
# against chl from CMEMS model

def image(title,m,axx,lons,lats,algae):
	x2, y2 = m(lons, lats)
	algae[algae>bounds[-1]]=bounds[-1]-0.01
	algae[algae<bounds[0]]=bounds[0]+0.01
	CS2 = m.contourf(x2,y2,algae,clevs=bounds,cmap=cmap2, norm=norm)
	m.plot(np.array([x2[0,0],x2[0,-1],x2[-1,-1],x2[-1,0],x2[0,0]]),np.array([y2[0,0],y2[0,-1],y2[-1,-1],y2[-1,0],y2[0,0]]),color='black',linewidth=1.0)

lo1 = lambda x: datetime(1900,1,1,0,0,0) + timedelta(days=x)
lo2 = lambda x: datetime(1993,1,1,0,0,0) + timedelta(seconds=x)
season_list = ['autumn-winter','spring','summer']
season = [[1,2,9,10,11,12],[3,4,5],[6,7,8]]

rr0 = Dataset(path+'/cmems_chl_sat_2006.nc', 'r', format='NETCDF4')

lon0,lat0 = np.meshgrid(rr0.variables['lon'][:],rr0.variables['lat'][:])
cmems_0 = np.nanmean(np.nanmean(rr0.variables['CHL'][:],axis=1),axis=1)

time_cmems = []
for i in range(12):
	time_cmems.append(lo1(int(rr0.variables['time'][i])))


year = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2011,2012,2013,2014,2015,2016,2017,2018,   2021,2022,2023]

thick = True
for y in range(len(year)): 

	nc0 = Dataset(path+'/Hindcast_CE2COAST_AVG_%s_2c_atm3.nc' %(year[y]), 'r', format='NETCDF4')
	units=nc0.variables['ocean_time'].units.split( )[2].split('-')
	lo2 = lambda x: datetime(int(units[0]),int(units[1]),int(units[2]),0,0,0) + timedelta(seconds=x)
	lats0,lons0,mask0 = nc0.variables['lat_rho'][:],nc0.variables['lon_rho'][:],nc0.variables['mask_rho'][:]

	light_koeffs = (1-np.exp(-0.1*np.cos(np.deg2rad(lats0))))
	light_base = (1-np.exp(-0.1*np.cos(np.deg2rad(52.0))))
	scale_light = light_koeffs / light_base

	chl_0 = nc0.variables['chlorophyll'][:,-1]
	c0 = np.nanmean(np.nanmean(chl_0,axis=1),axis=1)

	for s in range(3):

		if y == 0 and s == 0:
			chl_roms_arr = np.zeros((3,len(nc0.variables['chlorophyll'][0,-1]),len(nc0.variables['chlorophyll'][0,-1].T)))

		chl_sp0 = np.zeros((nc0.variables['chlorophyll'][0,-1].shape))
		
		k0 = 0

		for i in range(len(nc0.variables['ocean_time'])):
			if lo2(int(nc0.variables['ocean_time'][i])).month in season[s]:
				koeff = np.ones((nc0.variables['chlorophyll'][0,-1].shape)) * 0.05
				if thick == True:
					V = (nc0.variables['Cs_w'][-6:] - nc0.variables['Cs_w'][-7:-1])[:,np.newaxis,np.newaxis] * nc0.variables['h'][:][np.newaxis,:,:]
					nutr = np.nansum(V*nc0.variables['NO3'][i,-6:],axis=0) / np.nansum(V,axis=0)
					phyt = np.nansum(V*nc0.variables['phytoplankton'][i,-6:],axis=0) / np.nansum(V,axis=0)
				else:
					nutr = nc0.variables['NO3'][i,-1]
					phyt = nc0.variables['phytoplankton'][i,-1]
				koeff = koeff * nutr / (1.6 + nutr)						# !!! CHANGE HALF-SAT CONSTANT IF NECESSARY !!!
				chl_sp0 += phyt * 62 * koeff * scale_light
				k0 += 1

		chl_sp0 = chl_sp0/k0
		chl_sp0[chl_sp0>40] == np.nan

		chl_roms_arr[s] += chl_sp0
	nc0.close()
	print(y)



for s in range(3):
	chl_cm0 = np.zeros((rr0.variables['CHL'][0].shape))

	k = 0
	for i in range(12):
		if (i+1) in season[s]:
			chl_cm0 += rr0.variables['CHL'][i]
			k += 1

	chl_cm0 = chl_cm0/k
	chl_cm0[chl_cm0<0] = np.nan

	fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 10))

	m2 = Basemap(projection='merc',llcrnrlat=48.0,urcrnrlat=59.6,llcrnrlon=-3,urcrnrlon=10.8,lat_0=3,resolution='i', ax=ax[0])
	chl_cmems = chl_cm0
	m2.drawcoastlines()
	x2, y2 = m2(lon0, lat0)
	chl_cmems[chl_cmems>bounds[-1]]=bounds[-1]-0.01
	chl_cmems[chl_cmems<bounds[0]]=bounds[0]+0.01
	CS2 = m2.contourf(x2,y2,chl_cmems,clevs=bounds,cmap=cmap2, norm=norm, alpha=0.6)

	m1 = Basemap(projection='merc',llcrnrlat=48.0,urcrnrlat=59.6,llcrnrlon=-3,urcrnrlon=10.8,lat_0=3,resolution='i', ax=ax[1])
	cax = make_axes_locatable(ax[1]).append_axes("right", size=0.4, pad=0.15)
	chl_roms = chl_roms_arr[s] / len(year)
	m1.drawcoastlines()
	x2, y2 = m1(lons0, lats0)
	chl_roms[mask0==0]=np.nan
	chl_roms[chl_roms>bounds[-1]]=bounds[-1]-0.01
	chl_roms[chl_roms<bounds[0]]=bounds[0]+0.01
	CS2 = m1.contourf(x2,y2,chl_roms,clevs=bounds,cmap=cmap2, norm=norm, alpha=0.6)
	cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap2, norm=norm, boundaries = bounds, orientation='vertical')

	fig.savefig(path_fig + "/Map_CHL_%s.png" %(season_list[s]), dpi=300, bbox_inches='tight')
