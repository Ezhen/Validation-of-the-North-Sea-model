import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
import os
path_spm = os.environ.get('PATHSPM')
path = os.environ.get('PATH_OUTPUT')
sys.path.extend([path_spm])
from plotbcz import *

import matplotlib as mpl
import cmocean.cm as cm

# the script plots mud content over the NS sediments
# i distinct 10% and 30% as the perennial states (below - noncohesive, above - cohesive)
# in the early versions of ROMS modifications there was a distinction for processes 
# based on "cohesiveness" of the sediments. Anyway, we don't need to see the difference
# between 60% and 70%, only to locate the mud fields

year = 2009
file2 = path + '/Hindcast_CE2COAST_RST_2003_2c_atm3.nc' 

tt = Dataset(file2, 'r', format='NETCDF4')

lats,lons,h = tt.variables['lat_rho'][:], tt.variables['lon_rho'][:], tt.variables['h'][:]

av = tt.variables['mudfrac_01'][-1,0]*2*100

station = ['O100','O135','0190','D240','C300','C380','C450','C545','F640','F695']
lat = [54.149,54.419,54.875,55.173,55.623,56.069,56.587,57.360,58.200,58.839]
lon = [4.342,4.047,3.686,3.161,2.384,1.599,0.685,0.579,0.525,0.511]

plot_figure = True
mud_content = False
if plot_figure == True:

	if mud_content == True:
		cmap2 = cm.turbid
		bounds = [0,10,20,30,40,50,60,70,80,90,100]
		bounds = np.arange(0,102,2)
		norm = mpl.colors.BoundaryNorm(bounds, cmap2.N)

		fig = plt.figure(figsize=(12,8)) #10, 8))
		ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
		m2 = Basemap(projection='merc',llcrnrlat=50.5,urcrnrlat=59.6,llcrnrlon=-4,urcrnrlon=10.5,lat_0=3,resolution='i',ax=ax)
		cax = make_axes_locatable(ax).append_axes("right", size=0.4, pad=0.15)
		cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap2, boundaries = bounds, orientation='vertical',norm=norm)
		m2.drawcoastlines()
		x2, y2 = m2(lons, lats)

		varr = av
		varr[varr>bounds[-1]]=bounds[-1]-0.01
		varr[varr<bounds[0]]=bounds[0]+0.01
		CS2 = m2.contourf(x2,y2,varr,clevs=bounds,cmap=cmap2, norm=norm) #, alpha=0.8) #norm=norm, 

		for i in range(len(station)):
			x3, y3 = m2(lon[i],lat[i])
			m2.scatter(x3,y3,marker='*',c='m',s=150)
			a_shift,b_shift = 18000,18000
			ax.annotate('%s' %(station[i]),xy=(x3+a_shift, y3+b_shift), xycoords='data', xytext=(x3+a_shift, y3+b_shift), textcoords='data', fontsize=12)

		fig.savefig("./Figures/Mud_content.png", dpi=300, bbox_inches='tight')
	
	type_bed = True
	if type_bed == True:
		cmap2 = mpl.cm.inferno_r
		bounds = [0,1,2,3]
		av[av<10]=0.5
		av[av>30]=2.5
		av[av>2.5]=1.5
		norm = mpl.colors.BoundaryNorm(bounds, cmap2.N)
		fig = plt.figure(figsize=(12,8)) #10, 8))
		ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
		m2 = Basemap(projection='merc',llcrnrlat=51.0,urcrnrlat=59.2,llcrnrlon=-4,urcrnrlon=10.0,lat_0=3,resolution='i',ax=ax)
		cax = make_axes_locatable(ax).append_axes("right", size=0.4, pad=0.15)
		cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap2, boundaries = bounds, orientation='vertical',norm=norm)
		m2.drawcoastlines()
		x2, y2 = m2(lons, lats)
		varr = av
		varr[varr>bounds[-1]]=bounds[-1]-0.01
		varr[varr<bounds[0]]=bounds[0]+0.01
		CS2 = m2.contourf(x2,y2,varr,clevs=bounds,cmap=cmap2, norm=norm) #, alpha=0.8) #norm=norm, 
		CS3 = m2.contour(x2,y2,h,levels=[50,100,150],colors='blue') #, alpha=0.8) #norm=norm, 
		#plt.clabel(CS3, inline = True, fontsize=9, fmt='%d', colors = 'black')
		for i in range(len(station)):
			x3, y3 = m2(lon[i],lat[i])
			m2.scatter(x3,y3,marker='o',c='c',s=150)
			a_shift,b_shift = 18000,18000
			ax.annotate('%s' %(station[i]),xy=(x3+a_shift, y3+b_shift), xycoords='data', xytext=(x3+a_shift, y3+b_shift), textcoords='data', fontsize=10, backgroundcolor='w')
		fig.savefig("./Figures/Bed_type.png", dpi=300, bbox_inches='tight')

