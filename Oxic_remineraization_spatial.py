import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from scipy import spatial
from datetime import date, datetime, timedelta
import sys
sys.path.extend(['/home/ulg/mast/eivanov/Validation/SPM'])
import os
path_spm = os.environ.get('PATHSPM')
path = os.environ.get('PATH_OUTPUT')
path_fig = os.environ.get('PATH_FIG')
from plotbcz import *
import cmocean.cm as cm

# may be phyt is not bad in spring - we just captured diffrerent snippets of a situation
# Sdet 0.35 -> 0.3, nitrif 0.6 -> 0.5
# two variants of carbon remineralization map

station_1 = ['O100','O135','0190','D240','C300','C380','C450','C545','F640','F695','F770']
lat_1 = [54.149,54.419,54.875,55.173,55.623,56.069,56.587,57.360,58.200,58.839,59.415]
lon_1 = [4.342,4.047,3.686,3.161,2.384,1.599,0.685,0.579,0.525,0.511,0.509]
o2flux_1  = [10.27,17.52,16.03,11.45,11.82,14.11,15.02,8.24,12.15,12.77,12.27]
o2flux_std_1   = [1.32,7.45,4.35,4.25,2.94,1.43,1.08,2.52,4.22,5.41,8.23]

station_2 = ['11','20','30','38','45','52','56','59','62','65','71','80','88']
lat_2 = [53.2,54.4,55,56,57,57.5,58,58,58,58.5,59,60,61]
lon_2 = [2.5,8.1,5,2,5.251,7.5,-0.5,4.25,9.5,9.5,2.5,0.5,3.5]
o2flux_tou = [10.0,22.4,11.33,3.73,7.32,5.49,3.76,4.94,4.8,3.0,3.67,4.62,2.25]
o2_flux_dou = [0.2,2.07,6.95,1.06,1.67,0.83,1.47,4.33,2.49,1.88,1.76,1.25,1.87]

station_3 = ['120','130','330','780','BRN11','D6']
lat_3 = [51.185,51.271,51.433,51.471,51.308,51.552]; lon_3 = [2.702,2.905,2.808,3.058,2.602,2.923]
o2flux_3 = [24,48.8,28.9,78.9,36.6,23.3]; o2flux_std = [1,9.3,2.,3.9,2.4,5.5]

station_4 = ['NOAH-A','NOAH-B','NOAH-C','NOAH-D']
lat_4,lon_4 = [53.99,53.99,54.07,54.091],[6.24,6.87,8.02,7.358]
o2flux_4, o2flux_std = [15.0,25.8,31.9,27.8],[5.8,4.9,7.1,16.2]

station_5 = ['Up1','Up2','Up3','Up4','Up5','Up6','Rys1','Rys2','Rys3']
lat_5,lon_5 = [51.756,53.518,55.5,55.5,54.652,53.513,  58.23,58.08,57.81],[3.0,4.595,6.101,0.908,0.517,3.0,  9.53,10.06,10.16]
o2flux_5 = [8.83,12.55,12.03,10.33,7.76,12.1,  4.2, 12.5, 5.6]

latts = lat_1 #+ lat_2 + lat_3 + lat_4 + lat_5
lonns = lon_1 #+ lon_2 + lon_3 + lon_4 + lon_5
o2fxs = o2flux_1 #+ o2_flux_dou + o2flux_3 + o2flux_4 + o2flux_5

year = 2015
flag = 'atm3'
folder = '12dec_2025'



#ncdata2 = Dataset('/CECI/home/ulg/mast/eivanov/COAWST37_BGC/Projects/CE2COAST_BGC/Speed_up/NS_TEST_1_Jerlov2m.nc','r', format='NETCDF4')
#y_vert,x_vert = ncdata2.variables['lat_vert'][:],ncdata2.variables['lon_vert'][:]
#ncdata2.close()
ncdata2 = Dataset(path+'/Hindcast_CE2COAST_RST_2033_1c_mar_bcorr_owf.nc','r', format='NETCDF4')
z = ncdata2.variables['OWF_number'][:]
z[z>0]=1
lats = ncdata2.variables['lat_rho'][:]
lons = ncdata2.variables['lon_rho'][:]
mask = ncdata2.variables['mask_rho'][:]


logarythm = False

yrange = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2011,2012,2013,2014,2015,   2021,2022,2023]

for y in yrange:
	nc = Dataset(path+'/Hindcast_CE2COAST_HIS_%s_2c_%s.nc' %(y,flag), 'r', format='NETCDF4')
	if y == yrange[0]:
		var = np.zeros(np.shape(mask))

	var += np.nansum(np.nansum(nc.variables['carrem_01'][:,0:]+nc.variables['carrem_02'][:,0:],axis=0),axis=0) * 86400.0 * 30
	nc.close()

var = var / len(yrange)


#plt.imshow(boxcon[6],vmin=0,vmax=10)
#plt.show()

unit = 'mol' #'gram' #'mol' # 'gram'


llcrnrlon1=-3.5; urcrnrlon1=10.5; llcrnrlat1=48.3; urcrnrlat1=59.2
llcrnrlon1=-3.3; urcrnrlon1=10.2; llcrnrlat1=50.5; urcrnrlat1=59.3
fig, ax, cax, m1 = grid_instance(llcrnrlon=llcrnrlon1, urcrnrlon=urcrnrlon1, llcrnrlat=llcrnrlat1, urcrnrlat=urcrnrlat1, lat_ts=51.5, r='i', discr=4, caxx='vertical')
x1, y1 = m1(lons,lats)

cmap = cm.solar_r
if unit == 'mol':
	clevs = [0,1,2,3,4,5,6,7,8,9,10]
elif unit == 'gram':
	clevs = np.arange(10)*32
#clevs = [0,2,4,6,8,10,12,14,16,18,20]


import math

#var = np.nansum(boxcon[:],axis=0)
if logarythm == True:
	boxcon_log = np.copy(var)
	var = np.zeros(np.shape(boxcon_log))
	for i in range(len(var)):
		for j in range(len(var.T)):
			if mask[i,j]==1:
				#var[i,j]=math.log(boxcon_log[i,j],1.78)
				var[i,j]=math.log(boxcon_log[i,j],2.71)
	#var = np.log(boxcon_log)
	#clevs = np.array([0,1,2,3,4,5,6,7,8,9])
	clevs = np.array([0,1,2,3,4,5,6,7])
	print('lol')
else:
	clevs = [0,20,40,60,80,100,120,140]

norm = matplotlib.colors.Normalize(vmin=clevs[0], vmax=clevs[-1])
bounds = clevs
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,  boundaries = bounds, orientation='vertical',norm=norm)

xz,yz = bcz_bound()
xx, yy = m1(xz, yz)
m1.plot(xx,yy,color='black',linewidth=1.0)

#var = boxcon[5]
#if unit == 'gram':
#	var = var*32
var[var>bounds[-1]]=bounds[-1]
var[var<0]=0.01
var[mask==0]=np.nan

var[msk==1]=np.nan

test2 = False
if test2 != True:
	CS1 = m1.contourf(x1,y1,var,clevs=clevs,levels=bounds,cmap=cmap,norm=norm)

scatter = True
if scatter == True:
	for i in range(len(latts)-1):
		x,y = m1(lonns[i],latts[i])
		sc = ax.scatter(x, y, s=200, facecolor='none', edgecolors='k')
		ax.annotate(station_1[i],(x+17000,y), fontsize=14)
		#sc = ax.scatter(x, y, cmap=cmap, norm=norm, s=100, edgecolors='k')


if test2 == True:
	x_vert,y_vert = m1(x_vert,y_vert)
	for i in range(len(mask)):
		for j in range(len(mask.T)):
			if z[i,j] == 1:
				#xs,ys = plot_rec((x_vert[i,j],y_vert[i,j]), (x_vert[i,j+1],y_vert[i,j+1]), (x_vert[i+1,j],y_vert[i+1,j]), (x_vert[i+1,j+1],y_vert[i+1,j+1]))
				m1.plot([x_vert[i,j],x_vert[i,j+1]],[y_vert[i,j],y_vert[i,j+1]], '-k',linewidth=1)
				m1.plot([x_vert[i,j+1],x_vert[i+1,j+1]],[y_vert[i,j+1],y_vert[i+1,j+1]], '-k',linewidth=1)
				m1.plot([x_vert[i+1,j+1],x_vert[i+1,j]],[y_vert[i+1,j+1],y_vert[i+1,j]], '-k',linewidth=1)
				m1.plot([x_vert[i+1,j],x_vert[i,j]],[y_vert[i+1,j],y_vert[i,j]], '-k',linewidth=1)

if logarythm == False and test2 != True:
	fig.savefig(path_fig + "/Oxic_remineralization_spatial_%s_%s.png" %(flag,unit), dpi=500, bbox_inches='tight')
elif test2 == True:
	fig.savefig(path_fig + "/Oxic_remineralization_OWFs2.png", dpi=500, bbox_inches='tight')
else:
	fig.savefig(path_fig + "/Oxic_LOG_remineralization_spatial_%s_%s.png" %(flag,unit), dpi=300, bbox_inches='tight')
