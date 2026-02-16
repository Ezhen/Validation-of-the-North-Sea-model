import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from scipy import spatial

# The script plots simulation results for carbon remineralization in the upper 10cm of the marine sediments.
# Against field campaign results from De Borger 2021.

moment = 4  # once per month sample -> May sample has index 4 (normally 145,157 Julian range)

flag = 'atm3'         #'_extremeWilson
down = 25

# Data from De Borger's rapid remineralization 2021; F770 is outside of the range
station = ['O100','O135','0190','D240','C300','C380','C450','C545','F640','F695','F770']
lat = [54.149,54.419,54.875,55.173,55.623,56.069,56.587,57.360,58.200,58.839,59.415]
lon = [4.342,4.047,3.686,3.161,2.384,1.599,0.685,0.579,0.525,0.511,0.509]
depth = [45,43,41,26,72,79,224,84,148,135,135]
porosity = [0.57,0.46,0.45,0.39,0.48,0.43,0.57,0.46,0.71,0.49,0.61]
Oxic_Emil = np.array([7.02,12.14,10.66,6.42,4.61,6.95,9.01,2.28,6.21,6.36,4.35])
Denitr_Emil = np.array([0.29,0.1,0.04,0.25,0.18,0.48,0.18,0.32,1.06,0.46,0.54])
Anoxic_Emil = np.array([2.0,1.19,2.83,1.44,3.03,3.21,3.38,1.37,1.21,1.29,1.65])

# Acquire base data from ROMS simulations
nc = Dataset('/scratch/ulg/mast/eivanov/Output/CE2COAST_2006/Hindcast_CE2COAST_HIS_2010_2c_%s.nc' %(flag), 'r', format='NETCDF4')
nc2 = Dataset('/scratch/ulg/mast/eivanov/Output/CE2COAST_2006/Hindcast_CE2COAST_RST_2010_2c_%s.nc' %(flag), 'r', format='NETCDF4')
lats = nc.variables['lat_rho'][:]
lons = nc.variables['lon_rho'][:]
mask = nc.variables['mask_rho'][:]
aa=np.array((list(lons.flatten()), list(lats.flatten()))).T; width2=len(nc.variables['h'][:,:].T)
lt,ln,h = [],[],[]
nc.close()
nc2.close()


# years included into simulation span
years = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2021,2022,2023]
Oxic_Roms,Denitr_Roms,Anoxic_Roms = np.zeros((len(years),10)),np.zeros((len(years),10)),np.zeros((len(years),10))

for y in range(len(years)):
	nc = Dataset('/scratch/ulg/mast/eivanov/Output/CE2COAST_2006/Hindcast_CE2COAST_HIS_%s_2c_%s.nc' %(years[y],flag), 'r', format='NETCDF4')
	for i in range(len(lat)-1):
		# standard routine to locate the closest ROMS routine
		idxx,a,b = 0,0,-1
		while mask[a,b]==0:
			aa[idxx][0],aa[idxx][1]=0,0
			idxx=spatial.KDTree(aa).query([lon[i],lat[i]])[1]
			a=int(idxx/width2)
			b=int(idxx - a*width2)

		print('Station',station[i])
		#print(np.round(lats[a,b],2),np.round(lat[i],2),np.round(lons[a,b],2),np.round(lon[i],2),int(depth[i]),int(nc.variables['h'][a,b]),a,b)
		lt.append(float(lats[a,b]))
		ln.append(float(lons[a,b]))
		h.append(int(nc.variables['h'][a,b])*-1)


		denit = np.sum(nc.variables['carrem_04'][moment,0:down,a,b],axis=0)*86400
		anoxic = np.sum(nc.variables['carrem_06'][moment,0:down,a,b],axis=0)*86400
		oxic = np.sum(nc.variables['carrem_01'][moment,0:down,a,b]+nc.variables['carrem_02'][moment,0:down,a,b],axis=0)*86400 - denit - anoxic
		
		Oxic_Roms[y,i] = oxic
		Denitr_Roms[y,i] = denit
		Anoxic_Roms[y,i] = anoxic
		#print('Oxic:   ', 'Emil', Oxic_Emil[i], 'ROMS', Oxic_Roms[-1])
		#print('Denitr: ', 'Emil', Denitr_Emil[i], 'ROMS', Denitr_Roms[-1])
		#print('Anoxic: ', 'Emil', Anoxic_Emil[i], 'ROMS', Anoxic_Roms[-1])
		print('\n')
	nc.close()



Oxic_Roms = np.nanmean(Oxic_Roms,axis=0)
Denitr_Roms = np.nanmean(Denitr_Roms,axis=0)
Anoxic_Roms = np.nanmean(Anoxic_Roms,axis=0)

Oxic_yerr = np.nanstd(Oxic_Roms,axis=0)
Anoxic_yerr = np.nanstd(Anoxic_Roms,axis=0)

# Carbon remineralization
fig, ax = plt.subplots()

ind = np.arange(len(lat)-1)*1.3     # the x locations for the groups
width = 0.5         # the width of the bars
ax.bar(ind,  Oxic_Emil[0:-1], width*0.65, yerr = Oxic_yerr, bottom=0, color = 'lightblue',  label='Oxic')
ax.bar(ind, Denitr_Emil[0:-1], width*0.65, bottom = Oxic_Emil[0:-1], color = 'g',label='Denitrification')
ax.bar(ind, Anoxic_Emil[0:-1], width*0.65, yerr = Anoxic_yerr, bottom = Oxic_Emil[0:-1]+Denitr_Emil[0:-1], color='brown', label='Anoxic')

ax.bar(ind + width*0.7,  Oxic_Roms, width*0.65, bottom=0, color = 'lightblue')
ax.bar(ind + width*0.7, Denitr_Roms, width*0.65, bottom = Oxic_Roms, color = 'g')
ax.bar(ind + width*0.7, Anoxic_Roms,  width*0.65, bottom = Oxic_Roms+Denitr_Roms, color='brown')

ax.set_title('Carbon remineralization (De Borger, 2021 .vs. ROMS, 2023)')
ax.set_xticks(ind + width*0.7)
ax.set_xticklabels(station[0:-1])

ax.legend(loc='best')
ax.set_ylabel('mmol C m-2 d-1')
ax.autoscale_view()

fig.savefig("/home/users/e/i/eivanov/Validation/Climatological_range/Figures/Carbon_remineralization_%s_%s.png"  %(flag,down), dpi=300, bbox_inches='tight')


# Carbon remineralization (%)
fig, ax = plt.subplots()

Emil_TOC = (Oxic_Emil[0:-1] + Denitr_Emil[0:-1] + Anoxic_Emil[0:-1])
Roms_TOC = (Oxic_Roms + Denitr_Roms + Anoxic_Roms)

ind = np.arange(len(lat)-1)*1.3     # the x locations for the groups
width = 0.5         # the width of the bars
ax.bar(ind,  100*Oxic_Emil[0:-1]/Emil_TOC,width*0.65, bottom=0, color = 'lightblue',  label='Oxic')
ax.bar(ind, 100*Denitr_Emil[0:-1]/Emil_TOC, width*0.65, bottom = 100*Oxic_Emil[0:-1]/Emil_TOC, color = 'g',label='Denitrification')
ax.bar(ind, 100*Anoxic_Emil[0:-1]/Emil_TOC,  width*0.65, bottom = 100*(Oxic_Emil[0:-1]+Denitr_Emil[0:-1])/Emil_TOC, color='brown', label='Anoxic')

ax.bar(ind + width*0.7,  100*Oxic_Roms/Roms_TOC, width*0.65, bottom=0, color = 'lightblue') 
ax.bar(ind + width*0.7, 100*Denitr_Roms/Roms_TOC, width*0.65, bottom = 100*Oxic_Roms/Roms_TOC, color = 'g')
ax.bar(ind + width*0.7, 100*Anoxic_Roms/Roms_TOC,  width*0.65, bottom = 100*(Oxic_Roms+Denitr_Roms)/Roms_TOC, color='brown')

ax.set_title('Carbon remineralization % (De Borger, 2021 .vs. ROMS, 2023)')
ax.set_xticks(ind + width*0.7)
ax.set_xticklabels(station[0:-1])

ax.legend(loc='best')
ax.set_ylabel('%')
ax.autoscale_view()

fig.savefig("/home/users/e/i/eivanov/Validation/Climatological_range/Figures/Carbon_remineralization_percent_%s_%s.png" %(flag,down), dpi=300, bbox_inches='tight')


