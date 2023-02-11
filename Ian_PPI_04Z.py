#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:06:26 2023

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import pyart
from mpl_toolkits.basemap import Basemap
from tqdm import tqdm


file = 'KTBW_N0Q_20220929_035800.nc'

nc = Dataset(file, 'r')


#Extract attributes

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
z = nc.variables['bref'][0,:,:]

nc.close()


#Make grid for plotting

lat2,lon2 = np.meshgrid(lat,lon)


#Mask out bad Z values


z[z<5] = np.nan


#Plot the data


plt.figure(figsize=(20,20))

cmin = 0; cmax = 75; cint = 2; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-84,-79]); ylim = np.array([26,30])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawstates(); m.drawcountries(); m.drawcoastlines()
m.drawcounties()
cs = m.contourf(lon2,lat2,z.T, clevs, cmap = 'pyart_NWSRef', extend = 'both')

raxpol_lon = -81.94765
raxpol_lat =  28.75541


x2star,y2star = m(raxpol_lon,raxpol_lat)
m.plot(x2star,y2star,'bo',markersize=20)

cbar = plt.colorbar(fraction=0.046)
cbar.ax.tick_params(labelsize = 26)
cbar.set_label(label = '[dBZ]',size = 26)

plt.title(r'KTBW 9/29 0358UTC $0.5^{o} Z_{H}$', size = 40)

#plt.show()


#%%

#Do the same thing with Zdr; Then compute D_o and N_w

file_zdr = 'KTBW_N0X_20220929_035500.nc'


nc_zdr = Dataset(file_zdr, 'r')

zdr = nc_zdr.variables['zdr'][0,:,:] #Shape of (584,800)

#Let's plot Zdr first:

plt.figure(figsize=(20,20))

cmin = -1; cmax = 4; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-84,-79]); ylim = np.array([26,30])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawstates(); m.drawcountries(); m.drawcoastlines()
m.drawcounties()
cs = m.contourf(lon2,lat2,zdr.T, clevs, cmap = 'pyart_RefDiff', extend = 'both')

raxpol_lon = -81.94765
raxpol_lat =  28.75541


x2star,y2star = m(raxpol_lon,raxpol_lat)
m.plot(x2star,y2star,'bo',markersize=20)

cbar = plt.colorbar(fraction=0.046)
cbar.ax.tick_params(labelsize = 26)
cbar.set_label(label = '[dB]',size = 26)

plt.title(r'KTBW 9/29 0358 UTC $0.5^{o} Z_{DR}$', size = 40)

#plt.show()

#Now let's compute D_o
#First, reshape Zdr: 2D --> 1D

I,J = zdr.shape
zdr_reshape = zdr.reshape(I*J, order = 'F') #(467200)


d_o = []

for i in tqdm(range(len(zdr_reshape))):

    if zdr_reshape[i]>=1:


        d_o_compute = 0.0536*(zdr_reshape[i])**3 - 0.197*(zdr_reshape[i])**2 + 0.6261*(zdr_reshape[i]) + 1.0815
        d_o.append(d_o_compute)

    else:

        d_o_compute = 0.0424*(zdr_reshape[i])**4 - 0.4571*(zdr_reshape[i])**3 + 0.6215*(zdr_reshape[i])**2 + 0.457*(zdr_reshape[i]) + 0.8808

        d_o.append(d_o_compute)

#Change list to array


d_o = np.asarray(d_o)


#Alright now let's compute log10(Nw) as well.
#First reshape Z

X,Y = z.shape
z_reshape = z.reshape(X*Y, order = 'F')

log_nw = np.log10(19.76*(z_reshape/(d_o**7.46)))


#Let's mask out values of d_o and log(Nw) for any values of Zdr > 2.5:


d_o[zdr_reshape>2.5] == np.nan
log_nw[zdr_reshape>2.5] == np.nan

d_o[zdr_reshape<-0.5] == np.nan
log_nw[zdr_reshape<-0.5] == np.nan


#Now reshape D_o ^ Nw back into a grid:

log_nw_grid = log_nw.reshape(X,Y, order = 'F')
d_o_grid = d_o.reshape(I,J, order = 'F')

d_o_grid[d_o_grid>3] = np.nan

log_nw_grid[log_nw_grid<1] = np.nan



#Now let's plot?

plt.figure(figsize=(20,20))

cmin = 0; cmax = 4; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-84,-79]); ylim = np.array([26,30])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawstates(); m.drawcountries(); m.drawcoastlines()
m.drawcounties()
cs = m.contourf(lon2,lat2,d_o_grid.T, clevs, cmap = 'pyart_NWSRef', extend = 'both')

raxpol_lon = -81.94765
raxpol_lat =  28.75541


x2star,y2star = m(raxpol_lon,raxpol_lat)
m.plot(x2star,y2star,'bo',markersize=20)

cbar = plt.colorbar(fraction=0.046)
cbar.ax.tick_params(labelsize = 26)
cbar.set_label(label = '[mm]',size = 26)

plt.title(r'KTBW 9/29 0358 UTC $D_{M}$', size = 40)

#plt.show()






plt.figure(figsize=(20,20))

cmin = 1; cmax = 6; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-84,-79]); ylim = np.array([26,30])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawstates(); m.drawcountries(); m.drawcoastlines()
m.drawcounties()
cs = m.contourf(lon2,lat2,log_nw_grid.T, clevs, cmap = 'pyart_NWSRef', extend = 'both')

raxpol_lon = -81.94765
raxpol_lat =  28.75541


x2star,y2star = m(raxpol_lon,raxpol_lat)
m.plot(x2star,y2star,'bo',markersize=20)

cbar = plt.colorbar(fraction=0.046)
cbar.ax.tick_params(labelsize = 26)
cbar.set_label(label = '[$m^{-3}$ $ mm^{-1}$]',size = 26)

plt.title(r'KTBW 9/29 0358 UTC $log_{10}(N_{w})$', size = 40)

plt.show()
