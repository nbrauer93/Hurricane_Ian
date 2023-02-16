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
import glob


files_z = glob.glob('*N0Q*')
files_zdr = glob.glob('*N0X*')
files_rhohv = glob.glob('*N0C*')

#Order times

z_sort_index = np.argsort(files_z)
zdr_sort_index = np.argsort(files_zdr)
rhohv_sort_index = np.argsort(files_rhohv)

z_ordered = [files_z[k] for k in z_sort_index]
zdr_ordered = [files_zdr[k] for k in zdr_sort_index]
rhohv_ordered = [files_rhohv[k] for k in rhohv_sort_index]




for i in tqdm(range(len(z_ordered))):

    #First let's parse the time of each file:

    time_parse = z_ordered[i].split("_")[3] #Prints HHMMSS.nc
    time_parsed = time_parse.split(".")[0] #Only parse/keep the first portion of the string
    print(time_parsed) #Add to titles for figures to include timestamp



    nc_z = Dataset(z_ordered[i], 'r')
    nc_zdr = Dataset(zdr_ordered[i], 'r')
    nc_rhohv = Dataset(rhohv_ordered[i], 'r')



#Extract attributes

    lat = nc_z.variables['lat'][:]
    lon = nc_z.variables['lon'][:]
    z = nc_z.variables['bref'][0,:,:]
    zdr = nc_zdr.variables['zdr'][0,:,:]
    rhohv = nc_rhohv.variables['cc'][0,:,:]




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

    plt.title(r'KTBW 9/29 ' + time_parsed + ' UTC $0.5^{o} Z_{H}$', size = 40)

    #plt.show()





    #Let's plot Zdr first:

    plt.figure(figsize=(20,20))

    cmin = -1; cmax = 4; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
    xlim = np.array([-84,-79]); ylim = np.array([26,30])

    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
    m.drawstates(); m.drawcountries(); m.drawcoastlines()
    m.drawcounties()
    cs = m.contourf(lon2,lat2,zdr.T, clevs, cmap = 'pyart_RefDiff', extend = 'both')




    x2star,y2star = m(raxpol_lon,raxpol_lat)
    m.plot(x2star,y2star,'bo',markersize=20)

    cbar = plt.colorbar(fraction=0.046)
    cbar.ax.tick_params(labelsize = 26)
    cbar.set_label(label = '[dB]',size = 26)

    plt.title(r'KTBW 9/29 ' + time_parsed + ' UTC $0.5^{o} Z_{DR}$', size = 40)

    #plt.show()

    #Now let's compute D_o
    #First, reshape Zdr: 2D --> 1D

    I,J = zdr.shape
    zdr_reshape = zdr.reshape(I*J, order = 'F') #(467200)


    d_o = []

    for j in tqdm(range(len(zdr_reshape))):

        if zdr_reshape[j]>=1:


            d_o_compute = 0.0536*(zdr_reshape[j])**3 - 0.197*(zdr_reshape[j])**2 + 0.6261*(zdr_reshape[j]) + 1.0815
            d_o.append(d_o_compute)

        else:

            d_o_compute = 0.0424*(zdr_reshape[j])**4 - 0.4571*(zdr_reshape[j])**3 + 0.6215*(zdr_reshape[j])**2 + 0.457*(zdr_reshape[j]) + 0.8808

            d_o.append(d_o_compute)

            #Change list to array


    d_o = np.asarray(d_o)


    #Alright now let's compute log10(Nw) as well.
    #First reshape Z

    X,Y = z.shape
    z_reshape = z.reshape(X*Y, order = 'F')
    rhohv_reshape = rhohv.reshape(X*Y, order = 'F')

    log_nw = np.log10(19.76*(z_reshape/(d_o**7.46)))


    #Maks out values where rho_hv < 0.85

    d_o[rhohv_reshape<0.85] = np.nan
    log_nw[rhohv_reshape<0.85] = np.nan


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




    x2star,y2star = m(raxpol_lon,raxpol_lat)
    m.plot(x2star,y2star,'bo',markersize=20)

    cbar = plt.colorbar(fraction=0.046)
    cbar.ax.tick_params(labelsize = 26)
    cbar.set_label(label = '[mm]',size = 26)

    plt.title(r'KTBW 9/29 ' + time_parsed + ' UTC $D_{0}$', size = 40)
    plt.savefig('D_o ' + time_parsed + ' .png')

    #plt.show()






    plt.figure(figsize=(20,20))

    cmin = 1; cmax = 6; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
    xlim = np.array([-84,-79]); ylim = np.array([26,30])

    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
    m.drawstates(); m.drawcountries(); m.drawcoastlines()
    m.drawcounties()
    cs = m.contourf(lon2,lat2,log_nw_grid.T, clevs, cmap = 'pyart_NWSRef', extend = 'both')


    x2star,y2star = m(raxpol_lon,raxpol_lat)
    m.plot(x2star,y2star,'bo',markersize=20)

    cbar = plt.colorbar(fraction=0.046)
    cbar.ax.tick_params(labelsize = 26)
    cbar.set_label(label = '[$m^{-3}$ $ mm^{-1}$]',size = 26)

    plt.title(r'KTBW 9/29 ' + time_parsed + ' UTC $log_{10}(N_{w})$', size = 40)
    plt.savefig('log_nw ' + time_parsed + ' .png')

    #plt.show()
