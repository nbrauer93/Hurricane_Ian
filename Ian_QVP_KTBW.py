#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 20:15:25 2023

@author: noahbrauer
"""
import matplotlib.pyplot as plt
import numpy as np
import pyart
import cartopy.crs as ccrs
import glob
#from tqdm import tqdm

radar_file = 'KTBW20220928_215445_V06'


#import numpy as np

#from pyart.io import read
#from pyart.core import antenna_to_cartesian



def quasi_vertical_profile(filename, fields=None, gatefilter=None):

     radar = pyart.io.read(filename)

     if fields is None:
         fields = radar.fields

    #if gatefilter is not None:


     desired_angle = 20.0
     index = abs(radar.fixed_angle['data'] - desired_angle).argmin()
     #print(radar.fixed_angle['data'])
     #print(radar.elevation['data'][-1])

     qvp = {}

     for field in fields:
         this_field = radar.get_field(index, field).mean(axis = 0)
         qvp.update({field:this_field})

     qvp.update({'range': radar.range['data'], 'time': radar.time})
     x,y,z = pyart.core.antenna_to_cartesian(qvp['range']/1000.0, 0.0,
                                 radar.fixed_angle['data'][index])
     qvp.update({'height': z})
     #del radar
     return qvp

qvp = quasi_vertical_profile(radar_file)
#print(qvp['reflectivity'])

plt.plot(qvp['differential_reflectivity'],qvp['height']/1000.0)
plt.xlabel('Mean Differential Reflectivity (dB)')
plt.ylabel('Altitude')
plt.title('Quasi-Vertical Profile 2154 UTC')
#plt.show()


#Now do this for multiple files



files = glob.glob('KTBW*')

#Loop through file list to get times; Parse out of file radar_name
#Sort file names by time first:

order_index = np.argsort(files)
files_ordered = [files[i] for i in order_index]

times = []
qvps_z = np.ones((107,1832))*np.nan
qvps_zdr = np.ones((107,1832))*np.nan
qvps_rhohv = np.ones((107,1832))*np.nan

for i in range(len(files_ordered)):

    qvps = quasi_vertical_profile(files_ordered[i])

    qvps_z[i,:] = qvps['reflectivity']
    qvps_zdr[i,:] = qvps['differential_reflectivity']
    qvps_rhohv[i,:] = qvps['cross_correlation_ratio']

    parsed_time = files_ordered[i].split("_")[1]
    print(parsed_time)
    times.append(parsed_time)

#print(qvps_z.shape) # 107 x 1832
#print(np.nanmax(qvps_z))  # 44.495 dBZ

times = np.asarray(times)



#%%

qvps_z[qvps_z<5] = np.nan
qvps_zdr[qvps_zdr<=0] = np.nan
qvps_rhohv[qvps_rhohv<0.7] = np.nan
qvps_rhohv[qvps_rhohv>1] = np.nan


#Now plot:

font_size = 16    
title_size = 20

plt.figure(figsize=(10,10))
    
qvps_z[qvps_z<5] = np.nan
qvps_zdr[qvps_zdr<=0] = np.nan

plt.pcolormesh(times,qvp['height']/1000.0,qvps_z.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size= font_size)
plt.ylabel('Height (km)', size = font_size)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[dBZ]',size = font_size)
#plt.colorbar()
plt.ylim(0,12)

plt.title(r'KTBW Quasi-Vertical Profile 9/28 $Z_{H}$', size = title_size)
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()


plt.figure(figsize=(10,10))

plt.pcolormesh(times,qvp['height']/1000.0,qvps_zdr.T, cmap = 'pyart_RefDiff')
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('Height (km)', size = font_size)
plt.ylim(0,12)
plt.clim(0,3)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[dB]',size = font_size)
plt.title(r'KTBW Quasi-Vertical Profile 9/28 $Z_{DR}$', size = title_size)

plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()

plt.figure(figsize=(10,10))

plt.pcolormesh(times,qvp['height']/1000.0,qvps_rhohv.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('Height (km)', size = font_size)
plt.ylim(0,12)
plt.clim(0.9,1.)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = r'[$\rho_{HV}$]',size = font_size)
plt.title(r'KTBW Quasi-Vertical Profile 9/28 $\rho_{HV}$', size = title_size)

plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()


#%%

height = qvps['height']/1000
height_warm = height[4:39]
height_ice = height[64:87]


#Now compute the vertical slopes of Z and Zdr from the 1-4 km layer and the 6-8 km layer 

#Index to reflect each layer

z_warm = qvps_z[:,4:39]    
zdr_warm = qvps_zdr[:,4:39]  

z_ice = qvps_z[:,64:87]
zdr_ice = qvps_zdr[:,64:87]


slope_z_warm = []
slope_zdr_warm = []

slope_z_ice = []
slope_zdr_ice = []


#Compute slope via linear regression: Loop through all times

for i in range(z_warm.shape[0]):
    
    z_slope_warm = stats.linregress(z_warm[i,:],height_warm)[0]
    slope_z_warm.append(z_slope_warm)
    
    zdr_slope_warm = stats.linregress(zdr_warm[i,:],height_warm)[0]
    slope_zdr_warm.append(zdr_slope_warm)

    z_slope_ice = stats.linregress(z_ice[i,:],height_ice)[0]
    slope_z_ice.append(z_slope_ice)

    zdr_slope_ice = stats.linregress(zdr_ice[i,:],height_ice)[0]
    slope_zdr_ice.append(zdr_slope_ice)



#Change lists to arrays


slope_z_warm = np.asarray(slope_z_warm)
slope_zdr_warm = np.asarray(slope_zdr_warm)

slope_z_ice = np.asarray(slope_z_ice)
slope_zdr_ice = np.asarray(slope_zdr_ice)


#%%


#Now let's plot as time series
#Compute 3 point moving average to smooth time series

def moving_avg(x, n):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[n:] - cumsum[:-n]) / float(n)

n = 3

z_warm_avg = moving_avg(slope_z_warm,n)
z_ice_avg = moving_avg(slope_z_ice,n)

zdr_warm_avg = moving_avg(slope_zdr_warm,n)
zdr_ice_avg = moving_avg(slope_zdr_ice,n)







plt.figure(figsize=(10,10))
plt.plot(times[0:105],z_warm_avg, color = 'k', label = 'Slope in Liquid Phase')
plt.plot(times[0:105],z_ice_avg, color = 'b', label = 'Slope in Ice Phase')
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('dBZ/km', size = font_size)
plt.yticks(size= font_size)
plt.title(r'KTBW 9/28 Vertical Slope of $Z_{H}$ in Liquid Phase', size = title_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'r', linewidth =  3.0, linestyle = '--')
plt.legend()
plt.show()


    

   
plt.figure(figsize=(10,10))
plt.plot(times[0:105],zdr_warm_avg, color = 'k', label = 'Slope in Liquid Phase')
plt.plot(times[0:105],zdr_ice_avg, color = 'b', label = 'Slope in Ice Phase')
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('dB/km', size = font_size)
plt.yticks(size= font_size)
plt.title(r'KTBW 9/28 Vertical Slope of $Z_{DR}$ in Liquid Phase', size = title_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'r', linewidth =  3.0, linestyle = '--')
plt.legend()
plt.show() 



#Now let's estimate D_o and N_w in the liquid phase from z_warm and zdr_warm using DeHart and Bell (2021)

#z_warm has shape 107 x 35 (Time x Height)

d_o = np.ones((z_warm.shape))*np.nan


for j in range(z_warm.shape[0]):
    for k in range(z_warm.shape[1]):
        
        if zdr_warm[j,k]>=1:
            
            d_o[j,:] = 0.0536*(zdr_warm[j,:])**3 - 0.197*(zdr_warm[j,:])**2 + 0.6261*(zdr_warm[j,:]) + 1.0815
            
        else:
            
            d_o[j,:] = 0.0424*(zdr_warm[j,:])**4 - 0.4571*(zdr_warm[j,:])**3 + 0.6215*(zdr_warm[j,:])**2 + 0.457*(zdr_warm[j,:]) + 0.8808
            
            
        
        
#Now compute log10(N_w)

log_nw = np.log10(19.76*(z_warm/(d_o**7.46)))        
        
#%%    
    
#Now Plot


plt.figure(figsize=(10,10)) 
plt.pcolormesh(times,height_warm,d_o.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size= font_size)
plt.ylabel('Height (km)', size = font_size)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = '[mm]',size = font_size)


plt.title(r'KTBW Quasi-Vertical Profile 9/28 $D_{0}$', size = title_size)
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()


plt.figure(figsize=(10,10))
plt.pcolormesh(times,height_warm,log_nw.T, cmap = 'turbo')
plt.xlabel('Time (UTC)', size= font_size)
plt.ylabel('Height (km)', size = font_size)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = font_size)
cbar.set_label(label = r'[$m^{-3} $ $mm^{-1}$]',size = font_size)


plt.title(r'KTBW Quasi-Vertical Profile 9/28 $log_{10}(N_{W}$)', size = title_size)
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.yticks(size = font_size)
plt.show()




