#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 18:00:16 2023

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import pyart




file = 'RAXPOL-20220929-051136-A154.0-Z.nc'



#Parse out the time from the file name to make shit easier down the road

time = file.split("-")[2]


    #Read in file

nc = Dataset(file,'r')

#Extract attributes:

z = nc.variables['Intensity'][:]
    #zdr = nc.variables['Differential_Reflectivity'][:]
    #rhohv = nc.variables['RhoHV'][:]
#phidp = nc.variables['PhiDP'][:]
#w = nc.variables['Width'][:]
    #vel = nc.variables['Radial_Velocity'][:]
azimuth = nc.variables['Azimuth'][:]
elevation = nc.variables['Elevation'][:]


#Now calculate beam height and distance from the radar using Doviak and Zrnic



deg_to_rad = np.pi/180. #Convert from degrees to radians

range_to_first_cell = 0.
gate_width = nc.variables['GateWidth'][:]

range_gate = np.zeros_like(z[0,:])
range_gate = np.arange(0,len(range_gate))




range_gate_distance = range_gate.copy()
range_gate_distance = []

for j in range(len(range_gate)):

    dist = range_gate[j]*30
    range_gate_distance.append(dist)


range_gate_distance = np.asarray(range_gate_distance)



#Effective radius
a_e = (8494.66667)*1000

#Now compute height/distance from the radar

height = np.ones(z.shape)*np.nan
s = np.ones(z.shape)*np.nan

for k in range(0,len(elevation)):

    a = range_gate**2 + a_e**2 + (2.0*range_gate_distance*a_e*np.sin(elevation[k]*deg_to_rad))
    height[k,:] =((a**0.5)-a_e)+2
    s[k,:] = a_e*np.arcsin((range_gate_distance*np.cos(elevation[k]*deg_to_rad))/(a_e+height[k,:]))

#NaN out bad z values


z[z<0.1] = np.nan
    #zdr[zdr<0.1] = np.nan
    #rhohv[rhohv<0.8] = np.nan
#phidp[phidp<-10] = np.nan
#w[w<0] = np.nan
    #vel[vel<-100] = np.nan
#%%

#Let's rescale to match the DPR footprint (5 km); RaXPol is 30 m. 5km/0.03 km = 166.66667
#So average over every 166.7 range bins in the range dimension
#First we need to rescale the "s/1000", or range in km to match that of the DPR.
#In order to rescale, 2495 has a prime number of 155. Let's average over every 155 bins. 155 *0.03 km = 4.65 km resolution
    
zh_rescale = np.ones((z.shape[0],10))*np.nan
s_rescale = np.ones((z.shape[0],10))*np.nan
height_rescale = np.ones((z.shape[0],10))*np.nan




for i in range(zh_rescale.shape[1]):
    
    #if z.shape[1] = 2495:
    

    avg_interval = 166

    zh_rescale[:,i] = np.nanmean(z[:,(avg_interval*i):(avg_interval*i)+avg_interval], axis = 1)
    s_rescale[:,i] = s[:,avg_interval*i]
    height_rescale[:,i] = height[:,avg_interval*i]
        
        
        

        
  










#Now let's plot les data
#%%
label_size = 24
        





fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s_rescale/1000,height_rescale/1000, zh_rescale, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(12,60,step = 2))
    #plt.contourf(s/1000,height/1000, vel, extend = 'both', cmap = 'pyart_NWSVel', levels = np.arange(-16,16,step = 0.5))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[dBZ]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)

xlabels = np.arange(-50,60,10)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)
plt.xlim(-50,50)
plt.xticks(xlabels, size = label_size)
plt.yticks(ylabels, size = label_size)

plt.title(r'RaXPol $154^{o}$ Azimuth $Z_{H}$ 9/29 ' + time + ' UTC', size = label_size)
plt.savefig(time+ '_Z.png')
plt.show()