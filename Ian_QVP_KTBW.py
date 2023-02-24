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
from scipy import stats

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
#qvps_kdp = np.ones((107,1832))*np.nan

for i in range(len(files_ordered)):

    qvps = quasi_vertical_profile(files_ordered[i])

    qvps_z[i,:] = qvps['reflectivity']
    qvps_zdr[i,:] = qvps['differential_reflectivity']
    qvps_rhohv[i,:] = qvps['cross_correlation_ratio']
    #qvps_kdp[i,:] = pyart.retrieve.kdp_maesaka(qvps['differential_phase'])[0]

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

#Let's include LWC as well: 

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


#%%


#Let's compute the vertical slopes of D_o and N_w similar to what we did with Z and Zdr:
    
    
slope_d_o = []
slope_n_w = []


for i in range(d_o.shape[0]):
    
    d_o_slope = stats.linregress(d_o[i,:],height_warm)[0]
    slope_d_o.append(d_o_slope)
    
    n_w_slope = stats.linregress(log_nw[i,:],height_warm)[0]
    slope_n_w.append(n_w_slope)


#Convert lists to arrays

slope_d_o = np.asarray(slope_d_o)
slope_n_w = np.asarray(slope_n_w)

    
#Now apply a 3 point runnning average

n=10


d_o_avg = moving_avg(slope_d_o,n)    
n_w_avg = moving_avg(slope_n_w,n)    


    
plt.figure(figsize=(10,10))
plt.plot(times[0:98],d_o_avg, color = 'k')
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel('mm/km', size = font_size)
plt.yticks(size= font_size)
plt.title(r'KTBW 9/28 Vertical Slope of $D_{0}$ in Liquid Phase', size = title_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'r', linewidth =  3.0, linestyle = '--')
plt.show()


    

   
plt.figure(figsize=(10,10))
plt.plot(times[0:98],n_w_avg, color = 'k')
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel(r'($m^{-3}$ $mm^{-1}$)/km', size = font_size)
plt.yticks(size= font_size)
plt.title(r'KTBW 9/28 Vertical Slope of $log_{10}(N_{W})$ in Liquid Phase', size = title_size)
plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'r', linewidth =  3.0, linestyle = '--')

plt.show() 
    

#%%

#Let's now compute vertical averages of d_o and n_w


d_o_height_avg = []
nw_height_avg = []

for i in range(d_o.shape[0]):
    
    d_o_avg = np.nanmean(d_o[i,:], axis = 0)
    nw_avg = np.nanmean(log_nw[i,:], axis = 0)
    
    
    d_o_height_avg.append(d_o_avg)
    nw_height_avg.append(nw_avg)

d_o_height_avg = np.asarray(d_o_height_avg)
nw_height_avg = np.asarray(nw_height_avg)


#Apply 10 point running average

d_o_avg_vert = moving_avg(d_o_height_avg,n)    
n_w_avg_vert = moving_avg(nw_height_avg,n) 




    
plt.figure(figsize=(10,10))
plt.plot(times[0:98],n_w_avg_vert, color = 'k')
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel(r'($m^{-3}$ $mm^{-1}$)/km', size = font_size)
plt.yticks(size= font_size)
plt.title(r'KTBW 9/28 Vertical Average of $log_{10}(N_{W})$ in Liquid Phase', size = title_size)
#plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'r', linewidth =  3.0, linestyle = '--')

plt.show()     


plt.figure(figsize=(10,10))
plt.plot(times[0:98],d_o_avg_vert, color = 'k')
x = [0,18,36,54,72,90,108]
labels = np.array(['1200','1400','1600','1800','2000','2200','2400'])
plt.xticks(x,labels,size = font_size)
plt.xlabel('Time (UTC)', size = font_size)
plt.ylabel(r'mm/km', size = font_size)
plt.yticks(size= font_size)
plt.title(r'KTBW 9/28 Vertical Average of $D_{0}$ in Liquid Phase', size = title_size)
#plt.axhline(y = 0, xmin = -1.5, xmax = 4, color = 'r', linewidth =  3.0, linestyle = '--')

plt.show()  


#%%


#Plot on a scatter plot average d_0 versus vertical slope of d_0

text_size = 20


plt.figure(figsize=(10,10))
plt.scatter(d_o_avg_vert,d_o_avg)
plt.plot(np.unique(d_o_avg_vert), np.poly1d(np.polyfit(d_o_avg_vert,d_o_avg , 1))(np.unique(d_o_avg_vert)), color = 'r')
plt.xlabel(r'Mean 1-4 km AGL $D_{0}$ (mm)', size = text_size)
plt.ylabel(r'Slope of $D_{0}$ in Liquid Phase (mm/km)', size = text_size)
plt.xticks(size= text_size)
plt.yticks(size = text_size)
plt.title(r'KTBW $D_{0}$ 1200-2355 UTC 28 September', size = text_size)
plt.show()


#First remove last two digits of times:
    

times_cleaned = []

for i in range(len(times)):
    
    times_junk = times[i][0:4]
    times_cleaned.append(times_junk)


times_cleaned = np.asarray(times_cleaned, dtype = float)


#Apply 10 point running average

times_mean = moving_avg(times_cleaned,n)    




#Let's create a colormap and color each point by time:
    
    
    
colormap = ['black','dodgerblue', 'deepskyblue','lawngreen', 'lightgreen', 'green', 'yellow' ,'gold', 'darkorange', 'red', 'firebrick','maroon']

#Uncomment for raw, without time window mean

'''
time_color = []


for i in range(len(times_cleaned)):
    
    if times_cleaned[i]>1200 and times_cleaned[i]<1300:
        time_color.append(colormap[0])
        
    elif times_cleaned[i]>1300 and times_cleaned[i]<1400:
        time_color.append(colormap[1])
        
    elif times_cleaned[i]>1400 and times_cleaned[i]<=1500:
        time_color.append(colormap[2])
        
    elif times_cleaned[i]>1500 and times_cleaned[i]<1600:
        time_color.append(colormap[3])
        
    elif times_cleaned[i]>1600 and times_cleaned[i]<1700:
        time_color.append(colormap[4])
        
    elif times_cleaned[i]>1700 and times_cleaned[i]<1800:
        time_color.append(colormap[5])
        
    elif times_cleaned[i]>1800 and times_cleaned[i]<1900:
        time_color.append(colormap[6])
        
    elif times_cleaned[i]>1900 and times_cleaned[i]<2000:
        time_color.append(colormap[7])
        
    elif times_cleaned[i]>2000 and times_cleaned[i]<2100:
        time_color.append(colormap[8])
        
    elif times_cleaned[i]>2100 and times_cleaned[i]<2200:
        time_color.append(colormap[9])
        
    elif times_cleaned[i]>2200 and times_cleaned[i]<2300:
        time_color.append(colormap[10])
        
    elif times_cleaned[i]>2300 and times_cleaned[i]<2400:
        time_color.append(colormap[11])
'''       


#With 10 point running average: 

time_color = []


for i in range(len(times_mean)):
    
    if times_mean[i]>1200 and times_mean[i]<1300:
        time_color.append(colormap[0])
        
    elif times_mean[i]>=1300 and times_mean[i]<1400:
        time_color.append(colormap[1])
        
    elif times_mean[i]>1400 and times_mean[i]<=1500:
        time_color.append(colormap[2])
        
    elif times_mean[i]>1500 and times_mean[i]<1600:
        time_color.append(colormap[3])
        
    elif times_mean[i]>1600 and times_mean[i]<1700:
        time_color.append(colormap[4])
        
    elif times_mean[i]>1700 and times_mean[i]<1800:
        time_color.append(colormap[5])
        
    elif times_mean[i]>1800 and times_mean[i]<1900:
        time_color.append(colormap[6])
        
    elif times_mean[i]>1900 and times_mean[i]<2000:
        time_color.append(colormap[7])
        
    elif times_mean[i]>2000 and times_mean[i]<2100:
        time_color.append(colormap[8])
        
    elif times_mean[i]>2100 and times_mean[i]<2200:
        time_color.append(colormap[9])
        
    elif times_mean[i]>2200 and times_mean[i]<2300:
        time_color.append(colormap[10])
        
    elif times_mean[i]>2300 and times_mean[i]<2400:
        time_color.append(colormap[11])


from matplotlib.colors import ListedColormap
cmap_time = ListedColormap(colormap)    
    
    

    
#Setup plotting 
cmin = 1200; cmax = 2400; cint = 100; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_time,lut=nlevs)

import matplotlib

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
#scalar_mappable_object.set_array(time_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(d_o_avg_vert,d_o_avg, c = time_color, vmin = 1200, vmax = 2400, cmap = cmap, edgecolors = 'none')


    
plt.xlabel(r'Mean 1-4 km AGL $D_{0}$ (mm)', size = text_size)
plt.ylabel(r'Slope of $D_{0}$ in Liquid Phase (mm/km)', size = text_size)
plt.xticks(size = text_size)
plt.yticks(size = text_size)
plt.title(r'KTBW $D_{0}$ 1200-2355 UTC 28 September', size = text_size)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label('Time (UTC)', size = title_size)
color_bar_object.ax.tick_params(labelsize = text_size)

plt.show()
     
        
    
    


    
    


