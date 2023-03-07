#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  2 11:00:22 2022

@author: weihe
"""

open
#from pylab import *
import matplotlib.pyplot as plt
import glob
import netCDF4 as nc

import numpy as np
import datetime
#import matplotlib.ticker as ticker
from mpl_toolkits.basemap import maskoceans
#from mpl_toolkits.basemap import Basemap, shiftgrid,cm,maskoceans
from numpy import array
import os

from time import gmtime, strftime
import datetime as dt
import pandas as pd
from pandas import *
import matplotlib.dates as pltdt
import math
#import gdal
 
from scipy import signal

  


#################################################################
#   Plot for Figure 7:  regional anomalies 
#################################################################
 
    
'''    


AIMs (6) +  TBMs (9)
#2012:
[-0.09097514, -0.09631084, -0.05580726,  0.01444179, -0.06062256,
       -0.00916703]

[ -0.26348409, -0.14570927, -0.269025  , -0.09278803,
       -0.13071877,  0.04588336,  0.36026151, -0.03491898, -0.15559823]

#2015:
[-0.00166041,  0.01481946,  0.16060443,  0.07988457, -0.01066147,
       -0.01550863]

[ -0.17601254, -0.12713907, -0.03170502, -0.07483585,
       -0.00227846,  0.0263374 ,  0.19384429,  0.01836523, -0.09185225]



EUROCOM (7+mean)  +  satellite (Byrne2020, GCASv2, CCDAS)

#  2012
[-0.14973553, -0.16956765, -0.08046738, -0.15973426, -0.12046228,
       -0.13655968, -0.12292061, -0.14973553,
       -0.07526881, -0.11830681, -0.04695053]

[-0.14973553, -0.16956765, -0.08046738, -0.15973426, -0.12046228,
       -0.13655968, -0.12292061, 
       -0.07526881, -0.11830681, -0.04695053]


#  2015

[ 0.03969857,  0.02833359, -0.19756992, -0.12455267, -0.02905707,
        0.16666511,  0.01961452,  0.03969857,
        -0.0533832 , -0.09357872, -0.1669913 ]

[ 0.03969857,  0.02833359, -0.19756992, -0.12455267, -0.02905707,
        0.16666511,  0.01961452, \
        -0.0533832 , -0.09357872,-0.1669913 ]


'''

TBMs_2012 = [ -0.26348409, -0.14570927, -0.269025  , -0.09278803, \
       -0.13071877,  0.04588336,  0.36026151, -0.03491898, -0.15559823]

TBMs_2015 = [ -0.17601254, -0.12713907, -0.03170502, -0.07483585, \
       -0.00227846,  0.0263374 ,  0.19384429,  0.01836523, -0.09185225]



mean_2012 = np.nanmedian(TBMs_2012)    #-0.07623306   # -0.13071877
std_2012 = np.nanstd(TBMs_2012)      #0.18066144


mean_2015 = np.nanmedian(TBMs_2015)    #-0.02947514  # -0.03170502
std_2015 = np.nanstd(TBMs_2015)      #0.10164201




all_2012 = array([-0.09097514, -0.09631084, -0.05580726,  0.01444179, -0.06062256, -0.00916703, \
       -0.14973553, -0.16956765, -0.08046738, -0.15973426, -0.12046228,-0.13655968, -0.12292061, \
       -0.07526881, -0.06862784, -0.11830681, -0.04695053, \
       -0.13071877,\
       -0.04381018, -0.03289082] )  


all_2015 = array([ -0.00166041,  0.01481946,  0.16060443,  0.07988457, -0.01066147,  -0.01550863, \
        0.03969857,  0.02833359, -0.19756992, -0.12455267, -0.02905707, 0.16666511,  0.01961452, \
        -0.0533832 , -0.038638,-0.09357872,-0.1669913, \
        -0.03170502,\
        -0.02239878, -0.01200608])
        
 

all_2012_err = array([-0.0, -0.0, -0.0,  -0.0, -0.0, -0.0, \
        -0.0, -0.0, -0.0,  -0.0, -0.0, -0.0,  -0.0, \
        -0.0, -0.0,-0.0, -0.0, \
        0.18066144, \
        -0.0, -0.0] )  

all_2015_err = array([-0.0, -0.0, -0.0,  -0.0, -0.0, -0.0, \
        -0.0, -0.0, -0.0,  -0.0, -0.0, -0.0,  -0.0, \
        -0.0, -0.0, -0.0, -0.0, \
        0.10164201, \
        -0.0, -0.0] ) 
       
############################################ 

N = len(all_2012)
ind = np.arange(N)

#ind = np.arange(1)   # the x locations for the groups
width = 0.5*2  #5          # the width of the bars

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 6
fig_size[1] = 4.5
plt.rcParams["figure.figsize"] = fig_size

#colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["r", "limegreen", "orange", "royalblue", "gold", "magenta", "cyan", "k"]
          
colors = ["lightgreen", "skyblue", "r", "khaki", "gold", "magenta", "cyan", "k"]
      


fig = plt.figure(figsize=(8+1,4*2))
line =2
row=1

fz=18  #20
######################################################
# subplot R1
#from brokenaxes import brokenaxes


ax1=fig.add_subplot(line,row,1)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 

#
#data=np.zeros(shape=(1,17))
#data[0,:] = all_2012
#data[1,:] = data2[:8]

#for i in np.arange(0,17):
#   # clr= colors[i]
#    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
#    rects.append(ax1.bar( (i+1)*width, all_2012[i]*1e3, width,  hatch="//\\") )
#        
#        #color='white',ecolor=clr , edgecolor=clr , 
#    #clr= colors[5]
#

#rects.append(ax1.bar( ind[:6],  all_2012[:6]*1e3, width=width, color='white', ecolor='green', edgecolor='green',  hatch="//\\") )
#rects.append(ax1.bar( ind[6:13],  all_2012[6:13]*1e3, width=width, color='white', ecolor='orange', edgecolor='orange',  hatch="//\\") )
#rects.append(ax1.bar( ind[13:16],  all_2012[13:16]*1e3, width=width,color='white', edgecolor='red',    hatch="//\\") )
#rects.append(ax1.bar( ind[16],  all_2012[16]*1e3, yerr=all_2012_err[16]*1e3,width=width, color='white',ecolor='blue',edgecolor='blue',  hatch="\\\\") )
#rects.append(ax1.bar( ind[17:19],  all_2012[17:19]*1e3,width=width, color='white',ecolor='magenta',edgecolor='magenta',  hatch="//") )
 

rects.append(ax1.bar( ind[:6],  all_2012[:6]*1e3, width=width, color='lightgreen', ecolor='lightgreen', edgecolor='white') )
rects.append(ax1.bar( ind[6:13],  all_2012[6:13]*1e3, width=width, color='skyblue', ecolor='skyblue', edgecolor='white') )
rects.append(ax1.bar( ind[13:17],  all_2012[13:17]*1e3, width=width,color='red', edgecolor='white') )
rects.append(ax1.bar( ind[17],  all_2012[17]*1e3, yerr=all_2012_err[17]*1e3,width=width, color='khaki',ecolor='khaki',edgecolor='white') )
rects.append(ax1.bar( ind[18:20],  all_2012[18:20]*1e3,width=width, color='indianred',ecolor='indianred',edgecolor='white') )
      
            
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))

#ax.grid()
ax1.yaxis.grid(color='lightgray', linestyle='--')

ax1.set_ylabel('\u0394NEP [TgC yr$^{-1}$]',  fontsize=fz*0.9)
 
ax1.set_xticks(ind)
#ax.set_ylim([-170,0])
ax1.set_ylim([-200,200])
ax1.locator_params(axis='y',nticks=2)
 
ax1.set_xticklabels( ['CTE2018','CT2019B', 'CAMS_v18r2','Jena99_v43','JenaNEET_v43','JenaNEETW_v43', 'CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB', 'Byrne2020','CMS-Flux2020','GCASv2','CCDAS_SM+VOD', 'TRENDY_v6', 'FLUXCOM_RS+CRUJRA', 'FLUXCOM_RS+ERA5'],rotation=90, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')

#lgnames.append("CTE")
#lgnames.append("EnKF-RAMS")
#lgnames.append("FLEXINVERT")
#lgnames.append("CarboScope-Regional")
#lgnames.append("LUMIA")
#lgnames.append("PYVAR-CHIMERE")
#lgnames.append("NAME-HB")
#lgnames.append("Mean")

#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax1.legend( rects, lgnames, loc = 2, ncol=2, frameon=False, fontsize = fz*0.8)
plt.tick_params(labelsize=fz*0.9)  





ax2=fig.add_subplot(line,row,2)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 

#
#data=np.zeros(shape=(1,17))
#data[0,:] = all_2012
#data[1,:] = data2[:8]

#for i in np.arange(0,17):
#    #clr= colors[i]
#    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
#    rects.append(ax2.bar(  (i+1)*width, all_2015[i]*1e3, width,   hatch="//\\") )
#        #color='white', ecolor=clr , edgecolor=clr ,
#        
#    #clr= colors[5]

N = len(all_2015)
ind = np.arange(N)   # the x locations for the groups

 
#rects.append(ax2.bar( ind[:6],  all_2015[:6]*1e3, width=width, color='white', ecolor='green', edgecolor='green',  hatch="//\\") )
#rects.append(ax2.bar( ind[6:13],  all_2015[6:13]*1e3, width=width, color='white', ecolor='orange', edgecolor='orange',  hatch="//\\") )
#rects.append(ax2.bar( ind[13:16],  all_2015[13:16]*1e3, width=width,color='white', edgecolor='red',    hatch="//\\") )
#rects.append(ax2.bar( ind[16],  all_2015[16]*1e3, yerr=all_2015_err[16]*1e3, width=width, color='white',ecolor='blue',edgecolor='blue',  hatch="\\\\") )
#rects.append(ax2.bar( ind[17:19],  all_2015[17:19]*1e3,width=width, color='white',ecolor='magenta',edgecolor='magenta',  hatch="//") )
                 
rects.append(ax2.bar( ind[:6],  all_2015[:6]*1e3, width=width, color='lightgreen', ecolor='lightgreen', edgecolor='white' ) )
rects.append(ax2.bar( ind[6:13],  all_2015[6:13]*1e3, width=width, color='skyblue', ecolor='skyblue', edgecolor='white' ) )
rects.append(ax2.bar( ind[13:17],  all_2015[13:17]*1e3, width=width,color='r', edgecolor='white') )
rects.append(ax2.bar( ind[17],  all_2015[17]*1e3, yerr=all_2015_err[17]*1e3, width=width, color='khaki',ecolor='khaki',edgecolor='white') )
rects.append(ax2.bar( ind[18:20],  all_2015[18:20]*1e3,width=width, color='indianred',ecolor='indianred',edgecolor='white') )
           
 
 
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
 
#ax.grid()
ax2.yaxis.grid(color='lightgray', linestyle='--')

ax2.set_ylabel('\u0394NEP [TgC yr$^{-1}$]',  fontsize=fz*0.9)
 
ax2.set_xticks(ind)
#ax.set_ylim([-170,0])
ax2.set_ylim([-200,200])
ax2.locator_params(axis='y',nticks=2)
 
ax2.set_xticklabels( ['CTE2018','CT2019B', 'CAMS_v18r2','JenaStandard_v43','JenaNEET_v43','JenaNEETW_v43', 'CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB', 'Byrne2020','CMS-Flux2020','GCASv2','CCDAS_SM+VOD', 'TRENDY_v6', 'FLUXCOM_CRUJRA', 'FLUXCOM_ERA5'],rotation=90, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')

#lgnames.append("CTE")
#lgnames.append("EnKF-RAMS")
#lgnames.append("FLEXINVERT")
#lgnames.append("CarboScope-Regional")
#lgnames.append("LUMIA")
#lgnames.append("PYVAR-CHIMERE")
#lgnames.append("NAME-HB")
#lgnames.append("Mean")

#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax2.legend( rects, lgnames, loc = 2, ncol=2, frameon=False, fontsize = fz*0.8)
plt.tick_params(labelsize=fz*0.9)  


plt.text(0.75,0.85,'Drought 2012', horizontalalignment='left', verticalalignment='bottom', \
    transform = ax1.transAxes, fontsize=18, color='k')

plt.text(0.75,0.85,'Drought 2015', horizontalalignment='left', verticalalignment='bottom', \
    transform = ax2.transAxes, fontsize=18, color='k')    


yy = np.zeros(len(all_2012)) 
ax1.plot(yy, color='gray',linestyle='--',linewidth=0.5)
ax2.plot(yy, color='gray',linestyle='--',linewidth=0.5)
 

plt.setp(ax1.get_xticklabels(), visible=False) 

#--------------------------------------

fig.tight_layout() 

plt.subplots_adjust(left=0.15,wspace=0.10,hspace=0.10)
       

dir1='/Volumes/HEWEI_T5/Inversions/Byrne_inversions_GOSAT_surface_TCCON_2020/Byrneetal2020/monthly/'

fn='Bar_subcontinental_anomalies_drought_events_Europe_20221008_all'
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)


 