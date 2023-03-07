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

from scipy import signal
from scipy.stats import pearsonr

def linreg(X, Y):
    """
    return a,b in solution to y = ax + b such that root mean square distance between trend line and original points is minimized
    """
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in zip(X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    return (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det


def globarea(im=360, jm=180, silent=True):
    """ Function calculates the surface area according to TM5 definitions"""

    radius = 6.371e6  # the earth radius in meters
    deg2rad = np.pi / 180.
#    g = 9.80665 

    dxx = 360.0 / im * deg2rad 
    dyy = 180.0 / jm * deg2rad 
    lat = np.arange(-90 * deg2rad, 90 * deg2rad, dyy)
    dxy = dxx * (np.sin(lat + dyy) - np.sin(lat)) * radius ** 2
    area = np.resize(np.repeat(dxy, im, axis=0) , [jm, im])
    if not silent:
        print('total area of field = ', np.sum(area.flat))
        print('total earth area    = ', 4 * np.pi * radius ** 2)
    return area

  
   
file = '/Volumes/HEWEI_T5/codes_backup/20200610_final_paper/20200610_final_paper/Fig2_regional_analysis/EUROCOM_Europe_regions.nc'
   
regions_Europe=[]
 
ncf=nc.Dataset(file)
regions_Europe = ncf.variables['ID'][:]
ncf.close()

regions_Europe=np.array(regions_Europe)


mask = (regions_Europe==-9999)
regions_Europe[mask] = 0

regions_Europe_R1 = np.flipud(regions_Europe[0])
regions_Europe_R2 = np.flipud(regions_Europe[1])
regions_Europe_R3 = np.flipud(regions_Europe[2])
regions_Europe_R4 = np.flipud(regions_Europe[3])
regions_Europe_all = np.flipud(regions_Europe[4])

mask = (regions_Europe_all>=1)
regions_Europe_all[mask] = 1


def calc_seasonal_anomaly(data, M,P,N):
     
    a=[]
    
    for i in range(0,data.shape[0]): 
        me = np.nansum(data[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(M,M+N):
            baseline[i] = baseline[i] + np.nansum(data[12*j+i])/float(N)  #float(M+N)
    
    del data
 
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    if P>=1:
        A=a[12*M:-12*P] 
        B=a[12*M:-12*P]-a1
    
    if P==0:
        A=a[12*M:] 
        B=a[12*M:]-a1
    
    return A,B,baseline

 

#################################################################
#   TRENDY v6
#################################################################

dir1='/Volumes/HEWEI_T5/TRENDY_v6_s3/'

def flux_anomaly_TRENDYv6(item, varstr, M, N, latMin, latMax, lonMin, lonMax):
     
#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-15 #-10
#    lonMax=35 
#    

#    M=1
#    N=15
#    
#    item='NEP'
#    varstr='ensemble_mean'
    
    file= os.path.join(dir1,'TRENDY.v6.model.ensemble.'+item+'.0.5Deg.2000-2016.nc')
    
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = []
    
    for i in range(0,data1.shape[0]): 
        data2.append( np.transpose(data1[i]))   
     
    data2 = np.array(data2)
   
    
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
   # for i in range(0,data1.shape[0]):
   #     data2[i] = maskoceans(lon2, lat2, data1[i] )  #* area)
 
    data3 = data2[:-12,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  # Corrected July 3, 2017
    #print data2.shape, data3.shape
     
    #plt.imshow(data3[0])
    
 
    
    mask=(data3==9.96921e+36)
    data3[mask]=np.NaN
    
 
    plt.imshow(data3[0])
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]< 0 :
#                   data3[i,xx,yy]=np.NaN
#               else:      # kg C m-2 s-1
#                   data3[i,xx,yy]=data3[i,xx,yy]*365*24*3600*1e-15*1e3*1e3
 
    
    #data3=data3*area[2*latMin:2*latMax,2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
         
    del data1
    del data2
   
           
#    a=[]
#    
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(M,M+N):
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
#       
#    
#    del data3
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a[12*M:-12*1] 
#    B=a[12*M:-12*1]-a1
#    
#    return A,B,baseline

    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

 

 
#################################################################
#   CTE2018-flux1*1
#################################################################

dir11='/Volumes/HEWEI_T5/Inversions/'


def flux_anomaly_cte2018(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=1 
#    N=15  
#    varstr='bio_flux_opt'
    

    file= os.path.join(dir11,'CTE2018_monthly.nc')  #2000-2017  12*18=216
    
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    #data2 = data1
    area = globarea(360,180,True)
    
    lon=np.arange(0,360,1)
    lat=np.arange(0,180,1)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
 
    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2[(12*M):(-12*2),90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017


    #print data2.shape, data3.shape
     
    #plt.imshow(data3[0])

    mask=(data3==0)
    data3[mask]=np.NaN
    
    data3=data3*12.0*365*24*3600*1e-15
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3[0])
    
    data3=data3*area[90-latMax:90-latMin,180+lonMin:180+lonMax]  
      
    

#    a=[]
#          
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
   
    DATA = np.zeros(shape= (data3.shape[0],data3.shape[1]*2,data3.shape[2]*2) )
 
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*2):(j*2+2+1)] = np.repeat(data3[i,j], 2, axis=0)   
 
    M=0   
    P=0
    R1 = calc_seasonal_anomaly(DATA*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(DATA*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(DATA*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(DATA*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(DATA*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

#################################################################
#   CTE2019-flux1*1
#################################################################

#dir11b='/Volumes/HEWEI_WD/Inverions/CTE2019/'
#
#def flux_anomaly_cte2019(varstr, M, N, latMin, latMax, lonMin, lonMax):
#         
##    M=1 
##    N=17  
##    varstr='bio_flux_opt'
#    
#    file= os.path.join(dir11b,'flux1x1_all_months.nc')  #2000-2017  12*18=216
#    
#    data1=[]
#    
#     
#    ncf=nc.Dataset(file)
#    data1=ncf.variables[varstr][:]
#         
#    ncf.close()
#    
#    data1=np.array(data1)
#    #data2 = data1
#    area = globarea(360,180,True)
#    
#    lon=np.arange(0,360,1)
#    lat=np.arange(0,180,1)
#    lon2, lat2 = np.meshgrid(lon,lat) 
#    
#    #for i in range(0,data1.shape[0]):
#    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
#     
# 
#    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
#    data2=[]
#
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[(12*M):(-12*4),90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017
#
#
#    #print data2.shape, data3.shape
#     
#    #plt.imshow(data3[0])
#    
#
#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN
#    
#    data3=data3*12.0*365*24*3600*1e-15
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]==0 :
##                   data3[i,xx,yy]=np.NaN
##               else:
##                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
#    
#    #plt.imshow(data3[0])
#    
#    data3=data3*area[90-latMax:90-latMin,180+lonMin:180+lonMax]  
#      
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
   

 
#################################################################
#   CT2019B-flux1*1
#################################################################
 
dir12='/Volumes/HEWEI_T5/Inversions/'

def flux_anomaly_ct2019b(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-15   #-10
#    lonMax=35 
#
#     
#    M=1 
#    N=15  
#    varstr='bio_flux_opt'
 
        
    file= os.path.join(dir12,'CT2019B.flux1x1-monthly.nc')
    
    data1=[]
    #area=[]
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
    #area=ncf.variables['cell_area'][:]  #same as calculated
    
    ncf.close()
    
    data1=np.array(data1)
    #data2 = data1
    
    #area=np.array(area)
    area = globarea(360,180,True)
    

    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2[(12*M):(-12*3),90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017

    #print data2.shape, data3.shape
     
    #plt.imshow(data3[0])

    
#    mask=(data3==0)   #mask=(data3==-1.0E34)
#    data3[mask]=np.NaN
    data3=data3*12.0*365*24*3600*1e-15
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==-1.0E34 :  #0 for CTE2015
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3[0])
    
    data3=data3*area[90-latMax:90-latMin,180+lonMin:180+lonMax]  
        

#    a=[]      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a
#    B=a-a1
#    
#    return A,B
    
    DATA = np.zeros(shape= (data3.shape[0],data3.shape[1]*2,data3.shape[2]*2) )
 
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*2):(j*2+2+1)] = np.repeat(data3[i,j], 2, axis=0)   
    
    M=0  
    P=0
    R1 = calc_seasonal_anomaly(DATA*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(DATA*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(DATA*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(DATA*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(DATA*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT


#################################################################
#   CAMS-flux1*1 
#################################################################
 
dir13='/Volumes/HEWEI_T5/Inversions/CAMS_v18r2/1Deg/'

def flux_anomaly_cams(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=22 
#    N=15 
#    
#    latMin=35
#    latMax=70
#    lonMin=-10
#    lonMax=40
#    
#    varstr='flux_apos_bio'
    

    files = glob.glob(os.path.join(dir13,'Resampled_1Deg*.nc')) 
    files=sorted(files)
    
    data1=[]
    #area=[]
    lat=[]
    lon=[]
     
    ncf=nc.Dataset(files[0])
    lat=ncf.variables['lat'][:]
    lon=ncf.variables['lon'][:]
    ncf.close()   
        
    for file in files:
        ncf=nc.Dataset(file)
        data1.append(ncf.variables[varstr][:])
        #area=ncf.variables['cell_area'][:]  #same as calculated
     
        ncf.close()
    
    data1=np.array(data1)
    lat = np.array(lat)
    lon=np.array(lon)
    
        
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #area=np.array(area)
    area = globarea(360,180,True)
    
    data2=[]
    m=data1.shape[0]
    for i in range(0,m):
        xx = maskoceans(lon2, lat2, data1[i])
        data2.append(np.flipud(xx))
    
    data2=np.array(data2)
    
    #data problem, only valid from 1 to 147  all 147 months, 12 years
    
    data3=data2[(12*M):(-12*3),(90-latMax):(90-latMin),(180+lonMin):(180+lonMax)]   # Corrected on July 3, 2017
 
    #plt.imshow(data3[0])

    
    mask=(data3==0)
    data3[mask]=np.NaN
    
    # now kgC m-2 month-1
    data3=data3*12*1e3*1e-15   
   
    data3=data3*area[(90-latMax):(90-latMin),(180+lonMin):(180+lonMax)]  
    
    
#    a=[]      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a
#    B=a-a1
#    
#    return A,B

    DATA = np.zeros(shape= (data3.shape[0],data3.shape[1]*2,data3.shape[2]*2) )
 
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*2):(j*2+2+1)] = np.repeat(data3[i,j], 2, axis=0)   
                 
    M=0   
    P=0
    R1 = calc_seasonal_anomaly(DATA*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(DATA*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(DATA*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(DATA*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(DATA*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

#################################################################
#   MACTM-flux1*1
#################################################################

#dir15='/Users/weihe/Documents/Research_Projects/Drought_Carbon_Project/'
#
#def flux_anomaly_mactm(varstr, M, N, latMin, latMax, lonMin, lonMax):
#         
#    #M=5 
#    #N=15  
#    #varstr='flux_apos_land'
#    
#    file= os.path.join(dir15,'co2.MACTM-r84.JAMSTEC.v1.1996-2016_mask1x1.nc')
#    
#    data1=[]
#    lat=[]
#    lon=[]
#    #area=[]
#     
#    ncf=nc.Dataset(file)
#    data1=ncf.variables[varstr][:]
#    lat=ncf.variables['lat'][:]
#    lon=ncf.variables['lon'][:]
#    #area=ncf.variables['cell_area'][:]  #same as calculated
#    
#    ncf.close()
#    
#    data1=np.array(data1)
#    lat = np.array(lat)
#    lon=np.array(lon)
#        
#    lon2, lat2 = np.meshgrid(lon,lat) 
#    
#    #data2 = data1
#    
#    #area=np.array(area)
#    area = globarea(360,180,True)
#    
#    data2=[]
#
#    m=data1.shape[0]
#    for i in range(0,m):
#        #xx = maskoceans(lon2, lat2, data1[i])
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[(12*M):(-12*1),(90-latMax):(90-latMin),(180+lonMin):(180+lonMax)]   # Corrected on July 3, 2017
#
#    #print data2.shape, data3.shape
#     
#    plt.imshow(data3[0])
#    
#
#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN
#    
#    data3=data3*12*1e3*1e-15  #kg C/m2/month
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]==-1.0E34 :  #0 for CTE2015
##                   data3[i,xx,yy]=np.NaN
##               else:
##                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
#    
#    #plt.imshow(data3[0])
#    
#    data3=data3*area[90-latMax:90-latMin,180+lonMin:180+lonMax]  
#      
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a
#    B=a-a1
#    
#    return A,B

#
##################################################################
##   jenas99
##################################################################
#
#dir16='/Users/weihe/Documents/Research_Projects/Drought_Carbon_Project/Jena_CarboScope/monthlys99v3.8/2001-2015/'
#
#
#def flux_anomaly_jenas99(varstr, M, N, latMin, latMax, lonMin, lonMax):
#         
#    #M=0 
#    #N=15  
#    #varstr='co2flux_land'
#    
#    
#    files=glob.glob(os.path.join(dir16,'CarboScope_s99_v3.8_monthly_72x48_*.nc')) 
#    
#    
#    data1=[]
#    
#     
#    for i in files:       
#         ncf=nc.Dataset(i)
#         data1.append(ncf.variables[varstr][:])
#         ncf.close()
#    
#    data1=np.array(data1)
#    
#    
#    #data2 = data1
#    #area = globarea(72,48,True)
#    
#    #lon=np.arange(0,360,5)
#    #lat=np.arange(0,180,3.75)
#    #lon2, lat2 = np.meshgrid(lon,lat) 
#    
#    #for i in range(0,data1.shape[0]):
#    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
#     
# 
#    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
#    data2=[]
#
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[:,int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   # Corrected on July 3, 2017
#
#
#    #print data2.shape, data3.shape
#     
#    #plt.imshow(data3[0])
#    
#
#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN
#    
#    #data3=data3*12.0*12*1e12*1e-15
#    data3=data3  #*12.0*12*1e12*1e-15
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]==0 :
##                   data3[i,xx,yy]=np.NaN
##               else:
##                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
#    
#    #plt.imshow(data3[0])
#    
#    #data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
#      
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
   
#################################################################
#   jenas99
#################################################################

dir14='/Volumes/HEWEI_T5/Inversions/Jena_CarboScope/monthlys99v4.3/2001-2015/'

def flux_anomaly_jenas99v43(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=15  
#    varstr='co2flux_land'
    
    
    files=glob.glob(os.path.join(dir14,'CarboScope_s99_v4.3_monthly_72x48_*.nc')) 
    files=sorted(files)
    
    data1=[]
    
     
    for i in files:       
         ncf=nc.Dataset(i)
         data1.append(ncf.variables[varstr][:])
         ncf.close()
    
    data1=np.array(data1)
    
    
    #data2 = data1
#    area = globarea(72,48,True)
    
    #lon=np.arange(0,360,5)
    #lat=np.arange(0,180,3.75)
    #lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
 
    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2[:,int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   # Corrected on July 3, 2017


    #print data2.shape, data3.shape    
    #plt.imshow(data3[0])
    

    mask=(data3==0)
    data3[mask]=np.NaN
    
    #data3=data3*12.0*12*1e12*1e-15
    data3=data3  #*12.0*12*1e12*1e-15
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3[0])
    
#    data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
      

#    a=[]
          
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
   
    #  72 * 48 
  
    DATA = np.zeros(shape= (data3.shape[0],data3.shape[1]*15, data3.shape[2]*20) )
 
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*15):(j*15+15+1)] = np.repeat(data3[i,j], 20, axis=0)   
 
    DATA1 = np.zeros(shape= (data3.shape[0],int(data3.shape[1]*15/2), int(data3.shape[2]*20/2)) )
 
    for k in range(0,data3.shape[0]):
           for i in range(0,data3.shape[1]):
                  for j in range(0,DATA1.shape[2]):
                       DATA1[k,i:j] = np.nansum(DATA[k,(2*i):(2*i+2+1),(2*j):(2*j+2+1)])/300.0   
     
    
    P=0
    R1 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   

#################################################################
#   jenas04
#################################################################

#dir17='/Users/weihe/Documents/Research_Projects/Drought_Carbon_Project/Jena_CarboScope/monthlys04v4.1/'
#
#def flux_anomaly_jenas04(varstr, M, N, latMin, latMax, lonMin, lonMax):
#         
#         
#    #M=0 
#    #N=12 
#    #varstr='co2flux_land'
#    
#    
#    files=glob.glob(os.path.join(dir17,'CarboScope_s04_v4.1_monthly_72x48_*.nc')) 
#    
#    
#    data1=[]
#    
#     
#    for i in files:       
#         ncf=nc.Dataset(i)
#         data1.append(ncf.variables[varstr][:])
#         ncf.close()
#    
#    data1=np.array(data1)
#    
#    
#    #data2 = data1
#    #area = globarea(72,48,True)
#    
#    #lon=np.arange(0,360,5)
#    #lat=np.arange(0,180,3.75)
#    #lon2, lat2 = np.meshgrid(lon,lat) 
#    
#    #for i in range(0,data1.shape[0]):
#    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
#     
# 
#    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
#    data2=[]
#
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[:,int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   # Corrected on July 3, 2017
#
#
#    #print data2.shape, data3.shape
#     
#    #plt.imshow(data3[0])
#    
#
#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN
#    
#    #data3=data3*12.0*12*1e12*1e-15
#    data3=data3  #*12.0*12*1e12*1e-15
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]==0 :
##                   data3[i,xx,yy]=np.NaN
##               else:
##                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
#    
#    #plt.imshow(data3[0])
#    
#    #data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
#      
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,M+N):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
   

#################################################################
#   jenas04
#################################################################

#dir17a='/Users/weihe/Documents/Research_Projects/Drought_Carbon_Project/Jena_CarboScope/monthlys04v4.3/'
#
#def flux_anomaly_jenas04v43(varstr, M, N, latMin, latMax, lonMin, lonMax):
#         
#         
#    #M=0 
#    #N=12 
#    #varstr='co2flux_land'
#    
#    
#    files=glob.glob(os.path.join(dir17a,'CarboScope_s04_v4.3_monthly_72x48_*.nc')) 
#    
#    
#    data1=[]
#    
#     
#    for i in files:       
#         ncf=nc.Dataset(i)
#         data1.append(ncf.variables[varstr][:])
#         ncf.close()
#    
#    data1=np.array(data1)
#    
#    
#    #data2 = data1
#    #area = globarea(72,48,True)
#    
#    #lon=np.arange(0,360,5)
#    #lat=np.arange(0,180,3.75)
#    #lon2, lat2 = np.meshgrid(lon,lat) 
#    
#    #for i in range(0,data1.shape[0]):
#    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
#     
# 
#    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
#    data2=[]
#
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[:,int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   # Corrected on July 3, 2017
#
#
#    #print data2.shape, data3.shape
#     
#    #plt.imshow(data3[0])
#    
#
#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN
#    
#    #data3=data3*12.0*12*1e12*1e-15
#    data3=data3  #*12.0*12*1e12*1e-15
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]==0 :
##                   data3[i,xx,yy]=np.NaN
##               else:
##                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
#    
#    #plt.imshow(data3[0])
#    
#    #data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
#      
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,M+N):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
 

#################################################################
#   jenas NEET v4.3
#################################################################

dir15='/Volumes/HEWEI_T5/Inversions/Jena_CarboScope/monthlysEXTocNEETv4.3/2001-2015/'


def flux_anomaly_jenaNEETv43(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=15  
#    varstr='co2flux_land'
  
    files=glob.glob(os.path.join(dir15,'CarboScope_NEET_v4.3_monthly_72x48_*.nc')) 
    files=sorted(files)
    
    data1=[]
    
     
    for i in files:       
         ncf=nc.Dataset(i)
         data1.append(ncf.variables[varstr][:])
         ncf.close()
    
    data1=np.array(data1)
    
    
    #data2 = data1
#    area = globarea(72,48,True)
    
    #lon=np.arange(0,360,5)
    #lat=np.arange(0,180,3.75)
    #lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
 
    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2[:,int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   # Corrected on July 3, 2017


    #print data2.shape, data3.shape    
    #plt.imshow(data3[0])
    

    mask=(data3==0)
    data3[mask]=np.NaN
    
   # data3=data3*12.0*12*1e12*1e-15   
    data3=data3    #12.0*12*1e12*1e-15   
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3[0])
    
#    data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
      

#    a=[]
          
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,15):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,15):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B


    DATA = np.zeros(shape= (data3.shape[0],data3.shape[1]*15, data3.shape[2]*20) )
 
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*15):(j*15+15+1)] = np.repeat(data3[i,j], 20, axis=0)   
 
    DATA1 = np.zeros(shape= (data3.shape[0],int(data3.shape[1]*15/2), int(data3.shape[2]*20/2)) )
 
    for k in range(0,data3.shape[0]):
           for i in range(0,DATA1.shape[1]):
                  for j in range(0,DATA1.shape[2]):
                       DATA1[k,i,j] = np.nansum(DATA[k,(2*i):(2*i+2+1),(2*j):(2*j+2+1)])/300.0   
     
    P=0
    R1 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

 
    
   
#################################################################
#   jenas NEETW v4.3
#################################################################
dir16='/Volumes/HEWEI_T5/Inversions/Jena_v4.3_all/NEE-T-W-new/s99/monthly/2001-2015/'

def flux_anomaly_jenaNEETWv43(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    M=9 
#    N=6  
#    varstr='co2flux_land'
    
    
    files=glob.glob(os.path.join(dir16,'CarboScope_NEETW_v4.3_monthly_72x48_*.nc')) 
    files=sorted(files)
    
    data1=[]
    
     
    for i in files:       
         ncf=nc.Dataset(i)
         data1.append(ncf.variables[varstr][:])
         ncf.close()
    
    data1=np.array(data1)
    
    
    #data2 = data1
#    area = globarea(72,48,True)
    
    #lon=np.arange(0,360,5)
    #lat=np.arange(0,180,3.75)
    #lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
 
    #data3 = data2[:,180-latMax:180-latMin,180+lonMin:180+lonMax]
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2[12*M:,int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   # Corrected on July 3, 2017


    #print data2.shape, data3.shape    
    #plt.imshow(data3[0])
    

#    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN
    
   # data3=data3*12.0*12*1e12*1e-15   
    data3=data3    #12.0*12*1e12*1e-15   
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3[0])
    
#    data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
      
      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,N):    #baseline： 2001-2014
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
    
    DATA = np.zeros(shape= (data3.shape[0],data3.shape[1]*15, data3.shape[2]*20) )
 
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*15):(j*15+15+1)] = np.repeat(data3[i,j], 20, axis=0)   
 
    DATA1 = np.zeros(shape= (data3.shape[0],int(data3.shape[1]*15/2), int(data3.shape[2]*20/2)) )
 
    for k in range(0,data3.shape[0]):
           for i in range(0,DATA1.shape[1]):
                  for j in range(0,DATA1.shape[2]):
                       DATA1[k,i,j] = np.nansum(DATA[k,(2*i):(2*i+2+1),(2*j):(2*j+2+1)])/300.0   
     
    P=0
    R1 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(DATA1[:,5:-5,:]*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT
 

#################################################################
#   FLUXCOM
#################################################################

dir17='/Volumes/HEWEI_T5/FLUXCOM2020/RS+METEO/member/CRUJRA_v1/'
 
    
def flux_anomaly_fluxcom_RSMETEO1(varstr, ml, M, N, latMin, latMax, lonMin, lonMax):
       
#    M=0 
#    N=15  
#    varstr='NEE'
     
    files=glob.glob(os.path.join(dir17, varstr,'2001-2015',ml,'*.nc'))  
    files=sorted(files)
    
    data1=[]
    
    for file in files:
         ncf=nc.Dataset(file)
         data1.extend(ncf.variables[varstr][:])
         
         ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
      
    data3 = data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  # Corrected July 3, 2017
      
    #print data2.shape, data3.shape
    ###a=data3.mean(1).mean(1)
     
    #plt.imshow(data3[0])
    
#    a=[]
    
    mask=(data3==-9999)    #gC m-2 d-1
    data3[mask]=np.NaN
    
    data3=data3*365*1e-15
    

    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]
      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,N):   
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a[12*M:] 
#    B=a[12*M:]-a1
#    
#    return A,B,baseline

    P=0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

 
#################################################################
#   FLUXCOM
#################################################################

dir18='/Volumes/HEWEI_T5/FLUXCOM2020/RS+METEO/member/ERA5/'
    
def flux_anomaly_fluxcom_RSMETEO2(varstr,ml, M, N, latMin, latMax, lonMin, lonMax):
       
#    M=0 
#    N=15  
#    varstr='NEE'
    
  
    files=glob.glob(os.path.join(dir18,varstr,'2001-2015',ml,'*.nc'))  
    files=sorted(files)
    
    data1=[]
    
    for file in files:
         ncf=nc.Dataset(file)
         data1.extend(ncf.variables[varstr][:])
         
         ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
      
    data3 = data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  # Corrected July 3, 2017
      
    #print data2.shape, data3.shape
    ###a=data3.mean(1).mean(1)
     
    #plt.imshow(data3[0])
    
#    a=[]
    
    mask=(data3==-9999)    #gC m-2 d-1
    data3[mask]=np.NaN
    
    data3=data3*365*1e-15
    

    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]
      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(0,N):   
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
#   
#    
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a[12*M:] 
#    B=a[12*M:]-a1
#    
#    return A,B,baseline

    P=0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

 
#################################################################
#   EUROCOM - CTE, RAMS, FLEXPART,Jena, LUMIA, PYVAR_CHIMERE
#################################################################

dir20='/Volumes/HEWEI_T5/Inversions/EUROCOM-ICOS/'

def flux_anomaly_eurocom(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=10
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
#    
#    varstr='bio_posterior'
 
    file= os.path.join(dir20,dataset)
    
    data1=[]
    lon=[]
    lat=[]
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
    lat=ncf.variables['lat'][:]
    lon=ncf.variables['lon'][:]
         
    ncf.close()
    
    data1=np.array(data1)  #(120, 80, 100)
    lat=np.array(lat) 
    lon=np.array(lon) 
    area = globarea(720,360,True)
    
    
    plt.imshow(data1[0])
    
    
    #lon=np.arange(0,360,0.5)
    #lat=np.arange(0,180,1)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2  #[(12*M):,:,:]  #,90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017


    #print data2.shape, data3.shape
     
    #plt.imshow(data3[0])
    

#    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN   #kgC/m2/month
    
    data3=data3*12.0*1e-12
 
    data3=data3*area[2*(90-latMax):2*(90-latMin),int(2*(180+lonMin)-1):2*(180+lonMax)]   #(72, 80, 100)
      
      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):  
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
#   
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
     
    
    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT
 


#################################################################
#   Byrne GOSAT+insitu+TCCON CO2 inversion: 2009-2015
#################################################################

dir22='/Volumes/HEWEI_T5/Inversions/Byrne_inversions_GOSAT_surface_TCCON_2020/Byrneetal2020/monthly/'

def flux_anomaly_ByrneGOSAT(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=6
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
#  
#    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'      
#    varstr='ensemble_mean'

    #file= os.path.join(dir22,dataset)
    files = glob.glob(os.path.join(dir22,dataset+'*.nc')) 
    files=sorted(files)
    
    data1=[]
    lon=[]
    lat=[]
     
    ncf=nc.Dataset(files[1])

    lat=ncf.variables['Lat'][:]
    lon=ncf.variables['Lon'][:]
    area = ncf.variables['area'][:]
    ncf.close()
    
    for file in files:
        ncf=nc.Dataset(file)
        data1.append(ncf.variables[varstr][:])
        ncf.close()
          
    data1=np.array(data1)  #(120, 80, 100)
    lat=np.array(lat) 
    lon=np.array(lon) 
    area = np.array(area)  #globarea(720,360,True)
    
    
    #plt.imshow(data1[0])
 
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
 
       
    data3 = data2   #[12*M: ]  #, int((90-latMax)/3.75-1):int((90-latMin)/3.75-1), int((180+lonMin)/5.0):int((180+lonMax)/5.0) ]
 
    mask=(data3==0)
    data3[mask]=np.NaN
      
    data3 = data3  *12*1e-15
      
    data4 =  np.zeros(shape= (data3.shape[0],data3.shape[1]*16-16,data3.shape[2]*20) )
    for i in range(0,data3.shape[0]):
        for j in range(0,data3.shape[1]):
             data4[i,(j-1)*16:(j*16)] = np.repeat(data3[i,j], 20, axis=0)   
    
    #data4=data4[:,5:-5,:]/300.0 * mask_China2
    
    
    area = globarea(1440,720,True)
    data4 = data4*area[:,:]
    
    data4=data4[:,int(4*(90-latMax)):int(4*(90-latMin)), int(4*(180+lonMin)):int(4*(180+lonMax)) ]  


#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN   #gC/m2/month
#    
#    data3=data3*12.0*1e-15
# 
#    data3=data3*area[int(np.floor((90-latMax)/3.75-1)):int(np.ceil((90-latMin)/3.75-1)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   
#
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):  
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
#   
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
     
    data5 =  np.zeros(shape= (data4.shape[0],int(data4.shape[1]/2),int(data4.shape[2]/2)) )
    for k in range(0,data4.shape[0]):
        for i in range(0,int(data4.shape[1]/2)):
            for j in range(0,int(data4.shape[2]/2)):
                 data5[k,i,j] = np.nanmean(data4[k,i*2:(i+1)*2,j*2:(j+1)*2])   
    
    
    data3 = data5
    
    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT 

#################################################################
#   Liu et al.(2020) ESSD paper, CMS-Flux 2020, GOSAT inversions
#################################################################

dir23='/Volumes/HEWEI_T5/Inversions/NASA-CMS-Flux/ESSD_paper_v2/monthly_1Deg/'
  
def flux_anomaly_LiuGOSAT(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=6
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
#  
#    dataset='CMS-Flux-NBE-2020.monthly.grid'           
#    varstr = 'flux'
    
    #files=glob.glob(os.path.join(dir23,dataset+'*.nc'))  
    #files=sorted(files)
       
    files=glob.glob(os.path.join(dir23,'CMS-Flux*.nc'))  
    files=sorted(files)
    
     
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(files[0])
    lat.extend(ncf.variables['Lat'][:])
    lon.extend(ncf.variables['Lon'][:])
         
    for file in files:
        ncf=nc.Dataset(file)
        data1.append(ncf.variables[varstr][:])
        ncf.close()
    
    data1=np.array(data1)
    data2=data1
    
    area = globarea(360*2,180*2,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
    
    #data2a = np.empty( (12*N, data2.shape[1], data2.shape[2]), dtype=float, order='C')
 
  
    data3 = data2[:-12*3]   #[:, int(np.floor(90-latMax)):int(np.ceil(90-latMin)), int(np.floor(180+lonMin)):int(np.ceil(180+lonMax)) ]
     
  
    mask=(data3==0)
    data3[mask]=np.NaN
     
    
    data3 = data3*365*1e-15
    
    
    #data3=data3*area[ int(np.floor(90-latMax)):int(np.ceil(90-latMin)), int(np.floor(180+lonMin)):int(np.ceil(180+lonMax)) ] 
     
 
    DATA =  np.zeros(shape= (data3.shape[0],data3.shape[1]*2,data3.shape[2]*2) )
         
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*2):(j*2+2+1)] = np.repeat(data3[i,j], 2, axis=0)   
            
    data3 = DATA[:, (int(90-latMax)*2):(int(90-latMin)*2),(int(180+lonMin)*2):(int(180+lonMax)*2) ]   #12*M:

    data3 = data3*area[(int(90-latMax)*2):(int(90-latMin)*2),(int(180+lonMin)*2):(int(180+lonMax)*2)]  #*12*1e-15   # gC/m2/month    # kg C m-2 s-1 365*24*3600
  
    
    
    
#    a=[]
# 
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):  
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
#   
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
 
    
    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT 

#################################################################
#   Jiang et al.(2020) ACPD paper, GCASv2, GOSAT inversions
#################################################################

dir24='/Volumes/HEWEI_T5/Inversions/GCASv2/'
  
def flux_anomaly_JiangGOSAT(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
    M=0 
    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
        
    #varstr = 'bio_monthly_org'
    #varstr = 'bio_monthly_opt'
    
    file=os.path.join(dir24,'posterior.fluxes.GCASv2.nc')  
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(file)
    lat.extend(ncf.variables['lat'][:])
    lon.extend(ncf.variables['lon'][:])
    data1.append(ncf.variables[varstr][:]) 
    ncf.close()
    
    data1=np.array(data1[0])
    data2=data1
    
    area = globarea(360*2,180*2,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     

    data3 = np.zeros(shape= (data2.shape[0],data2.shape[1],data2.shape[2]) )
    mask = (data3==0)
    data3[mask] = np.nan
     
    for i in range(0,data2.shape[0]):
        data3[i, :, 0:180] = data2[i, :, 180:]
        data3[i, :, 180: ] = data2[i, :, 0:180]
     
        
    DATA =  np.zeros(shape= (data3.shape[0],data3.shape[1]*2,data3.shape[2]*2) )
         
    for i in range(0,data3.shape[0]):
            for j in range(0,data3.shape[1]):
                 DATA[i,(j*2):(j*2+2+1)] = np.repeat(data3[i,j], 2, axis=0)   
            
    data3 = DATA[:, (int(90-latMax)*2):(int(90-latMin)*2),(int(180+lonMin)*2):(int(180+lonMax)*2) ]   #12*M:

    data3 = data3*area[int(90-latMax)*2:int(90-latMin)*2,int(180+lonMin)*2:int(180+lonMax)*2]*12*1e-15   # gC/m2/month    # kg C m-2 s-1 365*24*3600
  
    
#    data3aa=data3aa*area[int(90-latMax):int(90-latMin),int(360+lonMin):]*12*1e-15   # gC/m2/month    # kg C m-2 s-1 365*24*3600
#    data3bb=data3bb*area[int(90-latMax):int(90-latMin),0:lonMax]*12*1e-15
  
#    a=[]
# 
#      
#    for i in range(0,data3aa.shape[0]): 
#        me = np.nansum(data3aa[i]) + np.nansum(data3bb[i])    
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):  
#            baseline[i] = baseline[i] + np.nansum(data3aa[12*j+i])/float(N) + np.nansum(data3bb[12*j+i])/float(N)
#   
#    del data3aa
#    del data3bb
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
  
    
    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

#################################################################
#   Wu et al.(2020) RSE paper, CCDAS-CO2+SM inversions
#################################################################

dir25='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_WuCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
    M=0 
    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
      
    dataset ='CCDAS_6yr_2deg_assim_results'    
    #varstr = 'rnep' 
    varstr = 'nep'
    file=os.path.join(dir25,dataset+'/'+'CCDAS_co2+sm_prior.nc')  
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(file)
    #lat.extend(ncf.variables['latitude'][:])
    #lon.extend(ncf.variables['longitude'][:])
    lat.extend(ncf.variables['lat'][:])
    lon.extend(ncf.variables['lon'][:])
    data1.append(ncf.variables[varstr][:]) 
    ncf.close()
    
    data1=np.array(data1[0])
    data2=data1
    
    area = globarea(180,90,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
        
    data3 = data2[12*M:, int((90-latMax)*0.5):int((90-latMin)*0.5),int((180+lonMin)*0.5):int((180+lonMax)*0.5)]
 
  
    mask=(data3==-9999)
    data3[mask]=np.NaN

 
    plt.imshow(data3[0])
 
    data3=data3*area[int((90-latMax)*0.5):int((90-latMin)*0.5),int((180+lonMin)*0.5):int((180+lonMax)*0.5)]*12*1e-15    #gC/m2/month
 
#    a=[]
# 
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])    
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):  
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  
#   
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
 
    
    data4 =  np.zeros(shape= (data3.shape[0],data3.shape[1]*4,data3.shape[2]*4) )
    for i in range(0,data3.shape[0]):
        for j in range(0,data3.shape[1]):
             #data4[i,(j-1)*4:(j*4+1)] = np.repeat(data3[i,j], 4, axis=0)/16.0   
             data4[i,j*4:(j*4+4)] = np.repeat(data3[i,j], 4, axis=0)/16.0   
    
    
   
    
    P = 0
    R1 = calc_seasonal_anomaly(data4*regions_Europe_R1, M, P, N)    #confirm later
    R2 = calc_seasonal_anomaly(data4*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data4*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data4*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data4*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT

  
#################################################################
#   Scholze et al.(2019) GRL paper, CCDAS-CO2+SM+VOD inversions
#################################################################

dir26='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_ScholzeCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
    M=0 
    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
 
    dataset ='ccdas_sm+vod'
    varstr = 'nep'
    
    file=os.path.join(dir26,dataset+'/'+'grid_nep_2010-2015.nc')  
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(file)
    lat.extend(ncf.variables['latitude'][:])
    lon.extend(ncf.variables['longitude'][:])
    data1.append(ncf.variables[varstr][:]) 
    ncf.close()
    
    data1=np.array(data1[0])
    data2=np.empty( (data1.shape[0], data1.shape[2],data1.shape[1]), dtype=float, order='C')
    #p.empty( (12*N, data2.shape[1], data2.shape[2]), dtype=float, order='C')
    
    area = globarea(1440,720,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.transpose(data1[i])
     
        
    data3 = data2[:, int((90-latMax)*4):int((90-latMin)*4),int((180+lonMin)*4):int((180+lonMax)*4)]   #12*M:
#    data3aa = data2[12*M:, int((90-latMax)*4):int((90-latMin)*4),int((360+lonMin)*4):]
#    data3bb = data2[12*M:, int((90-latMax)*4):int((90-latMin)*4),0:int(lonMax*4)]
#      
  
    mask=(data3==-9999)
    data3[mask]=np.NaN
    
 
    plt.imshow(data3[0])
 
    data3=data3*area[int((90-latMax)*4):int((90-latMin)*4),int((180+lonMin)*4):int((180+lonMax)*4)]*12*1e-15    #gC/m2/month
#    data3aa=data3aa*area[int((90-latMax)*4):int((90-latMin)*4),int((360+lonMin)*4):]*12*1e-15    #gC/m2/month
#    data3bb=data3bb*area[int((90-latMax)*4):int((90-latMin)*4),0:int(lonMax*4)]*12*1e-15    #gC/m2/month
      
#    a=[]
# 
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i]) 
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):  
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  
#   
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a 
#    B=a-a1
#    
#    return A,B
       
    data4 = data3
    data5 =  np.zeros(shape= (data4.shape[0],int(data4.shape[1]/2),int(data4.shape[2]/2)) )
    for k in range(0,data4.shape[0]):
        for i in range(0,int(data4.shape[1]/2)):
            for j in range(0,int(data4.shape[2]/2)):
                 data5[k,i,j] = np.nanmean(data4[k,i*2:(i+1)*2,j*2:(j+1)*2])   
    
    
    data3 = data5
    
    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1[:,:-1], M, P, N)    #confirm later
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2[:,:-1], M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3[:,:-1], M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4[:,:-1], M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all[:,:-1], M, P, N)

    return R1, R2, R3, R4, RT

 
    


#################################################################################
#                              CRUNCEP 
####################################################################################
dir2='/Volumes/HEWEI_T5/CRUNCEPv9_1981-2017/'

def anomaly_CRUNCEP(varstr, M, N,latMin, latMax, lonMin, lonMax):
 
    
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
#    M=0      
#    N=15     
#    varstr='tair'



    files=glob.glob(os.path.join(dir2, varstr,  'cruncep_v9*.nc'))
    
    files = sorted(files)
    
    data1=[]
    lon=[]
    lat=[]
    
    ncf=nc.Dataset(files[0])  
    lon = ncf.variables['Lon'][:] 
    lat = ncf.variables['Lat'][:] 
    ncf.close()    
    
    if varstr=='tair'  or varstr=='rain':
        for i in files[12*21:-12*2]:
            ncf=nc.Dataset(i)
            data1.append(ncf.variables[varstr][2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)])
            ncf.close()
    if varstr=='vpd':
         for i in files[:-12*2]:
            ncf=nc.Dataset(i)
            data1.append(ncf.variables[varstr][2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)])
            ncf.close()    
            
    data1=np.array(data1)
    lon=np.array(lon)
    lat=np.array(lat)
    
    #print data1.shape,lon.shape
    
    #data2 = data1
    data3 = data1  #data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]   #Corrected July 3, 2017
    del data1
    
    #mask=(data3<=-100)
    #data3[mask]=np.NaN
    
    plt.imshow(data3[0])
    
    
   
#    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
#    for i in range(0,N):
#        
#        if varstr=='tair':
#            #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
#            data4[i] = np.nanmean(data3[(12*i+0):(12*i+5)],axis=0)
#            #data4[i] = np.nanmean(data3[(12*i+2):(12*i+5)],axis=0)
#            #data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
#           
#        if varstr=='rain'  or varstr=='vpd':
#            data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
 
       
#    plt.imshow(data4[0])
#    
#    
#    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
#    
#    a=signal.detrend(a)    
#    mean = np.nanmean(a)
#    std = np.nanstd(a)
#    
#    b = (a-mean)/std
#
#    return a,b


    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   

 

####################################################################################
#                              SPEI index
####################################################################################

#dir2a='/Volumes/HEWEI_T5/SPEIbasev2.5_to2015/'
#
#def anomaly_spei(varstr, M, N, latMin,latMax,lonMin,lonMax):
#
#    
#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-15   #-10
#    lonMax=35    #40
#     
#    M=0      
#    N=15     
#    varstr='spei'
#
#    file=os.path.join(dir2a, 'spei06.nc') #09
#    
#    data1=[]
#    lon=[]
#    lat=[]
#     
#    ncf=nc.Dataset(file)
#    data1.append(ncf.variables[varstr][:])   
#    lon.append(ncf.variables['lon'][:]) 
#    lat.append(ncf.variables['lat'][:]) 
#         
#    ncf.close()
#    
#    data1=np.array(data1)
#    lon=np.array(lon)
#    lat=np.array(lat)
#     
#    #print data1.shape,lon.shape   #(1, 1368, 360, 720)  
# 
#    data2=[]
#    
#    for i in range(100*12,115*12):
#        data2.append(np.flipud(data1[0,i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]
#
#       
#    mask=(data3>=1.0e30)
#    data3[mask]=np.NaN
#       
#    plt.imshow(data3[0]) 
#  
#        
#
#    P = 0
#    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)
#    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
#    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
#    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
#    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)
#
#    return R1, R2, R3, R4, RT

 
####################################################################################
#                              GLEAM SM data
####################################################################################

dir3= '/Volumes/HEWEI_T5/GLEAMv3.3'

def anomaly_gleam(varstr, M, N, latMin,latMax,lonMin,lonMax):

    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
    M=0      
    N=15  
    varstr='SMroot'       
    
 
    data1=[]
    #lon=[]
    #lat=[]
    
    file = os.path.join(dir3,varstr+'_1980_2018_GLEAM_v3.3a_MO.nc')
 
        
    ncf=nc.Dataset(file)
    data1.extend(ncf.variables[varstr][:])    
        
    ncf.close()
    
    data1=np.array(data1)
    #lon=np.array(lon)  
    #lat=np.array(lat)
        
    #print data1.shape,lon.shape   #(1, 1368, 360, 720)
    
    data2=[]   
    for i in range(0,data1.shape[0]):
        data2.append(np.transpose(data1[i]))
    
    data2=np.array(data2)
 
     
    data3 = data2[12*21:-12*3,4*(90-latMax):4*(90-latMin),4*(180+lonMin):4*(180+lonMax)] 
    plt.imshow(data3[0])   
    
            
    mask=(data3==-999.0)    
    data3[mask]=np.NaN
    
    plt.imshow(data3[0]) 
 
 

#    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
#    for i in range(0,N):
#        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
#        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)


        
        
#    plt.imshow(data4[0])
#    
#    
#    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
#    
#    a=signal.detrend(a)
#    mean = np.nanmean(a)
#    std = np.nanstd(a)
#    
#    b = (a-mean)/std
#
#    return a,b
      
    data4 = data3
    
    data5 =  np.zeros(shape= (data4.shape[0],int(data4.shape[1]/2),int(data4.shape[2]/2)) )
    for k in range(0,data4.shape[0]):
        for i in range(0,int(data4.shape[1]/2)):
            for j in range(0,int(data4.shape[2]/2)):
                 data5[k,i,j] = np.nanmean(data4[k,i*2:(i+1)*2,j*2:(j+1)*2])   
     
    
    
    P = 0
    R1 = calc_seasonal_anomaly(data5*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data5*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data5*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data5*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data5*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   

  
####################################################################################
#                              CSIF  
####################################################################################
 
dir4='/Volumes/HEWEI_T5/CSIF_v2/monthly/'

def anomaly_CSIF(varstr, M, N, latMin, latMax, lonMin, lonMax):

#latMin=35
#latMax=70
#lonMin=-10
#lonMax=40
#varstr = 'clear_daily_sif'
     
#    M=0   
#    N=15   
#    varstr='SIF'

    data1=[]
    lon=[]
    lat=[]
     
    files=glob.glob(os.path.join(dir4,'OCO2.SIF.clear.*.nc'))
      
    ncf=nc.Dataset(files[0])    
    lon.append(ncf.variables['lon'][:])   
    lat.append(ncf.variables['lat'][:]) 
    
    for file in files[:-12*4]:
        ncf=nc.Dataset(file)
        #data1.extend(ncf.variables['clear_daily_SIF'][:])   #clear_daily_sif
        data1.append(ncf.variables['clear_daily_SIF'][:])
        ncf.close()
    
    data1=np.array(data1)
    lon=np.array(lon)[0]     
    lat=np.array(lat)[0]    
     
    #print data1.shape,lon.shape   #(96, 360, 720)  0.5*0.5
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2 = data1
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i])           # Added on July 3, 2017
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])           # Added on July 3, 2017
    
    
    #plt.imshow(data2[20])
    
    #data3 = data2[:,2*latMin:2*latMax,2*(180+lonMin):2*(180+lonMax)]
    #data3 = data2[:-92,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]   # Corrected on July 3, 2017
    data3 = data2[:,int(2*(90-latMax)):int(2*(90-latMin)),int(2*(180+lonMin)):int(2*(180+lonMax))]   # Corrected on July 3, 2017
     
#    a=[]
    
    mask=(data3==-999.9)
    data3[mask]=np.NaN
    
#plt.imshow(data3[20])
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0. :
#                   data3[i,xx,yy]=np.NaN

    
    
#    for i in range(0,data3.shape[0]): 
#        me = np.nanmean(data3[i])     
#        if math.isnan(np.nanmean(data3[i]))==True:
#            me = 0  # attention!!
#        if me < 0:
#            me = 0  # attention!!
#        a.append(me)
#        
#    a=np.array(a)
#    
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):
#            if math.isnan(np.nanmean(data3[12*j+i]))==True:
#                data3[12*j+i]=0    # attention!!
# 
#            baseline[i] = baseline[i] + np.nanmean(data3[12*j+i])/float(N)
#   
#    #print baseline
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#    
#    A=a 
#    B=a-a1
#    
#    return A,B
    

    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   


 
####################################   


dir4a='/Volumes/HEWEI_T5/MOD13C2/0.5Deg_NIRv1/'    
   
def anomaly_modisNIRv( M, N, latMin, latMax, lonMin, lonMax):
     
#    latMin=10
#    latMax=55
#    lonMin=70
#    lonMax=150
        
#    M=0 
#    N=15 
     
    files=glob.glob(os.path.join(dir4a,'*NIRv1.dat'))
     
    EVI=[]
    
    for i in files: 
    
        data = np.fromfile(i, dtype='float64', sep="")
        data = data.reshape(360,720)   
    
        EVI.append(data)
     
    EVI=np.array(EVI)
     
    data3 = EVI[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  
    
    #plt.imshow(data3[1])  
  
        
     
#    a=[]
#    
#    for i in range(0,data3.shape[0]): 
#        me = np.nanmean(data3[i])   
#        a.append(me)
#        
#    a=np.array(a)
#    
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        for j in range(0,N):     
#            baseline[i] = baseline[i] + np.nanmean(data3[12*j+i])/float(N)  
#     
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    
#       
#    A=a 
#    B=a-a1
#    
#    return A,B   

 
    P = 0
    R1 = calc_seasonal_anomaly(data3*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data3*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data3*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data3*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data3*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   



  
####################################################################################
#                              GLOBMAP LAI data
####################################################################################
 
dir5='/Volumes/HEWEI_T5/GLOBMAP_LAI_V3_0.5Deg_monthly_2001-2020/'

def anomaly_LAI(varstr, M, N, latMin,latMax,lonMin,lonMax):
  
        
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
 
   
    varstr = 'LAI'
    
    M=0
    N = 15  
           
    data1=[]
    lon=[]
    lat=[]
     
    files=glob.glob(os.path.join(dir5,'GLOBMAP_LAI*.nc'))    
    files=sorted(files)
    
    ncf=nc.Dataset(files[0])    
    lon.append(ncf.variables['lon'][:])   
    lat.append(ncf.variables['lat'][:]) 
    
    for file in files:
        ncf=nc.Dataset(file)
        data1.append(ncf.variables[varstr][:]) 
        
        ncf.close()
    
    
#    file = os.path.join(dir5c,'VODCA_X_1997_2018_yearly.nc')
#    
#    ncf=nc.Dataset(file)    
#    lon = ncf.variables['lon'][:] #longitude  
#    lat = ncf.variables['lat'][:]  #latitude
#    
#    ncf=nc.Dataset(file)
#    data1 = ncf.variables[varstr][:]
#    ncf.close()
#    
    
    data1=np.array(data1)    
    lon=np.array(lon)      
    lat=np.array(lat)     
    
    #print data1.shape,lon.shape   #(192, 720, 1440) 
    
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2 = np.zeros(shape=(data1.shape[0],data1.shape[2],data1.shape[1],) )
    for i in range(0,data1.shape[0]):
       data2[i] =  np.rot90(np.fliplr(data1[i]))         # Added on July 3, 2017
     
    #data2 = data1
    plt.imshow(data2[0])
     
    data3 = data2[M*12:(-4*12),2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  
    #data3 = data2[M*12:(-3*12),4*(90-latMax):4*(90-latMin),4*(180+lonMin):4*(180+lonMax)]  
     
    #del data1,data2
    
      
    mask=(data3>1000)
    data3[mask]=np.NaN
    
    plt.imshow(data3[0])
    
  
    
#    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
#    for i in range(0,N):
#        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
#        #data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
#        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
        
        
#    plt.imshow(data4[0])
#    
#    
#    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
#    
#    a=signal.detrend(a)
#    mean = np.nanmean(a)
#    std = np.nanstd(a)
#    
#    b = (a-mean)/std
#
#    return a,b

    data4 = data3

    P = 0
    R1 = calc_seasonal_anomaly(data4*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data4*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data4*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data4*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data4*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   


 
    

#################################################################
#   FluxSat GPP
#################################################################

#dir6='/Volumes/HEWEI_T5/FluxSat/0.5x0.5/2001-2015/'
dir6='/Volumes/HEWEI_T5/FluxSat_v2_partly/FluxSatv2_0.5Deg/2001-2015/' 
  
def flux_anomaly_fluxsat(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    M=0
#    N=15     
#    varstr = 'GPP'
 
        
    #files=glob.glob(os.path.join(dir7,'*/FluxSat_GPP_0.5_v1.1*.nc'))  
    files=glob.glob(os.path.join(dir6,'GPP_FluxSat_*.nc'))  
    files=sorted(files)
    
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(files[0])
#    lat.extend(ncf.variables['Latitude'][:])
#    lon.extend(ncf.variables['Longitude'][:])
    lat.extend(ncf.variables['lat'][:])
    lon.extend(ncf.variables['lon'][:])
         
    for file in files:
         ncf=nc.Dataset(file)
         data1.append(ncf.variables[varstr][:])
    
         ncf.close()
    
    data1=np.array(data1)
    data2=data1
    
    area = globarea(720,360,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
      
    data3 = data2[12*M:, 2*(90-latMax):2*(90-latMin), 2*(180+lonMin):2*(180+lonMax) ]
     
    
    
    mask=(data3== -9999.0)
    data3[mask]=np.NaN
     
    
    data3 = data3*365*1e-15/12.0
    
    data3=data3*area[2*(90-latMax):2*(90-latMin), 2*(180+lonMin):2*(180+lonMax) ]  
          
      
    
#    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
#    for i in range(0,N):
#        data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
#        
#        
#    plt.imshow(data4[0])
#    
#    
#    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
#    
#    a=signal.detrend(a)
#    mean = np.nanmean(a)
#    std = np.nanstd(a)
#    
#    b = (a-mean)/std
#
#    return a,b

    data4 = data3

    P = 0
    R1 = calc_seasonal_anomaly(data4*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data4*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data4*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data4*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data4*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   

 

#################################################################
#   GOSIF GPP v2 (Xiao Jingfeng)
#################################################################

dir7='/Volumes/HEWEI_T5/SIF_XiaoJF/GOSIF-GPP_v2/nc/GOSIF-GPP_v2/'
  
def flux_anomaly_gosifgpp(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    M=0
#    N=15     
#    varstr = 'GPP'
    
    files=glob.glob(os.path.join(dir7,'GOSIF_GPP_*.nc'))
    files=sorted(files)
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(files[0])
    lat.extend(ncf.variables['latitude'][:])
    lon.extend(ncf.variables['longitude'][:])
         
    for file in files:
         ncf=nc.Dataset(file)
         data1.append(ncf.variables[varstr][:])
    
         ncf.close()
    
    data1=np.array(data1)
    data2=data1
    
    area = globarea(720,360,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
      
    data3 = data2[10:-36, 2*(90-latMax):2*(90-latMin), 2*(180+lonMin):2*(180+lonMax)  ]
     
    mask=(data3==0)
    data3[mask]=np.NaN
     
    data3 = data3 *1e-15
    
    data3=data3*area[ 2*(90-latMax):2*(90-latMin), 2*(180+lonMin):2*(180+lonMax) ]  
          
      
   
    
#    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
#    for i in range(0,N):
#        data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
#        
#        
#    plt.imshow(data4[0])
#    
#    
#    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
#    
#    a=signal.detrend(a)
#    mean = np.nanmean(a)
#    std = np.nanstd(a)
#    
#    b = (a-mean)/std
#
#    return a,b
    
        
    data4 = data3

    P = 0
    R1 = calc_seasonal_anomaly(data4*regions_Europe_R1, M, P, N)    #(80, 100)
    R2 = calc_seasonal_anomaly(data4*regions_Europe_R2, M, P, N)
    R3 = calc_seasonal_anomaly(data4*regions_Europe_R3, M, P, N)
    R4 = calc_seasonal_anomaly(data4*regions_Europe_R4, M, P, N)
    RT = calc_seasonal_anomaly(data4*regions_Europe_all, M, P, N)

    return R1, R2, R3, R4, RT   



############################################
   
#def pipeline(latMin, latMax, lonMin, lonMax,flag, fn ) :    

fn = 'All_Flux_timeseries_Europe_AIMs_20221008'
 
    
latMin=33    #35
latMax=73    #70
lonMin=-15   #-10
lonMax=35    #40



N=15 
 
pydates1=[]
for i in range(0,N): 
    yr = 2001 + i
    tmp = datetime(yr, 1, 1)
    pydates1.append(tmp)  
    
#------------------------------


#M=1  
#N=15
#  
#iterm='NEP'
#varstr='ensemble_mean'
#trendy_nep = flux_anomaly_TRENDYv6(iterm, varstr,M, N, latMin, latMax, lonMin, lonMax) 

 
#------------------------------
 
M=1  
N=15 
  
varstr='bio_flux_opt'
cte2019nee =  flux_anomaly_cte2018(varstr,M, N, latMin, latMax, lonMin, lonMax)
     
varstr='fire_flux_imp'
cte2019fire = flux_anomaly_cte2018(varstr,M, N, latMin, latMax, lonMin, lonMax)
 
#------------------------------
     

M=1  
N=15 
  
varstr='bio_flux_opt'
ct2019nee =  flux_anomaly_ct2019b(varstr,M, N, latMin, latMax, lonMin, lonMax)
varstr='fire_flux_imp'
ct2019fire =  flux_anomaly_ct2019b(varstr,M, N, latMin, latMax, lonMin, lonMax)

 #------------------------------
       
M=22 
N=15 
  
varstr='flux_apos_bio'
camsnee =  flux_anomaly_cams(varstr,M, N, latMin, latMax, lonMin, lonMax)
  
    
#------------------------------
 
M=0 
N=15 
varstr='co2flux_land'
jena99v43nee =  flux_anomaly_jenas99v43(varstr,M, N, latMin, latMax, lonMin, lonMax)
 
 
M=0 
N=15 
varstr='co2flux_land'
jenaNEETv43nee =  flux_anomaly_jenaNEETv43(varstr,M, N, latMin, latMax, lonMin, lonMax)

     
#M=0 
#N=15 
#varstr='co2flux_land'
#jenaNEETWv43nee =  flux_anomaly_jenaNEETWv43(varstr,M, N, latMin, latMax, lonMin, lonMax)
 
 
  
dir15='/Volumes/HEWEI_T5/Inversions/Jena_v4.3_all/NEE-T-new/sEXT10/monthly/2001-2015/'
M=0 
N=15 
varstr='co2flux_land'
jenaNEETv43nee_sEXT10 =  flux_anomaly_jenaNEETv43(varstr,M, N, latMin, latMax, lonMin, lonMax)
 

  
dir16='/Volumes/HEWEI_T5/Inversions/Jena_v4.3_all/NEE-T-W-new/s99/monthly/2001-2015/'
M=0 
N=15 
varstr='co2flux_land'
jenaNEETWv43nee_s99 =  flux_anomaly_jenaNEETWv43(varstr,M, N, latMin, latMax, lonMin, lonMax)
   


dir16='/Volumes/HEWEI_T5/Inversions/Jena_v4.3_all/NEE-T-W-new/sEXT10/monthly/2001-2015/'
M=0 
N=15 
varstr='co2flux_land'
jenaNEETWv43nee_sEXT10 =  flux_anomaly_jenaNEETWv43(varstr,M, N, latMin, latMax, lonMin, lonMax)
   

#------------------------------

#M=0 
#N=15  
#varstr='NEE'
#ml = 'ANN' 
#fluxcomRSMETEOnee1 =  flux_anomaly_fluxcom_RSMETEO1(varstr,ml,M, N, latMin, latMax, lonMin, lonMax)
#
#M=0 
#N=15  
#varstr='NEE'
#ml = 'ANN' 
#fluxcomRSMETEOnee2 =  flux_anomaly_fluxcom_RSMETEO2(varstr,ml,M, N, latMin, latMax, lonMin, lonMax)
     


#------------------------------


M=0 
N=10
#M=4 
#N=6

latMin=33    #35
latMax=73    #70
lonMin=-14.5 #-10
lonMax=35    #40

   
dataset = 'EUROCOM_CarbonTrackerEurope_2006-2015.nc'
varstr='bio_posterior'
eurocom1nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
varstr='fire'
eurocom1fire=  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
  

varstr='bio_posterior'
dataset = 'EUROCOM_EnKF-RAMS_2006-2015.nc'
eurocom2nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)


varstr='bio_posterior' 
dataset = 'EUROCOM_FLEXINVERT_2006-2015.nc'
eurocom3nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
   
#M=4 
#N=6

varstr='bio_posterior'
dataset = 'EUROCOM_JenaCarboScopeRegional_2006-2015.nc'
eurocom4nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)

   
dataset = 'EUROCOM_LUMIA_2006-2015.nc'     
varstr='bio_posterior'
eurocom5nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
varstr='fire'
eurocom5fire =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)

   
varstr='bio_posterior'
    
dataset = 'EUROCOM_PYVAR-CHIMERE_2006-2015.nc'
eurocom6nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)



M=0 
N=5
varstr='bio_posterior'
        
dataset = 'EUROCOM_NAME-HB_2011-2015.nc'
eurocom7nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
    
 
#--------------------------------------   
latMin=33    #35
latMax=73    #70
lonMin=-14.5 #-10
lonMax=35    #40

lonMin=-15
 
M=0  
N=6 
  
dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
varstr='ensemble_mean'

ByrneGOSATnee_ens =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
    
#dataset='Byrne2020_NEE_GOSAT_only'
#varstr='ensemble_mean'
#
#ByrneGOSATnee_casa =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
  
#    dataset='Byrne2020_NEE_surface_only'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_fc =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#       
#    dataset='Byrne2020_NEE_TCCON_only'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_sib =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)



#------------------------------

M=0 
N=6

dataset='CMS-Flux-NBE-2020.monthly.grid'           
varstr = 'post-NBE' 
LiuGOSATnee =  flux_anomaly_LiuGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)


 
#------------------------------

M=0 
N=6
 
#dataset='CMS-Flux-NBE-2020.monthly.grid'        
varstr = 'bio_monthly_opt'
JiangGOSATnee =  flux_anomaly_JiangGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
  

 
#--------------------------------------   
 
M=0 
N=6
 
dataset='ccdas_sm+vod'        
varstr = 'nep'
ScholzeCCDASnee =  flux_anomaly_ScholzeCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
   
#dataset ='CCDAS_6yr_2deg_assim_results'
#varstr = 'nep'
#ScholzeCCDASnee =  flux_anomaly_WuCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
 
   
#--------------------------------------   

M=0      
N=15     
varstr='tair'
xx1_tair =  anomaly_CRUNCEP(varstr,M, N,latMin, latMax, lonMin, lonMax)


M=0      
N=15   
varstr='rain'
xx1_p =  anomaly_CRUNCEP(varstr,M, N,latMin, latMax, lonMin, lonMax)

 
 
M=0      
N=15   
varstr='vpd'
xx1_vpd =  anomaly_CRUNCEP(varstr,M, N,latMin, latMax, lonMin, lonMax)

 

M=0    
N=15   
varstr='SMroot'
xx1_smroot =  anomaly_gleam(varstr,M, N,latMin, latMax, lonMin, lonMax)  
   
#M=0   
#N=15    
#           
#varstr = 'LAI'
#xx1_lai = anomaly_LAI(varstr,M, N,latMin,latMax,lonMin,lonMax)
    
M=0   
N=15    
 
xx1_modisnirv = anomaly_modisNIRv(M, N,latMin,latMax,lonMin,lonMax)
 
    
  
M=0   
N=15   
varstr='SIF'
xx1_csif  =  anomaly_CSIF(varstr,M, N,latMin, latMax, lonMin, lonMax)
       

M=0
N=15     
varstr = 'GPP'
xx1_fluxsat =  flux_anomaly_fluxsat(varstr,M, N, latMin, latMax, lonMin, lonMax)

M=0
N=15     
varstr = 'GPP'
xx1_gosifgpp =  flux_anomaly_gosifgpp(varstr,M, N, latMin, latMax, lonMin, lonMax)
          
 
 

######################################################################

data=[]
        

for Reg in range(1,5+1):
   
    
    region = Reg-1   
    flag=1
          
    
    #######################################
    # raw fluxes
    #######################################
 
    cte2019neeann = []
    ct2019neeann = []
    camsneeann=[]
    jena99v43neeann=[]
    jenaNEETv43neeann=[]
    jenaNEETv43neeann_sEXT10=[]
    jenaNEETWv43neeann_s99=[]
    jenaNEETWv43neeann_sEXT10=[]
         

    
    eurocom1neeann=[]
    eurocom2neeann=[]
    eurocom3neeann=[]
    eurocom4neeann=[]
    eurocom5neeann=[]
    eurocom6neeann=[] 
    eurocom7neeann=[]
 
#    fluxcomRSMETEOneeann1=[]
#    fluxcomRSMETEOneeann2=[]
    
    ByrneGOSATneeann_ens=[]
    LiuGOSATneeann=[]
    JiangGOSATneeann=[]
    ScholzeCCDASneeann=[]
    
    tairann=[]
    prepann=[]
    vpdann=[]
    smrootann=[]
    nirvann=[]   
    csifann=[]  
    fluxsatann=[]
    gosifgppann=[]
    
    
    
    for i in range(0, 10):
#    for i in range(0, 6):
       eurocom1neeann.append(np.nanmean( eurocom1nee[region][flag][i*12:(i+1)*12] + eurocom1fire[region][flag][i*12:(i+1)*12] ) )
       eurocom2neeann.append(np.nanmean( eurocom2nee[region][flag][i*12:(i+1)*12] ) )
       eurocom3neeann.append(np.nanmean( eurocom3nee[region][flag][i*12:(i+1)*12] ) )
       eurocom4neeann.append(np.nanmean( eurocom4nee[region][flag][i*12:(i+1)*12] ) )
       eurocom5neeann.append(np.nanmean( eurocom5nee[region][flag][i*12:(i+1)*12] + eurocom5fire[region][flag][i*12:(i+1)*12] ) )
       eurocom6neeann.append(np.nanmean( eurocom6nee[region][flag][i*12:(i+1)*12] ) )
    
    for i in range(0, 15):
#       fluxcomRSMETEOneeann1.append(np.nanmean( fluxcomRSMETEOnee1[region][flag][i*12:(i+1)*12] ) )
#       fluxcomRSMETEOneeann2.append(np.nanmean( fluxcomRSMETEOnee2[region][flag][i*12:(i+1)*12] ) )
 
       tairann.append(np.nanmean( xx1_tair[region][flag][i*12:(i+1)*12] ) )
       prepann.append(np.nanmean( xx1_p[region][flag][i*12:(i+1)*12] ) )
       vpdann.append(np.nanmean( xx1_vpd[region][flag][i*12:(i+1)*12] ) )
       smrootann.append(np.nanmean( xx1_smroot[region][flag][i*12:(i+1)*12] ) )
       nirvann.append(np.nanmean( xx1_modisnirv[region][flag][i*12:(i+1)*12] ) )
       csifann.append(np.nanmean( xx1_csif[region][flag][i*12:(i+1)*12] ) )
       fluxsatann.append(np.nanmean( xx1_fluxsat[region][flag][i*12:(i+1)*12] ) )
       gosifgppann.append(np.nanmean( xx1_gosifgpp[region][flag][i*12:(i+1)*12] ) )
      
        
#    for i in range(0, 7):  
    for i in range(0, 5):
       eurocom7neeann.append(np.nanmean( eurocom7nee[region][flag][i*12:(i+1)*12] ) )
      
    for i in range(0, 6):
       ByrneGOSATneeann_ens.append(np.nanmean( ByrneGOSATnee_ens[region][flag][i*12:(i+1)*12] ) )
       LiuGOSATneeann.append(np.nanmean( LiuGOSATnee[region][flag][i*12:(i+1)*12] ) )
       JiangGOSATneeann.append(np.nanmean( JiangGOSATnee[region][flag][i*12:(i+1)*12] ) )
       ScholzeCCDASneeann.append(np.nanmean( ScholzeCCDASnee[region][flag][i*12:(i+1)*12] ) )
 
  
    eurocom1neeann=np.array(eurocom1neeann)
    eurocom2neeann=np.array(eurocom2neeann)
    eurocom3neeann=np.array(eurocom3neeann)
    eurocom4neeann=np.array(eurocom4neeann)
    eurocom5neeann=np.array(eurocom5neeann)
    eurocom6neeann=np.array(eurocom6neeann)
    eurocom7neeann=np.array(eurocom7neeann)
 
#    fluxcomRSMETEOneeann1=np.array(fluxcomRSMETEOneeann1)  
#    fluxcomRSMETEOneeann2=np.array(fluxcomRSMETEOneeann2)  
        
    ByrneGOSATneeann_ens=np.array(ByrneGOSATneeann_ens)
    LiuGOSATneeann=np.array(LiuGOSATneeann)
    JiangGOSATneeann=np.array(JiangGOSATneeann)
    ScholzeCCDASneeann=np.array(ScholzeCCDASneeann)
 
     
    tairann=np.array(tairann)
    prepann=np.array(prepann)
    vpdann=np.array(vpdann)
    smrootann=np.array(smrootann)
    nirvann=np.array(nirvann)
    csifann=np.array(csifann)
    fluxsatann=np.array(fluxsatann)
    gosifgppann=np.array(gosifgppann)
      
  
    
    
    for i in range(0, 15): 
    #       cte2016neeann.append(np.nanmean( cte2016nee[flag][i*12:(i+1)*12] + cte2016fire[flag][i*12:(i+1)*12] ) )   
       cte2019neeann.append(np.nanmean( cte2019nee[region][flag][i*12:(i+1)*12] + cte2019fire[region][flag][i*12:(i+1)*12]) )
    #       ct2016neeann.append(np.nanmean( ct2016nee[flag][i*12:(i+1)*12] + ct2016fire[flag][i*12:(i+1)*12]) )
       ct2019neeann.append(np.nanmean( ct2019nee[region][flag][i*12:(i+1)*12] + ct2019fire[region][flag][i*12:(i+1)*12]) )
       camsneeann.append(np.nanmean( camsnee[region][flag][i*12:(i+1)*12] ) )
    #       mactmneeann.append(np.nanmean( mactmnee[flag][i*12:(i+1)*12] ) )
    #       jena99neeann.append(np.nanmean( jena99nee[flag][i*12:(i+1)*12] ) )
       jena99v43neeann.append(np.nanmean( jena99v43nee[region][flag][i*12:(i+1)*12] ) )
       jenaNEETv43neeann.append(np.nanmean( jenaNEETv43nee[region][flag][i*12:(i+1)*12] ) )
#       jenaNEETWv43neeann.append(np.nanmean( jenaNEETWv43nee[region][flag][i*12:(i+1)*12] ) )  
       jenaNEETv43neeann_sEXT10.append(np.nanmean( jenaNEETv43nee_sEXT10[region][flag][i*12:(i+1)*12] ) )
       jenaNEETWv43neeann_s99.append(np.nanmean( jenaNEETWv43nee_s99[region][flag][i*12:(i+1)*12] ) )
       jenaNEETWv43neeann_sEXT10.append(np.nanmean( jenaNEETWv43nee_sEXT10[region][flag][i*12:(i+1)*12] ) )
  
 
    #    cte2016neeann=np.array(cte2016neeann)
    cte2019neeann=np.array(cte2019neeann)
    #    ct2016neeann=np.array(ct2016neeann)
    ct2019neeann=np.array(ct2019neeann)
    camsneeann=np.array(camsneeann)
    #mactmneeann=np.array(mactmneeann) 
    #jena99neeann=np.array(jena99neeann)*40  #/8000.0
    #jena04neeann=np.array(jena04neeann)*10  #/40000.0
    jena99v43neeann=np.array(jena99v43neeann)  #*5
    #jena04v43neeann=np.array(jena04v43neeann)*10 
    jenaNEETv43neeann=np.array(jenaNEETv43neeann)   #*20
    jenaNEETv43neeann_sEXT10=np.array(jenaNEETv43neeann_sEXT10)
    jenaNEETWv43neeann_s99=np.array(jenaNEETWv43neeann_s99) 
    jenaNEETWv43neeann_sEXT10=np.array(jenaNEETWv43neeann_sEXT10) 
               
    

    
#    trendy_nee =  -np.array(trendy_nep[region][flag])  
 
#    trendy_nee_ann=[]
   
#    N=15
#    for i in range(0, N):  
#       trendy_nee_ann.append(np.nanmean( trendy_nee[i*12:(i+1)*12] ) )
#  
#    trendy_nee_ann=np.array(trendy_nee_ann)
 
    #detrend
 
    #    cte2016neeann = signal.detrend(cte2016neeann)
    cte2019neeann = signal.detrend(cte2019neeann)
    #    ct2016neeann = signal.detrend(ct2016neeann)
    ct2019neeann = signal.detrend(ct2019neeann)
    camsneeann = signal.detrend(camsneeann)
    #    mactmneeann = signal.detrend(mactmneeann)
    #    jena99neeann = signal.detrend(jena99neeann)
    jena99v43neeann = signal.detrend(jena99v43neeann)
    #    jena04neeann = signal.detrend(jena04neeann)
    #    jena04v43neeann = signal.detrend(jena04v43neeann)                 
    jenaNEETv43neeann =signal.detrend(jenaNEETv43neeann) 
#    jenaNEETWv43neeann =signal.detrend(jenaNEETWv43neeann) 
     
    jenaNEETv43neeann_sEXT10 =signal.detrend(jenaNEETv43neeann_sEXT10) 
    jenaNEETWv43neeann_s99 =signal.detrend(jenaNEETWv43neeann_s99) 
    jenaNEETWv43neeann_sEXT10 =signal.detrend(jenaNEETWv43neeann_sEXT10) 
         
   
#    trendy_nee_ann =signal.detrend(trendy_nee_ann) 
    



     
     
    eurocom1neeann=signal.detrend(eurocom1neeann)
    eurocom2neeann=signal.detrend(eurocom2neeann)
    eurocom3neeann=signal.detrend(eurocom3neeann)
    eurocom4neeann=signal.detrend(eurocom4neeann)
    eurocom5neeann=signal.detrend(eurocom5neeann)
    eurocom6neeann=signal.detrend(eurocom6neeann)
    eurocom7neeann=signal.detrend(eurocom7neeann)
 
#    fluxcomRSMETEOneeann1=signal.detrend(fluxcomRSMETEOneeann1)  
#    fluxcomRSMETEOneeann2=signal.detrend(fluxcomRSMETEOneeann2)  
        
    ByrneGOSATneeann_ens=signal.detrend(ByrneGOSATneeann_ens)
    LiuGOSATneeann=signal.detrend(LiuGOSATneeann)
    JiangGOSATneeann=signal.detrend(JiangGOSATneeann)
    ScholzeCCDASneeann=signal.detrend(ScholzeCCDASneeann)
 
     
    tairann=signal.detrend(tairann)
    prepann=signal.detrend(prepann)
    vpdann=signal.detrend(vpdann)
    smrootann=signal.detrend(smrootann)
    nirvann=signal.detrend(nirvann)
    csifann=signal.detrend(csifann)
    fluxsatann=signal.detrend(fluxsatann)
    gosifgppann=signal.detrend(gosifgppann)
      
   
     
    #    cte2016neeann = cte2016neeann/np.std(cte2016neeann)
    #    cte2018neeann = cte2018neeann/np.std(cte2018neeann)
    #    ct2016neeann = ct2016neeann/np.std(ct2016neeann)
    #    ct2017neeann = ct2017neeann/np.std(ct2017neeann)
    #    camsneeann = camsneeann/np.std(camsneeann)
    #    jena99v43neeann = jena99v43neeann/np.std(jena99v43neeann)
    #    jenaNEETv43neeann = jenaNEETv43neeann/np.std(jenaNEETv43neeann)
    
    
    
    
    ##############################################################################################################################
     
    XXX=[]
     
     
    #tair
    X1=pearsonr(-cte2019neeann,tairann)
    X2=pearsonr(-ct2019neeann,tairann)
    X3=pearsonr(-camsneeann,tairann)
    X4=pearsonr(-jena99v43neeann,tairann)
    X5=pearsonr(-jenaNEETv43neeann,tairann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,tairann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,tairann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,tairann)
 
    X9=pearsonr(-eurocom1neeann,tairann[5:])
    X10=pearsonr(-eurocom2neeann,tairann[5:])
    X11=pearsonr(-eurocom3neeann,tairann[5:])
    X12=pearsonr(-eurocom4neeann,tairann[5:])
    X13=pearsonr(-eurocom5neeann,tairann[5:])
    X14=pearsonr(-eurocom6neeann,tairann[5:])   
    X15=pearsonr(-eurocom7neeann,tairann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,tairann[9:])
    X17=pearsonr(-LiuGOSATneeann,tairann[9:])
    X18=pearsonr(-JiangGOSATneeann,tairann[9:])
    X19=pearsonr(ScholzeCCDASneeann,tairann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0])
    
    
    
    
    #prep
    X1=pearsonr(-cte2019neeann,prepann)
    X2=pearsonr(-ct2019neeann,prepann)
    X3=pearsonr(-camsneeann,prepann)
    X4=pearsonr(-jena99v43neeann,prepann)
    X5=pearsonr(-jenaNEETv43neeann,prepann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,prepann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,prepann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,prepann)
 
    X9=pearsonr(-eurocom1neeann,prepann[5:])
    X10=pearsonr(-eurocom2neeann,prepann[5:])
    X11=pearsonr(-eurocom3neeann,prepann[5:])
    X12=pearsonr(-eurocom4neeann,prepann[5:])
    X13=pearsonr(-eurocom5neeann,prepann[5:])
    X14=pearsonr(-eurocom6neeann,prepann[5:])   
    X15=pearsonr(-eurocom7neeann,prepann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,prepann[9:])
    X17=pearsonr(-LiuGOSATneeann,prepann[9:])
    X18=pearsonr(-JiangGOSATneeann,prepann[9:])
    X19=pearsonr(ScholzeCCDASneeann,prepann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0])
    
    
    
    
    #vpd
    X1=pearsonr(-cte2019neeann,vpdann)
    X2=pearsonr(-ct2019neeann,vpdann)
    X3=pearsonr(-camsneeann,vpdann)
    X4=pearsonr(-jena99v43neeann,vpdann)
    X5=pearsonr(-jenaNEETv43neeann,vpdann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,vpdann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,vpdann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,vpdann)
 
    X9=pearsonr(-eurocom1neeann,vpdann[5:])
    X10=pearsonr(-eurocom2neeann,vpdann[5:])
    X11=pearsonr(-eurocom3neeann,vpdann[5:])
    X12=pearsonr(-eurocom4neeann,vpdann[5:])
    X13=pearsonr(-eurocom5neeann,vpdann[5:])
    X14=pearsonr(-eurocom6neeann,vpdann[5:])   
    X15=pearsonr(-eurocom7neeann,vpdann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,vpdann[9:])
    X17=pearsonr(-LiuGOSATneeann,vpdann[9:])
    X18=pearsonr(-JiangGOSATneeann,vpdann[9:])
    X19=pearsonr(ScholzeCCDASneeann,vpdann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0])
    
       
    
    #smroot
    X1=pearsonr(-cte2019neeann,smrootann)
    X2=pearsonr(-ct2019neeann,smrootann)
    X3=pearsonr(-camsneeann,smrootann)
    X4=pearsonr(-jena99v43neeann,smrootann)
    X5=pearsonr(-jenaNEETv43neeann,smrootann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,smrootann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,smrootann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,smrootann)
 
    X9=pearsonr(-eurocom1neeann,smrootann[5:])
    X10=pearsonr(-eurocom2neeann,smrootann[5:])
    X11=pearsonr(-eurocom3neeann,smrootann[5:])
    X12=pearsonr(-eurocom4neeann,smrootann[5:])
    X13=pearsonr(-eurocom5neeann,smrootann[5:])
    X14=pearsonr(-eurocom6neeann,smrootann[5:])   
    X15=pearsonr(-eurocom7neeann,smrootann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,smrootann[9:])
    X17=pearsonr(-LiuGOSATneeann,smrootann[9:])
    X18=pearsonr(-JiangGOSATneeann,smrootann[9:])
    X19=pearsonr(ScholzeCCDASneeann,smrootann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0]) 
    
    
    
        
       
    
    #lai
    X1=pearsonr(-cte2019neeann,nirvann)
    X2=pearsonr(-ct2019neeann,nirvann)
    X3=pearsonr(-camsneeann,nirvann)
    X4=pearsonr(-jena99v43neeann,nirvann)
    X5=pearsonr(-jenaNEETv43neeann,nirvann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,nirvann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,nirvann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,nirvann)
 
    X9=pearsonr(-eurocom1neeann,nirvann[5:])
    X10=pearsonr(-eurocom2neeann,nirvann[5:])
    X11=pearsonr(-eurocom3neeann,nirvann[5:])
    X12=pearsonr(-eurocom4neeann,nirvann[5:])
    X13=pearsonr(-eurocom5neeann,nirvann[5:])
    X14=pearsonr(-eurocom6neeann,nirvann[5:])   
    X15=pearsonr(-eurocom7neeann,nirvann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,nirvann[9:])
    X17=pearsonr(-LiuGOSATneeann,nirvann[9:])
    X18=pearsonr(-JiangGOSATneeann,nirvann[9:])
    X19=pearsonr(ScholzeCCDASneeann,nirvann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0]) 
    
    
    
        
       
    
    #csif
    X1=pearsonr(-cte2019neeann,csifann)
    X2=pearsonr(-ct2019neeann,csifann)
    X3=pearsonr(-camsneeann,csifann)
    X4=pearsonr(-jena99v43neeann,csifann)
    X5=pearsonr(-jenaNEETv43neeann,csifann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,csifann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,csifann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,csifann)
 
    X9=pearsonr(-eurocom1neeann,csifann[5:])
    X10=pearsonr(-eurocom2neeann,csifann[5:])
    X11=pearsonr(-eurocom3neeann,csifann[5:])
    X12=pearsonr(-eurocom4neeann,csifann[5:])
    X13=pearsonr(-eurocom5neeann,csifann[5:])
    X14=pearsonr(-eurocom6neeann,csifann[5:])   
    X15=pearsonr(-eurocom7neeann,csifann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,csifann[9:])
    X17=pearsonr(-LiuGOSATneeann,csifann[9:])
    X18=pearsonr(-JiangGOSATneeann,csifann[9:])
    X19=pearsonr(ScholzeCCDASneeann,csifann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0]) 
    
    
    
        
       
    
    #fluxsat
    X1=pearsonr(-cte2019neeann,fluxsatann)
    X2=pearsonr(-ct2019neeann,fluxsatann)
    X3=pearsonr(-camsneeann,fluxsatann)
    X4=pearsonr(-jena99v43neeann,fluxsatann)
    X5=pearsonr(-jenaNEETv43neeann,fluxsatann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,fluxsatann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,fluxsatann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,fluxsatann)
 
    X9=pearsonr(-eurocom1neeann,fluxsatann[5:])
    X10=pearsonr(-eurocom2neeann,fluxsatann[5:])
    X11=pearsonr(-eurocom3neeann,fluxsatann[5:])
    X12=pearsonr(-eurocom4neeann,fluxsatann[5:])
    X13=pearsonr(-eurocom5neeann,fluxsatann[5:])
    X14=pearsonr(-eurocom6neeann,fluxsatann[5:])   
    X15=pearsonr(-eurocom7neeann,fluxsatann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,fluxsatann[9:])
    X17=pearsonr(-LiuGOSATneeann,fluxsatann[9:])
    X18=pearsonr(-JiangGOSATneeann,fluxsatann[9:])
    X19=pearsonr(ScholzeCCDASneeann,fluxsatann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0]) 
    
    
    
        
       
    
    #gosifgpp
    X1=pearsonr(-cte2019neeann,gosifgppann)
    X2=pearsonr(-ct2019neeann,gosifgppann)
    X3=pearsonr(-camsneeann,gosifgppann)
    X4=pearsonr(-jena99v43neeann,gosifgppann)
    X5=pearsonr(-jenaNEETv43neeann,gosifgppann)
    X6=pearsonr(-jenaNEETv43neeann_sEXT10,gosifgppann)
    X7=pearsonr(-jenaNEETWv43neeann_s99,gosifgppann)
    X8=pearsonr(-jenaNEETWv43neeann_sEXT10,gosifgppann)
 
    X9=pearsonr(-eurocom1neeann,gosifgppann[5:])
    X10=pearsonr(-eurocom2neeann,gosifgppann[5:])
    X11=pearsonr(-eurocom3neeann,gosifgppann[5:])
    X12=pearsonr(-eurocom4neeann,gosifgppann[5:])
    X13=pearsonr(-eurocom5neeann,gosifgppann[5:])
    X14=pearsonr(-eurocom6neeann,gosifgppann[5:])   
    X15=pearsonr(-eurocom7neeann,gosifgppann[10:])
   
    X16=pearsonr(-ByrneGOSATneeann_ens,gosifgppann[9:])
    X17=pearsonr(-LiuGOSATneeann,gosifgppann[9:])
    X18=pearsonr(-JiangGOSATneeann,gosifgppann[9:])
    X19=pearsonr(ScholzeCCDASneeann,gosifgppann[9:])
  
    XXX.append(X1[0])
    XXX.append(X2[0])
    XXX.append(X3[0])
    XXX.append(X4[0])
    XXX.append(X5[0])
    XXX.append(X6[0])
    XXX.append(X7[0])
    XXX.append(X8[0])
    XXX.append(X9[0])
    XXX.append(X10[0])
    XXX.append(X11[0])
    XXX.append(X12[0])
    XXX.append(X13[0])
    XXX.append(X14[0])
    XXX.append(X15[0])
    XXX.append(X16[0])
    XXX.append(X17[0])
    XXX.append(X18[0])
    XXX.append(X19[0]) 
    
    
    #XXX = np.array(XXX)
    XXX=np.round(XXX, 2) 
    XXX=list(XXX) 
    print(XXX)
  
 
#    return data    
#
#   
#import timeit
#start = timeit.timeit()
#
#flag=1
#data = pipeline(35, 70, -10, 40,flag,'All_Flux_timeseries_Europe_AIMs_20201112_newextent_new_correct_S3_fire')
# 
# 
#end = timeit.timeit()
#print end - start
 
#R1= np.array([0.05, -0.03, 0.5, -0.09, -0.34, -0.04, -0.21, 0.17, 0.25, -0.24, -0.2, 0.61, -0.15, -0.15, 0.06, -0.21, -0.21, 0.16, 0.02, -0.53, -0.71, -0.35, -0.43, 0.3, 0.61, 0.35, 0.17, -0.59, -0.69, -0.61, 0.16, -0.45, -0.8, 0.84, -0.8, -0.82, -0.66, -0.67, 0.56, 0.37, 0.53, 0.2, -0.55, -0.56, -0.48, -0.15, 0.65, 0.25, 0.23, 0.5, 0.25, 0.18, -0.46, 0.31, 0.22, 0.49, 0.3, -0.39, -0.33, -0.42, -0.33, 0.32, 0.44, 0.4, -0.1, -0.62, -0.19, -0.29, -0.22, -0.11, -0.67, 0.49, -0.46, -0.49, -0.5, -0.47, 0.19, 0.0, 0.33, 0.21, -0.2, 0.01, -0.04, 0.34, 0.18, -0.3, -0.15, 0.69, 0.28, -0.2, 0.02, -0.31, -0.18, 0.3, 0.19, 0.4, 0.34, 0.32, 0.06, -0.13, -0.04, 0.01, 0.19, 0.29, 0.1, 0.1, 0.78, 0.32, -0.19, -0.46, 0.14, 0.11, 0.51, 0.3, 0.39, 0.13, 0.01, -0.23, 0.17, 0.26, 0.37, 0.33, 0.03, 0.09, 0.1, 0.76, 0.18, -0.46, -0.36, 0.13, -0.05, 0.33, 0.07, 0.64, 0.44, 0.49, 0.06, -0.28, -0.34, -0.21, 0.01, 0.45, 0.13, 0.15, 0.67, 0.56, 0.07, -0.66, 0.35, 0.37, 0.73, 0.55])
#R2= np.array([-0.45, 0.25, -0.12, 0.18, 0.19, 0.21, 0.28, 0.13, 0.27, -0.56, 0.12, 0.66, -0.04, 0.36, -0.68, -0.17, -0.23, 0.56, -0.42, -0.21, -0.39, -0.04, -0.51, 0.26, 0.23, 0.47, 0.45, -0.15, -0.15, 0.05, 0.18, 0.38, -0.0, 0.28, -0.41, 0.08, -0.23, -0.34, -0.38, 0.39, -0.46, 0.61, -0.55, -0.49, -0.64, -0.66, -0.29, -0.25, -0.19, -0.27, -0.25, -0.42, -0.63, 0.04, -0.24, 0.64, -0.15, 0.16, -0.37, 0.14, -0.6, 0.37, 0.36, 0.43, 0.5, 0.21, 0.11, 0.31, 0.31, 0.61, 0.33, 0.4, -0.04, 0.15, -0.52, -0.05, -0.08, -0.06, 0.23, -0.36, 0.6, 0.59, 0.73, 0.66, 0.57, -0.54, 0.23, 0.86, 0.12, 0.69, -0.61, -0.28, -0.3, 0.45, -0.53, -0.01, -0.02, 0.22, -0.39, 0.71, 0.71, 0.81, 0.76, 0.55, -0.33, 0.4, 0.89, 0.23, 0.74, -0.65, -0.23, -0.32, 0.5, -0.54, -0.09, 0.14, 0.22, -0.3, 0.68, 0.66, 0.76, 0.69, 0.55, -0.36, 0.26, 0.91, 0.14, 0.71, -0.54, -0.12, -0.23, 0.46, -0.44, -0.05, 0.05, 0.24, -0.23, 0.78, 0.76, 0.83, 0.76, 0.53, -0.2, 0.55, 0.87, 0.26, 0.73, -0.71, -0.06, -0.4, 0.6, -0.52])
#R3= np.array([0.26, 0.27, -0.1, -0.24, 0.32, 0.32, 0.38, 0.4, 0.6, -0.53, 0.09, 0.8, -0.06, -0.3, 0.19, -0.1, 0.08, 0.62, -0.38, -0.15, -0.45, -0.34, -0.24, 0.35, 0.38, 0.43, 0.21, -0.73, 0.53, 0.25, -0.41, 0.41, -0.28, 0.52, 0.45, -0.27, -0.6, 0.04, 0.19, 0.22, 0.07, -0.02, -0.14, -0.21, -0.41, -0.1, 0.64, -0.66, -0.32, 0.31, -0.09, 0.36, -0.89, -0.7, -0.16, 0.3, -0.17, 0.04, -0.27, -0.15, -0.15, 0.24, 0.33, 0.32, 0.13, -0.64, 0.77, 0.29, -0.26, 0.02, -0.34, 0.73, 0.58, 0.21, -0.19, 0.16, 0.27, 0.05, -0.26, -0.38, 0.4, 0.53, 0.74, 0.58, 0.27, -0.27, 0.29, 0.74, 0.07, -0.4, 0.6, 0.27, 0.19, 0.35, -0.26, 0.26, 0.03, -0.27, -0.29, 0.57, 0.68, 0.88, 0.69, 0.33, 0.02, 0.53, 0.69, 0.25, -0.38, 0.56, 0.23, -0.01, 0.57, -0.57, 0.32, 0.12, -0.23, -0.2, 0.69, 0.76, 0.9, 0.76, 0.3, 0.1, 0.56, 0.67, 0.25, -0.4, 0.56, 0.22, -0.02, 0.63, -0.61, 0.3, -0.01, -0.18, -0.08, 0.61, 0.7, 0.94, 0.73, 0.35, 0.19, 0.73, 0.69, 0.33, -0.28, 0.57, 0.23, 0.0, 0.67, -0.62])
#R4= np.array([0.33, -0.08, -0.37, -0.34, 0.08, -0.11, -0.03, 0.24, -0.0, 0.31, 0.38, 0.07, 0.34, 0.35, 0.11, 0.1, -0.33, 0.17, -0.39, -0.42, -0.37, 0.03, 0.35, 0.24, 0.01, 0.12, -0.09, -0.76, -0.63, 0.31, 0.26, 0.38, 0.36, -0.11, 0.9, 0.13, -0.18, 0.16, -0.05, 0.12, -0.5, -0.29, -0.6, -0.51, -0.45, -0.13, 0.48, 0.48, -0.14, -0.66, -0.36, -0.13, -0.59, -0.8, -0.26, -0.36, -0.37, -0.06, -0.04, 0.07, 0.55, 0.08, -0.18, -0.25, -0.48, -0.37, -0.4, 0.15, 0.62, 0.56, 0.32, 0.69, 0.6, 0.26, 0.36, 0.56, 0.17, 0.05, -0.01, 0.07, 0.57, 0.2, 0.27, 0.11, -0.28, -0.39, 0.07, 0.85, 0.41, 0.02, 0.59, 0.69, 0.48, 0.42, 0.51, 0.27, 0.14, 0.06, 0.05, 0.62, 0.27, 0.3, 0.18, -0.27, -0.28, 0.09, 0.8, 0.33, -0.01, 0.7, 0.65, 0.45, 0.53, 0.51, 0.22, 0.28, -0.1, 0.06, 0.49, 0.07, 0.11, 0.04, -0.19, -0.24, 0.05, 0.84, 0.32, -0.06, 0.64, 0.6, 0.56, 0.48, 0.56, 0.4, 0.21, 0.12, 0.1, 0.65, 0.35, 0.32, 0.22, -0.2, -0.1, 0.19, 0.8, 0.27, -0.02, 0.83, 0.52, 0.41, 0.71, 0.45])
#RT= np.array([0.28, 0.33, 0.18, -0.13, 0.22, 0.2, 0.24, 0.32, 0.67, -0.42, 0.06, 0.76, -0.08, 0.12, 0.86, -0.05, 0.03, 0.62, -0.19, -0.1, -0.32, -0.07, -0.25, 0.57, 0.38, 0.61, 0.38, -0.6, 0.31, 0.42, -0.1, 0.74, 0.2, -0.38, 0.6, -0.14, -0.52, -0.18, 0.01, 0.19, -0.34, 0.1, -0.36, -0.29, -0.45, -0.19, 0.48, -0.58, -0.14, 0.37, -0.4, 0.02, 0.44, -0.48, 0.05, 0.41, -0.14, -0.08, -0.12, -0.04, -0.06, 0.41, 0.26, 0.36, 0.16, -0.37, 0.58, 0.23, -0.02, 0.38, -0.1, 0.43, 0.72, 0.13, -0.13, 0.2, 0.15, 0.31, 0.14, -0.17, 0.51, 0.44, 0.58, 0.5, 0.47, -0.26, 0.16, 0.91, 0.09, -0.0, 0.74, 0.32, 0.2, 0.54, -0.09, 0.25, 0.25, 0.1, -0.22, 0.7, 0.57, 0.72, 0.62, 0.34, 0.03, 0.36, 0.83, 0.21, 0.01, 0.74, 0.49, 0.11, 0.61, -0.18, 0.24, 0.32, 0.05, -0.18, 0.71, 0.57, 0.69, 0.62, 0.31, 0.07, 0.35, 0.83, 0.21, -0.04, 0.75, 0.56, 0.15, 0.61, -0.15, 0.32, 0.23, 0.22, -0.07, 0.68, 0.53, 0.7, 0.59, 0.32, 0.11, 0.48, 0.82, 0.19, 0.08, 0.84, 0.47, 0.09, 0.75, -0.19])


#R1= np.array([0.05, -0.03, 0.5, -0.09, -0.34, -0.04, -0.21, 0.17, 0.25, -0.24, -0.2, 0.61, -0.15, -0.15, 0.06, -0.21, -0.21, 0.16, 0.14, -0.53, -0.71, -0.35, -0.43, 0.3, 0.61, 0.35, 0.17, -0.59, -0.69, -0.61, 0.16, -0.45, -0.8, 0.84, -0.8, -0.82, -0.57, 0.61, 0.56, 0.37, 0.53, 0.2, -0.55, -0.56, -0.48, -0.15, 0.65, 0.25, 0.23, 0.5, 0.25, 0.18, -0.46, 0.31, 0.22, 0.51, -0.21, -0.39, -0.33, -0.42, -0.33, 0.32, 0.44, 0.4, -0.1, -0.62, -0.19, -0.29, -0.22, -0.11, -0.67, 0.49, -0.46, -0.49, -0.39, 0.3, 0.19, 0.0, 0.33, 0.21, -0.2, 0.01, -0.04, 0.34, 0.18, -0.3, -0.15, 0.69, 0.28, -0.2, 0.02, -0.31, -0.18, 0.22, 0.15, 0.4, 0.34, 0.32, 0.06, -0.13, -0.04, 0.01, 0.19, 0.29, 0.1, 0.1, 0.78, 0.32, -0.19, -0.46, 0.14, 0.11, 0.54, -0.2, 0.39, 0.13, 0.01, -0.23, 0.17, 0.26, 0.37, 0.33, 0.03, 0.09, 0.1, 0.76, 0.18, -0.46, -0.36, 0.13, -0.05, 0.51, -0.19, 0.64, 0.44, 0.49, 0.06, -0.28, -0.34, -0.21, 0.01, 0.45, 0.13, 0.15, 0.67, 0.56, 0.07, -0.66, 0.35, 0.37, 0.73, 0.42])
#R2= np.array([-0.45, 0.25, -0.12, 0.18, 0.19, 0.21, 0.28, 0.13, 0.27, -0.56, 0.12, 0.66, -0.04, 0.36, -0.68, -0.17, -0.23, 0.66, -0.62, -0.21, -0.39, -0.04, -0.51, 0.26, 0.23, 0.47, 0.45, -0.15, -0.15, 0.05, 0.18, 0.38, -0.0, 0.28, -0.41, 0.08, 0.06, -0.39, -0.38, 0.39, -0.46, 0.61, -0.55, -0.49, -0.64, -0.66, -0.29, -0.25, -0.19, -0.27, -0.25, -0.42, -0.63, 0.04, -0.24, 0.43, -0.2, 0.16, -0.37, 0.14, -0.6, 0.37, 0.36, 0.43, 0.5, 0.21, 0.11, 0.31, 0.31, 0.61, 0.33, 0.4, -0.04, 0.15, -0.07, -0.22, -0.08, -0.06, 0.23, -0.36, 0.6, 0.59, 0.73, 0.66, 0.57, -0.54, 0.23, 0.86, 0.12, 0.69, -0.61, -0.28, -0.3, 0.43, -0.53, -0.01, -0.02, 0.22, -0.39, 0.71, 0.71, 0.81, 0.76, 0.55, -0.33, 0.4, 0.89, 0.23, 0.74, -0.65, -0.23, -0.32, 0.54, -0.62, -0.09, 0.14, 0.22, -0.3, 0.68, 0.66, 0.76, 0.69, 0.55, -0.36, 0.26, 0.91, 0.14, 0.71, -0.54, -0.12, -0.23, 0.71, -0.75, -0.05, 0.05, 0.24, -0.23, 0.78, 0.76, 0.83, 0.76, 0.53, -0.2, 0.55, 0.87, 0.26, 0.73, -0.71, -0.06, -0.4, 0.69, 0.72])
#R3= np.array([0.26, 0.27, -0.1, -0.24, 0.32, 0.32, 0.38, 0.4, 0.6, -0.53, 0.09, 0.8, -0.06, -0.3, 0.19, -0.1, 0.08, 0.03, -0.51, -0.15, -0.45, -0.34, -0.24, 0.35, 0.38, 0.43, 0.21, -0.73, 0.53, 0.25, -0.41, 0.41, -0.28, 0.52, 0.45, -0.27, 0.61, 0.21, 0.19, 0.22, 0.07, -0.02, -0.14, -0.21, -0.41, -0.1, 0.64, -0.66, -0.32, 0.31, -0.09, 0.36, -0.89, -0.7, -0.16, -0.76, 0.05, 0.04, -0.27, -0.15, -0.15, 0.24, 0.33, 0.32, 0.13, -0.64, 0.77, 0.29, -0.26, 0.02, -0.34, 0.73, 0.58, 0.21, 0.61, -0.13, 0.27, 0.05, -0.26, -0.38, 0.4, 0.53, 0.74, 0.58, 0.27, -0.27, 0.29, 0.74, 0.07, -0.4, 0.6, 0.27, 0.19, 0.42, -0.7, 0.26, 0.03, -0.27, -0.29, 0.57, 0.68, 0.88, 0.69, 0.33, 0.02, 0.53, 0.69, 0.25, -0.38, 0.56, 0.23, -0.01, 0.54, -0.76, 0.32, 0.12, -0.23, -0.2, 0.69, 0.76, 0.9, 0.76, 0.3, 0.1, 0.56, 0.67, 0.25, -0.4, 0.56, 0.22, -0.02, 0.55, -0.78, 0.3, -0.01, -0.18, -0.08, 0.61, 0.7, 0.94, 0.73, 0.35, 0.19, 0.73, 0.69, 0.33, -0.28, 0.57, 0.23, 0.0, 0.56, 0.81])
#R4= np.array([0.33, -0.08, -0.37, -0.34, 0.08, -0.11, -0.03, 0.24, -0.0, 0.31, 0.38, 0.07, 0.34, 0.35, 0.11, 0.1, -0.33, 0.06, -0.15, -0.42, -0.37, 0.03, 0.35, 0.24, 0.01, 0.12, -0.09, -0.76, -0.63, 0.31, 0.26, 0.38, 0.36, -0.11, 0.9, 0.13, 0.51, -0.53, -0.05, 0.12, -0.5, -0.29, -0.6, -0.51, -0.45, -0.13, 0.48, 0.48, -0.14, -0.66, -0.36, -0.13, -0.59, -0.8, -0.26, -0.8, 0.77, -0.06, -0.04, 0.07, 0.55, 0.08, -0.18, -0.25, -0.48, -0.37, -0.4, 0.15, 0.62, 0.56, 0.32, 0.69, 0.6, 0.26, 0.61, -0.54, 0.17, 0.05, -0.01, 0.07, 0.57, 0.2, 0.27, 0.11, -0.28, -0.39, 0.07, 0.85, 0.41, 0.02, 0.59, 0.69, 0.48, 0.88, -0.81, 0.27, 0.14, 0.06, 0.05, 0.62, 0.27, 0.3, 0.18, -0.27, -0.28, 0.09, 0.8, 0.33, -0.01, 0.7, 0.65, 0.45, 0.9, -0.82, 0.22, 0.28, -0.1, 0.06, 0.49, 0.07, 0.11, 0.04, -0.19, -0.24, 0.05, 0.84, 0.32, -0.06, 0.64, 0.6, 0.56, 0.89, -0.8, 0.4, 0.21, 0.12, 0.1, 0.65, 0.35, 0.32, 0.22, -0.2, -0.1, 0.19, 0.8, 0.27, -0.02, 0.83, 0.52, 0.41, 0.92, 0.85])
#RT= np.array([0.28, 0.33, 0.18, -0.13, 0.22, 0.2, 0.24, 0.32, 0.67, -0.42, 0.06, 0.76, -0.08, 0.12, 0.86, -0.05, 0.03, 0.28, -0.46, -0.1, -0.32, -0.07, -0.25, 0.57, 0.38, 0.61, 0.38, -0.6, 0.31, 0.42, -0.1, 0.74, 0.2, -0.38, 0.6, -0.14, 0.34, -0.27, 0.01, 0.19, -0.34, 0.1, -0.36, -0.29, -0.45, -0.19, 0.48, -0.58, -0.14, 0.37, -0.4, 0.02, 0.44, -0.48, 0.05, -0.13, -0.12, -0.08, -0.12, -0.04, -0.06, 0.41, 0.26, 0.36, 0.16, -0.37, 0.58, 0.23, -0.02, 0.38, -0.1, 0.43, 0.72, 0.13, 0.5, -0.22, 0.15, 0.31, 0.14, -0.17, 0.51, 0.44, 0.58, 0.5, 0.47, -0.26, 0.16, 0.91, 0.09, -0.0, 0.74, 0.32, 0.2, 0.6, -0.68, 0.25, 0.25, 0.1, -0.22, 0.7, 0.57, 0.72, 0.62, 0.34, 0.03, 0.36, 0.83, 0.21, 0.01, 0.74, 0.49, 0.11, 0.75, -0.83, 0.24, 0.32, 0.05, -0.18, 0.71, 0.57, 0.69, 0.62, 0.31, 0.07, 0.35, 0.83, 0.21, -0.04, 0.75, 0.56, 0.15, 0.82, -0.87, 0.32, 0.23, 0.22, -0.07, 0.68, 0.53, 0.7, 0.59, 0.32, 0.11, 0.48, 0.82, 0.19, 0.08, 0.84, 0.47, 0.09, 0.77, 0.85])


R1= np.array([0.05, -0.03, 0.5, -0.09, -0.34, -0.04, -0.21, 0.17, 0.25, -0.24, -0.2, 0.61, -0.15, -0.15, 0.06, -0.21, -0.21, 0.16, -0.14, -0.53, -0.71, -0.35, -0.43, 0.3, 0.61, 0.35, 0.17, -0.59, -0.69, -0.61, 0.16, -0.45, -0.8, 0.84, -0.8, -0.82, -0.57, -0.61, 0.56, 0.37, 0.53, 0.2, -0.55, -0.56, -0.48, -0.15, 0.65, 0.25, 0.23, 0.5, 0.25, 0.18, -0.46, 0.31, 0.22, 0.51, 0.21, -0.39, -0.33, -0.42, -0.33, 0.32, 0.44, 0.4, -0.1, -0.62, -0.19, -0.29, -0.22, -0.11, -0.67, 0.49, -0.46, -0.49, -0.39, -0.3, 0.19, 0.0, 0.33, 0.21, -0.2, 0.01, -0.04, 0.34, 0.18, -0.3, -0.15, 0.69, 0.28, -0.2, 0.02, -0.31, -0.18, 0.22, -0.15, 0.4, 0.34, 0.32, 0.06, -0.13, -0.04, 0.01, 0.19, 0.29, 0.1, 0.1, 0.78, 0.32, -0.19, -0.46, 0.14, 0.11, 0.54, 0.2, 0.39, 0.13, 0.01, -0.23, 0.17, 0.26, 0.37, 0.33, 0.03, 0.09, 0.1, 0.76, 0.18, -0.46, -0.36, 0.13, -0.05, 0.51, 0.19, 0.64, 0.44, 0.49, 0.06, -0.28, -0.34, -0.21, 0.01, 0.45, 0.13, 0.15, 0.67, 0.56, 0.07, -0.66, 0.35, 0.37, 0.73, 0.42])
R2= np.array([-0.45, 0.25, -0.12, 0.18, 0.19, 0.21, 0.28, 0.13, 0.27, -0.56, 0.12, 0.66, -0.04, 0.36, -0.68, -0.17, -0.23, 0.66, 0.62, -0.21, -0.39, -0.04, -0.51, 0.26, 0.23, 0.47, 0.45, -0.15, -0.15, 0.05, 0.18, 0.38, -0.0, 0.28, -0.41, 0.08, 0.06, 0.39, -0.38, 0.39, -0.46, 0.61, -0.55, -0.49, -0.64, -0.66, -0.29, -0.25, -0.19, -0.27, -0.25, -0.42, -0.63, 0.04, -0.24, 0.43, 0.2, 0.16, -0.37, 0.14, -0.6, 0.37, 0.36, 0.43, 0.5, 0.21, 0.11, 0.31, 0.31, 0.61, 0.33, 0.4, -0.04, 0.15, -0.07, 0.22, -0.08, -0.06, 0.23, -0.36, 0.6, 0.59, 0.73, 0.66, 0.57, -0.54, 0.23, 0.86, 0.12, 0.69, -0.61, -0.28, -0.3, 0.43, 0.53, -0.01, -0.02, 0.22, -0.39, 0.71, 0.71, 0.81, 0.76, 0.55, -0.33, 0.4, 0.89, 0.23, 0.74, -0.65, -0.23, -0.32, 0.54, 0.62, -0.09, 0.14, 0.22, -0.3, 0.68, 0.66, 0.76, 0.69, 0.55, -0.36, 0.26, 0.91, 0.14, 0.71, -0.54, -0.12, -0.23, 0.71, 0.75, -0.05, 0.05, 0.24, -0.23, 0.78, 0.76, 0.83, 0.76, 0.53, -0.2, 0.55, 0.87, 0.26, 0.73, -0.71, -0.06, -0.4, 0.69, 0.72])
R3= np.array([0.26, 0.27, -0.1, -0.24, 0.32, 0.32, 0.38, 0.4, 0.6, -0.53, 0.09, 0.8, -0.06, -0.3, 0.19, -0.1, 0.08, 0.03, 0.51, -0.15, -0.45, -0.34, -0.24, 0.35, 0.38, 0.43, 0.21, -0.73, 0.53, 0.25, -0.41, 0.41, -0.28, 0.52, 0.45, -0.27, 0.61, -0.21, 0.19, 0.22, 0.07, -0.02, -0.14, -0.21, -0.41, -0.1, 0.64, -0.66, -0.32, 0.31, -0.09, 0.36, -0.89, -0.7, -0.16, -0.76, -0.05, 0.04, -0.27, -0.15, -0.15, 0.24, 0.33, 0.32, 0.13, -0.64, 0.77, 0.29, -0.26, 0.02, -0.34, 0.73, 0.58, 0.21, 0.61, 0.13, 0.27, 0.05, -0.26, -0.38, 0.4, 0.53, 0.74, 0.58, 0.27, -0.27, 0.29, 0.74, 0.07, -0.4, 0.6, 0.27, 0.19, 0.42, 0.7, 0.26, 0.03, -0.27, -0.29, 0.57, 0.68, 0.88, 0.69, 0.33, 0.02, 0.53, 0.69, 0.25, -0.38, 0.56, 0.23, -0.01, 0.54, 0.76, 0.32, 0.12, -0.23, -0.2, 0.69, 0.76, 0.9, 0.76, 0.3, 0.1, 0.56, 0.67, 0.25, -0.4, 0.56, 0.22, -0.02, 0.55, 0.78, 0.3, -0.01, -0.18, -0.08, 0.61, 0.7, 0.94, 0.73, 0.35, 0.19, 0.73, 0.69, 0.33, -0.28, 0.57, 0.23, 0.0, 0.56, 0.81])
R4= np.array([0.33, -0.08, -0.37, -0.34, 0.08, -0.11, -0.03, 0.24, -0.0, 0.31, 0.38, 0.07, 0.34, 0.35, 0.11, 0.1, -0.33, 0.06, 0.15, -0.42, -0.37, 0.03, 0.35, 0.24, 0.01, 0.12, -0.09, -0.76, -0.63, 0.31, 0.26, 0.38, 0.36, -0.11, 0.9, 0.13, 0.51, 0.53, -0.05, 0.12, -0.5, -0.29, -0.6, -0.51, -0.45, -0.13, 0.48, 0.48, -0.14, -0.66, -0.36, -0.13, -0.59, -0.8, -0.26, -0.8, -0.77, -0.06, -0.04, 0.07, 0.55, 0.08, -0.18, -0.25, -0.48, -0.37, -0.4, 0.15, 0.62, 0.56, 0.32, 0.69, 0.6, 0.26, 0.61, 0.54, 0.17, 0.05, -0.01, 0.07, 0.57, 0.2, 0.27, 0.11, -0.28, -0.39, 0.07, 0.85, 0.41, 0.02, 0.59, 0.69, 0.48, 0.88, 0.81, 0.27, 0.14, 0.06, 0.05, 0.62, 0.27, 0.3, 0.18, -0.27, -0.28, 0.09, 0.8, 0.33, -0.01, 0.7, 0.65, 0.45, 0.9, 0.82, 0.22, 0.28, -0.1, 0.06, 0.49, 0.07, 0.11, 0.04, -0.19, -0.24, 0.05, 0.84, 0.32, -0.06, 0.64, 0.6, 0.56, 0.89, 0.8, 0.4, 0.21, 0.12, 0.1, 0.65, 0.35, 0.32, 0.22, -0.2, -0.1, 0.19, 0.8, 0.27, -0.02, 0.83, 0.52, 0.41, 0.92, 0.85])
RT= np.array([0.28, 0.33, 0.18, -0.13, 0.22, 0.2, 0.24, 0.32, 0.67, -0.42, 0.06, 0.76, -0.08, 0.12, 0.86, -0.05, 0.03, 0.28, 0.46, -0.1, -0.32, -0.07, -0.25, 0.57, 0.38, 0.61, 0.38, -0.6, 0.31, 0.42, -0.1, 0.74, 0.2, -0.38, 0.6, -0.14, 0.34, 0.27, 0.01, 0.19, -0.34, 0.1, -0.36, -0.29, -0.45, -0.19, 0.48, -0.58, -0.14, 0.37, -0.4, 0.02, 0.44, -0.48, 0.05, -0.13, 0.12, -0.08, -0.12, -0.04, -0.06, 0.41, 0.26, 0.36, 0.16, -0.37, 0.58, 0.23, -0.02, 0.38, -0.1, 0.43, 0.72, 0.13, 0.5, 0.22, 0.15, 0.31, 0.14, -0.17, 0.51, 0.44, 0.58, 0.5, 0.47, -0.26, 0.16, 0.91, 0.09, -0.0, 0.74, 0.32, 0.2, 0.6, 0.68, 0.25, 0.25, 0.1, -0.22, 0.7, 0.57, 0.72, 0.62, 0.34, 0.03, 0.36, 0.83, 0.21, 0.01, 0.74, 0.49, 0.11, 0.75, 0.83, 0.24, 0.32, 0.05, -0.18, 0.71, 0.57, 0.69, 0.62, 0.31, 0.07, 0.35, 0.83, 0.21, -0.04, 0.75, 0.56, 0.15, 0.82, 0.87, 0.32, 0.23, 0.22, -0.07, 0.68, 0.53, 0.7, 0.59, 0.32, 0.11, 0.48, 0.82, 0.19, 0.08, 0.84, 0.47, 0.09, 0.77, 0.85])

R1_matrix=R1.reshape(8,19)
R2_matrix=R2.reshape(8,19)
R3_matrix=R3.reshape(8,19)
R4_matrix=R4.reshape(8,19)
RT_matrix=RT.reshape(8,19)



RR_arr=np.zeros(shape=[8,5,19])

for i in np.arange(0,8):
    
      RR_arr[i,0] = R1_matrix[i]
      RR_arr[i,1] = R2_matrix[i]
      RR_arr[i,2] = R3_matrix[i]
      RR_arr[i,3] = R4_matrix[i]
      RR_arr[i,4] = RT_matrix[i]
    

RR_arr= np.array(RR_arr)  #(40, 19)




'''
RR_arr
Out[141]: 
  
array([[[ 0.05, -0.03,  0.5 , -0.09, -0.34, -0.04, -0.21,  0.17,  0.25,
         -0.24, -0.2 ,  0.61, -0.15, -0.15,  0.06, -0.21, -0.21,  0.16,
          0.02],
        [-0.45,  0.25, -0.12,  0.18,  0.19,  0.21,  0.28,  0.13,  0.27,
         -0.56,  0.12,  0.66, -0.04,  0.36, -0.68, -0.17, -0.23,  0.56,
         -0.42],
        [ 0.26,  0.27, -0.1 , -0.24,  0.32,  0.32,  0.38,  0.4 ,  0.6 ,
         -0.53,  0.09,  0.8 , -0.06, -0.3 ,  0.19, -0.1 ,  0.08,  0.62,
         -0.38],
        [ 0.33, -0.08, -0.37, -0.34,  0.08, -0.11, -0.03,  0.24, -0.  ,
          0.31,  0.38,  0.07,  0.34,  0.35,  0.11,  0.1 , -0.33,  0.17,
         -0.39],
        [ 0.28,  0.33,  0.18, -0.13,  0.22,  0.2 ,  0.24,  0.32,  0.67,
         -0.42,  0.06,  0.76, -0.08,  0.12,  0.86, -0.05,  0.03,  0.62,
         -0.19]],

       [[-0.53, -0.71, -0.35, -0.43,  0.3 ,  0.61,  0.35,  0.17, -0.59,
         -0.69, -0.61,  0.16, -0.45, -0.8 ,  0.84, -0.8 , -0.82, -0.66,
         -0.67],
        [-0.21, -0.39, -0.04, -0.51,  0.26,  0.23,  0.47,  0.45, -0.15,
         -0.15,  0.05,  0.18,  0.38, -0.  ,  0.28, -0.41,  0.08, -0.23,
         -0.34],
        [-0.15, -0.45, -0.34, -0.24,  0.35,  0.38,  0.43,  0.21, -0.73,
          0.53,  0.25, -0.41,  0.41, -0.28,  0.52,  0.45, -0.27, -0.6 ,
          0.04],
        [-0.42, -0.37,  0.03,  0.35,  0.24,  0.01,  0.12, -0.09, -0.76,
         -0.63,  0.31,  0.26,  0.38,  0.36, -0.11,  0.9 ,  0.13, -0.18,
          0.16],
        [-0.1 , -0.32, -0.07, -0.25,  0.57,  0.38,  0.61,  0.38, -0.6 ,
          0.31,  0.42, -0.1 ,  0.74,  0.2 , -0.38,  0.6 , -0.14, -0.52,
         -0.18]],

       [[ 0.56,  0.37,  0.53,  0.2 , -0.55, -0.56, -0.48, -0.15,  0.65,
          0.25,  0.23,  0.5 ,  0.25,  0.18, -0.46,  0.31,  0.22,  0.49,
          0.3 ],
        [-0.38,  0.39, -0.46,  0.61, -0.55, -0.49, -0.64, -0.66, -0.29,
         -0.25, -0.19, -0.27, -0.25, -0.42, -0.63,  0.04, -0.24,  0.64,
         -0.15],
        [ 0.19,  0.22,  0.07, -0.02, -0.14, -0.21, -0.41, -0.1 ,  0.64,
         -0.66, -0.32,  0.31, -0.09,  0.36, -0.89, -0.7 , -0.16,  0.3 ,
         -0.17],
        [-0.05,  0.12, -0.5 , -0.29, -0.6 , -0.51, -0.45, -0.13,  0.48,
          0.48, -0.14, -0.66, -0.36, -0.13, -0.59, -0.8 , -0.26, -0.36,
         -0.37],
        [ 0.01,  0.19, -0.34,  0.1 , -0.36, -0.29, -0.45, -0.19,  0.48,
         -0.58, -0.14,  0.37, -0.4 ,  0.02,  0.44, -0.48,  0.05,  0.41,
         -0.14]],

       [[-0.39, -0.33, -0.42, -0.33,  0.32,  0.44,  0.4 , -0.1 , -0.62,
         -0.19, -0.29, -0.22, -0.11, -0.67,  0.49, -0.46, -0.49, -0.5 ,
         -0.47],
        [ 0.16, -0.37,  0.14, -0.6 ,  0.37,  0.36,  0.43,  0.5 ,  0.21,
          0.11,  0.31,  0.31,  0.61,  0.33,  0.4 , -0.04,  0.15, -0.52,
         -0.05],
        [ 0.04, -0.27, -0.15, -0.15,  0.24,  0.33,  0.32,  0.13, -0.64,
          0.77,  0.29, -0.26,  0.02, -0.34,  0.73,  0.58,  0.21, -0.19,
          0.16],
        [-0.06, -0.04,  0.07,  0.55,  0.08, -0.18, -0.25, -0.48, -0.37,
         -0.4 ,  0.15,  0.62,  0.56,  0.32,  0.69,  0.6 ,  0.26,  0.36,
          0.56],
        [-0.08, -0.12, -0.04, -0.06,  0.41,  0.26,  0.36,  0.16, -0.37,
          0.58,  0.23, -0.02,  0.38, -0.1 ,  0.43,  0.72,  0.13, -0.13,
          0.2 ]],

       [[ 0.19,  0.  ,  0.33,  0.21, -0.2 ,  0.01, -0.04,  0.34,  0.18,
         -0.3 , -0.15,  0.69,  0.28, -0.2 ,  0.02, -0.31, -0.18,  0.3 ,
          0.19],
        [-0.08, -0.06,  0.23, -0.36,  0.6 ,  0.59,  0.73,  0.66,  0.57,
         -0.54,  0.23,  0.86,  0.12,  0.69, -0.61, -0.28, -0.3 ,  0.45,
         -0.53],
        [ 0.27,  0.05, -0.26, -0.38,  0.4 ,  0.53,  0.74,  0.58,  0.27,
         -0.27,  0.29,  0.74,  0.07, -0.4 ,  0.6 ,  0.27,  0.19,  0.35,
         -0.26],
        [ 0.17,  0.05, -0.01,  0.07,  0.57,  0.2 ,  0.27,  0.11, -0.28,
         -0.39,  0.07,  0.85,  0.41,  0.02,  0.59,  0.69,  0.48,  0.42,
          0.51],
        [ 0.15,  0.31,  0.14, -0.17,  0.51,  0.44,  0.58,  0.5 ,  0.47,
         -0.26,  0.16,  0.91,  0.09, -0.  ,  0.74,  0.32,  0.2 ,  0.54,
         -0.09]],

       [[ 0.4 ,  0.34,  0.32,  0.06, -0.13, -0.04,  0.01,  0.19,  0.29,
          0.1 ,  0.1 ,  0.78,  0.32, -0.19, -0.46,  0.14,  0.11,  0.51,
          0.3 ],
        [-0.01, -0.02,  0.22, -0.39,  0.71,  0.71,  0.81,  0.76,  0.55,
         -0.33,  0.4 ,  0.89,  0.23,  0.74, -0.65, -0.23, -0.32,  0.5 ,
         -0.54],
        [ 0.26,  0.03, -0.27, -0.29,  0.57,  0.68,  0.88,  0.69,  0.33,
          0.02,  0.53,  0.69,  0.25, -0.38,  0.56,  0.23, -0.01,  0.57,
         -0.57],
        [ 0.27,  0.14,  0.06,  0.05,  0.62,  0.27,  0.3 ,  0.18, -0.27,
         -0.28,  0.09,  0.8 ,  0.33, -0.01,  0.7 ,  0.65,  0.45,  0.53,
          0.51],
        [ 0.25,  0.25,  0.1 , -0.22,  0.7 ,  0.57,  0.72,  0.62,  0.34,
          0.03,  0.36,  0.83,  0.21,  0.01,  0.74,  0.49,  0.11,  0.61,
         -0.18]],

       [[ 0.39,  0.13,  0.01, -0.23,  0.17,  0.26,  0.37,  0.33,  0.03,
          0.09,  0.1 ,  0.76,  0.18, -0.46, -0.36,  0.13, -0.05,  0.33,
          0.07],
        [-0.09,  0.14,  0.22, -0.3 ,  0.68,  0.66,  0.76,  0.69,  0.55,
         -0.36,  0.26,  0.91,  0.14,  0.71, -0.54, -0.12, -0.23,  0.46,
         -0.44],
        [ 0.32,  0.12, -0.23, -0.2 ,  0.69,  0.76,  0.9 ,  0.76,  0.3 ,
          0.1 ,  0.56,  0.67,  0.25, -0.4 ,  0.56,  0.22, -0.02,  0.63,
         -0.61],
        [ 0.22,  0.28, -0.1 ,  0.06,  0.49,  0.07,  0.11,  0.04, -0.19,
         -0.24,  0.05,  0.84,  0.32, -0.06,  0.64,  0.6 ,  0.56,  0.48,
          0.56],
        [ 0.24,  0.32,  0.05, -0.18,  0.71,  0.57,  0.69,  0.62,  0.31,
          0.07,  0.35,  0.83,  0.21, -0.04,  0.75,  0.56,  0.15,  0.61,
         -0.15]],

       [[ 0.64,  0.44,  0.49,  0.06, -0.28, -0.34, -0.21,  0.01,  0.45,
          0.13,  0.15,  0.67,  0.56,  0.07, -0.66,  0.35,  0.37,  0.73,
          0.55],
        [-0.05,  0.05,  0.24, -0.23,  0.78,  0.76,  0.83,  0.76,  0.53,
         -0.2 ,  0.55,  0.87,  0.26,  0.73, -0.71, -0.06, -0.4 ,  0.6 ,
         -0.52],
        [ 0.3 , -0.01, -0.18, -0.08,  0.61,  0.7 ,  0.94,  0.73,  0.35,
          0.19,  0.73,  0.69,  0.33, -0.28,  0.57,  0.23,  0.  ,  0.67,
         -0.62],
        [ 0.4 ,  0.21,  0.12,  0.1 ,  0.65,  0.35,  0.32,  0.22, -0.2 ,
         -0.1 ,  0.19,  0.8 ,  0.27, -0.02,  0.83,  0.52,  0.41,  0.71,
          0.45],
        [ 0.32,  0.23,  0.22, -0.07,  0.68,  0.53,  0.7 ,  0.59,  0.32,
          0.11,  0.48,  0.82,  0.19,  0.08,  0.84,  0.47,  0.09,  0.75,
         -0.19]]])
        
    
    '''

