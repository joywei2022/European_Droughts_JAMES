# -*-coding:UTF-8 -*-

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
import pandas as pdprin
from pandas import *
import matplotlib.dates as pltdt
import math

from scipy import signal


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
    
    data3 = data2[(12*M):,:,:]  #,90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017


    #print data2.shape, data3.shape
     
    #plt.imshow(data3[0])
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN   #kgC/m2/month
    
    data3=data3*12.0*1e-12
 
    data3=data3*area[2*(90-latMax):2*(90-latMin),int(2*(180+lonMin)-1):2*(180+lonMax)]  
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
 
 
    
 
#################################################################
#   EUROCOM for Drought 2018  
#################################################################

dir21='/Volumes/HEWEI_T5/Inversions/Drought-2018_v2/'

def flux_anomaly_eurocom_Drought2018(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=7
#    
#    latMin=30   #33 #35
#    latMax=75   #73 #70
#    lonMin=-14.5
#    lonMax=35
    
#    dataset='co2flux_monthly_flexinvert_2009_2018.nc'
#    varstr='fnee_post'
 
    file= os.path.join(dir21,dataset)
    
    data1=[]
    lon=[]
    lat=[]
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
    lat=ncf.variables['latitude'][:]
    lon=ncf.variables['longitude'][:]
         
    ncf.close()
    
    data1=np.array(data1)  #(120, 80, 100)
    lat=np.array(lat) 
    lon=np.array(lon) 
    area = globarea(720,360,True)
    
    
   # plt.imshow(data1[0])
    
    
    #lon=np.arange(0,360,0.5)
    #lat=np.arange(0,180,1)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
    
    data3 = data2[(12*M):,:,:]  #,90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017
 
    #plt.imshow(data3[0])
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN   #kgC/m2/month
    
    data3=data3*24*365*1e-12
 
     
    if  dataset=='co2flux_monthly_carboscope_clim_2009_2018.nc' or dataset=='co2flux_monthly_carboscope_core_2009_2018.nc' or dataset=='co2flux_monthly_carboscope_all_2009_2018.nc' :
        area = globarea(720*2,360*2,True)
        data3=data3*area[int(4*(90-latMax)):int(4*(90-latMin)),int(4*(180+lonMin)-2):int(4*(180+lonMax))]  
    else:
        area = globarea(720,360,True)
        data3=data3*area[int(2*(90-latMax)):int(2*(90-latMin)),int(2*(180+lonMin)-1):int(2*(180+lonMax))]  
   
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
 
 
##################################################################
##   EUROCOM for Drought 2018 (select 2009-2015, 7 years)
##################################################################
#
#dir21='/Volumes/HEWEI_T5/Inversions/Drought-2018_v2/'
#
#def flux_anomaly_eurocom_Drought2018_cs(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
#         
##    M=0 
##    N=7
##    
##    latMin=30   #33 #35
##    latMax=75   #73 #70
##    lonMin=-14.5
##    lonMax=35
#    
##    dataset='co2flux_monthly_flexinvert_2009_2018.nc'
##    varstr='fnee_post'
# 
#    file= os.path.join(dir21,dataset)
#    
#    data1=[]
#    lon=[]
#    lat=[]
#     
#    ncf=nc.Dataset(file)
#    data1=ncf.variables[varstr][:]
#    lat=ncf.variables['latitude'][:]
#    lon=ncf.variables['longitude'][:]
#         
#    ncf.close()
#    
#    data1=np.array(data1)  #(120, 80, 100)
#    lat=np.array(lat) 
#    lon=np.array(lon) 
#    area = globarea(720*2,360*2,True)
#    
#    
#   # plt.imshow(data1[0])
#    
#    
#    #lon=np.arange(0,360,0.5)
#    #lat=np.arange(0,180,1)
#    lon2, lat2 = np.meshgrid(lon,lat) 
#    
#    data2=[]
#
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#    
#    data3 = data2[(12*M):(-12*3),:,:]  #,90-latMax:90-latMin,180+lonMin:180+lonMax]   # Corrected on July 3, 2017
# 
#    #plt.imshow(data3[0])
#    
#
#    a=[]
#    
#    mask=(data3==0)
#    data3[mask]=np.NaN   #kgC/m2/month
#    
#    data3=data3*24*365*1e-12
# 
#    data3=data3*area[int(4*(90-latMax)):int(4*(90-latMin)),int(4*(180+lonMin)-2):int(4*(180+lonMax))]  
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
#    


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
    
    data3 = data2[(12*M):,int(np.floor((90-latMax)/4-1)):int(np.ceil((90-latMin)/4-1)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   

    plt.imshow(data3[0])
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN   #gC/m2/month
    
    data3=data3*12.0*1e-15
 
    data3=data3*area[int(np.floor((90-latMax)/4-1)):int(np.ceil((90-latMin)/4-1)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   

      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
 

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
    
    area = globarea(360,180,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
    
    #data2a = np.empty( (12*N, data2.shape[1], data2.shape[2]), dtype=float, order='C')
 
  
    data3 = data2[:, int(np.floor(90-latMax)):int(np.ceil(90-latMin)), int(np.floor(180+lonMin)):int(np.ceil(180+lonMax)) ]
     
  
    mask=(data3==0)
    data3[mask]=np.NaN
     
    
    data3 = data3*365*1e-15
    
    
    data3=data3*area[ int(np.floor(90-latMax)):int(np.ceil(90-latMin)), int(np.floor(180+lonMin)):int(np.ceil(180+lonMax)) ] 
     
    plt.imshow(data3[0])
    
    

    a=[]
 
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
  


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
    
    area = globarea(360,180,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
        
    data3aa = data2[:, int(90-latMax):int(90-latMin),int(360+lonMin):]
    data3bb = data2[:, int(90-latMax):int(90-latMin),0:lonMax]
    
 
  
    mask=(data3aa==0)
    data3aa[mask]=np.NaN

    mask=(data3bb==0)
    data3bb[mask]=np.NaN
          
    plt.imshow(data3aa[0])
    plt.imshow(data3bb[0])
    
    data3aa=data3aa*area[int(90-latMax):int(90-latMin),int(360+lonMin):]*12*1e-15   # gC/m2/month    # kg C m-2 s-1 365*24*3600
    data3bb=data3bb*area[int(90-latMax):int(90-latMin),0:lonMax]*12*1e-15
  
    a=[]
 
      
    for i in range(0,data3aa.shape[0]): 
        me = np.nansum(data3aa[i]) + np.nansum(data3bb[i])    
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3aa[12*j+i])/float(N) + np.nansum(data3bb[12*j+i])/float(N)
   
    del data3aa
    del data3bb
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
  
#################################################################
#   Wu et al.(2020) RSE paper, CCDAS-CO2+SM inversions
#################################################################

dir4='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_WuCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=1 
#    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
      
    #dataset ='ccdas_sm+fapar'
    
    #varstr = 'rnep'
    varstr = 'nep'
    
    if dataset =='CCDAS_6yr_2deg_assim_results':
         file=os.path.join(dir4,dataset+'/'+'CCDAS_co2+sm.nc') 
    
    if dataset =='ccdas_sm+fapar':
         file=os.path.join(dir4,dataset+'/'+'CCDAS_co2+sm+fapar.nc')  
    
    data1=[]
    lat=[]
    lon=[]
    #time=[]
    
    ncf=nc.Dataset(file)
#    lat.extend(ncf.variables['latitude'][:])
#    lon.extend(ncf.variables['longitude'][:])
    lat.extend(ncf.variables['lat'][:])
    lon.extend(ncf.variables['lon'][:])
    
    data1.append(ncf.variables[varstr][:]) 
    #time = ncf.variables['time'][:]
    
    
    
    ncf.close()
    
    data1=np.array(data1[0])
    data2=data1
    
    #area = globarea(180,90,True)
    area = globarea(180*2,90*2,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
        
    #data3 = data2[12*M:, int((90-latMax)*0.5):int((90-latMin)*0.5),int((180+lonMin)*0.5):int((180+lonMax)*0.5)]
    data3 = data2[12*M:, int((90-latMax)*1):int((90-latMin)*1),int((180+lonMin)*1):int((180+lonMax)*1)]
  
  
    mask=(data3==-9999)
    data3[mask]=np.NaN

 
    plt.imshow(data3[0])
 
    #data3=data3*area[int((90-latMax)*0.5):int((90-latMin)*0.5),int((180+lonMin)*0.5):int((180+lonMax)*0.5)]*12*1e-15    #gC/m2/month
    data3=data3*area[int((90-latMax)*1):int((90-latMin)*1),int((180+lonMin)*1):int((180+lonMax)*1)]*12*1e-15    #gC/m2/month
 
    a=[]
 
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])    
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
  

#################################################################
#   Wu et al.(2020) RSE paper, CCDAS-CO2+SM inversions
#################################################################

dir4a='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_WuCCDASv2022(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=1 
#    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
      
    #dataset ='ccdas_sm+fapar'
    
    #varstr = 'rnep'
    varstr = 'nep'
   
    file=os.path.join(dir4a,dataset+'/CCDAS_v2022/'+'CCDAS_NEP_2010_2020_05deg_gridcorrection_new.nc')  
    
    data1=[]
    lat=[]
    lon=[]
    #time=[]
    
    ncf=nc.Dataset(file)
#    lat.extend(ncf.variables['latitude'][:])
#    lon.extend(ncf.variables['longitude'][:])
    lat.extend(ncf.variables['lat'][:])
    lon.extend(ncf.variables['lon'][:])
    
    data1.append(ncf.variables[varstr][:]) 
    #time = ncf.variables['time'][:]
    
    
    
    ncf.close()
    
    data1=np.array(data1[0])
    data2=data1
    
    #area = globarea(180,90,True)
    area = globarea(180*4,90*4,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
        
    #data3 = data2[12*M:, int((90-latMax)*0.5):int((90-latMin)*0.5),int((180+lonMin)*0.5):int((180+lonMax)*0.5)]
    data3 = data2[12*M:-12*5, int((90-latMax)*2):int((90-latMin)*2),int((180+lonMin)*2):int((180+lonMax)*2)]
  
  
    mask=(data3==-9999)
    data3[mask]=np.NaN

 
    plt.imshow(data3[0])
 
    #data3=data3*area[int((90-latMax)*0.5):int((90-latMin)*0.5),int((180+lonMin)*0.5):int((180+lonMax)*0.5)]*12*1e-15    #gC/m2/month
    data3=data3*area[int((90-latMax)*2):int((90-latMin)*2),int((180+lonMin)*2):int((180+lonMax)*2)]*12*1e-15    #gC/m2/month
 
    a=[]
 
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])    
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
  
  
#################################################################
#   Scholze et al.(2019) GRL paper, CCDAS-CO2+SM+VOD inversions
#################################################################

dir5='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_ScholzeCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
    M=0 
    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
 
    dataset ='ccdas_sm+vod'
    varstr = 'nep'
    
    file=os.path.join(dir5,dataset+'/'+'grid_nep_2010-2015.nc')  
    
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
     
        
    data3 = data2[12*M:, int((90-latMax)*4):int((90-latMin)*4),int((180+lonMin)*4):int((180+lonMax)*4)]
#    data3aa = data2[12*M:, int((90-latMax)*4):int((90-latMin)*4),int((360+lonMin)*4):]
#    data3bb = data2[12*M:, int((90-latMax)*4):int((90-latMin)*4),0:int(lonMax*4)]
#      
  
    mask=(data3==-9999)
    data3[mask]=np.NaN
    
 
    plt.imshow(data3[0])
 
    data3=data3*area[int((90-latMax)*4):int((90-latMin)*4),int((180+lonMin)*4):int((180+lonMax)*4)]*12*1e-15    #gC/m2/month
#    data3aa=data3aa*area[int((90-latMax)*4):int((90-latMin)*4),int((360+lonMin)*4):]*12*1e-15    #gC/m2/month
#    data3bb=data3bb*area[int((90-latMax)*4):int((90-latMin)*4),0:int(lonMax*4)]*12*1e-15    #gC/m2/month
      
    a=[]
 
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i]) 
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B




#################################################################
#   CAMS-flux1*1 
#################################################################

#dir6='/Volumes/HEWEI_T5/Inversions/CAMS_satellite_GOSAT+OCO2_FT21r2_to202112/1Deg/2010-2015/'
dir6='/Volumes/HEWEI_T5/Inversions/CAMS_satellite_GOSAT+OCO2_FT21r2_to202112/'
        
def flux_anomaly_cams_oco2_surface(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
    M=0 
    N=6 
 
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35 
    
    varstr='flux_apos_bio'
    
    #files = glob.glob(os.path.join(dir6,'Resampled_1Deg*.nc')) 
    files = glob.glob(os.path.join(dir6,'cams73_v21r2*.nc')) 
    
    data1=[]
    area=[]
    lat=[]
    lon=[]
     
    ncf=nc.Dataset(files[0])
#    lat=ncf.variables['lat'][:]
#    lon=ncf.variables['lon'][:]
    
    lat.extend(ncf.variables['latitude'][:])
    lon.extend(ncf.variables['longitude'][:])
    
    ncf.close()   
        
    #for file in files[7:]:
    for file in files[7:-12*6]:
        ncf=nc.Dataset(file)
        data1.append(ncf.variables[varstr][:])  
        area=ncf.variables['area'][:]  #same as calculated
     
        ncf.close()
    
    data1=np.array(data1)
    lat = np.array(lat)
    lon=np.array(lon)
    
        
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    area=np.array(area)
    #area = globarea(360,180,True)
    
    data2=[]
    m=data1.shape[0]
    for i in range(0,m):
        xx = maskoceans(lon2, lat2, data1[i])
        plt.imshow(xx)
        data2.append(np.flipud(xx))
    
    data2=np.array(data2)
    plt.imshow(data2[0])
    
    #data problem, only valid from 1 to 147  all 147 months, 12 years
    
    #data3=data2[(12*M):(-12),(90-latMax):(90-latMin),(180+lonMin):(180+lonMax)]   # Corrected on July 3, 2017
    data3=data2[(12*M):]   # Corrected on July 3, 2017
 

    

    mask=(data3==0)
    #mask=(data3==-9999.0)
    data3[mask]=np.NaN

    #gC/m2/yr 
    
    # now kgC m-2 month-1
    data3=data3*12*1e3*1e-15   
   
    data3=data3* np.flipud(area)   
    #data3 = data3[:,(90-latMax):(90-latMin),(180+lonMin):(180+lonMax)]  
    data3 = data3[:,int((90-latMax)*96/180.0):int((90-latMin)*96/180.0),int((180+lonMin)*96/360.0):int((180+lonMax)*96/360.0)]  
        
      
    plt.imshow(data3[0])
        
#    data3a = np.zeros(shape=[data3.shape[0],data3.shape[1]*2,  data3.shape[2]*2])*np.NaN
#
#    for k in range(0,data3.shape[0]): 
#        for i in range(0,data3.shape[1]):
#            for j in range(0,data3.shape[2]):
#                 data3a[k, i*2:(i+1)*2, j*2:(j+1)*2] = data3[k,i,j]
#
#
#    data3 = data3a/4.0

 
    a=[]
 

    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,N):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a
    B=a-a1
   
 
    return A,B 



#################################################################
#   Jiang et al.(2020) ACPD paper, GCASv2, GOSAT inversions
#################################################################

dir7='/Volumes/HEWEI_T5/Inversions/GCASv2/v20210630/'
  
def flux_anomaly_JiangGCAS2021(dataset,varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=10
#    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
#    dataset = "posterior.fluxes.GOSAT.nc"        
#    varstr = 'bio_monthly_opt'
    
    file=os.path.join(dir7,dataset)  
    
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
    
    area = globarea(360,180,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    for i in range(0,data1.shape[0]):
        data2[i] = np.flipud(data1[i])
     
        
        
    data3 = np.zeros(shape= (data2.shape[0],data2.shape[1],data2.shape[2]) )
    
    if dataset == "posterior.fluxes.GCAS2020.nc":
            for i in range(0,data2.shape[0]):
                data3[i, :, 0:180] = data2[i, :, 180:]
                data3[i, :, 180: ] = data2[i, :, 0:180]
                
    if dataset == "posterior.fluxes.OCO2.nc":
            for i in range(0,data2.shape[0]):
                data3[i, :, 0:180] = data2[i, :, 180:]
                data3[i, :, 180: ] = data2[i, :, 0:180] 
                
    if dataset == "posterior.fluxes.GOSAT.nc":
            for i in range(0,data2.shape[0]):
                data3[i, 0:90, :] = data2[i, 90:180, :]
                data3[i, 90:180, : ] = data2[i, 0:90, :]

      
            data3a = np.zeros(shape= (data2.shape[0],data2.shape[1],data2.shape[2]) )
             
            for i in range(0,data3.shape[0]):
                data3a[i, :, 0:300] = data3[i, :, 60:360]
                data3a[i, :, 300:360] = data3[i, :, 0:60]
                 
            data3 = data3a 
    
                
    mask=(data3==0)
    data3[mask]=np.NaN
 
    plt.imshow(data3[0])
   
    if dataset == "posterior.fluxes.GCAS2020.nc":
          data3 = data3[:-1, int(90-latMax):int(90-latMin),int(180+lonMin):int(180+lonMax) ]
    if dataset == "posterior.fluxes.OCO2.nc":
          data3 = data3[:-1, int(90-latMax):int(90-latMin),int(180+lonMin):int(180+lonMax) ]
    if dataset == "posterior.fluxes.GOSAT.nc":
          data3 = data3[7:-2, int(90-latMax):int(90-latMin),int(180+lonMin):int(180+lonMax) ]
          
    data3=data3*area[int(90-latMax):int(90-latMin),int(180+lonMin):int(180+lonMax)] *12*1e-15     # gC/m2/month    # kg C m-2 s-1 365*24*3600
  
    #plt.imshow(data3[11])
    
    
    a=[]
 
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N) 
  
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
  


#################################################################
#   OCO-2 v10 MIP CO2 inversion: 2015-2020
#################################################################

dir8='/Volumes/HEWEI_T5/OCO-2_v10_MIP/LNLG/'      #    IS

def flux_anomaly_oco2v10mip(dir8,dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
     
#    M=0
#    N=4    
#    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
#  
#    dataset='UT_gridded_fluxes_LNLG.nc4'      
#    varstr='land'

    file= os.path.join(dir8,dataset) 
    
    data1=[]
    lon=[]
    lat=[]
     
    ncf=nc.Dataset(file)
    lat=ncf.variables['latitude'][:]
    lon=ncf.variables['longitude'][:]
    data1 = ncf.variables[varstr][:]
    ncf.close()  
   
    data1=np.array(data1)    #( 48, 180, 360)
    lat=np.array(lat) 
    lon=np.array(lon) 
  
    #plt.imshow(data1[0])
  
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)
 
  
    #plt.imshow(data2[0])
  
    # gC/m2/yr 
      
    area = globarea(360,180,True)
 
    data3 = data2 * area
    #
    data3 = data3[:,int(90-latMax):int(90-latMin), int(180+lonMin):int(180+lonMax) ]  
    data3 = data3  * 1e-15  # * mask_China1   #*12
    
    #mask=(data3==0)
    #data3[mask]=np.NaN   #gC/m2/month
    
    plt.imshow( np.nansum(data3,axis=0)) 
    
   

    a=[]
    
   
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i] )/float(N)   #* mask_China1
   
    
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
       
    #XX = np.nansum( data3, axis = 0 )/float(N)  #* mask_China1
    
    return A,B    #,baseline,XX
    

############################################
   
def pipeline(latMin, latMax, lonMin, lonMax,flag, fn ) :    
        
#    latMin=35
#    latMax=70
#    lonMin=-10
#    lonMax=40
        
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
    
    
    flag=1
#    fn='NA' 
     
    N=15 
     
    pydates1=[]
    for i in range(0,N): 
        yr = 2001 + i
        tmp = datetime(yr, 1, 1)
        pydates1.append(tmp)  
        
        
    N=20    
    pydates2=[]
    for i in range(0,N): 
        yr = 2001 + i
        tmp = datetime(yr, 1, 1)
        pydates2.append(tmp)  
        
   #--------------------------------------   
    
      
    M=4 
    N=6
    
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
   
    M=4 
    N=6
    
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
    
      
    
    #------------------------------    
#    M=1 
#    N=9
#
#    
#    latMin=30   #33 #35
#    latMax=75   #73 #70
#    lonMin=-14.5
#    lonMax=35
#    
#    varstr='fnee_post' 
#    dataset = 'co2flux_monthly_flexinvert_clim_2009_2018.nc'
#    eurocom8nee =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
# 
#    varstr='fnee_post' 
#    dataset = 'co2flux_monthly_flexinvert_core_2009_2018.nc'
#    eurocom8nee_select =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#
#          
#    latMin=33  #35
#    latMax=73  #70
#    lonMin=-14.5
#    lonMax=35
#    
#    dataset = 'co2flux_monthly_lumia_clim_2009_2018.nc'     
#    varstr='fnee_post'
#    eurocom9nee =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
# 
#    dataset = 'co2flux_monthly_lumia_all_2009_2018.nc'     
#    varstr='fnee_post'
#    eurocom9nee_all =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset = 'co2flux_monthly_lumia_core_2009_2018.nc'     
#    varstr='fnee_post'
#    eurocom9nee_select =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
# 
#         
#    latMin=32.5  #35
#    latMax=73.5  #70
#    lonMin=-15
#    lonMax=35.5
#    
#    dataset = 'co2flux_monthly_name_clim_2009_2018.nc'     
#    varstr='fnee_post'
#    eurocom10nee =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
# 
#         
#    dataset = 'co2flux_monthly_name_core_2009_2018.nc'     
#    varstr='fnee_post'
#    eurocom10nee_select =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
# 
#     
#    latMin=31.5   #33 #35
#    latMax=74     #73 #70
#    lonMin=-15
#    lonMax=35
#    
#    dataset = 'co2flux_monthly_pyvar_all_2009_2018.nc'
#    varstr='fnee_post'
#    eurocom11nee_all =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#         
#    M=1 
#    N=9
#    
#    latMin=33  #35
#    latMax=73  #70
#    lonMin=-14.5
#    lonMax=35
#    dataset = 'co2flux_monthly_carboscope_clim_2009_2018.nc'
#    varstr='fnee_post'
#    eurocom12nee =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#        
#    
#    dataset = 'co2flux_monthly_carboscope_all_2009_2018.nc'
#    varstr='fnee_post'
#    eurocom12nee_all =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#        
#    
#    dataset = 'co2flux_monthly_carboscope_core_2009_2018.nc'
#    varstr='fnee_post'
#    eurocom12nee_select =  flux_anomaly_eurocom_Drought2018(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
        
     
  
    #------------------------------
    latMin=33    #35
    latMax=73    #70
    lonMin=-14.5 #-10
    lonMax=35    #40
    
    lonMin=-15
     
    M=0  
    N=6 
      
#    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_ens =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
        
#    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
#    varstr='flux_CASA'
#    
#    ByrneGOSATnee_casa =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#      
#    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
#    varstr='flux_FLUXCOM'
#    
#    ByrneGOSATnee_fc =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#       
#    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
#    varstr='flux_SiB'
#    
#    ByrneGOSATnee_sib =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
        
    
    
    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
    varstr='ensemble_mean'
    
    ByrneGOSATnee_all =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
        
    dataset='Byrne2020_NEE_GOSAT_only'
    varstr='ensemble_mean'
    
    ByrneGOSATnee_gosat =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
      
    dataset='Byrne2020_NEE_surface_only'
    varstr='ensemble_mean'
    
    ByrneGOSATnee_surface =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
       
    dataset='Byrne2020_NEE_TCCON_only'
    varstr='ensemble_mean'
    
    ByrneGOSATnee_tccon =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    
    
    
    #------------------------------

    M=0 
    N=9
    
    dataset='CMS-Flux-NBE-2020.monthly.grid'           
    varstr = 'post-NBE'  #'flux'
    LiuGOSATnee =  flux_anomaly_LiuGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)

      
    #------------------------------

    M=0 
    N=6
     
    #dataset='CMS-Flux-NBE-2020.monthly.grid'        
    varstr = 'bio_monthly_opt'
    JiangGOSATnee =  flux_anomaly_JiangGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
  
    
    
    #------------------------------

    #M=1 
    M=0
    N=6
     
    dataset='CCDAS_6yr_2deg_assim_results'        
    varstr = 'rnep'
    WuCCDASnee1 =  flux_anomaly_WuCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
  
    
    #------------------------------
    
    #M=1 
    M=0
    N=6
     
    dataset='ccdas_sm+fapar'        
    varstr = 'rnep'
    #WuCCDASnee2 =  flux_anomaly_WuCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
  
    WuCCDASnee2 =  flux_anomaly_WuCCDASv2022(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
   
    #------------------------------
    
    M=0 
    N=6
     
    dataset='ccdas_sm+vod'        
    varstr = 'nep'
    ScholzeCCDASnee =  flux_anomaly_ScholzeCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
   
    #------------------------------
      
    M=0 
    N=6 
      
    varstr='flux_apos_bio'
    camsnee =  flux_anomaly_cams_oco2_surface(varstr,M, N, latMin, latMax, lonMin, lonMax)
  
    
    
    #------------------------------
    
    M=0 
    N=5
    dataset = 'posterior.fluxes.OCO2.nc'       
    varstr = 'bio_monthly_opt'
    JiangOCO2nee =  flux_anomaly_JiangGCAS2021(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
   
    #------------------------------
 
#    M=0
#    N=6    
#    
#    
#    varstr='land'
#    
#    dir8='/Volumes/HEWEI_T5/OCO-2_v10_MIP/IS/'    
#    
#    dataset='EnsMean_gridded_fluxes_IS.nc4'    
#    oco2mipnee_v10mean_IS =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='EnsStd_gridded_fluxes_IS.nc4'    
#    oco2mipnee_v10std_IS =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#
#
#    dir8='/Volumes/HEWEI_T5/OCO-2_v10_MIP/LNLG/'    
#    
#    dataset='EnsMean_gridded_fluxes_LNLG.nc4'    
#    oco2mipnee_v10mean_LNLG =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='EnsStd_gridded_fluxes_LNLG.nc4'    
#    oco2mipnee_v10std_LNLG =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#
#
#    dir8='/Volumes/HEWEI_T5/OCO-2_v10_MIP/LNLGIS/'    
#    
#    dataset='EnsMean_gridded_fluxes_LNLGIS.nc4'    
#    oco2mipnee_v10mean_LNLGIS =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='EnsStd_gridded_fluxes_LNLGIS.nc4'    
#    oco2mipnee_v10std_LNLGIS =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)


#    dir8='/Volumes/HEWEI_T5/OCO-2_v10_MIP/Prior/'    
 

 
   
#    #------------------------------
#    M=0
#    N=6    
#   
#     
#    dir8='/Volumes/HEWEI_T5/OCO-2_v10_MIP/LNLG/'   
#    varstr='land'
#    
#    dataset='Ames_gridded_fluxes_LNLG.nc4'
#    oco2mip1nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='Baker_gridded_fluxes_LNLG.nc4' 
#    oco2mip2nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='CAMS_gridded_fluxes_LNLG.nc4' 
#    oco2mip3nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='CMS-Flux_gridded_fluxes_LNLG.nc4'     
#    oco2mip4nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='COLA_gridded_fluxes_LNLG.nc4' 
#    oco2mip5nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='CSU_gridded_fluxes_LNLG.nc4'     
#    oco2mip6nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#      
#    dataset='CT_gridded_fluxes_LNLG.nc4'     
#    oco2mip7nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#      
#    dataset='JHU_gridded_fluxes_LNLG.nc4' 
#    oco2mip8nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='LoFI_gridded_fluxes_LNLG.nc4'     
#    oco2mip9nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='NIES_gridded_fluxes_LNLG.nc4'     
#    oco2mip10nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='OU_gridded_fluxes_LNLG.nc4' 
#    oco2mip11nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='TM5-4DVAR_gridded_fluxes_LNLG.nc4' 
#    oco2mip12nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='UT_gridded_fluxes_LNLG.nc4' 
#    oco2mip13nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset='WOMBAT_gridded_fluxes_LNLG.nc4' 
#    oco2mip14nee =  flux_anomaly_oco2v10mip(dir8,dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)


    
    #######################################
    # raw fluxes
    #######################################

    
    eurocom1neeann=[]
    eurocom2neeann=[]
    eurocom3neeann=[]
    eurocom4neeann=[]
    eurocom5neeann=[]
    eurocom6neeann=[]
  
    eurocom7neeann=[]

#    eurocom8neeann=[]
#    eurocom9neeann=[]  
#    eurocom9allneeann=[] 
#    eurocom10neeann=[]
#    eurocom11allneeann=[]    
#    eurocom12neeann=[]
#    eurocom12allneeann=[]
#    
#    eurocom8selectneeann=[]
#    eurocom9selectneeann=[]   
#    eurocom10selectneeann=[]
#    eurocom12selectneeann=[]
   
    ByrneGOSATneeann_all=[]
    ByrneGOSATneeann_gosat=[]
    ByrneGOSATneeann_surface=[]
    ByrneGOSATneeann_tccon=[]
    
    LiuGOSATneeann=[]
    JiangGOSATneeann=[]
    WuCCDASneeann1=[]
    WuCCDASneeann2=[]
    ScholzeCCDASneeann=[]
    camsneeann=[]
    
    JiangOCO2neeann=[]
#    oco2mipneeann_v10mean_IS=[]
#    oco2mipneeann_v10std_IS=[]
#    oco2mipneeann_v10mean_LNLG=[]
#    oco2mipneeann_v10std_LNLG=[]
#    oco2mipneeann_v10mean_LNLGIS=[]
#    oco2mipneeann_v10std_LNLGIS=[]
#    
#    
#    oco2mip1neeann=[]
#    oco2mip2neeann=[]
#    oco2mip3neeann=[]
#    oco2mip4neeann=[]
#    oco2mip5neeann=[]
#    oco2mip6neeann=[]
#    oco2mip7neeann=[]
#    oco2mip8neeann=[]
#    oco2mip9neeann=[] 
#    oco2mip10neeann=[] 
#    oco2mip11neeann=[] 
#    oco2mip12neeann=[] 
#    oco2mip13neeann=[] 
#    oco2mip14neeann=[] 


    for i in range(0, 6):
       eurocom1neeann.append(np.nanmean( eurocom1nee[flag][i*12:(i+1)*12] + eurocom1fire[flag][i*12:(i+1)*12] ) )
       eurocom2neeann.append(np.nanmean( eurocom2nee[flag][i*12:(i+1)*12] ) )
       eurocom3neeann.append(np.nanmean( eurocom3nee[flag][i*12:(i+1)*12] ) )
       eurocom4neeann.append(np.nanmean( eurocom4nee[flag][i*12:(i+1)*12] ) )
       eurocom5neeann.append(np.nanmean( eurocom5nee[flag][i*12:(i+1)*12]  + eurocom5fire[flag][i*12:(i+1)*12] ) )
       eurocom6neeann.append(np.nanmean( eurocom6nee[flag][i*12:(i+1)*12] ) )
    
    for i in range(0, 5):
       eurocom7neeann.append(np.nanmean( eurocom7nee[flag][i*12:(i+1)*12] ) )


#    for i in range(0, 9):
#       eurocom8neeann.append(np.nanmean( eurocom8nee[flag][i*12:(i+1)*12] ) )
#       eurocom9neeann.append(np.nanmean( eurocom9nee[flag][i*12:(i+1)*12] ) )
#       eurocom9allneeann.append(np.nanmean( eurocom9nee_all[flag][i*12:(i+1)*12] ) ) 
#       eurocom10neeann.append(np.nanmean( eurocom10nee[flag][i*12:(i+1)*12] ) )
#       eurocom11allneeann.append(np.nanmean( eurocom11nee_all[flag][i*12:(i+1)*12] ) )
#       eurocom12neeann.append(np.nanmean( eurocom12nee[flag][i*12:(i+1)*12] ) )
#       eurocom12allneeann.append(np.nanmean( eurocom12nee_all[flag][i*12:(i+1)*12] ) )
#      
#        
#       eurocom8selectneeann.append(np.nanmean( eurocom8nee_select[flag][i*12:(i+1)*12] ) )
#       eurocom9selectneeann.append(np.nanmean( eurocom9nee_select[flag][i*12:(i+1)*12] ) )
#       eurocom10selectneeann.append(np.nanmean( eurocom10nee_select[flag][i*12:(i+1)*12] ) )
#       eurocom12selectneeann.append(np.nanmean( eurocom12nee_select[flag][i*12:(i+1)*12] ) )
     
    
    for i in range(0, 5):
       JiangOCO2neeann.append(np.nanmean( JiangOCO2nee[flag][i*12:(i+1)*12]) )   # + JiangGOSATnee3_fire[flag][i*12:(i+1)*12] ) )
    
       
    for i in range(0, 6):
       ByrneGOSATneeann_all.append(np.nanmean( ByrneGOSATnee_all[flag][i*12:(i+1)*12] ) )
       ByrneGOSATneeann_gosat.append(np.nanmean( ByrneGOSATnee_gosat[flag][i*12:(i+1)*12] ) )
       ByrneGOSATneeann_surface.append(np.nanmean( ByrneGOSATnee_surface[flag][i*12:(i+1)*12] ) )
       ByrneGOSATneeann_tccon.append(np.nanmean( ByrneGOSATnee_tccon[flag][i*12:(i+1)*12] ) )
       LiuGOSATneeann.append(np.nanmean( LiuGOSATnee[flag][i*12:(i+1)*12] ) )
       JiangGOSATneeann.append(np.nanmean( JiangGOSATnee[flag][i*12:(i+1)*12] ) )
       WuCCDASneeann1.append(np.nanmean( WuCCDASnee1[flag][i*12:(i+1)*12] ) )
       WuCCDASneeann2.append(np.nanmean( WuCCDASnee2[flag][i*12:(i+1)*12] ) )
       ScholzeCCDASneeann.append(np.nanmean( ScholzeCCDASnee[flag][i*12:(i+1)*12] ) )
       camsneeann.append(np.nanmean( camsnee[flag][i*12:(i+1)*12] ) )
     
#    for i in range(0, 6):
#       oco2mipneeann_v10mean_IS.append(np.nanmean( oco2mipnee_v10mean_IS[flag][i*12:(i+1)*12] ) )
#       oco2mipneeann_v10std_IS.append(np.nanmean( oco2mipnee_v10std_IS[flag][i*12:(i+1)*12] ) )
#       oco2mipneeann_v10mean_LNLG.append(np.nanmean( oco2mipnee_v10mean_LNLG[flag][i*12:(i+1)*12] ) )
#       oco2mipneeann_v10std_LNLG.append(np.nanmean( oco2mipnee_v10std_LNLG[flag][i*12:(i+1)*12] ) )
#       oco2mipneeann_v10mean_LNLGIS.append(np.nanmean( oco2mipnee_v10mean_LNLGIS[flag][i*12:(i+1)*12] ) )
#       oco2mipneeann_v10std_LNLGIS.append(np.nanmean( oco2mipnee_v10std_LNLGIS[flag][i*12:(i+1)*12] ) )
#  
#    for i in range(0, 6):
#       oco2mip1neeann.append(np.nanmean( oco2mip1nee[flag][i*12:(i+1)*12] ) )
#       oco2mip2neeann.append(np.nanmean( oco2mip2nee[flag][i*12:(i+1)*12] ) )
#       oco2mip3neeann.append(np.nanmean( oco2mip3nee[flag][i*12:(i+1)*12] ) )
#       oco2mip4neeann.append(np.nanmean( oco2mip4nee[flag][i*12:(i+1)*12] ) )
#       oco2mip5neeann.append(np.nanmean( oco2mip5nee[flag][i*12:(i+1)*12] ) )
#       oco2mip6neeann.append(np.nanmean( oco2mip6nee[flag][i*12:(i+1)*12] ) )
#       oco2mip7neeann.append(np.nanmean( oco2mip7nee[flag][i*12:(i+1)*12] ) )
#       oco2mip8neeann.append(np.nanmean( oco2mip8nee[flag][i*12:(i+1)*12] ) )
#       oco2mip9neeann.append(np.nanmean( oco2mip9nee[flag][i*12:(i+1)*12] ) )
#       oco2mip10neeann.append(np.nanmean( oco2mip10nee[flag][i*12:(i+1)*12] ) )
#       oco2mip11neeann.append(np.nanmean( oco2mip11nee[flag][i*12:(i+1)*12] ) )
#       oco2mip12neeann.append(np.nanmean( oco2mip12nee[flag][i*12:(i+1)*12] ) )
#       oco2mip13neeann.append(np.nanmean( oco2mip13nee[flag][i*12:(i+1)*12] ) )
#       oco2mip14neeann.append(np.nanmean( oco2mip14nee[flag][i*12:(i+1)*12] ) )                                                                    

    
    eurocom1neeann=np.array(eurocom1neeann)
    eurocom2neeann=np.array(eurocom2neeann)
    eurocom3neeann=np.array(eurocom3neeann)
    eurocom4neeann=np.array(eurocom4neeann)
    eurocom5neeann=np.array(eurocom5neeann)
    eurocom6neeann=np.array(eurocom6neeann)
    eurocom7neeann=np.array(eurocom7neeann)

#    eurocom8neeann=np.array(eurocom8neeann)
#    eurocom9neeann=np.array(eurocom9neeann)
#    eurocom9allneeann=np.array(eurocom9allneeann)
#    eurocom10neeann=np.array(eurocom10neeann)
#    eurocom11allneeann=np.array(eurocom11allneeann)
#    eurocom12neeann=np.array(eurocom12neeann)
#    eurocom12allneeann=np.array(eurocom12allneeann)
#     
#    
#    eurocom8selectneeann=np.array(eurocom8selectneeann)
#    eurocom9selectneeann=np.array(eurocom9selectneeann)
#    eurocom10selectneeann=np.array(eurocom10selectneeann)   
#    eurocom12selectneeann=np.array(eurocom12selectneeann)
                
              
    ByrneGOSATneeann_all=np.array(ByrneGOSATneeann_all)
    ByrneGOSATneeann_gosat=np.array(ByrneGOSATneeann_gosat)
    ByrneGOSATneeann_surface=np.array(ByrneGOSATneeann_surface)
    ByrneGOSATneeann_tccon=np.array(ByrneGOSATneeann_tccon)
    LiuGOSATneeann=np.array(LiuGOSATneeann)
    JiangGOSATneeann=np.array(JiangGOSATneeann)
    WuCCDASneeann1=np.array(WuCCDASneeann1)
    WuCCDASneeann2=np.array(WuCCDASneeann2)
    ScholzeCCDASneeann=np.array(ScholzeCCDASneeann)
    camsneeann=np.array(camsneeann)
    JiangOCO2neeann=np.array(JiangOCO2neeann)

#    oco2mipneeann_v10mean_IS=np.array(oco2mipneeann_v10mean_IS)
#    oco2mipneeann_v10std_IS=np.array(oco2mipneeann_v10std_IS)
#    oco2mipneeann_v10mean_LNLG=np.array(oco2mipneeann_v10mean_LNLG)
#    oco2mipneeann_v10std_LNLG=np.array(oco2mipneeann_v10std_LNLG)
#    oco2mipneeann_v10mean_LNLGIS=np.array(oco2mipneeann_v10mean_LNLGIS)
#    oco2mipneeann_v10std_LNLGIS=np.array(oco2mipneeann_v10std_LNLGIS)
              
#    oco2mip1neeann=np.array(oco2mip1neeann)
#    oco2mip2neeann=np.array(oco2mip2neeann)
#    oco2mip3neeann=np.array(oco2mip3neeann)
#    oco2mip4neeann=np.array(oco2mip4neeann)
#    oco2mip5neeann=np.array(oco2mip5neeann)
#    oco2mip6neeann=np.array(oco2mip6neeann)
#    oco2mip7neeann=np.array(oco2mip7neeann)
#    oco2mip8neeann=np.array(oco2mip8neeann)
#    oco2mip9neeann=np.array(oco2mip9neeann)
#    oco2mip10neeann=np.array(oco2mip10neeann)
#    oco2mip11neeann=np.array(oco2mip11neeann)
#    oco2mip12neeann=np.array(oco2mip12neeann)
#    oco2mip13neeann=np.array(oco2mip13neeann)
#    oco2mip14neeann=np.array(oco2mip14neeann)


    eurocom1neeann=signal.detrend(eurocom1neeann)
    eurocom2neeann=signal.detrend(eurocom2neeann)
    eurocom3neeann=signal.detrend(eurocom3neeann)
    eurocom4neeann=signal.detrend(eurocom4neeann)
    eurocom5neeann=signal.detrend(eurocom5neeann)
    eurocom6neeann=signal.detrend(eurocom6neeann)
    eurocom7neeann=signal.detrend(eurocom7neeann)

#    eurocom8neeann=signal.detrend(eurocom8neeann)
#    eurocom9neeann=signal.detrend(eurocom9neeann)
#    eurocom9allneeann=signal.detrend(eurocom9allneeann)
#    eurocom10neeann=signal.detrend(eurocom10neeann)
#    eurocom11allneeann=signal.detrend(eurocom11allneeann)
#    eurocom12neeann=signal.detrend(eurocom12neeann)
#    eurocom12allneeann=signal.detrend(eurocom12allneeann)
#    
#    eurocom8selectneeann=signal.detrend(eurocom8selectneeann)
#    eurocom9selectneeann=signal.detrend(eurocom9selectneeann)
#    eurocom10selectneeann=signal.detrend(eurocom10selectneeann)
#    eurocom12selectneeann=signal.detrend(eurocom12selectneeann)
      
 

       
    ByrneGOSATneeann_all=signal.detrend(ByrneGOSATneeann_all)
    ByrneGOSATneeann_gosat=signal.detrend(ByrneGOSATneeann_gosat)
    ByrneGOSATneeann_surface=signal.detrend(ByrneGOSATneeann_surface)
    ByrneGOSATneeann_tccon=signal.detrend(ByrneGOSATneeann_tccon)
    LiuGOSATneeann=signal.detrend(LiuGOSATneeann)
    JiangGOSATneeann=signal.detrend(JiangGOSATneeann)
    WuCCDASneeann1=signal.detrend(WuCCDASneeann1)
    WuCCDASneeann2=signal.detrend(WuCCDASneeann2)
    ScholzeCCDASneeann=signal.detrend(ScholzeCCDASneeann)
    camsneeann=signal.detrend(camsneeann)
    #JiangOCO2neeann=signal.detrend(JiangOCO2neeann)
    #oco2mipneeann_v10mean=signal.detrend(oco2mipneeann_v10mean)
    
     
#    oco2mip1neeann=signal.detrend(oco2mip1neeann)
#    oco2mip2neeann=signal.detrend(oco2mip2neeann)
#    oco2mip3neeann=signal.detrend(oco2mip3neeann)
#    oco2mip4neeann=signal.detrend(oco2mip4neeann)
#    oco2mip5neeann=signal.detrend(oco2mip5neeann)
#    oco2mip6neeann=signal.detrend(oco2mip6neeann)
#    oco2mip7neeann=signal.detrend(oco2mip7neeann)
#    oco2mip8neeann=signal.detrend(oco2mip8neeann)
#    oco2mip9neeann=signal.detrend(oco2mip9neeann)
#    oco2mip10neeann=signal.detrend(oco2mip10neeann)
#    oco2mip11neeann=signal.detrend(oco2mip11neeann)
#    oco2mip12neeann=signal.detrend(oco2mip12neeann)
#    oco2mip13neeann=signal.detrend(oco2mip13neeann)
#    oco2mip14neeann=signal.detrend(oco2mip14neeann)    
      
    #######################    
    
    eurocom7neeannnew = []
    eurocom7neeannnew.extend([0])
    eurocom7neeannnew.extend(eurocom7neeann)
    eurocom7neeannnew = np.array(eurocom7neeannnew)
    eurocom7neeannnew[0] = np.nan
      
    
    #####################
    #  plot out
    #####################
    fig = plt.figure(figsize=(6, 12))  #8, 10, 12,5*3-2
    
    fontsize=40  #32  #40 #32  #26
    mlw=1  #0.8
 
      
     
    ax1=fig.add_subplot(3,1,1)
 
#    mean = np.nanmedian([-eurocom1neeann,-eurocom2neeann,-eurocom3neeann,-eurocom4neeann,-eurocom5neeann,-eurocom6neeann,-eurocom7neeannnew], axis=0)  
#    ax4.plot(pydates1[9:],mean,linewidth=mlw*1.2,color='b',linestyle='-', marker='.', markersize = 5)   
    
#    ax4.plot(pydates1[9:],-eurocom1neeann,linewidth=mlw*1.0,color='r',linestyle='-', marker='.', markersize = 5) 
#    ax4.plot(pydates1[9:],-eurocom2neeann,linewidth=mlw*1.0,color='limegreen',linestyle='-', marker='.', markersize = 5) 
#    ax4.plot(pydates1[9:],-eurocom3neeann,linewidth=mlw*1.0,color='orange',linestyle='-', marker='.', markersize = 5) 
#    ax4.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='RoyalBlue',linestyle='-', marker='.', markersize = 5)     
#    ax4.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)     
#    ax4.plot(pydates1[9:],-eurocom6neeann,linewidth=mlw*1.0,color='magenta',linestyle='-', marker='.', markersize = 5) 
#    ax4.plot(pydates1[(9+1):],-eurocom7neeann,linewidth=mlw*2.0,color='cyan',linestyle='-', marker='.', markersize = 5) 

    mean = np.nanmean([-eurocom1neeann,-eurocom2neeann,-eurocom3neeann,-eurocom4neeann,-eurocom5neeann,-eurocom6neeann,-eurocom7neeannnew], axis=0)  
    ax1.plot(pydates1[9:],mean,linewidth=mlw*2,color='black',linestyle='-', marker='o', markersize = 3)   
 

    ax1.plot(pydates1[9:],-eurocom1neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax1.plot(pydates1[9:],-eurocom2neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax1.plot(pydates1[9:],-eurocom3neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax1.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)     
    ax1.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)     
    ax1.plot(pydates1[9:],-eurocom6neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax1.plot(pydates1[(9+1):],-eurocom7neeann,linewidth=mlw*2.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 

 
  
    #a,b = linreg(range(len(nepmean)),nepmean)

    #trendline=[a*index + b for index in range(15)]
    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
     
    #ax4.yaxis.grid(color='whitesmoke', linestyle='--')
    ax1.set_ylabel('\u0394NEP [PgC yr$^{-1}$]', fontsize=fontsize*0.35)
      
    ax1.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.25)
    ax1.locator_params(nbins=8)
    
    #ax4.legend(['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB', 'mean'], ncol = 3,  frameon=False, loc = 2, fontsize = fontsize*0.2)
    ax1.legend(['EUROCOM ensemble','individual'], ncol = 3,  frameon=False, loc = 4, fontsize = fontsize*0.3)
 
    plt.setp(ax1.get_xticklabels(), visible=True) 
    #MIN=min( np.nanmin(-cte2016neeann),np.nanmin(-cte2018neeann),np.nanmin(-ct2016neeann),np.nanmin(-ct2017neeann),np.nanmin(-camsneeann),np.nanmin(-jena99v43neeann),np.nanmin(-jena04v43neeann),np.nanmin(-jenaNEETv43neeann) ) 
    #MAX=max( np.nanmax(-cte2016neeann),np.nanmax(-cte2018neeann),np.nanmax(-ct2016neeann),np.nanmax(-ct2017neeann),np.nanmax(-camsneeann),np.nanmax(-jena99v43neeann),np.nanmax(-jena04v43neeann),np.nanmax(-jenaNEETv43neeann)  ) 
 
    #DEV = max(MAX, abs(MIN))
    #ax4.set_ylim([-DEV*1.2,DEV*1.2])
    ax1.set_ylim([-0.55,0.55])   
  
    
    ####################### 
      
    
#    ax2=fig.add_subplot(3,2,2)
# 
#
#    mean = np.nanmean([-eurocom9allneeann,-eurocom11allneeann,-eurocom12allneeann], axis=0)  
#    ax2.plot(pydates2[9:-2],mean,linewidth=mlw*2,color='k',linestyle='-', marker='o', markersize = 4)   
# 
#    mean = np.nanmean([-eurocom8selectneeann,-eurocom9selectneeann,-eurocom10selectneeann,-eurocom12selectneeann], axis=0)  
#    ax2.plot(pydates2[9:-2],mean,linewidth=mlw*2,color='k',linestyle='--', marker='o', markersize = 4)   
#    
#    mean = np.nanmean([-eurocom8neeann,-eurocom9neeann,-eurocom10neeann,-eurocom12neeann], axis=0)  
#    ax2.plot(pydates2[9:-2],mean,linewidth=mlw*2,color='gray',linestyle='--', marker='o', markersize = 4)   
#        
#    
#    ax2.plot(pydates2[9:-2],-eurocom9allneeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)     
#    ax2.plot(pydates2[9:-2],-eurocom11allneeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)         
#    ax2.plot(pydates2[9:-2],-eurocom12allneeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)     
#  
#    ax2.plot(pydates2[9:-2],-eurocom8selectneeann,linewidth=mlw*1.0,color='palegreen',linestyle='-', marker='.', markersize = 1.5) 
#    ax2.plot(pydates2[9:-2],-eurocom9selectneeann,linewidth=mlw*1.0,color='palegreen',linestyle='-', marker='.', markersize = 1.5) 
#    ax2.plot(pydates2[9:-2],-eurocom10selectneeann,linewidth=mlw*1.0,color='palegreen',linestyle='-', marker='.', markersize = 1.5) 
#    ax2.plot(pydates2[9:-2],-eurocom12selectneeann,linewidth=mlw*1.0,color='palegreen',linestyle='-', marker='.', markersize = 1.5)     
#   
#    ax2.plot(pydates2[9:-2],-eurocom8neeann,linewidth=mlw*1.0,color='mistyrose',linestyle='-', marker='.', markersize = 1.5) 
#    ax2.plot(pydates2[9:-2],-eurocom9neeann,linewidth=mlw*1.0,color='mistyrose',linestyle='-', marker='.', markersize = 1.5) 
#    ax2.plot(pydates2[9:-2],-eurocom10neeann,linewidth=mlw*1.0,color='mistyrose',linestyle='-', marker='.', markersize = 1.5) 
#    ax2.plot(pydates2[9:-2],-eurocom12neeann,linewidth=mlw*1.0,color='mistyrose',linestyle='-', marker='.', markersize = 1.5)     
 
  
#    #a,b = linreg(range(len(nepmean)),nepmean)
#
#    #trendline=[a*index + b for index in range(15)]
#    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
#     
#    #ax4.yaxis.grid(color='whitesmoke', linestyle='--')
#    ax2.set_ylabel('NEP anomaly [PgC/yr]', fontsize=fontsize*0.35)
#      
#    ax2.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.25)
#    ax2.locator_params(nbins=8)
#    
#    #ax4.legend(['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB', 'mean'], ncol = 3,  frameon=False, loc = 2, fontsize = fontsize*0.2)
#    ax2.legend(['EUROCOM v2 ensemble','individual'], ncol = 3,  frameon=False, loc = 4, fontsize = fontsize*0.3)
# 
#    plt.setp(ax2.get_xticklabels(), visible=True) 
#    #MIN=min( np.nanmin(-cte2016neeann),np.nanmin(-cte2018neeann),np.nanmin(-ct2016neeann),np.nanmin(-ct2017neeann),np.nanmin(-camsneeann),np.nanmin(-jena99v43neeann),np.nanmin(-jena04v43neeann),np.nanmin(-jenaNEETv43neeann) ) 
#    #MAX=max( np.nanmax(-cte2016neeann),np.nanmax(-cte2018neeann),np.nanmax(-ct2016neeann),np.nanmax(-ct2017neeann),np.nanmax(-camsneeann),np.nanmax(-jena99v43neeann),np.nanmax(-jena04v43neeann),np.nanmax(-jenaNEETv43neeann)  ) 
# 
#    #DEV = max(MAX, abs(MIN))
#    #ax4.set_ylim([-DEV*1.2,DEV*1.2])
#    ax2.set_ylim([-0.55,0.55])   
    
    ####################### 
    
    ax3=fig.add_subplot(3,1,2)
       
     
#    mean = np.nanmean([-eurocom1neeann,-eurocom2neeann,-eurocom3neeann,-eurocom4neeann,-eurocom5neeann,-eurocom6neeann,-eurocom7neeannnew], axis=0)  
#    std = np.nanstd([-eurocom1neeann,-eurocom2neeann,-eurocom3neeann,-eurocom4neeann,-eurocom5neeann,-eurocom6neeann,-eurocom7neeannnew], axis=0)  
#    
#    ax5.plot(pydates1[9:],mean,linewidth=mlw*1.2,color='blue',linestyle='-', marker='^', markersize = 4)   
#    ax5.fill_between(pydates1[9:], mean-std, mean+std,facecolor = "skyblue")
   
#    mean = np.nanmedian([-eurocom1neeann,-eurocom2neeann,-eurocom3neeann,-eurocom4neeann,-eurocom5neeann,-eurocom6neeann,-eurocom7neeannnew], axis=0)  
#    ax5.plot(pydates1[9:],mean,linewidth=mlw*1.2,color='k',linestyle='-', marker='.', markersize = 5)   
 
#    mean = np.nanmean([-ByrneGOSATneeann_casa,-ByrneGOSATneeann_fc,-ByrneGOSATneeann_sib], axis=0)  
#    std = np.nanstd([-ByrneGOSATneeann_casa,-ByrneGOSATneeann_fc,-ByrneGOSATneeann_sib], axis=0)  
    
#    ax5.plot(pydates1[9:],mean,linewidth=mlw*1.2,color='r',linestyle='-', marker='o', markersize = 4)   
    ax3.plot(pydates1[9:],-ByrneGOSATneeann_all,linewidth=mlw*1.2,color='r',linestyle='-', marker='o', markersize = 4)   
#    ax5.fill_between(pydates1[9:], mean-std, mean+std,facecolor = "salmon")
   
    
    
#    ax5.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='RoyalBlue',linestyle='-', marker='.', markersize = 5)     
#    ax5.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)     
#    ax5.plot(pydates1[9:],-ByrneGOSATneeann_ens,linewidth=mlw*1.0,color='r',linestyle='-', marker='.', markersize = 5)   
#    ax5.plot(pydates1[9:],-ByrneGOSATneeann_casa,linewidth=mlw*1.0,color='b',linestyle='-', marker='.', markersize = 5)   
    ax3.plot(pydates1[9:],-JiangGOSATneeann,linewidth=mlw*1.0,color='g',linestyle='-', marker='*', markersize = 5)   
    ax3.plot(pydates1[9:],-LiuGOSATneeann,linewidth=mlw*1.0,color='royalblue',linestyle='-', marker='s', markersize = 5)    #'orange'
#    ax5.plot(pydates1[9:],-camsneeann,linewidth=mlw*1.0,color='gray',linestyle='-', marker='+', markersize = 5)   
    
#    ax5.plot(pydates1[9:],WuCCDASneeann1,linewidth=mlw*1.0,color='pink',linestyle='-', marker='.', markersize = 5)   
#    ax5.plot(pydates1[9:],WuCCDASneeann2,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)   

    ax3.plot(pydates1[9:],ScholzeCCDASneeann,linewidth=mlw*1.0,color='tan',linestyle='-', marker='^', markersize = 5)   
    
  
    print("Byrne2020:", mean)
    print("CMS-Flux:", -LiuGOSATneeann)
    print("GCASv2:", -JiangGOSATneeann)
    print("CCDAS:", ScholzeCCDASneeann)

    #a,b = linreg(range(len(nepmean)),nepmean)

    #trendline=[a*index + b for index in range(15)]
    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
     
    #ax5.yaxis.grid(color='whitesmoke', linestyle='--')
    ax3.set_ylabel('\u0394NEP [PgC yr$^{-1}$]', fontsize=fontsize*0.35)
      
    ax3.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.3)
    ax3.locator_params(nbins=8)
    
    #ax4.legend(['CarboScope-Regional','LUMIA','Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.23)
#    ax5.legend(['EUROCOM mean','CarboScope-Regional','LUMIA','Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.2)
    #ax5.legend(['EUROCOM','Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.3)
    #ax5.legend(['Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.3)
    #ax5.legend(['Byrne2020','CMS-Flux2020','GCASv2','CCDAS_SM+FAPAR','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 4, fontsize = fontsize*0.3)
    ax3.legend(['Byrne2020','GCASv2','CMS-Flux2020','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 4, fontsize = fontsize*0.3)
     
    plt.setp(ax3.get_xticklabels(), visible=True) 
#    MIN=min( np.nanmin(-ByrneGOSATneeann_fc),np.nanmin(-LiuGOSATneeann),np.nanmin(-JiangGOSATneeann)) 
#    MAX=max( np.nanmax(-ByrneGOSATneeann_fc),np.nanmax(-LiuGOSATneeann),np.nanmax(-JiangGOSATneeann)) 
 
#    DEV = max(MAX, abs(MIN))
#    ax3.set_ylim([-DEV*1.2,DEV*1.2])
    ax3.set_ylim([-0.55,0.55])
  
  
    
    ####################### 
    
    ax4=fig.add_subplot(3,1,3)
 
#    mean = np.nanmean([-ByrneGOSATneeann_casa,-ByrneGOSATneeann_fc,-ByrneGOSATneeann_sib], axis=0)  
#    std = np.nanstd([-ByrneGOSATneeann_casa,-ByrneGOSATneeann_fc,-ByrneGOSATneeann_sib], axis=0)  
    
#    ax5.plot(pydates1[9:],mean,linewidth=mlw*1.2,color='r',linestyle='-', marker='o', markersize = 4)   
    ax4.plot(pydates1[9:],-ByrneGOSATneeann_all,linewidth=mlw*1.2,color='r',linestyle='-', marker='o', markersize = 4)   
#    ax5.fill_between(pydates1[9:], mean-std, mean+std,facecolor = "salmon")
   
    
    
#    ax5.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='RoyalBlue',linestyle='-', marker='.', markersize = 5)     
#    ax5.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)     
#    ax5.plot(pydates1[9:],-ByrneGOSATneeann_ens,linewidth=mlw*1.0,color='r',linestyle='-', marker='.', markersize = 5)   
#    ax5.plot(pydates1[9:],-ByrneGOSATneeann_casa,linewidth=mlw*1.0,color='b',linestyle='-', marker='.', markersize = 5)   
    ax4.plot(pydates1[9:],-ByrneGOSATneeann_surface,linewidth=mlw*1.0,color='g',linestyle='--', marker='o', markersize = 5)    #'orange'
    ax4.plot(pydates1[9:],-ByrneGOSATneeann_gosat,linewidth=mlw*1.0,color='royalblue',linestyle='--', marker='o', markersize = 5)   
    ax4.plot(pydates1[9:],-ByrneGOSATneeann_tccon,linewidth=mlw*1.0,color='tan',linestyle='--', marker='o', markersize = 5)   
#    ax5.plot(pydates1[9:],-camsneeann,linewidth=mlw*1.0,color='gray',linestyle='-', marker='+', markersize = 5)   
    
#    ax5.plot(pydates1[9:],WuCCDASneeann1,linewidth=mlw*1.0,color='pink',linestyle='-', marker='.', markersize = 5)   
#    ax5.plot(pydates1[9:],WuCCDASneeann2,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)   

#    ax5.plot(pydates1[9:],ScholzeCCDASneeann,linewidth=mlw*1.0,color='tan',linestyle='-', marker='^', markersize = 5)   
    
  
    print("Byrne2020_surface:", -ByrneGOSATneeann_surface)
    print("Byrne2020_GOSAT:", -ByrneGOSATneeann_gosat)
    print("Byrne2020_TCCON:", -ByrneGOSATneeann_tccon)
#    print("CAMS FT21r2:", -camsneeann)

    #a,b = linreg(range(len(nepmean)),nepmean)

    #trendline=[a*index + b for index in range(15)]
    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
     
    #ax5.yaxis.grid(color='whitesmoke', linestyle='--')
    ax4.set_ylabel('\u0394NEP [PgC yr$^{-1}$]', fontsize=fontsize*0.35)
      
    ax4.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.3)
    ax4.locator_params(nbins=8)
 
    #ax5.legend(['Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.3)
    #ax5.legend(['Byrne2020','CMS-Flux2020','GCASv2','CCDAS_SM+FAPAR','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 4, fontsize = fontsize*0.3)
    ax4.legend(['Byrne2020','Byrne2020_Surface','Byrne2020_GOSAT','Byrne2020_TCCON'], ncol = 2,  frameon=False, loc = 4, fontsize = fontsize*0.3)
     
    plt.setp(ax4.get_xticklabels(), visible=True) 
#    MIN=min( np.nanmin(-ByrneGOSATneeann_fc),np.nanmin(-LiuGOSATneeann),np.nanmin(-JiangGOSATneeann)) 
#    MAX=max( np.nanmax(-ByrneGOSATneeann_fc),np.nanmax(-LiuGOSATneeann),np.nanmax(-JiangGOSATneeann)) 
 
#    DEV = max(MAX, abs(MIN))
#    ax4.set_ylim([-DEV*1.2,DEV*1.2])
    ax4.set_ylim([-0.55,0.55])
  
    
    
    
   ####################### 
    
#    ax5=fig.add_subplot(3,2,5)
#       
#  
##    ax5.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='RoyalBlue',linestyle='-', marker='.', markersize = 5)     
##    ax5.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)     
##    ax5.plot(pydates1[9:],-ByrneGOSATneeann_ens,linewidth=mlw*1.0,color='r',linestyle='-', marker='.', markersize = 5)   
##    ax5.plot(pydates1[9:],-ByrneGOSATneeann_casa,linewidth=mlw*1.0,color='b',linestyle='-', marker='.', markersize = 5)   
##    ax5.plot(pydates1[9:],-LiuGOSATneeann,linewidth=mlw*1.0,color='royalblue',linestyle='-', marker='*', markersize = 5)    #'orange'
#    ax5.plot(pydates2[9+5:-1],-JiangOCO2neeann,linewidth=mlw*1.0,color='g',linestyle='-', marker='s', markersize = 5)   
##    ax5.plot(pydates2[9+5:],-oco2mipneeann_v10mean_IS,linewidth=mlw*1.0,color='gray',linestyle='-', marker='o', markersize = 5)   
#    ax5.plot(pydates2[9+5:],-oco2mipneeann_v10mean_LNLG,linewidth=mlw*1.0,color='r',linestyle='-', marker='^', markersize = 5)   
##    ax5.plot(pydates2[9+5:],-oco2mipneeann_v10mean_LNLGIS,linewidth=mlw*1.0,color='pink',linestyle='-', marker='^', markersize = 5)   
#     
#    ax5.plot(pydates2[9+5:],-oco2mip1neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip2neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip3neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip4neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip5neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip6neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip7neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip8neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#   # ax5.plot(pydates2[9+5:],-oco2mip9neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip10neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip11neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip12neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip13neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
#    ax5.plot(pydates2[9+5:],-oco2mip14neeann,linewidth=mlw*1.0,color='skyblue',linestyle='-', marker='.', markersize = 1.5)   
# 
#  
##    median = np.nanmedian([-oco2mip1neeann,-oco2mip2neeann,-oco2mip3neeann,-oco2mip4neeann,-oco2mip5neeann,-oco2mip6neeann,-oco2mip7neeann, \
##                           -oco2mip8neeann,-oco2mip9neeann,-oco2mip10neeann,-oco2mip11neeann,-oco2mip13neeann,-oco2mip13neeann,-oco2mip14neeann], axis=0)  
##
##    ax5.plot(pydates2[9+5:],median,linewidth=mlw*1.0,color='blue',linestyle='-', marker='^', markersize = 5)   
# 
#
#    print("GCASv2 OCO-2:", -JiangOCO2neeann)
#    print("OCO-2 v10 MIP:", -oco2mipneeann_v10mean_LNLG)
# 
#    #a,b = linreg(range(len(nepmean)),nepmean)
#
#    #trendline=[a*index + b for index in range(15)]
#    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
#     
#    #ax5.yaxis.grid(color='whitesmoke', linestyle='--')
#    ax5.set_ylabel('NEP anomaly [PgC/yr]', fontsize=fontsize*0.35)
#      
#    ax5.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.3)
#    ax5.locator_params(nbins=8)
#    
#    ax5.legend(['GCASv2 OCO-2', 'OCO2MIP_LNLG ensemble','individual'], ncol = 2,  frameon=False, loc = 4, fontsize = fontsize*0.3)
#   
#    plt.setp(ax5.get_xticklabels(), visible=True) 
##    MIN=min( np.nanmin(-ByrneGOSATneeann_fc),np.nanmin(-LiuGOSATneeann),np.nanmin(-JiangGOSATneeann)) 
##    MAX=max( np.nanmax(-ByrneGOSATneeann_fc),np.nanmax(-LiuGOSATneeann),np.nanmax(-JiangGOSATneeann)) 
# 
##    DEV = max(MAX, abs(MIN))
##    ax6.set_ylim([-DEV*1.2,DEV*1.2])
#    ax5.set_ylim([-0.25,0.25])
#  
# 
#    
#    ####################### 
#    
#    ax6=fig.add_subplot(3,2,6)
#       
#  
##    ax5.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='RoyalBlue',linestyle='-', marker='.', markersize = 5)     
##    ax5.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)     
##    ax5.plot(pydates1[9:],-ByrneGOSATneeann_ens,linewidth=mlw*1.0,color='r',linestyle='-', marker='.', markersize = 5)   
##    ax5.plot(pydates1[9:],-ByrneGOSATneeann_casa,linewidth=mlw*1.0,color='b',linestyle='-', marker='.', markersize = 5)   
##    ax5.plot(pydates1[9:],-LiuGOSATneeann,linewidth=mlw*1.0,color='royalblue',linestyle='-', marker='*', markersize = 5)    #'orange'
##    ax6.plot(pydates2[9+5:-1],-JiangOCO2neeann,linewidth=mlw*1.0,color='g',linestyle='-', marker='s', markersize = 5)   
#    ax6.plot(pydates2[9+5:],-oco2mipneeann_v10mean_IS,linewidth=mlw*1.0,color='royalblue',linestyle='-', marker='o', markersize = 5)   
#    ax6.plot(pydates2[9+5:],-oco2mipneeann_v10mean_LNLG,linewidth=mlw*1.0,color='r',linestyle='-', marker='^', markersize = 5)   
#    ax6.plot(pydates2[9+5:],-oco2mipneeann_v10mean_LNLGIS,linewidth=mlw*1.0,color='g',linestyle='-', marker='s', markersize = 5)   
#  
#  
#    print("OCO-2 v10 MIP IS:", -oco2mipneeann_v10mean_IS)
#    print("OCO-2 v10 MIP LNLG:", -oco2mipneeann_v10mean_LNLG)
#    print("OCO-2 v10 MIP LNLGIS:", -oco2mipneeann_v10mean_LNLGIS)
#    #a,b = linreg(range(len(nepmean)),nepmean)
#
#    #trendline=[a*index + b for index in range(15)]
#    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
#     
#    #ax5.yaxis.grid(color='whitesmoke', linestyle='--')
#    ax6.set_ylabel('NEP anomaly [PgC/yr]', fontsize=fontsize*0.35)
#      
#    ax6.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.3)
#    ax6.locator_params(nbins=8)
#    
#    ax6.legend(['OCO2MIP_IS','OCO2MIP_LNLG','OCO2MIP_LNLGIS'], ncol = 2,  frameon=False, loc = 4, fontsize = fontsize*0.3)
#   
#    plt.setp(ax6.get_xticklabels(), visible=True) 
##    MIN=min( np.nanmin(-ByrneGOSATneeann_fc),np.nanmin(-LiuGOSATneeann),np.nanmin(-JiangGOSATneeann)) 
##    MAX=max( np.nanmax(-ByrneGOSATneeann_fc),np.nanmax(-LiuGOSATneeann),np.nanmax(-JiangGOSATneeann)) 
# 
##    DEV = max(MAX, abs(MIN))
##    ax6.set_ylim([-DEV*1.2,DEV*1.2])
#    ax6.set_ylim([-0.25,0.25])
    

 ####################### 
    
#    ax7=fig.add_subplot(4,1,4)
#   
##    ax5.plot(pydates1[9:],-eurocom4neeann,linewidth=mlw*1.0,color='RoyalBlue',linestyle='-', marker='.', markersize = 5)     
##    ax5.plot(pydates1[9:],-eurocom5neeann,linewidth=mlw*1.0,color='gold',linestyle='-', marker='.', markersize = 5)     
##    ax5.plot(pydates1[9:],-ByrneGOSATneeann_ens,linewidth=mlw*1.0,color='r',linestyle='-', marker='.', markersize = 5)   
##    ax5.plot(pydates1[9:],-ByrneGOSATneeann_casa,linewidth=mlw*1.0,color='b',linestyle='-', marker='.', markersize = 5)   
##    ax5.plot(pydates1[9:],-LiuGOSATneeann,linewidth=mlw*1.0,color='royalblue',linestyle='-', marker='*', markersize = 5)    #'orange'
##    ax5.plot(pydates1[9:],-JiangGOSATneeann,linewidth=mlw*1.0,color='g',linestyle='-', marker='s', markersize = 5)   
##    ax5.plot(pydates1[9:],-camsneeann,linewidth=mlw*1.0,color='gray',linestyle='-', marker='+', markersize = 5)   
#    
##    ax6.plot(pydates1[9:],WuCCDASneeann1,linewidth=mlw*1.0,color='pink',linestyle='-', marker='^', markersize = 5)   
#    ax7.plot(pydates1[9:],WuCCDASneeann2,linewidth=mlw*1.0,color='magenta',linestyle='-', marker='^', markersize = 5)   
#    ax7.plot(pydates1[9:],ScholzeCCDASneeann,linewidth=mlw*1.0,color='tan',linestyle='-', marker='^', markersize = 5)   
#    
#  
# 
#    print("CCDAS_CO2+SM+FAPAR", WuCCDASneeann2)
#    print("CCDAS_CO2+SM+VOD:", ScholzeCCDASneeann)
# 
#
#    #a,b = linreg(range(len(nepmean)),nepmean)
#
#    #trendline=[a*index + b for index in range(15)]
#    #ax4.plot(pydates1,trendline,linewidth=mlw,color='blue',linestyle='--') 
#     
#    #ax5.yaxis.grid(color='whitesmoke', linestyle='--')
#    ax7.set_ylabel('NEP anomaly [PgC/yr]', fontsize=fontsize*0.35)
#      
#    ax7.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.3)
#    ax7.locator_params(nbins=8)
#    
#    #ax4.legend(['CarboScope-Regional','LUMIA','Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.23)
##    ax5.legend(['EUROCOM mean','CarboScope-Regional','LUMIA','Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.2)
#    #ax5.legend(['EUROCOM','Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.3)
#    #ax5.legend(['Byrne2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.3)
#    ax7.legend(['CCDAS_CO2+SM+FAPAR','CCDAS_CO2+SM+VOD'], ncol = 1,  frameon=False, loc = 4, fontsize = fontsize*0.3)
#   
#    plt.setp(ax7.get_xticklabels(), visible=True) 
#    MIN=min( np.nanmin(-ByrneGOSATneeann_fc),np.nanmin(-LiuGOSATneeann),np.nanmin(-JiangGOSATneeann)) 
#    MAX=max( np.nanmax(-ByrneGOSATneeann_fc),np.nanmax(-LiuGOSATneeann),np.nanmax(-JiangGOSATneeann)) 
# 
#    DEV = max(MAX, abs(MIN))
#    ax7.set_ylim([-DEV*1.2,DEV*1.2])
#    ax7.set_ylim([-0.55,0.55])
  
       
    
    yy = np.zeros(len(pydates1))
    yy2 = np.zeros(len(pydates2))
    
    if flag == 1:
      ax1.plot(pydates1[9:],yy[9:], color='lightgray',linestyle='--',linewidth=0.5)
#      ax2.plot(pydates2[9:-2],yy2[9:-2], color='lightgray',linestyle='--',linewidth=0.5)
      ax3.plot(pydates1[9:],yy[9:], color='lightgray',linestyle='--',linewidth=0.5)
      #ax4.plot(pydates1[9:],yy[9:], color='lightgray',linestyle='--',linewidth=0.5)
      ax4.plot(pydates1[9:],yy[9:], color='lightgray',linestyle='--',linewidth=0.5)
#      ax5.plot(pydates2[9+5:],yy2[9+5:], color='lightgray',linestyle='--',linewidth=0.5)
#      ax6.plot(pydates2[9+5:],yy2[9+5:], color='lightgray',linestyle='--',linewidth=0.5)
#      ax7.plot(pydates1[9:],yy[9:], color='lightgray',linestyle='--',linewidth=0.5)
        
    ax1.xaxis.set_major_locator(pltdt.YearLocator())
    ax1.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #
#    ax2.xaxis.set_major_locator(pltdt.YearLocator())
#    ax2.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #
    ax3.xaxis.set_major_locator(pltdt.YearLocator())
    ax3.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #
    ax4.xaxis.set_major_locator(pltdt.YearLocator())
    ax4.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #
#    ax5.xaxis.set_major_locator(pltdt.YearLocator())
#    ax5.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #
#    ax6.xaxis.set_major_locator(pltdt.YearLocator())
#    ax6.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #
#    ax7.xaxis.set_major_locator(pltdt.YearLocator())
#    ax7.xaxis.set_major_formatter(pltdt.DateFormatter('%Y')) #

    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax1.get_xticklabels()]
    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax1.get_yticklabels()]
    
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax2.get_xticklabels()]
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax2.get_yticklabels()]
        
    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax3.get_xticklabels()]
    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax3.get_yticklabels()]
    
    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax4.get_xticklabels()]
    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax4.get_yticklabels()]
     
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax5.get_xticklabels()]
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax5.get_yticklabels()]
#      
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax6.get_xticklabels()]
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax6.get_yticklabels()]
       
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax7.get_xticklabels()]
#    dummy = [lab.set_fontsize(0.3*fontsize) for lab in ax7.get_yticklabels()]
         
    
    fig = plt.gcf()
    fig.set_size_inches(6, 3*3) # (6, 4.5)  #(6,10)
    
    
    fig.tight_layout() 
    plt.subplots_adjust(hspace = 0.15)
          
    #fn='Multiple_flux_timeseries_monthly'
    outp=os.path.join(dir22,fn+'.png')
    print(outp)
    plt.savefig(outp,dpi=300)
    
     
    data=[]
    
    data.append(-eurocom1neeann)
    data.append(-eurocom2neeann) 
    data.append(-eurocom3neeann) 
    data.append(-eurocom4neeann) 
    data.append(-eurocom5neeann) 
    data.append(-eurocom6neeann) 
    data.append(-eurocom7neeann) 

#    data.append(-eurocom8neeann) 
#    data.append(-eurocom9neeann) 
#    data.append(-eurocom9allneeann) 
#    data.append(-eurocom10neeann) 
#    data.append(-eurocom11allneeann) 
#    data.append(-eurocom12neeann) 
#    data.append(-eurocom12allneeann) 
     
    data.append(-ByrneGOSATneeann_all)  
#    data.append(-ByrneGOSATneeann_surface)  
#    data.append(-ByrneGOSATneeann_gosat)  
#    data.append(-ByrneGOSATneeann_tccon)  
    data.append(-LiuGOSATneeann)  
    data.append(-JiangGOSATneeann)  
#    data.append(WuCCDASneeann1)    
#    data.append(WuCCDASneeann2)   
    data.append(ScholzeCCDASneeann) 
#    data.append(camsneeann) 
      
    data=np.array(data)
    
    return data    

   
import timeit
start = timeit.timeit()

flag=1
data = pipeline(35, 70, -10, 40,flag,'All_Flux_timeseries_Europe_AIMs_20211201_EUROCOM+GOSAT2+GCASv2+CCDAS2_final')
 
 
end = timeit.timeit()
print(end - start)

print(data)     


#data

#np.array([-0.07796869,  0.16029364, -0.14973553,  0.17016346, -0.14245146,
#        0.03969857]),
#       np.array([ 0.03635836,  0.09542384, -0.16956765, -0.0642361 ,  0.07368795,
#        0.02833359]),
#       np.array([ 0.05698552,  0.02557602, -0.08046738, -0.3413053 ,  0.53678107,
#       -0.19756992]),
#       np.array([-0.10844791,  0.18667693, -0.15973426,  0.06867671,  0.1373812 ,
#       -0.12455267]),
#       np.array([ 0.13075221, -0.10266092, -0.12046228, -0.00315861,  0.12458666,
#       -0.02905707]),
#       np.array([ 0.23153015, -0.11316518, -0.13655968, -0.14684059, -0.00162981,
#        0.16666511]),
#       np.array([ 0.1477639 , -0.12292061, -0.17783597,  0.13337816,  0.01961452]),
#      
#        
#       np.array([ 0.07309177, -0.08063309, -0.13566975,  0.12901775,  0.10604726,
#       -0.09185394]),
#       np.array([ 0.01930987, -0.07348228, -0.06161834,  0.18161103,  0.01901275,
#       -0.08483303]),
#       np.array([-0., -0., -0., -0., -0., -0.]),
#       
#       
#       np.array([-0.01296051,  0.04737569, -0.17032468,  0.11297656,  0.16032072,
#       -0.13738777]),
#       np.array([-0.12920947,  0.17686924, -0.11287834, -0.0060121 ,  0.2892296 ,
#       -0.21799893])



#2022-10-08
#[array([-0.07796869,  0.16029364, -0.14973553,  0.17016346, -0.14245146,
#        0.03969857])
# array([ 0.03635836,  0.09542384, -0.16956765, -0.0642361 ,  0.07368795,
#        0.02833359])
# array([ 0.05698552,  0.02557602, -0.08046738, -0.3413053 ,  0.53678107,
#       -0.19756992])
# array([-0.10844791,  0.18667693, -0.15973426,  0.06867671,  0.1373812 ,
#       -0.12455267])
# array([ 0.13075221, -0.10266092, -0.12046228, -0.00315861,  0.12458666,
#       -0.02905707])
# array([ 0.23153015, -0.11316518, -0.13655968, -0.14684059, -0.00162981,
#        0.16666511])
# array([ 0.1477639 , -0.12292061, -0.17783597,  0.13337816,  0.01961452])
# array([ 0.07188867, -0.0792685 , -0.13504031,  0.12854501,  0.10566152,
#       -0.09178641])
# array([-0.04233668,  0.04380452, -0.13279528,  0.23493794, -0.03502472,
#       -0.06858578])
# array([-0.01296051,  0.04737569, -0.17032468,  0.11297656,  0.16032072,
#       -0.13738777])
# array([-0.12920947,  0.17686924, -0.11287834, -0.0060121 ,  0.2892296 ,
#       -0.21799893])]
    