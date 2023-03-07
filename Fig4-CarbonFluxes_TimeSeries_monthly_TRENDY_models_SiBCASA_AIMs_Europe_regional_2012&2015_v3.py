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

from scipy import signal

 
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
#   TRENDYv6
#################################################################
 
#dir1='/Volumes/HEWEI_T5/TRENDY_v6_s3/'

#def flux_anomaly_TRENDYv6(item, varstr, M, N, latMin, latMax, lonMin, lonMax):
#     
##    latMin=33    #35
##    latMax=73    #70
##    lonMin=-15 #-10
##    lonMax=35 
#    
##    M=1
##    N=15
##    
##    item='NEP'
##    varstr='ensemble_mean'
#    
#    file= os.path.join(dir1,'TRENDY.v6.model.ensemble.'+item+'.0.5Deg.2000-2016.nc')
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
#    data2 = []
#    
#    for i in range(0,data1.shape[0]): 
#        data2.append( np.transpose(data1[i]))   
#     
#    data2 = np.array(data2)
#   
#    
#    area = globarea(720,360,True)
#    
#    lon=np.arange(0,360,0.5)
#    lat=np.arange(0,180,0.5)
#    lon2, lat2 = np.meshgrid(lon,lat) 
#    
#   # for i in range(0,data1.shape[0]):
#   #     data2[i] = maskoceans(lon2, lat2, data1[i] )  #* area)
# 
#    data3 = data2[12*M:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  # Corrected July 3, 2017
#    #print data2.shape, data3.shape
#     
#    #plt.imshow(data3[0])
#    
#    
#    a=[]
#    
#    mask=(data3==9.96921e+36)
#    data3[mask]=np.NaN
#    
#    plt.imshow(data3[0])
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]< 0 :
##                   data3[i,xx,yy]=np.NaN
##               else:      # kg C m-2 s-1
##                   data3[i,xx,yy]=data3[i,xx,yy]*365*24*3600*1e-15*1e3*1e3
# 
#    
#    #data3=data3*area[2*latMin:2*latMax,2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
#    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
#         
#    del data1
#    del data2
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
#    A=a[:-12*1] 
#    B=a[:-12*1]-a1
#    
#    return A,B,baseline

 
#################################################################
#   ORCHIDEE
#################################################################
 
dir1='/Volumes/HEWEI_T5/TRENDY_v6_s3/'

def flux_anomaly_orchidee(varstr, M, N, latMin, latMax, lonMin, lonMax):
     
#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-15 #-10
#    lonMax=35 
#    
#    M=0
#    N=15
#    
#    varstr='gpp'
    
    file= os.path.join(dir1,'ORCHIDEE_S3_'+varstr+'.nc')
    
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
   # for i in range(0,data1.shape[0]):
   #     data2[i] = maskoceans(lon2, lat2, data1[i] )  #* area)
 
    data3 = data2[140*12:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  # Corrected July 3, 2017
    #print data2.shape, data3.shape
     
    #plt.imshow(data3[0])
    
    
    a=[]
    
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
    
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(M,M+N):
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
       
    
    del data3
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline

#################################################################
#   ORCHIDEE-MICT
#################################################################
 
 
def flux_anomaly_orchideemict(varstr, M, N, latMin, latMax, lonMin, lonMax):
      
#    M=6
#    N=8  
#    varstr='gpp'
    
    file= os.path.join(dir1,'ORCHIDEE-MICT_S3_'+varstr+'.nc')
    
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(360,180,True)
    
   # lon=np.arange(0.5,360,1)
   # lat=np.arange(0.5,180,1)
   # lon2, lat2 = np.meshgrid(lon,lat) 
    
   # for i in range(0,data1.shape[0]):
   #     data2[i] = maskoceans(lon2, lat2, data1[i] )  #* area)
 
    data3 = data2[141*12:,1*(90-latMax):1*(90-latMin),1*(180+lonMin):1*(180+lonMax)]  # Corrected July 3, 2017
    #print data2.shape, data3.shape
     
    plt.imshow(data3[0])
    
    
    a=[]
    
    mask=(data3==-99999)
    data3[mask]=np.NaN
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]< 0 :
#                   data3[i,xx,yy]=np.NaN
#               else:      # kg C m-2 s-1
#                   data3[i,xx,yy]=data3[i,xx,yy]*365*24*3600*1e-15*1e3*1e3
    
    #plt.imshow(data3[0])
    
    data3=data3*area[1*(90-latMax):1*(90-latMin),1*(180+lonMin):1*(180+lonMax)]*365*24*3600*1e-15*1e3
         
    del data1
    del data2
    
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(M,M+N):
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
       
    
    del data3
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline

#################################################################
#   CABLE
#################################################################

dir2=dir1


def flux_anomaly_cable(varstr, M, N, latMin, latMax, lonMin, lonMax):

    #latMin=10
    #latMax=80
    #lonMin=-175
    #lonMax=-55
    #        
   # M=6
   # N=8  
   # varstr='gpp'
            
    file=os.path.join(dir2,'CABLE_S3_'+varstr+'.nc')
 
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
    data3 = data2[141*12:, 2*(90-latMax):2*(90-latMin) ,2*(180+lonMin):2*(180+lonMax)]   #Corrected July 3, 2017
    #print data2.shape, data3.shape
    
    del data1
    del data2
    #plt.imshow(data3[0])
    
    
    a=[]
        
    mask=(data3==-1.e+33)
    data3[mask]=np.NaN
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]< 0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    plt.imshow(data3[0])
    
    #data3=data3*area[2*latMin:2*latMax ,2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
    data3=data3*area[2*(90-latMax):2*(90-latMin) ,2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(M,M+N):  
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N) #float(M+N)
       
    del data3
        
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline

    
#################################################################
#   DLEM
#################################################################
 

dir3=dir1


def flux_anomaly_dlem(varstr, M, N, latMin, latMax, lonMin, lonMax):

    #latMin=10
    #latMax=80
    #lonMin=-175
    #lonMax=-55
    #        
    #M=6
    #N=8  
    #varstr='gpp'
            
    file=os.path.join(dir3,'DLEM_S3_'+varstr+'.nc')
 
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    #data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
    data2=[]

    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)  
    
    data3 = data2[100*12:, 2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]

    
    del data1
    del data2
     
    #plt.imshow(data3[0])
    
    
    a=[]
        
    mask=(data3==-99999.0)
    data3[mask]=np.NaN

   
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==-99999.0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*365*24*3600*1e-15*1e3*1e3
#    
    plt.imshow(data3[0])
    
    #data3=data3*area[2*(180-latMax):2*(180-latMin),2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(M,M+N): 
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
       
    del data3
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline
  

#################################################################
#   ISAM
#################################################################
 
dir4=dir1
 

def flux_anomaly_isam(varstr, M, N, latMin, latMax, lonMin, lonMax):  #speacial geographic extent

#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-15 #-10
#    lonMax=35 
#           
#    M=0
#    N=15  
#    varstr='gpp'
            
    file=os.path.join(dir4,'ISAM_S3_'+varstr+'.nc')
     
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    #lon=np.arange(0,360,0.5)
    #lat=np.arange(0,180,0.5)
    #lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
    data2=[]
    
    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)  
    
    #data3 = data2[139:,:, 2*(180-latMax):2*(180-latMin),2*(360+lonMin):2*(360+lonMax)]
    data3aa = data2[141*12:, 2*(90-latMax):2*(90-latMin),2*(360+lonMin):]
    data3bb = data2[141*12:, 2*(90-latMax):2*(90-latMin),0:2*lonMax]
    
    del data1
    del data2
     
    #plt.imshow(data3[0])
    
    
    a=[]
    
    mask=(data3aa==-9999.0)
    data3aa[mask]=np.NaN
    
    mask=(data3bb==-9999.0)
    data3bb[mask]=np.NaN
   
    plt.imshow(data3aa[0])
    
    #data3=data3*area[2*(180-latMax):2*(180-latMin),2*(360+lonMin):2*(360+lonMax)]*12.0*1e-15*1e3   #kgC/m2/month
    #data3=data3*area[2*(90-latMax):2*(90-latMin),2*(360+lonMin):2*(360+lonMax)]*12.0*1e-15*1e3   #kgC/m2/month
    data3aa=data3aa*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*180]*365*24*3600*1e-15*1e3   #   #kg C m-2 s-1 365*24*3600
    data3bb=data3bb*area[2*(90-latMax):2*(90-latMin),2*180:2*(180+lonMax)]*365*24*3600*1e-15*1e3
        
    
    for i in range(0,data3aa.shape[0]): 
        #for j in range(0,data3.shape[1]):
        me = np.nansum(data3aa[i]) +  np.nansum(data3bb[i])
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,M+N):  
            baseline[i] = baseline[i] + np.nansum(data3aa[12*j+i])/float(N) + np.nansum(data3bb[12*j+i])/float(N)  #float(M+N)
       
    del data3aa
    del data3bb
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline
 
    
#################################################################
#   LPX-Bern
#################################################################
 
#
#dir5='/Users/weihe/Documents/LPX-Bern/S2/'
#
#
#def flux_anomaly_lpxbern(varstr, M, N, latMin, latMax, lonMin, lonMax):
#
#    #latMin=10
#    #latMax=80
#    #lonMin=-175
#    #lonMax=-55
#    #        
#    #M=6
#    #N=8  
#    #varstr='gpp'
#            
#    file=os.path.join(dir5,'LPX_S2_'+varstr+'.nc')
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
#    data2=[]
#    
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
#    
#    data2=np.array(data2)
#
#    #plt.imshow(data2[0])
#    #data3 = data2[139*12:, 180-latMax:180-latMin,180+lonMin:180+lonMax]
#    data3 = data2[139*12:, 90-latMax:90-latMin,180+lonMin:180+lonMax]
#    
#    #print data2.shape, data3.shape
#    
#    del data1
#    del data2
#     
#    #plt.imshow(data3[0])
#    
#    
#    
#    a=[]
#        
#    mask=(data3==-99999.0)
#    data3[mask]=np.NaN
#    
##    for i in range(0,data3.shape[0]): 
##        for xx in range(0,data3[i].shape[0]) :
##            for yy in range(0,data3[i].shape[1]) :
##               if data3[i,xx,yy]== -99999.0 :
##                   data3[i,xx,yy]=np.NaN
##               else:
##                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
#    
#    plt.imshow(data3[0])
#    
#    #data3=data3*area[180-latMax:180-latMin,180+lonMin:180+lonMax]*365*24*3600*1e-15*1e3
#    data3=data3*area[90-latMax:90-latMin,180+lonMin:180+lonMax]*365*24*3600*1e-15*1e3
#      
#      
#    for i in range(0,data3.shape[0]): 
#        me = np.nansum(data3[i])   
#        a.append(me)
#      
#    baseline=np.zeros(12)
#    for i in range(0,12):
#        #for j in range(0,M+N):  
#        for j in range(9-1,9-1+N): 
#            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N) #float(M+N)
#       
#    del data3
#    #a1=np.repeat(baseline,N)
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=np.array(a1)
#     
#    A=a[12*M:-12] 
#    B=a[12*M:-12]-a1
#    
#    return A,B,baseline
#    
#################################################################
#   VEGAS
#################################################################
 

dir6=dir1


def flux_anomaly_vegas(varstr, M, N, latMin, latMax, lonMin, lonMax):   #speacial geographic extent

#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-15 #-10
#    lonMax=35 
#             
#    M=0
#    N=15  
#    varstr='gpp'
            
    file=os.path.join(dir6,'VEGAS_S3_'+varstr+'.nc')
 
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
#    data2=[]
    
#    m=data1.shape[0]
#    for i in range(0,m):
#        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)  
        
    #data3 = data2[139*12:, 2*(180-latMax):2*(180-latMin),2*(360+lonMin):2*(360+lonMax)]
    data3aa = data2[141*12:, 2*(90-latMax):2*(90-latMin),2*(360+lonMin):]
    data3bb = data2[141*12:, 2*(90-latMax):2*(90-latMin),0:2*lonMax]
    
    
    del data1
    del data2
    #plt.imshow(data3[0])
    
    a=[]
           
    mask=(data3aa==-99999.0)
    data3aa[mask]=np.NaN
    
    mask=(data3bb==-99999.0)
    data3bb[mask]=np.NaN

#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==-99999.0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3aa[0])
    
    #data3=data3*area[2*(180-latMax):2*(180-latMin),2*(360+lonMin):2*(360+lonMax)]*365*24*3600*1e-15*1e3
    data3aa=data3aa*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*180]*365*24*3600*1e-15*1e3
    data3bb=data3bb*area[2*(90-latMax):2*(90-latMin),2*180:2*(180+lonMax)]*365*24*3600*1e-15*1e3
        
      
    for i in range(0,data3aa.shape[0]): 
        me = np.nansum(data3aa[i]) +  np.nansum(data3bb[i])
        a.append(me)
       
    baseline=np.zeros(12)
    for i in range(0,12):
        for j in range(0,M+N):  
            baseline[i] = baseline[i] + np.nansum(data3aa[12*j+i])/float(N)  + np.nansum(data3bb[12*j+i])/float(N) #float(M+N)
       
    del data3aa
    del data3bb
 
 
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline
    
#################################################################
#   VISIT
#################################################################
 
dir7=dir1


def flux_anomaly_visit(varstr, M, N, latMin, latMax, lonMin, lonMax):

    latMin=33    #35
    latMax=73    #70
    lonMin=-15 #-10
    lonMax=35 
             
    M=0
    N=15  
    varstr='gpp'
            
    file=os.path.join(dir7,'VISIT_S3_'+varstr+'.nc')
 
    data1=[]
    
     
    ncf=nc.Dataset(file)
    data1=ncf.variables[varstr][:]
         
    ncf.close()
    
    data1=np.array(data1)
    data2 = data1
    area = globarea(720,360,True)
    
    lon=np.arange(0,360,0.5)
    lat=np.arange(0,180,0.5)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
     
    #data3 = data2[139*12:, 2*latMin:2*latMax,2*(180+lonMin):2*(180+lonMax)]
    data3 = data2[141*12:, 2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]
  
    del data1
    del data2
    #plt.imshow(data3[0])
    
    
    a=[]
        
    mask=(data3==-99999.0)
    data3[mask]=np.NaN

   
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==-99999.0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*365*24*3600*1e-15*1e3*1e3
#    
    plt.imshow(data3[0])
    
    #data3=data3*area[2*latMin:2*latMax,2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]*365*24*3600*1e-15*1e3
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(M,M+N): 
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
       
    del data3   
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:-12*1] 
    B=a[12*M:-12*1]-a1
    
    return A,B,baseline
  

#################################################################
#   SIBCASA
#################################################################
dir8='/Volumes/HEWEI_T5/forXiaocui/sibcasa_fluxes_monthly1×1/2001-2015/'  
 
def flux_anomaly_sibcasa(varstr, M, N, latMin, latMax, lonMin, lonMax):

    #M=0 
    #N=15     
    #latMin=35
    #latMax=70
    #lonMin=-10
    #lonMax=40
    #varstr = 'GPP'

    files=glob.glob(os.path.join(dir8,'SiBCASA.monthly.flux1x1.global*.nc')) 
    
    data1=[]
        
    for file in files:
         ncf=nc.Dataset(file)
         data1.append(ncf.variables[varstr][:])
         ncf.close()
    
    data1=np.array(data1)
    #data2 = data1
    area = globarea(360,180,True)
     
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
    data2=[]
    
    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)  
    
    #data3 = data2[:, (180-latMax):(180-latMin), (180+lonMin):(180+lonMax) ] 
    data3 = data2[:, (90-latMax):(90-latMin), (180+lonMin):(180+lonMax) ]    # Corrected July 3, 2017
           
    #plt.imshow(data3[0])
    
    a=[]
    
    #mask=(data3==0)
    #data3[mask]=np.NaN
  
    #plt.imshow(data3[0])
    
    #data3=data3*area[180-latMax:180-latMin, (180+lonMin):(180+lonMax)]*12*365.0*24*3600*1e-6*1e-15 
    data3=data3*area[(90-latMax):(90-latMin), (180+lonMin):(180+lonMax)]*12*365.0*24*3600*1e-6*1e-15 
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,14):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/14.0  #float(M+N)
    
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:] 
    B=a[12*M:]-a1
    
    return A,B


#################################################################
#   BEPS
#################################################################
#dir9='/Volumes/WEIHE/Global_0.5CarbonFlux/GPP/nc/'
  
def flux_anomaly_beps(varstr, M, N, latMin, latMax, lonMin, lonMax):
        
#    M=0
#    N=15     
#    latMin=10
#    latMax=80
#    lonMin=-175
#    lonMax=-55
#    varstr = 'GPP'
    
    dir9 = os.path.join('/Volumes/HEWEI_WD2/Global_0.5CarbonFlux',varstr,"2001-2015/")
    files=glob.glob(os.path.join(dir9,'BEPS.fluxes.halfDeg.global.*.nc'))  
    
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
    #data2 = data1
    area = globarea(720,360,True)
    
    lat=np.array(lat)
    lon=np.array(lon)
    
    #lon=np.arange(0,360,1)
    #lat=np.arange(0,180,1)
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    #for i in range(0,data1.shape[0]):
    #    data2[i] = maskoceans(lon2, lat2, data1[i] * area)
     
    data2=[]
    
    m=data1.shape[0]
    for i in range(0,m):
        data2.append(np.flipud(data1[i]))
    
    data2=np.array(data2)  
    
    #data3 = data2[:, (180-latMax):(180-latMin), (180+lonMin):(180+lonMax) ] 
    data3 = data2[:, 2*(90-latMax):2*(90-latMin), 2*(180+lonMin):2*(180+lonMax) ]    # Corrected July 3, 2017
           
#    plt.imshow(data3[0])
    
    a=[]
    
    #mask=(data3==0)
    #data3[mask]=np.NaN
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0 :
#                   data3[i,xx,yy]=np.NaN
#               else:       #umol/m^2/s
#                   data3[i,xx,yy]=data3[i,xx,yy]*12*365.0*24*3600*1e-6*1e-15*1e3
#    
    plt.imshow(data3[0])
    
    #data3=data3*area[180-latMax:180-latMin, (180+lonMin):(180+lonMax)]*12*365.0*24*3600*1e-6*1e-15 
    data3=data3*area[2*(90-latMax):2*(90-latMin), 2*(180+lonMin):2*(180+lonMax)]*12*1e-15 
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,14):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/14.0  #float(M+N)
    
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:] 
    B=a[12*M:]-a1
    
    return A,B



 

#################################################################
#   CTE2018-flux1*1
#################################################################

dir11='/Volumes/HEWEI_T5/Inversions/'

def flux_anomaly_cte2018(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=1 
#    N=17  
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
    

    a=[]
    
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
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,15):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,15):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
   
 


#################################################################
#   CT2019-flux1*1
#################################################################

dir13b='/Volumes/HEWEI_T5/Inversions/'
 
def flux_anomaly_ct2019b(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=1 
#    N=15  
#    varstr='bio_flux_opt'
    
    file= os.path.join(dir13b,'CT2019B.flux1x1-monthly.nc')
    
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
    

    a=[]
    
    mask=(data3==-1.0E34)
    data3[mask]=np.NaN
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
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,15):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,15):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a
    B=a-a1
    
    return A,B


#################################################################
#   CAMS-flux1*1 
#################################################################

#dir14='/Users/weihe/Documents/Research_Projects/Drought_Carbon_Project/'
dir14='/Volumes/HEWEI_T5/Inversions/CAMS_v18r2/1Deg/'

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
    
    files = glob.glob(os.path.join(dir14,'Resampled_1Deg*.nc')) 
    
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
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN
    
    # now kgC m-2 month-1
    data3=data3*12*1e3*1e-15   
   
    data3=data3*area[(90-latMax):(90-latMin),(180+lonMin):(180+lonMax)]  
    
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,15):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,15):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a
    B=a-a1
    
    return A,B



##################################################################
##   MACTM-flux1*1
##################################################################
#
#dir15='/Volumes/HEWEI_WD/Inverions/'
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


##################################################################
##   jenas99
##################################################################
#
#dir16='/Volumes/HEWEI_WD/Inverions/Jena_CarboScope/monthlys99v3.8/2001-2015/'
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

dir16a='/Volumes/HEWEI_T5/Inversions/Jena_CarboScope/monthlys99v4.3/2001-2015/'
 
def flux_anomaly_jenas99v43(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=15  
#    varstr='co2flux_land'
    
    
    files=glob.glob(os.path.join(dir16a,'CarboScope_s99_v4.3_monthly_72x48_*.nc')) 
    
    
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
    

    a=[]
    
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
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,15):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,15):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
   
  

#################################################################
#   jenas04
#################################################################

#dir17='/Volumes/HEWEI_WD/Inverions/Jena_CarboScope/monthlys04v4.1/'
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
#   

#################################################################
#   jenas04
#################################################################

#dir17a='/Volumes/HEWEI_WD/Inverions/Jena_CarboScope/monthlys04v4.3/'
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

dir18='/Volumes/HEWEI_T5/Inversions/Jena_CarboScope/monthlysEXTocNEETv4.3/2001-2015/'


def flux_anomaly_jenaNEETv43(varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=15  
#    varstr='co2flux_land'
    
    
    files=glob.glob(os.path.join(dir18,'CarboScope_NEET_v4.3_monthly_72x48_*.nc')) 
    
    
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
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN
    
    #data3=data3*12.0*12*1e12*1e-15  
    data3=data3
    
#    for i in range(0,data3.shape[0]): 
#        for xx in range(0,data3[i].shape[0]) :
#            for yy in range(0,data3[i].shape[1]) :
#               if data3[i,xx,yy]==0 :
#                   data3[i,xx,yy]=np.NaN
#               else:
#                   data3[i,xx,yy]=data3[i,xx,yy]*12.0*365*24*3600*1e-15*1e3
    
    #plt.imshow(data3[0])
    
#    data3=data3*area[int(np.floor((90-latMax)/3.75)):int(np.ceil((90-latMin)/3.75)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]  
      
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,15):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/15.0  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,15):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
   
  
#################################################################
#   jenas NEETW v4.3
#################################################################

#dir19='/Volumes/HEWEI_WD/Inverions/Jena_v4.3_all/NEE-T-W/monthly/2001-2015/'
dir19='/Volumes/HEWEI_T5/Inversions/Jena_v4.3_all/NEE-T-W-new/s99/monthly/2001-2015/'
  

def flux_anomaly_jenaNEETWv43(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    M=9 
#    N=6  
#    varstr='co2flux_land'
    
    
    files=glob.glob(os.path.join(dir19,'CarboScope_NEETW_v4.3_monthly_72x48_*.nc')) 
    
    
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
      
     
    a=[]
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,N):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B
 
    


#################################################################
#   FluxSat GPP
#################################################################

dir20='/Volumes/HEWEI_T5/FluxSat/0.5x0.5/2001-2015/'
 
  
def flux_anomaly_fluxsat(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    M=0
#    N=15     
#    varstr = 'GPP'
    
    files=glob.glob(os.path.join(dir20,'*/FluxSat_GPP_0.5_v1.1*.nc'))  
    files=sorted(files)
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(files[0])
    lat.extend(ncf.variables['Latitude'][:])
    lon.extend(ncf.variables['Longitude'][:])
         
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
          
      
    a=[]
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,N):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B

 

#################################################################
#   GOSIF GPP v2 (Xiao Jingfeng)
#################################################################

dir21='/Volumes/HEWEI_T5/SIF_XiaoJF/nc/GOSIF-GPP_v2/'
 
  
def flux_anomaly_gosifgpp(varstr, M, N, latMin, latMax, lonMin, lonMax):
    
#    M=0
#    N=15     
#    varstr = 'GPP'
    
    files=glob.glob(os.path.join(dir21,'GOSIF_GPP_*.nc'))
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
          
      
    a=[]
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,N):    #baseline： 2001-2014
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)  #float(M+N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a 
    B=a-a1
    
    return A,B



 
#################################################################
#   FLUXCOM
#################################################################

dir22='/Volumes/HEWEI_T5/FLUXCOM2020/RS+METEO/member/CRUJRA_v1/NEE/2001-2015/ANN/'

    
def flux_anomaly_fluxcom_RSMETEO1(varstr, M, N, latMin, latMax, lonMin, lonMax):
       
#    M=0 
#    N=15  
#    varstr='NEE'
    
    files=glob.glob(os.path.join(dir22,'*.nc'))  
    files = sorted(files)
     
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
    
    a=[]
    
    mask=(data3==-9999)    #gC m-2 d-1
    data3[mask]=np.NaN
    
    data3=data3*365*1e-15
    

    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,N):   
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:] 
    B=a[12*M:]-a1
    
    return A,B,baseline

 
#################################################################
#   FLUXCOM
#################################################################

dir23='/Volumes/HEWEI_T5/FLUXCOM2020/RS+METEO/member/ERA5/NEE/2001-2015/ANN/'
    
def flux_anomaly_fluxcom_RSMETEO2(varstr, M, N, latMin, latMax, lonMin, lonMax):
       
#    M=0 
#    N=15  
#    varstr='NEE'
    
    files=glob.glob(os.path.join(dir23,'*.nc'))  
    files = sorted(files)
     
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
    
    a=[]
    
    mask=(data3==-9999)    #gC m-2 d-1
    data3[mask]=np.NaN
    
    data3=data3*365*1e-15
    

    data3=data3*area[2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]
      
    for i in range(0,data3.shape[0]): 
        me = np.nansum(data3[i])   
        a.append(me)
      
    baseline=np.zeros(12)
    for i in range(0,12):
        #for j in range(0,M+N):  
        for j in range(0,N):   
            baseline[i] = baseline[i] + np.nansum(data3[12*j+i])/float(N)
   
    
    #a1=np.repeat(baseline,N)
    a1=[]
    for k in range(0,N):     
       a1.extend(baseline)  
    
    a1=np.array(a1)
     
    A=a[12*M:] 
    B=a[12*M:]-a1
    
    return A,B,baseline


#################################################################
#   Byrne GOSAT+insitu+TCCON CO2 inversion: 2009-2015
#################################################################

dir24='/Volumes/HEWEI_T5/Inversions/Byrne_inversions_GOSAT_surface_TCCON_2020/Byrneetal2020/monthly/'

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

    files = glob.glob(os.path.join(dir24,dataset+'*.nc')) 
    
    
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
    
    data3 = data2[(12*M):,int(np.floor((90-latMax)/3.75-1)):int(np.ceil((90-latMin)/3.75-1)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   

    plt.imshow(data3[0])
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN   #gC/m2/month
    
    data3=data3*12.0*1e-15
 
    data3=data3*area[int(np.floor((90-latMax)/3.75-1)):int(np.ceil((90-latMin)/3.75-1)),int(np.floor((180+lonMin)/5)):int(np.ceil((180+lonMax)/5))]   

      
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

dir25='/Volumes/HEWEI_T5/Inversions/NASA-CMS-Flux/ESSD_paper_v2/monthly_1Deg/'
  
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
    
    files=glob.glob(os.path.join(dir25,'CMS-Flux*.nc'))  
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

dir26='/Volumes/HEWEI_T5/Inversions/GCASv2/'
  
def flux_anomaly_JiangGOSAT(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=6
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
        
    varstr = 'bio_monthly_opt'
    
    file=os.path.join(dir26,'posterior.fluxes.GCASv2.nc')  
    
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
            baseline[i] = baseline[i] + np.nansum(data3aa[12*j+i])/float(N) + + np.nansum(data3bb[12*j+i])/float(N)
   
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

dir27='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_WuCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=1 
#    N=6
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
      
    #dataset ='CCDAS_6yr_2deg_assim_results'
    varstr = 'rnep'
    
    file=os.path.join(dir27,dataset+'/'+'rnep.nc')  
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(file)
    lat.extend(ncf.variables['latitude'][:])
    lon.extend(ncf.variables['longitude'][:])
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

dir28='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_ScholzeCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=6
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
 
    dataset ='ccdas_sm+vod'
    varstr = 'nep'
    
    file=os.path.join(dir28,dataset+'/'+'grid_nep_2010-2015.nc')  
    
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
#   EUROCOM - CTE, RAMS, FLEXPART,Jena, LUMIA, PYVAR_CHIMERE
#################################################################

dir29='/Volumes/HEWEI_T5/Inversions/EUROCOM-ICOS/'

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
 
    file= os.path.join(dir29,dataset)
    
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
    
    #print(data2.shape)
    
    if lonMax>35:
        lonMax=35
     
    if lonMin<-14.5:
        lonMin=-14.5
        
    #data3 = data2[(12*M):,2*(73-latMax):2*(73-latMin),int(2*(lonMin+14.5)-1):int(2*(lonMax+14.5))]   
    
    data3 = data2[(12*M):,:,:] 
 
    #plt.imshow(data3[0])
    

    a=[]
    
    mask=(data3==0)
    data3[mask]=np.NaN   #kgC/m2/month
    
    data3=data3*12.0*1e-12
 
    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35 
    
#    latMin=35  
#    latMax=50  
#    lonMin=10
#    lonMax=40
 
 
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
  
############################################

def pipeline(latMin, latMax, lonMin, lonMax,flag, Year, fn ) :    
 
          
#    latMin=35
#    latMax=70
#    lonMin=-10
#    lonMax=40
    flag=1
    
  
    N=15 
     
    pydates=[]
    for i in range(1,12*N+1):   #+ (Year-2003)
        yr = 2001-1  + int(math.ceil(i/float(12)))
        if i%12 == 0 :
            mn = 12
        else:
            mn = i-int(math.floor(i/float(12)))*12
        tmp = datetime(yr, mn, 1)
        pydates.append(tmp)   
        
    #------------------------------
     
#    M=1  
#    N=15 
#      
#    iterm='GPP'
#    varstr='ensemble_mean'
##    trendy_gpp = flux_anomaly_TRENDYv6(iterm, varstr,M, N, latMin, latMax, lonMin, lonMax) 
#     
#    iterm='NEP'
#    varstr='ensemble_mean'
#    trendy_nep = flux_anomaly_TRENDYv6(iterm, varstr,M, N, latMin, latMax, lonMin, lonMax) 
        
    # individual models
    
            
    #------------------------------
     
    M=0  
    N=15 
      
    varstr='gpp'
    orchidee_gpp = flux_anomaly_orchidee(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    orchidee_ra = flux_anomaly_orchidee(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    orchidee_rh = flux_anomaly_orchidee(varstr,M, N, latMin, latMax, lonMin, lonMax) 
    #------------------------------
     
    M=0  
    N=15 
      
    varstr='gpp'
    orchideemict_gpp = flux_anomaly_orchideemict(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    orchideemict_ra = flux_anomaly_orchideemict(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    orchideemict_rh = flux_anomaly_orchideemict(varstr,M, N, latMin, latMax, lonMin, lonMax) 

    #------------------------------   
    M=0  
    N=15 
      
    varstr='gpp'
    cable_gpp = flux_anomaly_cable(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    cable_ra = flux_anomaly_cable(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    cable_rh = flux_anomaly_cable(varstr,M, N, latMin, latMax, lonMin, lonMax) 

    #------------------------------   
    M=0  
    N=15 
      
    varstr='gpp'
    dlem_gpp = flux_anomaly_dlem(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    dlem_ra = flux_anomaly_dlem(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    dlem_rh = flux_anomaly_dlem(varstr,M, N, latMin, latMax, lonMin, lonMax) 


    #------------------------------   
    M=0  
    N=15 
      
    varstr='gpp'
    isam_gpp = flux_anomaly_isam(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    isam_ra = flux_anomaly_isam(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    isam_rh = flux_anomaly_isam(varstr,M, N, latMin, latMax, lonMin, lonMax) 

    #------------------------------   
    M=0  
    N=15 
      
    varstr='gpp'
    vegas_gpp = flux_anomaly_vegas(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    vegas_ra = flux_anomaly_vegas(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    vegas_rh = flux_anomaly_vegas(varstr,M, N, latMin, latMax, lonMin, lonMax) 


    #------------------------------   
    M=0 
    N=15 
      
    varstr='gpp'
    visit_gpp = flux_anomaly_visit(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='ra'
    visit_ra = flux_anomaly_visit(varstr,M, N, latMin, latMax, lonMin, lonMax) 
     
    varstr='rh'
    visit_rh = flux_anomaly_visit(varstr,M, N, latMin, latMax, lonMin, lonMax) 


    #------------------------------
#    M=0  #1
#    N=15
#     
#    varstr='GPP'
#    sibcasa_gpp =  flux_anomaly_sibcasa(varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    varstr='Rtotal'   
#    sibcasa_reco =  flux_anomaly_sibcasa(varstr,M, N, latMin, latMax, lonMin, lonMax) 
    
       
    #------------------------------ 
 
      
#    M=0
#    N=15
#     
#    varstr='GPP'
#    beps_gpp =  flux_anomaly_beps(varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    varstr='NEP'   
#    beps_nep =  flux_anomaly_beps(varstr,M, N, latMin, latMax, lonMin, lonMax) 
    
    
     
    #------------------------------
    #    M=1  
    #    N=15 
    #      
    #    varstr='bio_flux_opt'
    #    cte2016nee =  flux_anomaly_cte2016(varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    M=1  
    N=15 
      
    varstr='bio_flux_opt'
    cte2018nee =  flux_anomaly_cte2018(varstr,M, N, latMin, latMax, lonMin, lonMax)
    varstr='fire_flux_imp'
    cte2018fire = flux_anomaly_cte2018(varstr,M, N, latMin, latMax, lonMin, lonMax)
 
    
    #------------------------------
    #    M=1  
    #    N=15 
    #      
    #    varstr='bio_flux_opt'
    #    ct2016nee =  flux_anomaly_ct2016(varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    
    M=1  
    N=15 
      
    varstr='bio_flux_opt'
    ct2019nee =  flux_anomaly_ct2019b(varstr,M, N, latMin, latMax, lonMin, lonMax)
    varstr='fire_flux_imp'
    ct2019fire =  flux_anomaly_ct2019b(varstr,M, N, latMin, latMax, lonMin, lonMax)

      
    
    
    M=22 
    N=15 
      
    varstr='flux_apos_bio'
    camsnee =  flux_anomaly_cams(varstr,M, N, latMin, latMax, lonMin, lonMax)
      
        
    
    
    #    M=5
    #    N=15 
    #      
    #    varstr='flux_apos_land'
    #    mactmnee =  flux_anomaly_mactm(varstr,M, N, latMin, latMax, lonMin, lonMax)
      
        
    
    #    M=0 
    #    N=15 
    #    varstr='co2flux_land'
    #    jena99nee =  flux_anomaly_jenas99(varstr,M, N, latMin, latMax, lonMin, lonMax)
     
    M=0 
    N=15 
    varstr='co2flux_land'
    jena99v43nee =  flux_anomaly_jenas99v43(varstr,M, N, latMin, latMax, lonMin, lonMax)
     
    
 
    M=0 
    N=15 
    varstr='co2flux_land'
    jenaNEETv43nee =  flux_anomaly_jenaNEETv43(varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    
  
    #dir19='/Volumes/HEWEI_T5/Inversions/Jena_v4.3_all/NEE-T-W-new/s99/monthly/2001-2015/'
    M=0 
    N=15 
    varstr='co2flux_land'
    jenaNEETWv43nee =  flux_anomaly_jenaNEETWv43(varstr,M, N, latMin, latMax, lonMin, lonMax)
   
 
#    M=0
#    N=15     
#    varstr = 'GPP'
#    fluxsat =  flux_anomaly_fluxsat(varstr,M, N, latMin, latMax, lonMin, lonMax)
#
#    M=0
#    N=15     
#    varstr = 'GPP'
#    gosifgpp =  flux_anomaly_gosifgpp(varstr,M, N, latMin, latMax, lonMin, lonMax)
        
    
    
        
    M=0 
    N=15  
    varstr='NEE'
    fluxcomRSMETEOnee1 =  flux_anomaly_fluxcom_RSMETEO1(varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    M=0 
    N=15  
    varstr='NEE'
    fluxcomRSMETEOnee2 =  flux_anomaly_fluxcom_RSMETEO2(varstr,M, N, latMin, latMax, lonMin, lonMax)
         




    M=0  
    N=6 
      
    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
    varstr='ensemble_mean'
    
    ByrneGOSATnee_ens =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
        
    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
    varstr='flux_CASA'
    
    ByrneGOSATnee_casa =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
      
    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
    varstr='flux_FLUXCOM'
    
    ByrneGOSATnee_fc =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
       
    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
    varstr='flux_SiB'
    
    ByrneGOSATnee_sib =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
        
    #------------------------------


    M=0 
    N=9
    
    dataset='CMS-Flux-NBE-2020.monthly.grid'           
    varstr = 'post-NBE' 
    LiuGOSATnee =  flux_anomaly_LiuGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)

    #------------------------------

    M=0 
    N=6
     
    #dataset='CMS-Flux-NBE-2020.monthly.grid'        
    varstr = 'bio_monthly_opt'
    JiangGOSATnee =  flux_anomaly_JiangGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    
    #------------------------------
#
#    M=1 
#    N=6
#     
#    dataset='CCDAS_6yr_2deg_assim_results'        
#    varstr = 'rnep'
#    WuCCDASnee1 =  flux_anomaly_WuCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#  
#    
#    #------------------------------
#    
#    M=1 
#    N=6
#     
#    dataset='ccdas_sm+fapar'        
#    varstr = 'rnep'
#    WuCCDASnee2 =  flux_anomaly_WuCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
 
        
    #------------------------------
    
    M=0 
    N=6
     
    dataset='ccdas_sm+vod'        
    varstr = 'nep'
    ScholzeCCDASnee =  flux_anomaly_ScholzeCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
   
    
    
    #------------------------------
        
#    M=4 
#    N=6
#    
##    latMin=33    #35
##    latMax=73    #70
##    lonMin=-14.5 #-10
##    lonMax=35    #40
#    
#    varstr='bio_posterior'
# 
#    dataset = 'EUROCOM_JenaCarboScopeRegional_2006-2015.nc'
#    JenaRegionalnee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    dataset = 'EUROCOM_LUMIA_2006-2015.nc'
#    LUMIAnee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
     
    
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
    
         
    
    #######################################
    # raw fluxes
    #######################################
  
    cte2018neeann = []
    ct2019neeann = []
    camsneeann=[]
    jena99v43neeann=[]
    jenaNEETv43neeann=[]
    jenaNEETWv43neeann=[]
    
    
    for i in range(0, 15): 
       cte2018neeann.append(np.nanmean( cte2018nee[flag][i*12:(i+1)*12] + cte2018fire[flag][i*12:(i+1)*12]) )
       ct2019neeann.append(np.nanmean( ct2019nee[flag][i*12:(i+1)*12] + ct2019fire[flag][i*12:(i+1)*12]) )
       camsneeann.append(np.nanmean( camsnee[flag][i*12:(i+1)*12] ) )
       jena99v43neeann.append(np.nanmean( jena99v43nee[flag][i*12:(i+1)*12] ) )
       jenaNEETv43neeann.append(np.nanmean( jenaNEETv43nee[flag][i*12:(i+1)*12] ) )
       jenaNEETWv43neeann.append(np.nanmean( jenaNEETWv43nee[flag][i*12:(i+1)*12] ) )
      
 
    cte2018neeann=array(cte2018neeann)
    ct2019neeann=array(ct2019neeann)
    camsneeann=array(camsneeann)

    jena99v43neeann=array(jena99v43neeann)   
    jenaNEETv43neeann=array(jenaNEETv43neeann)   
    jenaNEETWv43neeann=array(jenaNEETWv43neeann) 
     
    
    
        
#    trendy_gpp2 = array(trendy_gpp[flag])
#    sibcasa_gpp2 = array(sibcasa_gpp[flag])
#    beps_gpp2 = array(beps_gpp[flag])
#       
#    sibcasa_reco2 = array(sibcasa_reco[flag])  
#    
#    #beps_reco = array(beps_gpp[flag])-array(beps_nep[flag])
#       
#    trendy_nee2 =  -trendy_nep[flag] 
#    sibcasa_nee2 =  sibcasa_reco2-sibcasa_gpp2 
#    beps_nee2 = -array(beps_nep[flag])
#    
#    beps_reco2 = beps_gpp1+beps_nee1 
    
#    trendy_gpp_ann=[]
#    sibcasa_gpp_ann=[]
#    beps_gpp_ann=[]
#    
##    trendy_reco_ann=[]
#    sibcasa_reco_ann=[]
#    beps_reco_ann=[]
#    
#    trendy_nee_ann=[]
#    sibcasa_nee_ann=[]
#    beps_nee_ann=[]
    
#    N=15
#    for i in range(0, N):  
#       trendy_gpp_ann.append(np.nanmean( trendy_gpp2[i*12:(i+1)*12] ) )
#       sibcasa_gpp_ann.append(np.nanmean( sibcasa_gpp2[i*12:(i+1)*12] ) )
#       beps_gpp_ann.append(np.nanmean( beps_gpp2[i*12:(i+1)*12] ) )
#       
#       #trendy_reco_ann.append(np.nanmean( trendy_reco2[i*12:(i+1)*12] ) )
#       sibcasa_reco_ann.append(np.nanmean( sibcasa_reco2[i*12:(i+1)*12] ) )
#       beps_reco_ann.append(np.nanmean( beps_reco2[i*12:(i+1)*12] ) )
#       
#       trendy_nee_ann.append(np.nanmean( trendy_nee2[i*12:(i+1)*12] ) )
#       sibcasa_nee_ann.append(np.nanmean( sibcasa_nee2[i*12:(i+1)*12] ) )
#       beps_nee_ann.append(np.nanmean( beps_nee2[i*12:(i+1)*12] ) )
     
#    trendy_gpp_ann=array(trendy_gpp_ann)
#    sibcasa_gpp_ann=array(sibcasa_gpp_ann)
#    beps_gpp_ann=array(beps_gpp_ann)
#     
##    trendy_reco_ann=array(trendy_reco_ann)
#    sibcasa_reco_ann=array(sibcasa_reco_ann)
#    beps_reco_ann=array(beps_reco_ann)
#    
#    trendy_nee_ann=array(trendy_nee_ann)
#    sibcasa_nee_ann=array(sibcasa_nee_ann)
#    beps_nee_ann=array(beps_nee_ann)
 
    #detrend
#    trendy_gpp_ann=signal.detrend(trendy_gpp_ann)
#    sibcasa_gpp_ann=signal.detrend(sibcasa_gpp_ann)
#    beps_gpp_ann=signal.detrend(beps_gpp_ann)
#     
##    trendy_reco_ann=signal.detrend(trendy_reco_ann)
#    sibcasa_reco_ann=signal.detrend(sibcasa_reco_ann)
##    beps_reco_ann=signal.detrend(beps_reco_ann)
#    
#    trendy_nee_ann=signal.detrend(trendy_nee_ann)
#    sibcasa_nee_ann=signal.detrend(sibcasa_nee_ann)
#    beps_nee_ann=signal.detrend(beps_nee_ann)
     
    N=15   
    for i in range(0, N): 
        cte2018neeann = signal.detrend(cte2018neeann)
        ct2019neeann = signal.detrend(ct2019neeann)
        camsneeann = signal.detrend(camsneeann)
        jena99v43neeann = signal.detrend(jena99v43neeann)
        jenaNEETv43neeann =signal.detrend(jenaNEETv43neeann) 
        jenaNEETWv43neeann =signal.detrend(jenaNEETWv43neeann) 
    
    #IX = Year - 2001   
    #pydates1=pydates[12*IX:12*(IX+1)] 
     
    
    # consider detrend, similar to the yearly analysis
#    fluxsat_gpp = array(fluxsat[flag])
#    gosif_gpp = array(gosifgpp[flag])
   
#    fluxsat_gpp_ann=[]
#    gosif_gpp_ann=[]
#    
#    N=15
#    for i in range(0, N):  
#       fluxsat_gpp_ann.append(np.nanmean( fluxsat_gpp[i*12:(i+1)*12] ) )
#       gosif_gpp_ann.append(np.nanmean( gosif_gpp[i*12:(i+1)*12] ) )
#   
#    fluxsat_gpp_ann=signal.detrend(fluxsat_gpp_ann)
#    gosif_gpp_ann=signal.detrend(gosif_gpp_ann)
    
    
  
    
    eurocom1neeann=[]
    eurocom2neeann=[]
    eurocom3neeann=[]
    eurocom4neeann=[]
    eurocom5neeann=[]
    eurocom6neeann=[]
    eurocom7neeann=[]
    #JenaRegionalneeann = []
    #LUMIAneeann = []
    ByrneGOSATneeann=[]
    ByrneGOSATneeann=[]
    ByrneGOSATneeann=[]
    ByrneGOSATneeann=[]
    JiangGOSATneeann=[]
    LiuGOSATneeann=[]
    ScholzeCCDASneeann=[]
    
 
    for i in range(0, 6):
        eurocom1neeann.append(np.nanmean( eurocom1nee[flag][i*12:(i+1)*12] + eurocom1fire[flag][i*12:(i+1)*12] ) )
        eurocom2neeann.append(np.nanmean( eurocom2nee[flag][i*12:(i+1)*12] ) )
        eurocom3neeann.append(np.nanmean( eurocom3nee[flag][i*12:(i+1)*12] ) )
        eurocom4neeann.append(np.nanmean( eurocom4nee[flag][i*12:(i+1)*12] ) )
        eurocom5neeann.append(np.nanmean( eurocom5nee[flag][i*12:(i+1)*12]  + eurocom5fire[flag][i*12:(i+1)*12] ) )
        eurocom6neeann.append(np.nanmean( eurocom6nee[flag][i*12:(i+1)*12] ) )
        #JenaRegionalneeann.append(np.nanmean( JenaRegionalnee[flag][i*12:(i+1)*12]  ) ) 
        #LUMIAneeann.append(np.nanmean( LUMIAnee[flag][i*12:(i+1)*12]  ) )
        ByrneGOSATneeann.append(np.nanmean( ByrneGOSATnee_ens[flag][i*12:(i+1)*12]  ) )
        JiangGOSATneeann.append(np.nanmean( JiangGOSATnee[flag][i*12:(i+1)*12]  ) ) 
        LiuGOSATneeann.append(np.nanmean( LiuGOSATnee[flag][i*12:(i+1)*12]  ) ) 
        ScholzeCCDASneeann.append(np.nanmean( ScholzeCCDASnee[flag][i*12:(i+1)*12]  ) ) 
       
#    for i in range(0, 7):  
    for i in range(0, 5):
        eurocom7neeann.append(np.nanmean( eurocom7nee[flag][i*12:(i+1)*12] ) )
    
    
    eurocom1neeann=np.array(eurocom1neeann)
    eurocom2neeann=np.array(eurocom2neeann)
    eurocom3neeann=np.array(eurocom3neeann)
    eurocom4neeann=np.array(eurocom4neeann)
    eurocom5neeann=np.array(eurocom5neeann)
    eurocom6neeann=np.array(eurocom6neeann)
    eurocom7neeann=np.array(eurocom7neeann)
    #JenaRegionalneeann=np.array(JenaRegionalneeann)     
    #LUMIAneeann=np.array(LUMIAneeann) 
    ByrneGOSATneeann=np.array(ByrneGOSATneeann) 
    JiangGOSATneeann=np.array(JiangGOSATneeann)
    LiuGOSATneeann=np.array(LiuGOSATneeann)
    ScholzeCCDASneeann=np.array(ScholzeCCDASneeann)  

    eurocom1neeann=signal.detrend(eurocom1neeann)
    eurocom2neeann=signal.detrend(eurocom2neeann)
    eurocom3neeann=signal.detrend(eurocom3neeann)
    eurocom4neeann=signal.detrend(eurocom4neeann)
    eurocom5neeann=signal.detrend(eurocom5neeann)
    eurocom6neeann=signal.detrend(eurocom6neeann)
    eurocom7neeann=signal.detrend(eurocom7neeann)
    #JenaRegionalneeann=signal.detrend(JenaRegionalneeann)
    #LUMIAneeann=signal.detrend(LUMIAneeann)
    ByrneGOSATneeann=signal.detrend(ByrneGOSATneeann)
    JiangGOSATneeann=signal.detrend(JiangGOSATneeann)
    LiuGOSATneeann=signal.detrend(LiuGOSATneeann)
    ScholzeCCDASneeann=signal.detrend(ScholzeCCDASneeann)
  
    
    
     
       
    eurocom7neeannnew = []
    eurocom7neeannnew.extend([0])
    eurocom7neeannnew.extend(eurocom7neeann)
    eurocom7neeannnew = np.array(eurocom7neeannnew)
    eurocom7neeannnew[0] = np.nan 
  
    eurocom7neenew = []
    eurocom7neenew.extend([0,0,0,0,0,0,0,0,0,0,0,0])
    eurocom7neenew.extend(eurocom7nee[flag])
    eurocom7neenew = np.array(eurocom7neenew)
    eurocom7neenew[:12] = np.repeat(np.nan,12)
    
    
    
    
    IX = Year - 2001
#    trendy_gpp1 = np.array(trendy_gpp[flag])[12*IX:12*(IX+1)] 
#    sibcasa_gpp1 = np.array(sibcasa_gpp[flag])[12*IX:12*(IX+1)]  
#    beps_gpp1 = np.array(beps_gpp[flag])[12*IX:12*(IX+1)]  
       
     
#    sibcasa_reco1 = np.array(sibcasa_reco[flag])[12*IX:12*(IX+1)]   
#    beps_reco1 = np.array(beps_gpp[flag])[12*IX:12*(IX+1)] - np.array(beps_nep[flag])[12*IX:12*(IX+1)] 
          
    
#    trendy_nee1 = -np.array(trendy_nep[flag])[12*IX:12*(IX+1)] 
#    sibcasa_nee1 =  sibcasa_reco1-sibcasa_gpp1 
#    beps_nee1 = -np.array(beps_nep[flag])[12*IX:12*(IX+1)] 
#    beps_reco = beps_gpp + beps_nee
        
    orchidee_gpp = np.array(orchidee_gpp[flag])
    orchideemict_gpp = np.array(orchideemict_gpp[flag])
    
    cable_gpp = np.array(cable_gpp[flag]) 
    dlem_gpp = np.array(dlem_gpp[flag]) 
    isam_gpp = np.array(isam_gpp[flag]) 
    #lpxbern_gpp = np.array(lpxbern_gpp[flag]) 
    vegas_gpp = np.array(vegas_gpp[flag]) 
    visit_gpp = np.array(visit_gpp[flag]) 
#    sibcasa_gpp = np.array(sibcasa_gpp[flag])
#    beps_gpp = np.array(beps_gpp[flag])
   
     
    orchidee_reco = np.array(orchidee_ra[flag]) + np.array(orchidee_rh[flag])  
    orchideemict_reco = np.array(orchideemict_ra[flag]) + np.array(orchideemict_rh[flag]) 
    cable_reco = np.array(cable_ra[flag]) + np.array(cable_rh[flag]) 
    dlem_reco = np.array(dlem_ra[flag]) + np.array(dlem_rh[flag]) 
    isam_reco = np.array(isam_ra[flag]) + np.array(isam_rh[flag])
    #lpxbern_reco = np.array(lpxbern_ra[flag])[12*4:12*6] + np.array(lpxbern_rh[flag])
    vegas_reco = np.array(vegas_ra[flag]) + np.array(vegas_rh[flag]) 
    visit_reco = np.array(visit_ra[flag]) + np.array(visit_rh[flag]) 
#    sibcasa_reco = np.array(sibcasa_reco[flag])      
#    beps_reco = np.array(beps_gpp[flag])-np.array(beps_nep[flag])
   
    orchidee_nee =  orchidee_reco-orchidee_gpp  
    orchideemict_nee =  orchideemict_reco-orchideemict_gpp  
    cable_nee =  cable_reco-cable_gpp  
    dlem_nee =  dlem_reco-dlem_gpp  
    isam_nee =  isam_reco-isam_gpp  
    #lpxbern_nee =  lpxbern_reco-lpxbern_gpp  
    vegas_nee =  vegas_reco-vegas_gpp  
    visit_nee =  visit_reco-visit_gpp  
     
    
    orchidee_nee1 = orchidee_nee[12*IX:12*(IX+1)]
    orchideemict_nee1 = orchideemict_nee[12*IX:12*(IX+1)]
    cable_nee1 = cable_nee[12*IX:12*(IX+1)]
    dlem_nee1 = dlem_nee[12*IX:12*(IX+1)]
    isam_nee1 = isam_nee[12*IX:12*(IX+1)]
    vegas_nee1 = vegas_nee[12*IX:12*(IX+1)]
    visit_nee1 = visit_nee[12*IX:12*(IX+1)]
    
    
    
    cte2018nee1=cte2018nee[flag][12*IX:12*(IX+1)] + cte2018fire[flag][12*IX:12*(IX+1)]  
    ct2019nee1=ct2019nee[flag][12*IX:12*(IX+1)] + ct2019fire[flag][12*IX:12*(IX+1)]  
    camsnee1=camsnee[flag][12*IX:12*(IX+1)] 
    jena99v43nee1=jena99v43nee[flag][12*IX:12*(IX+1)] 
    jenaNEETv43nee1=jenaNEETv43nee[flag][12*IX:12*(IX+1)] 
    jenaNEETWv43nee1=jenaNEETWv43nee[flag][12*IX:12*(IX+1)] 
    
    #fluxsat=fluxsat[flag][12*IX:12*(IX+1)] 
    #gosifgpp=gosifgpp[flag][12*IX:12*(IX+1)] 
     
    fluxcomRSMETEOneeC = fluxcomRSMETEOnee1[flag][12*IX:12*(IX+1)] 
    fluxcomRSMETEOneeE = fluxcomRSMETEOnee2[flag][12*IX:12*(IX+1)]    
 
    
    IX = Year - 2001   
    pydates1=pydates[12*IX:12*(IX+1)] 
     
  
    
    #####################
    #  plot out
    #####################
    fig = plt.figure(figsize=(4.5, 7.5))  #12,5*3-2
    
    fontsize=24
    mlw=1.2 
    
    
    #######################
    
    ax1=fig.add_subplot(5,1,1)
        
    from matplotlib.ticker import FormatStrFormatter
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
 
    
    mean = np.mean([-cte2018nee1/12.0*1e3, -ct2019nee1/12.0*1e3, -camsnee1/12.0*1e3, -jena99v43nee1/12.0*1e3, -jenaNEETv43nee1/12.0*1e3, -jenaNEETWv43nee1/12.0*1e3], axis=0)
    std = np.std([-cte2018nee1/12.0*1e3, -ct2019nee1/12.0*1e3, -camsnee1/12.0*1e3, -jena99v43nee1/12.0*1e3, -jenaNEETv43nee1/12.0*1e3, -jenaNEETWv43nee1/12.0*1e3], axis=0)
    ax1.plot(pydates1,mean, linewidth=mlw*1.5,color='black',linestyle='-', marker='o', markersize = 3) 
#    ax1.fill_between(pydates1, mean-std, mean+std,facecolor = "lightgray")
    
    
    ax1.plot(pydates1,-cte2018nee1/12.0*1e3,linewidth=mlw,color='lightgreen',linestyle='-', marker='o', markersize = 1.2) 
    ax1.plot(pydates1,-ct2019nee1/12.0*1e3,linewidth=mlw,color='lightgreen',linestyle='-', marker='o', markersize = 1.2) 
    ax1.plot(pydates1,-camsnee1/12.0*1e3,linewidth=mlw,color='lightgreen',linestyle='-', marker='o', markersize = 1.2) 
    ax1.plot(pydates1,-jena99v43nee1/12.0*1e3,linewidth=mlw,color='lightgreen',linestyle='-', marker='o', markersize = 1.2) 
    ax1.plot(pydates1,-jenaNEETv43nee1/12.0*1e3,linewidth=mlw,color='blue',linestyle='-', marker='o', markersize = 1.2) 
    ax1.plot(pydates1,-jenaNEETWv43nee1/12.0*1e3,linewidth=mlw,color='deeppink',linestyle='-', marker='o', markersize = 1.2) 
                

    ax1.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    #ax1.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.5)
    
     
    ax1.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax1.locator_params(nbins=8)
    
    if Year ==2012:
         #ax1.legend(['CTE2018','CT2019', 'CAMS_v18r2','Jena99_v43','JenaNEET_v43','JenaNEETW_v43'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.3)
         ax1.legend(['global AIMs ensemble', 'individual',' ',' ',' ','JenaNEET_v4.3','JenaNEETW_v4.3'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.33)
  
    
    #MIN=abs(np.nanmin(gppmean2))
    #MAX=abs(np.nanmax(gppmean2))  
    
    #ax1.set_ylim([MIN*1.1,MAX*1.1])
       # ax2.set_ylim([MIN*1.1,MAX*1.1])
    
  #  MIN=min( np.nanmin(-cte2018nee1/12.0*1e3),np.nanmin(-ct2019nee1/12.0*1e3), np.nanmin(-camsnee1/12.0*1e3),np.nanmin(-jena99v43nee1/12.0*1e3), np.nanmin(-jenaNEETv43nee1/12.0*1e3), np.nanmin(-jenaNEETWv43nee1/12.0*1e3) ) 
  #  MAX=max( np.nanmax(-cte2018nee1/12.0*1e3),np.nanmax(-ct2019nee1/12.0*1e3), np.nanmax(-camsnee1/12.0*1e3),np.nanmax(-jena99v43nee1/12.0*1e3), np.nanmax(-jenaNEETv43nee1/12.0*1e3), np.nanmax(-jenaNEETWv43nee1/12.0*1e3) )     
    #DEV = max(abs(MIN),abs(MAX)) 
    #ax1.set_ylim([MIN*1.1,MAX*1.1])
    ax1.set_ylim([-110.0,110.0])
     
    yy = np.zeros(len(pydates1))
    
    if flag == 1:
    #  ax1.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
      #ax2.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
      ax1.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
    
    ax1.tick_params(axis='both', which='major', labelsize=fontsize*0.5)      
    
    #ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    #ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    
    ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax1.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax1.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax1.get_yticklabels()]
    
    plt.setp(ax1.get_xticklabels(), visible=False) 
  

    
    #################
    #ax2=fig.add_subplot(3,1,2)
    ax2=fig.add_subplot(5,1,2)
    
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
 
    IX = Year - 2010   
    pydates2=pydates[12*IX:12*(IX+1)] 
     
    
    mean = np.nanmean([-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,  \
                       -np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-eurocom7neenew[12*IX:12*(IX+1)]/12.0*1e3  ], axis=0)  
    
    std = np.nanstd([-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,  \
                       -np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-eurocom7neenew[12*IX:12*(IX+1)]/12.0*1e3  ], axis=0)  
    
     
    
    ax2.plot(pydates2,mean,linewidth=mlw*1.5,color='k',linestyle='-', marker='o', markersize = 3)     
#    ax2.fill_between(pydates2, mean-std, mean+std,facecolor = "lightgray")
   
#    ax2.plot(pydates1,-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='r',linestyle='-', marker='o', markersize = 4)  
#    ax2.plot(pydates1,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='limegreen',linestyle='-', marker='o', markersize = 4)         
#    ax2.plot(pydates1,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='orange',linestyle='-', marker='o', markersize = 4) 
#    ax2.plot(pydates1,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='royalblue',linestyle='-', marker='^', markersize = 4) 
#    ax2.plot(pydates1,-np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='gold',linestyle='-', marker='o', markersize = 4) 
#    ax2.plot(pydates1,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='magenta',linestyle='-', marker='s', markersize = 4) 
#    ax2.plot(pydates1,-np.array(eurocom7neenew)[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='cyan',linestyle='-', marker='^', markersize = 4) 
      
    
    ax2.plot(pydates2,-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5)  
    ax2.plot(pydates2,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5)         
    ax2.plot(pydates2,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax2.plot(pydates2,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax2.plot(pydates2,-np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax2.plot(pydates2,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5) 
    ax2.plot(pydates2,-np.array(eurocom7neenew)[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 1.5)     
       
    
    ax2.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    #ax2.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.5)
    
     
    ax2.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax2.locator_params(nbins=8)
    
#    if Year ==2015:
#        ax2.legend(['Byrne2020_ensemble','Liu2020','Jiang2020','CCDAS_SM','CCDAS_SM+FAPAR','CCDAS_SM+VOD'], ncol = 3,  frameon=True, loc = 3, fontsize = fontsize*0.32)
 
    if Year ==2012:
        #ax2.legend(['CTE','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Mean'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.32)
        ax2.legend(['EUROCOM ensemble','individual'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.35)
  
    ax2.set_ylim([-110,110])   #[-90,90]
     
    yy = np.zeros(len(pydates2))
    
    if flag == 1:
        ax2.plot(pydates2,yy, color='gray',linestyle='--',linewidth=1)
    
    ax2.tick_params(axis='both', which='major', labelsize=fontsize*0.5)    
    
    #ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    #ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    
    ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax2.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax2.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax2.get_yticklabels()]
    
    plt.setp(ax2.get_xticklabels(), visible=False) 
        
    
    #################        
    #ax3=fig.add_subplot(3,1,3)
    ax3=fig.add_subplot(5,1,3)
    
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
#    mean = np.nanmean([-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,  \
#                       -np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-eurocom7neenew[12*IX:12*(IX+1)]/12.0*1e3  ], axis=0)  
#    
#    std = np.nanstd([-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,  \
#                       -np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-eurocom7neenew[12*IX:12*(IX+1)]/12.0*1e3  ], axis=0)  
    
     
#    ax1.plot(pydates1,mean,linewidth=mlw,color='b',linestyle='-', marker='^', markersize = 4)     
#    ax1.fill_between(pydates1, mean-std, mean+std,facecolor = "skyblue")
   
  
    
    mean = np.nanmean([-np.array(ByrneGOSATnee_casa[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(ByrneGOSATnee_fc[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(ByrneGOSATnee_sib[flag])[12*IX:12*(IX+1)]/12.0*1e3 ], axis=0)    
    std = np.nanstd([-np.array(ByrneGOSATnee_casa[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(ByrneGOSATnee_fc[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(ByrneGOSATnee_sib[flag])[12*IX:12*(IX+1)]/12.0*1e3 ], axis=0)  
    
    ax3.plot(pydates2,mean,linewidth=mlw,color='r',linestyle='-', marker='^', markersize = 3)   
#    ax1.fill_between(pydates1, mean-std, mean+std,facecolor = "salmon")
   
    
#    ax2.plot(pydates1,-np.array(JenaRegionalnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='royalblue',linestyle='-', marker='o', markersize = 4)  
#    ax2.plot(pydates1,-np.array(LUMIAnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='gold',linestyle='-', marker='o', markersize = 4)         
#    ax1.plot(pydates1,-np.array(ByrneGOSATnee_ens[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='r',linestyle='-', marker='o', markersize = 4) 
    ax3.plot(pydates2,-np.array(JiangGOSATnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='g',linestyle='-', marker='o', markersize = 3) 
    ax3.plot(pydates2,-np.array(LiuGOSATnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='royalblue',linestyle='-', marker='s', markersize = 3) 
  #    ax2.plot(pydates1,WuCCDASnee1,linewidth=mlw,color='pink',linestyle='-', marker='o', markersize = 4) 
#    ax2.plot(pydates1,WuCCDASnee2,linewidth=mlw,color='royalblue',linestyle='-', marker='s', markersize = 4) 
    ax3.plot(pydates2,np.array(ScholzeCCDASnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='tan',linestyle='-', marker='o', markersize = 3) 
       
    ax3.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    ax3.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.55)
    
     
    ax3.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax3.locator_params(nbins=8)
    
#    if Year ==2015:
#        ax2.legend(['Byrne2020_ensemble','Liu2020','Jiang2020','CCDAS_SM','CCDAS_SM+FAPAR','CCDAS_SM+VOD'], ncol = 3,  frameon=True, loc = 3, fontsize = fontsize*0.32)
 
    if Year ==2012:
        ax3.legend(['Byrne2020','GCASv2','CMS-Flux2020','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.35)
  
    ax3.set_ylim([-110,110])
     
    yy = np.zeros(len(pydates2))
    
    if flag == 1:
        ax3.plot(pydates2,yy, color='gray',linestyle='--',linewidth=1)
    
    ax3.tick_params(axis='both', which='major', labelsize=fontsize*0.5)    
    
    #ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    #ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    
    ax3.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax3.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax3.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax3.get_yticklabels()]
    
    plt.setp(ax3.get_xticklabels(), visible=False) 
    
 
    
    
   #################
    ax4=fig.add_subplot(5,1,4)
     
    nepmean2 = np.mean([-orchidee_nee1,-orchideemict_nee1, -cable_nee1,-dlem_nee1,-isam_nee1,-vegas_nee1,-visit_nee1], axis=0)
    nepstd2 = np.std([-orchidee_nee1,-orchideemict_nee1, -cable_nee1,-dlem_nee1,-isam_nee1,-vegas_nee1,-visit_nee1], axis=0)
    ax4.plot(pydates1,nepmean2/12.0*1e3, linewidth=mlw*1.5,color='black',linestyle='-', marker='o', markersize = 3) 
#    ax4.fill_between(pydates1, nepmean2/12.0*1e3-nepstd2/12.0*1e3, nepmean2/12.0*1e3+nepstd2/12.0*1e3,facecolor = "lightgray")
     
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax4.plot(pydates1,-orchidee_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2) 
    ax4.plot(pydates1,-orchideemict_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2) 
    ax4.plot(pydates1,-cable_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2)     
    ax4.plot(pydates1,-dlem_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2) 
    ax4.plot(pydates1,-isam_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2)   
    ax4.plot(pydates1,-vegas_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2)  
    ax4.plot(pydates1,-visit_nee1/12.0*1e3,linewidth=mlw,color='khaki',linestyle='-', marker='o', markersize = 1.2) 
#    ax2.plot(pydates1,sibcasa_nee1,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 3) 
#    ax2.plot(pydates1,beps_nee1,linewidth=mlw,color='skyblue',linestyle='-', marker='.', markersize = 3)
     
  
    ax4.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    #ax4.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.5)
    
     
    ax4.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax4.locator_params(nbins=8)
    
    if Year ==2012:
         ax4.legend(['TRENDY_v6 ensemble','individual'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.35)
       
        
    ax4.set_ylim([-110.0,110.0])
   
    yy = np.zeros(len(pydates1))
    
    if flag == 1:
         ax4.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
    
    ax4.tick_params(axis='both', which='major', labelsize=fontsize*0.5)      
        
    ax4.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax4.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax4.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax4.get_yticklabels()]
    
    plt.setp(ax4.get_xticklabels(), visible=False) 
    
    
    #################       
    ax5=fig.add_subplot(5,1,5)
     
    
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax5.plot(pydates1,-fluxcomRSMETEOneeC/12.0*1e3,linewidth=mlw,color='indianred',linestyle='-', marker='o', markersize = 1.2) 
    ax5.plot(pydates1,-fluxcomRSMETEOneeE/12.0*1e3,linewidth=mlw,color='seagreen',linestyle='-', marker='o', markersize = 1.2)     
  
    
    ax5.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    #ax5.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.5)
    
     
    ax5.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax5.locator_params(nbins=4)
    
    if Year ==2012:
         ax5.legend(['FLUXCOM_CRUJRA','FLUXCOM_ERA5'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.35)
    
        
    ax5.set_ylim([-20.0,20.0])
    
    yy = np.zeros(len(pydates1))
    
    if flag == 1:
         ax5.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
    
    ax5.tick_params(axis='both', which='major', labelsize=fontsize*0.5)      
        
    ax5.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax5.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax5.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax5.get_yticklabels()]
 
 
    
        
    fig = plt.gcf()
    #fig.set_size_inches(4.5, 7.5)  #(6,10)
    fig.set_size_inches(4.5, 7.5+0.75) 
     
    
    fig.tight_layout() 
    #plt.subplots_adjust(hspace = 0.15)  #0.10
    plt.subplots_adjust(hspace = 0.05)
          
       
    fn='All_Flux_timeseries_Europe_20221008_'+str(Year)
    outp=os.path.join(dir1,fn+'.png')
    print(outp)
    plt.savefig(outp,dpi=300)
    
    outp=os.path.join(dir1,str(Year)+fn+'.eps')
    print(outp)
    plt.savefig(outp,dpi=600,format='eps')    
    
    
    data=[]
    
 
    data.append(-cte2018neeann[IX])
    data.append(-ct2019neeann[IX])
    data.append(-camsneeann[IX])
    data.append(-jena99v43neeann[IX])
    data.append(-jenaNEETv43neeann[IX])
    data.append(-jenaNEETWv43neeann[IX])
     
    #data.append(fluxsat_gpp[IX])
    #data.append(gosif_gpp[IX])      
   
    data.append(-eurocom1neeann[IX])
    data.append(-eurocom2neeann[IX])
    data.append(-eurocom3neeann[IX])
    data.append(-eurocom4neeann[IX])
    data.append(-eurocom5neeann[IX])
    data.append(-eurocom6neeann[IX])
    data.append(-ByrneGOSATneeann[IX])
    data.append(-JiangGOSATneeann[IX])
    data.append(-LiuGOSATneeann[IX])
    data.append(-ScholzeCCDASneeann[IX])
    
  
    data=np.array(data)
    
    print("CMS-Flux2020 for year %s"%Year)
    print(-LiuGOSATneeann[IX])
    
    
    
    return data    
    
       
import timeit
start = timeit.timeit()

flag=1
 
#Year = 2003 
#data = pipeline(35, 70, -10, 40, flag, Year,'All_Flux_timeseries_Europe_20200428')
#
#Year = 2006
#data = pipeline(35, 70, -10, 40, flag, Year, 'All_Flux_timeseries_Europe_20200428')

Year = 2012
#data = pipeline(35, 70, -10, 40, flag, Year, 'All_Flux_timeseries_Europe_20200428')
data11 = pipeline(35, 50, 10, 40, flag, Year, 'All_Flux_timeseries_Europe_20221008')


print(data11)
 
Year = 2015
#data = pipeline(35, 70, -10, 40, flag, Year, 'All_Flux_timeseries_Europe_20200428')
data22 = pipeline(45, 65, 0, 35, flag, Year, 'All_Flux_timeseries_Europe_20221008')



end = timeit.timeit()
print(end - start)
     





##data11
#[-0.09097514, -0.09631084, -0.05580726,  0.01444179, -0.06062256,
#       -0.00916703, -0.26348409, -0.14570927, -0.269025  , -0.09278803,
#       -0.13071877,  0.04588336,  0.36026151, -0.03491898, -0.15559823,
#        0.36837873,  0.17191247,  0.46642819,  0.12501272,  0.08572502,
#       -0.03019259,  0.36026151,  0.13276507,  0.30414449,  0.10489465,
#        0.0262032 ,  0.19740319,  0.03222469, -0.04499374,  0.01569078,
#        0.72052303,  0.09784609,         nan, -0.02677872, -0.01876761]
#
##data22
#[-0.00166041,  0.01481946,  0.16060443,  0.07988457, -0.01066147,
#       -0.01550863, -0.17601254, -0.12713907, -0.03170502, -0.07483585,
#       -0.00227846,  0.0263374 ,  0.19384429,  0.01836523, -0.09185225,
#        0.24474932,  0.1729141 , -0.00261306,  0.05826731,  0.07526981,
#       -0.01860145,  0.19384429,  0.02888189,  0.29529432,  0.06873678,
#        0.04577503, -0.03431808, -0.01656854,  0.07299135,  0.00773596,
#        0.38768858,  0.04724712,         nan,  0.03340919,  0.01939947]


'''
 
############################################ 
    
ind = np.arange(2)   # the x locations for the groups
width = 0.125          # the width of the bars

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 6
fig_size[1] = 4.5
plt.rcParams["figure.figsize"] = fig_size

#colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["royalblue", "gold", "r", "g", "tan"]

#colors = ["orange","r","g", "royalblue", "magenta","yellow"]
colors = ["royalblue", "gold", "r", "g", "tan", "magenta"]

fig = plt.figure(figsize=(8,4))
line =1
row=1

fz=18  #20
######################################################
# subplot R1
#from brokenaxes import brokenaxes


ax=fig.add_subplot(line,row,1)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 


data=np.zeros(shape=(2,6))
data[0,:] = data11[:6]
data[1,:] = data22[:6]

for i in np.arange(0,6):
    clr= colors[i]
    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
    rects.append(ax.bar(ind + (i+1)*width, data[:,i]*1e3, width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
    
    #clr= colors[5]
    
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
 
#ax.grid()
ax.set_ylabel('NEP anomaly [TgC/yr]',  fontsize=fz)
 
ax.set_xticks(ind+3*width)
ax.set_ylim([-170,170])
#ax.set_ylim([-210,210])
ax.locator_params(axis='y',nticks=2)
 
ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')

lgnames.append("CTE2018")
lgnames.append("CT2019")
lgnames.append("CAMS_v18r2")
lgnames.append("JenaStandard_v43")
lgnames.append("JenaNEET_v43")
lgnames.append("JenaNEETW_v43") 

#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax.legend( rects, lgnames, loc = 2, ncol=2, frameon=False, fontsize = fz*0.8)
plt.tick_params(labelsize=fz*0.9)  

#--------------------------------------
#fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
       
fig.tight_layout() 

fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_AIMs'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)

     
     



############################################ 
    
ind = np.arange(2)   # the x locations for the groups
width = 0.75    #0.10      # the width of the bars

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 6
fig_size[1] = 4.5
plt.rcParams["figure.figsize"] = fig_size

colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["royalblue", "gold", "r", "g", "tan", "chocolate" ,"crimson", "IndianRed", "forestgreen"]
          
fig = plt.figure(figsize=(8,4))
line =1
row=1

fz=18  #20
######################################################
# subplot R1
#from brokenaxes import brokenaxes


ax=fig.add_subplot(line,row,1)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 


data=np.zeros(shape=(2,1))
err=np.zeros(shape=(2,1))
#data[0,:] = data11[5:]
#data[1,:] = data22[5:]

data[0,:] = np.nanmedian(data11[6:15])
data[1,:] = np.nanmedian(data22[6:15])


err[0,:] = np.nanstd(data11[6:15])
err[1,:] = np.nanstd(data22[6:15])


#for i in np.arange(0,9):
#    clr= colors[i]
#    rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )

clr= colors[0]
#rects.append(ax.errorbar(ind + (0+1)*width, data[:,0]/12.0*1e3, width, color=clr, yerr= err[:,0]/12.0*1e3,  ecolor=clr) )
#rects.append(ax.bar(ind + (2+1)*width, data[:,0]/12.0*1e3, width, color=clr) )
rects.append(ax.bar(ind + (2+1)*width, data[:,0]*1e3, yerr=err[:,0]*1e3, width=width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
    
    
    #clr= colors[5]
    
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
 
#ax.grid()
ax.set_ylabel('NEP anomaly [TgC/yr]',  fontsize=fz)
 
ax.set_xticks(ind+3*width)
#ax.set_ylim([-170,0])
ax.set_ylim([-170,170])
#ax.set_ylim([-210,0])
ax.locator_params(axis='y',nticks=2)
 
ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')

#lgnames.append("ORC")
#lgnames.append("ORC-MICT")
#lgnames.append("CABLE")
#lgnames.append("DLEM")
#lgnames.append("ISAM")
#lgnames.append("VEGAS")
#lgnames.append("VISIT")
#lgnames.append("SiBCASA")
#lgnames.append("BEPS")

lgnames.append("TBMs")  #NEP
#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax.legend( rects, lgnames, loc = 4, ncol=2, frameon=False, fontsize = fz*0.75)
plt.tick_params(labelsize=fz*0.9)  

#--------------------------------------
#fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
       
fig.tight_layout() 

fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_TBMsv2-NEP'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)







############################################ 
    
ind = np.arange(2)   # the x locations for the groups
width = 0.75    #0.10      # the width of the bars

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 6
fig_size[1] = 4.5
plt.rcParams["figure.figsize"] = fig_size

colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["royalblue", "gold", "r", "g", "tan", "chocolate" ,"crimson", "IndianRed", "forestgreen"]
          
fig = plt.figure(figsize=(8,4))
line =1
row=1

fz=18  #20
######################################################
# subplot R1
#from brokenaxes import brokenaxes


ax=fig.add_subplot(line,row,1)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 


data=np.zeros(shape=(2,1))
err=np.zeros(shape=(2,1))
#data[0,:] = data11[5:]
#data[1,:] = data22[5:]

data[0,:] = np.nanmedian(data11[15:24])
data[1,:] = np.nanmedian(data22[15:24])


err[0,:] = np.nanstd(data11[15:24])
err[1,:] = np.nanstd(data22[15:24])


#for i in np.arange(0,9):
#    clr= colors[i]
#    rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )

clr= colors[0]
#rects.append(ax.errorbar(ind + (0+1)*width, data[:,0]/12.0*1e3, width, color=clr, yerr= err[:,0]/12.0*1e3,  ecolor=clr) )
#rects.append(ax.bar(ind + (2+1)*width, data[:,0]/12.0*1e3, width, color=clr) )
rects.append(ax.bar(ind + (2+1)*width, data[:,0]*1e3, yerr=err[:,0]*1e3, width=width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
    
    
    #clr= colors[5]
    
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
 
#ax.grid()
ax.set_ylabel('GPP anomaly [TgC/yr]',  fontsize=fz)
 
ax.set_xticks(ind+3*width)
#ax.set_ylim([-170,0])
ax.locator_params(axis='y',nticks=2)
 
ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')

#lgnames.append("ORC")
#lgnames.append("ORC-MICT")
#lgnames.append("CABLE")
#lgnames.append("DLEM")
#lgnames.append("ISAM")
#lgnames.append("VEGAS")
#lgnames.append("VISIT")
#lgnames.append("SiBCASA")
#lgnames.append("BEPS")

lgnames.append("TBMs GPP")
#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax.legend( rects, lgnames, loc = 1, ncol=2, frameon=False, fontsize = fz*0.75)
plt.tick_params(labelsize=fz*0.9)  

#--------------------------------------      
#fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
       
fig.tight_layout() 

fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_TBMsv2-GPP'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)






   

############################################ 
    
ind = np.arange(2)   # the x locations for the groups
width = 0.75    #0.10      # the width of the bars

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 6
fig_size[1] = 4.5
plt.rcParams["figure.figsize"] = fig_size

colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["royalblue", "gold", "r", "g", "tan", "chocolate" ,"crimson", "IndianRed", "forestgreen"]
          
fig = plt.figure(figsize=(8,4))
line =1
row=1

fz=18  #20
######################################################
# subplot R1
#from brokenaxes import brokenaxes


ax=fig.add_subplot(line,row,1)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 


data=np.zeros(shape=(2,1))
err=np.zeros(shape=(2,1))
#data[0,:] = data11[5:]
#data[1,:] = data22[5:]

data[0,:] = np.nanmedian(data11[24:33])
data[1,:] = np.nanmedian(data22[24:33])


err[0,:] = np.nanstd(data11[24:33])
err[1,:] = np.nanstd(data22[24:33])


#for i in np.arange(0,9):
#    clr= colors[i]
#    rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )

clr= colors[0]
#rects.append(ax.errorbar(ind + (0+1)*width, data[:,0]/12.0*1e3, width, color=clr, yerr= err[:,0]/12.0*1e3,  ecolor=clr) )
#rects.append(ax.bar(ind + (2+1)*width, data[:,0]/12.0*1e3, width, color=clr) )
rects.append(ax.bar(ind + (2+1)*width, data[:,0]*1e3, yerr=err[:,0]*1e3, width=width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
    
    
    #clr= colors[5]
    
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
 
#ax.grid()
ax.set_ylabel('Reco anomaly [TgC/yr]',  fontsize=fz)
 
ax.set_xticks(ind+3*width)
#ax.set_ylim([-170,0])
ax.locator_params(axis='y',nticks=2)
 
ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')

#lgnames.append("ORC")
#lgnames.append("ORC-MICT")
#lgnames.append("CABLE")
#lgnames.append("DLEM")
#lgnames.append("ISAM")
#lgnames.append("VEGAS")
#lgnames.append("VISIT")
#lgnames.append("SiBCASA")
#lgnames.append("BEPS")

lgnames.append("TBMs Reco")
#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax.legend( rects, lgnames, loc = 1, ncol=2, frameon=False, fontsize = fz*0.75)
plt.tick_params(labelsize=fz*0.9)  

#--------------------------------------
#fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
       
fig.tight_layout() 

fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_TBMsv2-Reco'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)




 
############################################ 
    
ind = np.arange(2)   # the x locations for the groups
width = 0.375    #0.10      # the width of the bars

# Get current size
fig_size = plt.rcParams["figure.figsize"]

#print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 6
fig_size[1] = 4.5
plt.rcParams["figure.figsize"] = fig_size

colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["royalblue", "gold", "r", "g", "tan", "chocolate" ,"crimson", "IndianRed", "forestgreen"]
colors = ["r","b"]
          
fig = plt.figure(figsize=(8,4))
line =1
row=1

fz=18  #20
######################################################
# subplot R1
#from brokenaxes import brokenaxes


ax=fig.add_subplot(line,row,1)
 
rects=[]
lgnames=[]


#ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 


data=np.zeros(shape=(2,2))
err=np.zeros(shape=(2,2))
#data[0,:] = data11[5:]
#data[1,:] = data22[5:]

data[0,:] = data11[33:]
data[1,:] = data22[33:]


#err[0,:] = data11[15:]
#err[1,:] = data22[15:]

 
for i in np.arange(0,2):
    clr= colors[i]
    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
    rects.append(ax.bar(ind + (i+1)*width, data[:,i]*1e3, width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
            
#clr= colors[0]
##rects.append(ax.errorbar(ind + (0+1)*width, data[:,0]/12.0*1e3, width, color=clr, yerr= err[:,0]/12.0*1e3,  ecolor=clr) )
#rects.append(ax.bar(ind + (2+1)*width, data[:,0]/12.0*1e3, width, color=clr) )
 
    
    #clr= colors[5]
    
#for i in np.arange(0,5):
#    clr= colors[1+i]
#    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
 
#ax.grid()
ax.set_ylabel('GPP anomaly [TgC/yr]',  fontsize=fz)
 
ax.set_xticks(ind+1.5*width)
#ax.set_ylim([-170,0])
ax.locator_params(axis='y',nticks=2)
 
ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
 
#ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')
 
lgnames.append("FluxSat")
lgnames.append("GOSIF GPP")

#ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
ax.legend( rects, lgnames, loc = 4, ncol=2, frameon=False, fontsize = fz*0.75)
plt.tick_params(labelsize=fz*0.9)  

#--------------------------------------
#fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
       
fig.tight_layout() 

fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_GOSIFGPP'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
'''