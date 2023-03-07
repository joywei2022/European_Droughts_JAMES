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
#   Byrne GOSAT+insitu+TCCON CO2 inversion: 2009-2015
#################################################################

dir1='/Volumes/HEWEI_T5/Inversions/Byrne_inversions_GOSAT_surface_TCCON_2020/Byrneetal2020/monthly/'

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

    files = glob.glob(os.path.join(dir1,dataset+'*.nc')) 
    
    
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

dir2='/Volumes/HEWEI_T5/Inversions/NASA-CMS-Flux/ESSD_paper_v2/monthly_1Deg/'
  
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
    
    files=glob.glob(os.path.join(dir2,'CMS-Flux*.nc'))  
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
 
  
    data3 = data2[:-12*3, int(np.floor(90-latMax)):int(np.ceil(90-latMin)), int(np.floor(180+lonMin)):int(np.ceil(180+lonMax)) ]
     
  
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


dir3='/Volumes/HEWEI_T5/Inversions/GCASv2/'
  
def flux_anomaly_JiangGOSAT(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
#    M=0 
#    N=6
#    
#    latMin=33 #35
#    latMax=73 #70
#    lonMin=-14.5
#    lonMax=35
        
#    varstr = 'bio_monthly_opt'
    
    file=os.path.join(dir3,'posterior.fluxes.GCASv2.nc')  
    
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

dir4='/Volumes/HEWEI_T5/Inversions/'
  
def flux_anomaly_WuCCDAS(dataset, varstr, M, N, latMin, latMax, lonMin, lonMax):
         
    M=0 
    N=6
    
    latMin=33 #35
    latMax=73 #70
    lonMin=-14.5
    lonMax=35
      
    ##dataset ='CCDAS_6yr_2deg_assim_results'
    #varstr = 'rnep'
    #file=os.path.join(dir4,dataset+'/'+'rnep.nc')  
       
    dataset ='CCDAS_6yr_2deg_assim_results'
    varstr = 'nep'
    file=os.path.join(dir4,dataset+'/'+'CCDAS_co2+sm_prior.nc')  
    
    
    data1=[]
    lat=[]
    lon=[]
    
    ncf=nc.Dataset(file) 
    lat.extend(ncf.variables['lat'][:])  #latitude
    lon.extend(ncf.variables['lon'][:])   #longitude
    data1.append(ncf.variables[varstr][:]) 
    ncf.close()
    
    data1=np.array(data1[0])
    data2=data1
    
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
#   Scholze et al.(2019) GRL paper, CCDAS-CO2+SM+VOD inversions
#################################################################

dir5='/Volumes/HEWEI_T5/Inversions/'
  
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
#   EUROCOM - CTE, RAMS, FLEXPART,Jena, LUMIA, PYVAR_CHIMERE
#################################################################

dir6='/Volumes/HEWEI_T5/Inversions/EUROCOM-ICOS/'

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
 
    file= os.path.join(dir6,dataset)
    
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
    
    fn='All_Flux_timeseries_Europe_20201012+CCDAS'
         
     
    N=6 
     
    pydates=[]
    for i in range(1,12*N+1): 
        yr = 2010-1 + int(math.ceil(i/float(12)))
        if i%12 == 0 :
            mn = 12
        else:
            mn = i-int(math.floor(i/float(12)))*12
        tmp = datetime(yr, mn, 1)
        pydates.append(tmp)   
        
 
    
       
    # when varstr == 'ndvi' 
    N=15
    pydates2=[]
    for i in range(1,12*2*N+1): 
        yr = 2001-1 + int(math.ceil(i/float(24)))
        if i%2==0:  
            dd = 15
        else:
            dd = 1
        if i%24 == 0 :
            mn = 12
        else:
            mn = int(math.ceil( (i - int(math.floor(i/24.0))*24)/2.0  )  )
      
        tmp = dt.datetime(yr, mn, dd)
        pydates2.append(tmp)  
    
 
    N=15
    pydates3=[]
    for i in range(0,N): 
        for j in range(0,92): 
            t = dt.datetime(2001+i, 1, 1) + dt.timedelta(days=4*j)
            pydates3.append(t)             


    
  
    
    #------------------------------
#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-14.5 #-10
#    lonMax=35    #40
    
 
       
    #lonMin=-15
     
#    M=0  
#    N=6 
#      
#    dataset='Byrne2020_NEE_GOSAT_surface_TCCON'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_ens =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#        
#    dataset='Byrne2020_NEE_GOSAT_only'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_casa =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#      
#    dataset='Byrne2020_NEE_surface_only'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_fc =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#       
#    dataset='Byrne2020_NEE_TCCON_only'
#    varstr='ensemble_mean'
#    
#    ByrneGOSATnee_sib =  flux_anomaly_ByrneGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
#       
#    #------------------------------
#
 
    M=0 
    N=6
    
    dataset='CMS-Flux-NBE-2020.monthly.grid'           
    varstr = 'prior-NBE' 
    LiuGOSATnee =  flux_anomaly_LiuGOSAT(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)


   #------------------------------

    M=0 
    N=6
     
    dataset=' '        
    varstr = 'bio_monthly_org'
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
     
    dataset='ccdas_sm+vod'         #only prior, same as WuCCDAS 
    varstr = 'nep'
    #ScholzeCCDASnee =  flux_anomaly_ScholzeCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
    ScholzeCCDASnee =  flux_anomaly_WuCCDAS(dataset,varstr,M, N, latMin, latMax, lonMin, lonMax)
 
    
    
    #------------------------------
        
#    M=4 
#    N=6
    
#    latMin=33    #35
#    latMax=73    #70
#    lonMin=-14.5 #-10
#    lonMax=35    #40
    
#    varstr='bio_prior'
#    dataset = 'EUROCOM_JenaCarboScopeRegional_2006-2015.nc'
#    JenaRegionalnee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    varstr='bio_prior'
#    dataset = 'EUROCOM_LUMIA_2006-2015.nc'
#    LUMIAnee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
#    
#    varstr='fire'
#    LUMIAfire =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
 
    M=4 
    N=6
    
    latMin=33    #35
    latMax=73    #70
    lonMin=-14.5 #-10
    lonMax=35    #40
    
   
    dataset = 'EUROCOM_CarbonTrackerEurope_2006-2015.nc'
    varstr='bio_prior'
    eurocom1nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
    varstr='fire'
    eurocom1fire=  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
      
    
    varstr='bio_prior'
    dataset = 'EUROCOM_EnKF-RAMS_2006-2015.nc'
    eurocom2nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)


    varstr='bio_prior' 
    dataset = 'EUROCOM_FLEXINVERT_2006-2015.nc'
    eurocom3nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
   
    M=4 
    N=6
    
    varstr='bio_prior'
    dataset = 'EUROCOM_JenaCarboScopeRegional_2006-2015.nc'
    eurocom4nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
    
   
    dataset = 'EUROCOM_LUMIA_2006-2015.nc'     
    varstr='bio_prior'
    eurocom5nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)

    varstr='fire'
    eurocom5fire =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)


    varstr='bio_prior'
        
    dataset = 'EUROCOM_PYVAR-CHIMERE_2006-2015.nc'
    eurocom6nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    
    M=0 
    N=5
    varstr='bio_prior'
        
    dataset = 'EUROCOM_NAME-HB_2011-2015.nc'
    eurocom7nee =  flux_anomaly_eurocom(dataset, varstr,M, N, latMin, latMax, lonMin, lonMax)
    
    
    
    
    
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
    #JenaRegionalneeann = []
    #LUMIAneeann = []
    #ByrneGOSATneeann=[]
    LiuGOSATneeann=[]
    JiangGOSATneeann=[]
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
        #ByrneGOSATneeann.append(np.nanmean( ByrneGOSATnee_ens[flag][i*12:(i+1)*12]  ) )
        LiuGOSATneeann.append(np.nanmean( LiuGOSATnee[flag][i*12:(i+1)*12]  ) )
        JiangGOSATneeann.append(np.nanmean( JiangGOSATnee[flag][i*12:(i+1)*12]  ) ) 
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
    #ByrneGOSATneeann=np.array(ByrneGOSATneeann) 
    LiuGOSATneeann=np.array(LiuGOSATneeann) 
    JiangGOSATneeann=np.array(JiangGOSATneeann)   
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
    #ByrneGOSATneeann=signal.detrend(ByrneGOSATneeann)
    LiuGOSATneeann=signal.detrend(LiuGOSATneeann)
    JiangGOSATneeann=signal.detrend(JiangGOSATneeann)
    ScholzeCCDASneeann=signal.detrend(ScholzeCCDASneeann)
  
    IX = Year - 2010
 
#    ByrneGOSATnee_ens = signal.detrend(np.array(ByrneGOSATnee_ens[flag]))[12*IX:12*(IX+1)] 
#    ByrneGOSATnee_casa = signal.detrend(np.array(ByrneGOSATnee_casa[flag]))[12*IX:12*(IX+1)] 
#    ByrneGOSATnee_fc = signal.detrend(np.array(ByrneGOSATnee_fc[flag]))[12*IX:12*(IX+1)] 
#    ByrneGOSATnee_sib = signal.detrend(np.array(ByrneGOSATnee_sib[flag]))[12*IX:12*(IX+1)] 
#    LiuGOSATnee = signal.detrend(np.array(LiuGOSATnee[flag]))[12*IX:12*(IX+1)] 
#    JiangGOSATnee = signal.detrend(np.array(JiangGOSATnee[flag]))[12*IX:12*(IX+1)] 
##    WuCCDASnee1 = signal.detrend(np.array(WuCCDASnee1[flag]))[12*IX:12*(IX+1)] 
##    WuCCDASnee2 = signal.detrend(np.array(WuCCDASnee2[flag]))[12*IX:12*(IX+1)] 
#    ScholzeCCDASnee = signal.detrend(np.array(ScholzeCCDASnee[flag]))[12*IX:12*(IX+1)] 
#    JenaRegionalnee = signal.detrend(np.array(JenaRegionalnee[flag]))[12*IX:12*(IX+1)] 
                            
#    ByrneGOSATnee_ens =  np.array(ByrneGOSATnee_ens[flag])[12*IX:12*(IX+1)] 
#    ByrneGOSATnee_casa =  np.array(ByrneGOSATnee_casa[flag])[12*IX:12*(IX+1)] 
#    ByrneGOSATnee_fc =  np.array(ByrneGOSATnee_fc[flag])[12*IX:12*(IX+1)] 
#    ByrneGOSATnee_sib = np.array(ByrneGOSATnee_sib[flag])[12*IX:12*(IX+1)] 
#    LiuGOSATnee =  np.array(LiuGOSATnee[flag])[12*IX:12*(IX+1)] 
#    JiangGOSATnee =  np.array(JiangGOSATnee[flag])[12*IX:12*(IX+1)] 
##    WuCCDASnee1 =  np.array(WuCCDASnee1[flag])[12*IX:12*(IX+1)] 
##    WuCCDASnee2 =  np.array(WuCCDASnee2[flag])[12*IX:12*(IX+1)] 
#    ScholzeCCDASnee =  np.array(ScholzeCCDASnee[flag])[12*IX:12*(IX+1)] 
#    JenaRegionalnee =  np.array(JenaRegionalnee[flag])[12*IX:12*(IX+1)] 
#    LUMIAnee =  np.array(LUMIAnee[flag])[12*IX:12*(IX+1)]
    
    pydates1=pydates[12*IX:12*(IX+1)] 
    #pydates2=pydates2[24*IX:24*(IX+1)] 
    #pydates3=pydates3[92*IX:92*(IX+1)] 
    
   
    #####################
    #  plot out
    #####################
    fig = plt.figure(figsize=(6.75, 5))   #4.5, 5 #12,5*3-2
  
    fontsize=24
    mlw=1.2 
    
 
    ax3=fig.add_subplot(3,1,1)
     
    from matplotlib.ticker import FormatStrFormatter
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#     
#    ax3.plot(pydates1,-np.array(ByrneGOSATnee_ens[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='r',linestyle='-', marker='o', markersize = 4) 
#    ax3.plot(pydates1,-np.array(ByrneGOSATnee_casa[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='limegreen',linestyle='-', marker='s', markersize = 4) 
#    ax3.plot(pydates1,-np.array(ByrneGOSATnee_fc[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='royalblue',linestyle='-', marker='^', markersize = 4) 
#    ax3.plot(pydates1,-np.array(ByrneGOSATnee_sib[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='tan',linestyle='-', marker='x', markersize = 4) 
#  
#    ax3.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
#    ax3.set_ylabel('NEP anomaly [TgC/month]', fontsize=fontsize*0.5)
#    
#     
#    ax3.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
#    ax3.locator_params(nbins=8)
#    
#    #if Year ==2015:
#    #    ax3.legend(['Byrne2020_ensemble','Byrne2020_CASA','Byrne2020_FLUXCOM','Byrne2020_SiB3'], ncol = 2,  frameon=True, loc = 3, fontsize = fontsize*0.32)
#    
#    if Year ==2012:
#        ax3.legend(['Byrne2020_ensemble','Byrne2020_CASA','Byrne2020_FLUXCOM','Byrne2020_SiB3'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.32)
#        
#    ax3.set_ylim([-90,90])
#    
#    yy = np.zeros(len(pydates1))
#    
#    if flag == 1:
#    #  ax1.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
#      #ax2.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
#      ax3.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
#    
#    ax3.tick_params(axis='both', which='major', labelsize=fontsize*0.5)    
#    
#    #ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
#    #ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
#    
#    ax3.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
#    ax3.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
#    
#    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax3.get_xticklabels()]
#    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax3.get_yticklabels()]
# 
#    plt.setp(ax3.get_xticklabels(), visible=False) 
    
    
 
    
       
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
    
    ax2=fig.add_subplot(3,1,2)
    
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    
    mean = np.nanmean([-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,  \
                       -np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,-eurocom7neenew[12*IX:12*(IX+1)]/12.0*1e3  ], axis=0)  
    ax2.plot(pydates1,mean,linewidth=mlw*1.2,color='k',linestyle='-', marker='o', markersize = 4)   
  
    
    ax2.plot(pydates1,-np.array(eurocom1nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-')  
    ax2.plot(pydates1,-np.array(eurocom2nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-', )         
    ax2.plot(pydates1,-np.array(eurocom3nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-') 
    ax2.plot(pydates1,-np.array(eurocom4nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-') 
    ax2.plot(pydates1,-np.array(eurocom5nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-') 
    ax2.plot(pydates1,-np.array(eurocom6nee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-') 
    ax2.plot(pydates1,-np.array(eurocom7neenew)[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='skyblue',linestyle='-') 
       
    
   
    ax2.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    ax2.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.5)
    
     
    ax2.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax2.locator_params(nbins=8)
    
#    if Year ==2015:
#        ax2.legend(['Byrne2020_ensemble','Liu2020','Jiang2020','CCDAS_SM','CCDAS_SM+FAPAR','CCDAS_SM+VOD'], ncol = 3,  frameon=True, loc = 3, fontsize = fontsize*0.32)
 
    if Year ==2012:
        #ax2.legend(['CTE','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Mean'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.32)
        ax2.legend(['EUROCOM ensemble','individual'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.35)
  
    ax2.set_ylim([-110,110])   #[-90,90]
     
    yy = np.zeros(len(pydates1))
    
    if flag == 1:
        ax2.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
    
    ax2.tick_params(axis='both', which='major', labelsize=fontsize*0.5)    
    
    #ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    #ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    
    ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax2.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax2.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax2.get_yticklabels()]
    
    
    
        
    ax1=fig.add_subplot(3,1,3)
    
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
#    ax2.plot(pydates1,-np.array(JenaRegionalnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='royalblue',linestyle='-', marker='o', markersize = 4)  
#    ax2.plot(pydates1,-np.array(LUMIAnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='gold',linestyle='-', marker='o', markersize = 4)         
#    ax1.plot(pydates1,-np.array(ByrneGOSATnee_ens[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='r',linestyle='-', marker='o', markersize = 4) 
    ax1.plot(pydates1,-np.array(LiuGOSATnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='royalblue',linestyle='-', marker='s', markersize = 4) 
    ax1.plot(pydates1,-np.array(JiangGOSATnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='g',linestyle='-', marker='o', markersize = 4) 
#    ax2.plot(pydates1,WuCCDASnee1,linewidth=mlw,color='pink',linestyle='-', marker='o', markersize = 4) 
#    ax2.plot(pydates1,WuCCDASnee2,linewidth=mlw,color='royalblue',linestyle='-', marker='s', markersize = 4) 
    ax1.plot(pydates1,np.array(ScholzeCCDASnee[flag])[12*IX:12*(IX+1)]/12.0*1e3,linewidth=mlw,color='tan',linestyle='-', marker='o', markersize = 4) 
       
    ax1.xaxis.grid(color='lightgray', linestyle='--', linewidth=1)
    ax1.set_ylabel('\u0394NEP [TgC mon$^{-1}$]', fontsize=fontsize*0.5)
    
     
    ax1.tick_params(axis='both', which='major',  top='off', labelsize=fontsize*0.4)
    ax1.locator_params(nbins=8)
    
#    if Year ==2015:
#        ax2.legend(['Byrne2020_ensemble','Liu2020','Jiang2020','CCDAS_SM','CCDAS_SM+FAPAR','CCDAS_SM+VOD'], ncol = 3,  frameon=True, loc = 3, fontsize = fontsize*0.32)
 
    if Year ==2012:
        ax1.legend(['CMS-Flux2020','GCASv2','CCDAS_SM+VOD'], ncol = 2,  frameon=False, loc = 2, fontsize = fontsize*0.35)
  
    ax1.set_ylim([-110,110])
     
    yy = np.zeros(len(pydates1))
    
    if flag == 1:
        ax1.plot(pydates1,yy, color='gray',linestyle='--',linewidth=1)
    
    ax1.tick_params(axis='both', which='major', labelsize=fontsize*0.5)    
    
    #ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    #ax2.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    
    ax1.xaxis.set_major_locator(pltdt.MonthLocator([3,6,9,12]))
    ax1.xaxis.set_major_formatter(pltdt.DateFormatter('%Y%b')) #
    
    dummy = [lab.set_fontsize(0.2*2*fontsize) for lab in ax1.get_xticklabels()]
    dummy = [lab.set_fontsize(0.23*2*fontsize) for lab in ax1.get_yticklabels()]
    
    
    
        
    fig = plt.gcf()
    fig.set_size_inches(4.5, 7.5)  #(6,10)
     
    
    fig.tight_layout() 
    plt.subplots_adjust(hspace = 0.15)   #0.1
 
    
    fn='All_Flux_timeseries_Europe_GOSATInv+CCDAS+Jenaregional_20221008_prior_'+str(Year)
    outp=os.path.join(dir1,fn+'.png')
    print(outp)
    plt.savefig(outp,dpi=300)
    
        
    data=[]
    
#    data.append(np.nansum(-JenaRegionalnee))
#    data.append(np.nansum(-LUMIAnee))
#    data.append( np.nansum(-ByrneGOSATnee_ens))
##    data.append(ByrneGOSATnee_casa)
##    data.append(ByrneGOSATnee_fc)
##    data.append(ByrneGOSATnee_sib)
##    data.append(LiuGOSATnee)
#    data.append(np.nansum(-JiangGOSATnee))
##    data.append(WuCCDASnee1)    
##    data.append(WuCCDASnee2)    
#    data.append(np.nansum(ScholzeCCDASnee))   

   # data.append(-JenaRegionalneeann[IX])
   # data.append(-LUMIAneeann[IX])
    data.append(-eurocom1neeann[IX])
    data.append(-eurocom2neeann[IX])
    data.append(-eurocom3neeann[IX])
    data.append(-eurocom4neeann[IX])
    data.append(-eurocom5neeann[IX])
    data.append(-eurocom6neeann[IX])
    data.append(-eurocom7neeannnew[IX])
    
    #annualmean = np.nanmean([-eurocom1neeann[IX],-eurocom1neeann[IX],-eurocom1neeann[IX],-eurocom1neeann[IX], \
    #                    -eurocom1neeann[IX], -eurocom1neeann[IX], -eurocom1neeann[IX] ], axis=0)  
 
    #data.append(annualmean)
    
    #data.append(-ByrneGOSATneeann[IX])
    data.append(-LiuGOSATneeann[IX])
    data.append(-JiangGOSATneeann[IX])
    data.append(ScholzeCCDASneeann[IX])
     
    data=np.array(data)
    
    print("CMS-Flux prior for year%s"%Year)
    print(-LiuGOSATneeann[IX])
    
    return data    
    
       
import timeit
start = timeit.timeit()

flag=1

#data = pipeline(35, 70, -10, 40, flag,Year,'All_Flux_timeseries_Europe_20180801')

#Year = 2003 
##data = pipeline(40, 50, 0, 25, flag, Year,'All_Flux_timeseries_Europe_20190609')
#data = pipeline(40, 55, 0, 25, flag, Year,'All_Flux_timeseries_Europe_20200428')
#
#Year = 2006
#data = pipeline(53, 68, 10, 30, flag, Year, 'All_Flux_timeseries_Europe_20200428')
#
#Year = 2012
#data = pipeline(42, 52, 15, 40, flag, Year, 'All_Flux_timeseries_Europe_20200428')
#
#Year = 2015
##data = pipeline(47, 55, 20, 30, flag, Year, 'All_Flux_timeseries_Europe_20190609')
#data = pipeline(45, 55, 15, 30, flag, Year, 'All_Flux_timeseries_Europe_20200428')

 
Year = 2012
#data = pipeline(35, 70, -10, 40, flag, Year, 'All_Flux_timeseries_Europe_20201012_benchmark+CCDAS')
data1 = pipeline(35, 50, 10, 40, flag, Year, 'All_Flux_timeseries_Europe_20210601_benchmark+CCDAS+JenaRegional')

Year = 2015
data2 = pipeline(45, 65, 0, 35, flag, Year, 'All_Flux_timeseries_Europe_20210601_benchmark+CCDAS+JenaRegional')



#2012
np.array([-0.08654909, -0.00105949, -0.16559372, -0.16678975, -0.2540477 ,
       -0.25102952, -0.20727155, 0, -0.00097638, -0.06956293, -0.00703166])


#2015

np.array([ 0.00852351,  0.00016917, -0.11863217, -0.15051872, -0.07982911,
       -0.10496698, -0.14281232,  0, -0.0028067, 0.01553903, -0.02054539])



#
############################################# 
#    
#ind = np.arange(2)   # the x locations for the groups
#width = 0.11  #5          # the width of the bars
#
## Get current size
#fig_size = plt.rcParams["figure.figsize"]
#
##print "Current size:", fig_size
#
## Set figure width to 12 and height to 9
#fig_size[0] = 6
#fig_size[1] = 4.5
#plt.rcParams["figure.figsize"] = fig_size
#
##colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = ["r", "limegreen", "orange", "royalblue", "gold", "magenta", "cyan", "k"]
#          
#fig = plt.figure(figsize=(8,4))
#line =1
#row=1
#
#fz=18  #20
#######################################################
## subplot R1
##from brokenaxes import brokenaxes
#
#
#ax=fig.add_subplot(line,row,1)
# 
#rects=[]
#lgnames=[]
#
#
##ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 
#
#
#data=np.zeros(shape=(2,8))
#data[0,:] = data1[:8]
#data[1,:] = data2[:8]
#
#for i in np.arange(0,8):
#    clr= colors[i]
#    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
#    rects.append(ax.bar(ind + (i+1)*width, data[:,i]*1e3, width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
#        
#        
#    #clr= colors[5]
#    
##for i in np.arange(0,5):
##    clr= colors[1+i]
##    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
# 
##ax.grid()
#ax.set_ylabel('NEP anomaly [TgC/yr]',  fontsize=fz)
# 
#ax.set_xticks(ind+4.5*width)
##ax.set_ylim([-170,0])
#ax.set_ylim([-300,300])
#ax.locator_params(axis='y',nticks=2)
# 
#ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
# 
##ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')
#
#lgnames.append("CTE")
#lgnames.append("EnKF-RAMS")
#lgnames.append("FLEXINVERT")
#lgnames.append("CarboScope-Regional")
#lgnames.append("LUMIA")
#lgnames.append("PYVAR-CHIMERE")
#lgnames.append("NAME-HB")
#lgnames.append("Mean")
#
##ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
#ax.legend( rects, lgnames, loc = 2, ncol=2, frameon=False, fontsize = fz*0.8)
#plt.tick_params(labelsize=fz*0.9)  
#
##--------------------------------------
##fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
#plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
#       
#fig.tight_layout() 
#
#fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_EUROCOM_prior'
#outp=os.path.join(dir1,fn+'.png')
#print outp
#plt.savefig(outp,dpi=300)
#
#
#  
# 
#
#
#######################################################
#   
#ind = np.arange(2)   # the x locations for the groups
#width = 0.15          # the width of the bars
#
## Get current size
#fig_size = plt.rcParams["figure.figsize"]
#
##print "Current size:", fig_size
#
## Set figure width to 12 and height to 9
#fig_size[0] = 6
#fig_size[1] = 4.5
#plt.rcParams["figure.figsize"] = fig_size
#
##colors = ["gray","crimson", "IndianRed", "forestgreen", "limegreen", "blue","royalblue","chocolate", "sandybrown"]
#colors = [ "g", "tan"]
#          
#fig = plt.figure(figsize=(8,4))
#line =1
#row=1
#
#fz=18  #20
#######################################################
## subplot R1
##from brokenaxes import brokenaxes
#
#
#ax=fig.add_subplot(line,row,1)
# 
#rects=[]
#lgnames=[]
#
#
##ax = brokenaxes( ylims=((0, 150), (400, 450)), hspace=0.05, despine=False) 
#
#
#data=np.zeros(shape=(2,3-1))
#data[0,:] = data1[8:10]
#data[1,:] = data2[8:10]
#
#for i in np.arange(0,3-1):
#    clr= colors[i]
#    #rects.append(ax.bar(ind + (i+1)*width, data[:,i]/12.0*1e3, width, color=clr) )
#    rects.append(ax.bar(ind + (i+1)*width, data[:,i]*1e3, width, color='white', ecolor=clr , edgecolor=clr , hatch="//\\") )
#        
#        
#    #clr= colors[5]
#    
##for i in np.arange(0,5):
##    clr= colors[1+i]
##    rects.append(ax.bar(2*ind+ (6+(i+1))*width, data2[i], width, color=clr))
# 
##ax.grid()
#ax.set_ylabel('NEP anomaly [TgC/yr]',  fontsize=fz)
# 
#ax.set_xticks(ind+1.5*width)
#ax.set_ylim([-170,170])
#ax.locator_params(axis='y',nticks=2)
# 
#ax.set_xticklabels( ['2012','2015'],rotation=0, fontsize=fz*0.9)
# 
##ax.legend( (rects1[0], rects2[0]), ('Prior', 'Optimized'), fontsize = 'x-small')
# 
##lgnames.append("Byrne2020")
#lgnames.append("GCASv2")
#lgnames.append("CCDAS_SM+VOD")
# 
#
##ax.legend( (rectspri[0], rectsopt[0],rectspri[1], rectsopt[1], rectspri[2], rectsopt[2]), ('Prior', 'Optimized'), fontsize = 'x-small')
#ax.legend( rects, lgnames, loc = 2, ncol=2, frameon=False, fontsize = fz*0.8)
#plt.tick_params(labelsize=fz*0.9)  
#
##--------------------------------------
##fig.subplots_adjust(bottom=bottom,left=left,right=right,top=top,wspace=wspace,hspace=hspace)
#plt.subplots_adjust(left=0.002,wspace=0.00,hspace=0.02)
#       
#fig.tight_layout() 
#
#fn='Bar_subcontinental_anomalies_drought_events_Europe_20210601_GOSATInv+CCDAS_prior'
#outp=os.path.join(dir1,fn+'.png')
#print outp
#plt.savefig(outp,dpi=300)
#
#
#        