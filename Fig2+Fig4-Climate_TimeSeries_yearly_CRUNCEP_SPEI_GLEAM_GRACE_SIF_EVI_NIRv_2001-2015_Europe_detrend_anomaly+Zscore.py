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
import matplotlib.dates as pltdt
#import matplotlib.ticker as mtick
import math
#import numpy.ma as ma
#from calendar import monthrange

 
    
import scipy

from scipy import signal
from scipy.stats import pearsonr

 
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
        print( 'total earth area    = ', 4 * np.pi * radius ** 2)
    return area
 


#################################################################################
#                              CRUNCEP 
####################################################################################
dir1='/Volumes/HEWEI_T5/CRUNCEPv9_1981-2017/'

def anomaly_CRUNCEP(varstr, M, N,latMin, latMax, lonMin, lonMax):
 
    
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
#    M=0      
#    N=15     
#    varstr='tair'



    files=glob.glob(os.path.join(dir1, varstr,  'cruncep_v9*.nc'))
    
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
    
    
   
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        
        if varstr=='tair':
            #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
            data4[i] = np.nanmean(data3[(12*i+0):(12*i+5)],axis=0)
            #data4[i] = np.nanmean(data3[(12*i+2):(12*i+5)],axis=0)
            #data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
           
        if varstr=='rain'  or varstr=='vpd':
            data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
            
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)    
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b


 
####################################################################################
#                              SPEI index
####################################################################################

dir2='/Volumes/HEWEI_T5/SPEIbasev2.5_to2015/'

def anomaly_spei(varstr, M, N, latMin,latMax,lonMin,lonMax):

    
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
    M=0      
    N=15     
    varstr='spei'

    file=os.path.join(dir2, 'spei03.nc') #09
    
    data1=[]
    lon=[]
    lat=[]
     
    ncf=nc.Dataset(file)
    data1.append(ncf.variables[varstr][:])   
    lon.append(ncf.variables['lon'][:]) 
    lat.append(ncf.variables['lat'][:]) 
         
    ncf.close()
    
    data1=np.array(data1)
    lon=np.array(lon)
    lat=np.array(lat)
     
    #print data1.shape,lon.shape   #(1, 1368, 360, 720)  
 
    data2=[]
    
    for i in range(100*12,115*12):
        data2.append(np.flipud(data1[0,i]))
    
    data2=np.array(data2)
    
    data3 = data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]

       
    mask=(data3>=1.0e30)
    data3[mask]=np.NaN
       
    plt.imshow(data3[0]) 
    
    
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        #data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
        
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    #b = (a-mean)/std
    b = a

    return a,b


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
 
 

    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b



####################################################################################
#                              GRACE-EWT index
####################################################################################

dir4='/Volumes/HEWEI_T5/GRACE_REC/'

def anomaly_ewt(varstr, M, N, latMin,latMax,lonMin,lonMax):

    file=os.path.join(dir4, 'GRACE_REC_v02.nc') 
 
    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
    M=0      
    N=15  
 
    
    varstr='rec_ensemble_mean'
    
    data1=[]
    lon=[]
    lat=[]
     
    ncf=nc.Dataset(file)
    data1.append(ncf.variables[varstr][:])   
    lon.append(ncf.variables['lon'][:]) 
    lat.append(ncf.variables['lat'][:]) 
     
    ncf.close()
    
    data1=np.array(data1)
    lon=np.array(lon)
    lat=np.array(lat)

    #print data1.shape,lon.shape   #(1, 456, 360, 720)

    data2=[]
    
    for i in range(22*12,37*12):
        data2.append(np.flipud(data1[0,i]))
    
    data2=np.array(data2)
    
    data3 = data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  # corrected on July 3, 2017


    a=[]
    
    mask=(data3==-9999.0)
    data3[mask]=np.NaN
       
    plt.imshow(data3[0]) 
 
    
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        #data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b


 
####################################################################################
#                              CSIF  
####################################################################################

dir5='/Volumes/HEWEI_T5/CSIF_v2/monthly/'

def anomaly_CSIF(varstr, M, N, latMin, latMax, lonMin, lonMax):

    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
    M=0      
    N=15  
    
    varstr = 'clear_daily_SIF'
        
    data1=[]
    lon=[]
    lat=[]
     
    files=glob.glob(os.path.join(dir5,'OCO2.SIF.clear.*.nc'))
      
    ncf=nc.Dataset(files[0])    
    lon.append(ncf.variables['lon'][:])   
    lat.append(ncf.variables['lat'][:]) 
    
    for file in files[:-12*3]:
        ncf=nc.Dataset(file)
        data1.append(ncf.variables[varstr][:]) 
        
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
    data3 = data2[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]

    del data2
    
 
    mask=(data3==-999.9)
    data3[mask]=np.NaN
    
    plt.imshow(data3[0])
   
 
    
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        #data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b


 
####################################################################################
#                MODIS EVI index  (MOD13C2 original)  - NIRv
####################################################################################
dir6a='/Volumes/HEWEI_T5/MOD13C2/0.5Deg_NIRv1/'  

def anomaly_NIRv(M, N, VAR, latMin,latMax,lonMin,lonMax):
     

    latMin=33    #35
    latMax=73    #70
    lonMin=-15   #-10
    lonMax=35    #40
     
    M=0      
    N=15  
    
    
    VAR = 'NIRv1'
    
    files=glob.glob(os.path.join(dir6a,'*'+VAR+'.dat'))
    files1 = sorted(files)
    
    NIRv=[]
 
    for i in files1: 
            
        data = np.fromfile(i,dtype='float64', sep="")
        data = data.reshape(360, 720)
        NIRv.append(data)
 
    NIRv= np.array(NIRv)
 
    data3 = NIRv[:,2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  
   
    del NIRv
  

   
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        #data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b

 
    
 
####################################################################################
#                              GLOBMAP LAI data
####################################################################################
 
dir6='/Volumes/HEWEI_T5/GLOBMAP_LAI_V3_0.5Deg_monthly_2001-2020/'

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
     
    files=glob.glob(os.path.join(dir6,'GLOBMAP_LAI*.nc'))    
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
    
    data1=array(data1)    
    lon=array(lon)      
    lat=array(lat)     
    
    #print data1.shape,lon.shape   #(192, 720, 1440) 
    
    
    lon2, lat2 = np.meshgrid(lon,lat) 
    
    data2 = np.zeros(shape=(data1.shape[0],data1.shape[2],data1.shape[1],) )
    for i in range(0,data1.shape[0]):
       data2[i] =  np.rot90(np.fliplr(data1[i]))         # Added on July 3, 2017
     
    #data2 = data1
    plt.imshow(data2[0])
     
    data3 = data2[M:(-4),2*(90-latMax):2*(90-latMin),2*(180+lonMin):2*(180+lonMax)]  
    #data3 = data2[M*12:(-3*12),4*(90-latMax):4*(90-latMin),4*(180+lonMin):4*(180+lonMax)]  
     
    #del data1,data2
    
      
    mask=(data3>1000)
    data3[mask]=np.NaN
    
    plt.imshow(data3[0])
    
  
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        #data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        #data4[i] = np.nanmean(data3[(12*i+2):(12*i+8)],axis=0)
        data4[i] = np.nanmean(data3[(12*i+3):(12*i+9)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b




#    a=[]
#    for i in range(0,data3.shape[0]): 
#        me = np.nanmean(data3[i]) 
#        a.append(me)
#        
#    a=array(a)
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
#    
#    a1=[]
#    for k in range(0,N):     
#       a1.extend(baseline)  
#    
#    a1=array(a1)
#    
#    
#    a1 = np.nanmean(a)
#    A=a 
#    B=a-a1
# 
#    return A,B
 

 
    

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
          
      
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b

 

#################################################################
#   GOSIF GPP v2 (Xiao Jingfeng)
#################################################################

dir21='/Volumes/HEWEI_T5/SIF_XiaoJF/GOSIF-GPP_v2/nc/GOSIF-GPP_v2/'
 
  
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
          
      
   
    
    data4 = np.zeros(shape=[N, data3.shape[1],data3.shape[2] ])*np.nan
    for i in range(0,N):
        data4[i] = np.nanmean(data3[12*i:12*(i+1)],axis=0)
        
        
    plt.imshow(data4[0])
    
    
    a = np.nansum(np.nansum(data4,axis=1),axis=1)/(data3.shape[1]*data3.shape[2] )
    
    a=signal.detrend(a)
    mean = np.nanmean(a)
    std = np.nanstd(a)
    
    b = (a-mean)/std

    return a,b

    
####################################   
  
#def pipeline(latMin=40,latMax=50,lonMin=-110,lonMax=-100,flag=1,fn='R1'):

latMin=33    #35
latMax=73    #70
lonMin=-15   #-10
lonMax=35    #40
 
flag=1

fn='climate_anomanly_yearly_2001-2015_Europe_20221004'
 


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

 

#M=0      
#N=15     
#varstr='spei'
#xx1_spei =  anomaly_spei(varstr,M, N,latMin, latMax, lonMin, lonMax)

      
M=0    
N=15   
varstr='SMroot'
xx1_smroot =  anomaly_gleam(varstr,M, N,latMin, latMax, lonMin, lonMax)  
   
         
#M=0    
#N=15   
#varstr='rec_ensemble_mean'
#xx1_ewt =  anomaly_ewt(varstr,M, N,latMin, latMax, lonMin, lonMax)  
 

M=0   
N=15   
varstr='SIF'
xx1_csif =  anomaly_CSIF(varstr,M, N,latMin, latMax, lonMin, lonMax)
   
    
M=0   
N=15    
  
VAR = 'NIRv1'
xx1_nirv = anomaly_NIRv(M, N,VAR,latMin,latMax,lonMin,lonMax)



#M=0   
#N=15    
#           
#varstr = 'LAI'
#xx1_lai = anomaly_LAI(varstr,M, N,latMin,latMax,lonMin,lonMax)

  
#M=0
#N=15     
#varstr = 'GPP'
#fluxsat =  flux_anomaly_fluxsat(varstr,M, N, latMin, latMax, lonMin, lonMax)
#
#M=0
#N=15     
#varstr = 'GPP'
#gosifgpp =  flux_anomaly_gosifgpp(varstr,M, N, latMin, latMax, lonMin, lonMax)
          
  
#####################
#  plot out
#####################      
plt.clf()

#from mpl_toolkits.axes_grid1 import host_subplot
#import mpl_toolkits.axisartist as AA
  
fontsize= 32 #28  #22 #24
mlw=2  #1.5


N=15  
pydates1=[]
for i in range(1,N+1): 
    yr = 2001-1 + i    #2001-1
    tmp = dt.datetime(yr, 1, 1)
    pydates1.append(tmp) 
    
 



xx2_tair= np.array(xx1_tair[flag])
xx2_p= np.array(xx1_p[flag])
xx2_vpd=np.array(xx1_vpd[flag])
xx2_smroot= np.array(xx1_smroot[flag])
#xx2_lai= np.array(xx1_lai[flag])
xx2_csif= np.array(xx1_csif[flag])
xx2_nirv= np.array(xx1_nirv[flag])
  
 
 
 
   
    
    
######################################################################
#   plot figure
######################################################################
fig = plt.figure(figsize=(18,9))


N = 15
ind = np.arange(N)   # the x locations for the groups
width = 0.5  #0.08  #0.10         # the width of the bars





######################################################################
ax1=fig.add_subplot(6,1,1) 
#ax2 = host_subplot(312, axes_class=AA.Axes)

#xs=pydates1
#ys=signal.detrend(np.array(xx1_tair[0]))
#d = scipy.zeros(len(ys))
#ax2.fill_between(xs, ys, where=ys>=d, interpolate=True, color='limegreen')
#ax2.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')

#ax2.bar(ind+1*width,xx1_spei[flag],width,color='red') 
ax1.bar(ind+1*width,xx2_tair,width,color='red') 

 
#p2, =ax2.plot(pydates1,signal.detrend(np.array(xx1_spei[0])) ,linewidth=mlw,color='tan')  #magenta
#units='[$^\circ$C]'
units='[Z-score]'
ax1.set_ylabel('Ta %s'%units, fontsize= fontsize*0.4, color='k')

#ax2.tick_params(axis='y', which='major', labelsize=fontsize*0.3)  #'both'
#ax2.locator_params(tight=False, axis='y',nticks=8)
#ax2.locator_params(tight=False, axis='x',nticks=15)
plt.setp(ax1.get_xticklabels(), visible=False)
 
#DEV=max( abs(np.nanmin(np.array(signal.detrend(np.array(xx1_spei[flag])) ))), abs(np.nanmax(np.array(signal.detrend(np.array(xx1_spei[flag])) )))  )
DEV=max( abs(np.nanmin(np.array(xx2_tair))), abs(np.nanmax(np.array(xx2_tair)))  )
ax1.set_ylim([-DEV*1.1,DEV*1.1])
#ax2.set_ylim([-2.5, 2.5])
plt.xticks(ind+1*width,['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'],fontsize=fontsize*0.9)  #,'Mean'

 
######################################################################
ax2=fig.add_subplot(6,1,2) 
 
#xs=pydates1
#ys=signal.detrend(np.array(xx1_ewt[flag]))
#d = scipy.zeros(len(ys))
#ax4.fill_between(xs, ys, where=ys>=d, interpolate=True, color='limegreen')
#ax4.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')


#ax4.bar(ind+1*width,xx1_csif[flag], width, color='limegreen') 
ax2.bar(ind+1*width,xx2_p,width,color='blue') 
#ax4.bar(ind+1*width,xx1_csif[flag],width,color='limegreen') 
#ax4.bar(ind+4*width,xx1_smroot[flag],width,color='skyblue') 
#ax4.bar(ind+1*width,xx1_ewt[flag],width,color='royalblue') 



#p3, =ax3.plot(pydates1,xx1_smroot[flag],linewidth=mlw,color='blue')  #navy   
#units='[mm]'
units='[Z-score]'
#ax4.set_ylabel('TWS %s'%units, fontsize= fontsize*0.6, color='k')
ax2.set_ylabel('Prep %s'%units, fontsize= fontsize*0.4, color='k')

#ax1.set_xticks(ind+1*width);
plt.xticks(ind+1*width,['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'],fontsize=fontsize*0.45)  #,'Mean'

#ax1.set_xticks(ind+1*width) 
#ax1.set_xticklabels( ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'] , rotation=30, fontsize= fontsize*0.7, color='k')

#ax3.tick_params(axis='y', which='major', labelsize=fontsize*0.3)  #'both'
#ax4.locator_params(tight=False, axis='y',nticks=8)
#ax4.locator_params(tight=False, axis='x',nticks=15)
plt.setp(ax2.get_xticklabels(), visible=False)
 

 

DEV=max( abs(np.nanmin(np.array(xx2_p))), abs(np.nanmax(np.array(xx2_p)))  )
 
#DEV=max( abs(np.nanmin(np.array(xx1_nirv[flag]))), abs(np.nanmax(np.array(xx1_nirv[flag])))  )
ax2.set_ylim([-DEV*1.2,DEV*1.2])
    



  
######################################################################
ax3=fig.add_subplot(6,1,3) 
 
#xs=pydates1
#ys=signal.detrend(np.array(xx1_ewt[flag]))
#d = scipy.zeros(len(ys))
#ax4.fill_between(xs, ys, where=ys>=d, interpolate=True, color='limegreen')
#ax4.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')


#ax4.bar(ind+1*width,xx1_csif[flag], width, color='limegreen') 
ax3.bar(ind+1*width,xx2_vpd,width,color='magenta') 
#ax4.bar(ind+1*width,xx1_csif[flag],width,color='limegreen') 
#ax4.bar(ind+4*width,xx1_smroot[flag],width,color='skyblue') 
#ax4.bar(ind+1*width,xx1_ewt[flag],width,color='royalblue') 



#p3, =ax3.plot(pydates1,xx1_smroot[flag],linewidth=mlw,color='blue')  #navy   
#units='[mm]'
units='[Z-score]'
#ax4.set_ylabel('TWS %s'%units, fontsize= fontsize*0.6, color='k')
ax3.set_ylabel('VPD %s'%units, fontsize= fontsize*0.4, color='k')

#ax1.set_xticks(ind+1*width)
plt.xticks(ind+1*width,['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'],fontsize=fontsize*0.45)  #,'Mean'

#ax1.set_xticks(ind+1*width) 
#ax1.set_xticklabels( ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'] , rotation=30, fontsize= fontsize*0.7, color='k')

#ax3.tick_params(axis='y', which='major', labelsize=fontsize*0.3)  #'both'
#ax4.locator_params(tight=False, axis='y',nticks=8)
#ax4.locator_params(tight=False, axis='x',nticks=15)
plt.setp(ax3.get_xticklabels(), visible=False)
 

 

DEV=max( abs(np.nanmin(np.array(xx2_vpd))), abs(np.nanmax(np.array(xx2_vpd)))  )
 
#DEV=max( abs(np.nanmin(np.array(xx1_nirv[flag]))), abs(np.nanmax(np.array(xx1_nirv[flag])))  )
ax3.set_ylim([-DEV*1.2,DEV*1.2])
    

######################################################################
ax4=fig.add_subplot(6,1,4) 
 
#xs=pydates1
#ys=xx1_smroot[flag]
#d = scipy.zeros(len(ys))
#ax3.fill_between(xs, ys, where=ys>=d, interpolate=True, color='limegreen')
#ax3.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')

ax4.bar(ind+1*width,xx2_smroot,width,color='royalblue') 

#p3, =ax3.plot(pydates1,xx1_smroot[flag],linewidth=mlw,color='blue')  #navy   
#units='[m$^{3}$m$^{-3}$]'
units='[Z-score]'
ax4.set_ylabel('SM %s'%units, fontsize= fontsize*0.4, color='k')

#ax3.tick_params(axis='y', which='major', labelsize=fontsize*0.3)  #'both'
#ax3.locator_params(tight=False, axis='y',nticks=8)
#ax3.locator_params(tight=False, axis='x',nticks=15)
plt.setp(ax4.get_xticklabels(), visible=False)
  
DEV=max( abs(np.nanmin(np.array(xx2_smroot))), abs(np.nanmax(np.array(xx2_smroot)))  )
ax4.set_ylim([-DEV*1.1,DEV*1.1])
    
plt.xticks(ind+1*width,['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'],fontsize=fontsize*0.9)  #,'Mean'

  
######################################################################
ax5=fig.add_subplot(6,1,5) 
 
#xs=pydates1
#ys=signal.detrend(np.array(xx1_ewt[flag]))
#d = scipy.zeros(len(ys))
#ax4.fill_between(xs, ys, where=ys>=d, interpolate=True, color='limegreen')
#ax4.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')


#ax4.bar(ind+1*width,xx1_csif[flag], width, color='limegreen') 
ax5.bar(ind+1*width,xx2_nirv,width,color='gold') 
#ax4.bar(ind+1*width,xx1_csif[flag],width,color='limegreen') 
#ax4.bar(ind+4*width,xx1_smroot[flag],width,color='skyblue') 
#ax4.bar(ind+1*width,xx1_ewt[flag],width,color='royalblue') 



#p3, =ax3.plot(pydates1,xx1_smroot[flag],linewidth=mlw,color='blue')  #navy   
#units='[mm]'
units='[Z-score]'
#ax4.set_ylabel('TWS %s'%units, fontsize= fontsize*0.6, color='k')
ax5.set_ylabel('NIRv %s'%units, fontsize= fontsize*0.4, color='k')

#ax1.set_xticks(ind+1*width)
plt.xticks(ind+1*width,['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'],fontsize=fontsize*0.45)  #,'Mean'

#ax1.set_xticks(ind+1*width) 
#ax1.set_xticklabels( ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'] , rotation=30, fontsize= fontsize*0.7, color='k')

#ax3.tick_params(axis='y', which='major', labelsize=fontsize*0.3)  #'both'
#ax4.locator_params(tight=False, axis='y',nticks=8)
#ax4.locator_params(tight=False, axis='x',nticks=15)
plt.setp(ax5.get_xticklabels(), visible=False)
 

 

#DEV=max( abs(np.nanmin(np.array(xx2_lai))), abs(np.nanmax(np.array(xx2_lai)))  )
 
DEV=max( abs(np.nanmin(np.array(xx1_nirv[flag]))), abs(np.nanmax(np.array(xx1_nirv[flag])))  )
ax5.set_ylim([-DEV*1.2,DEV*1.2])
    

######################################################################
ax6=fig.add_subplot(6,1,6)
#         
#ax1 = host_subplot(311, axes_class=AA.Axes)
#plt.subplots_adjust(right=0.75)
 
#xs=pydates1
#ys=xx1_csif[flag]
#d = scipy.zeros(len(ys))
#ax1.fill_between(xs, ys, where=ys>=d, interpolate=True, color='limegreen')
#ax1.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')
 
ax6.bar(ind+1*width,xx2_csif,width,color='limegreen') 

plt.xticks(ind+1*width,['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'],fontsize=fontsize*0.5)  #,'Mean'

 

#p1, =ax1.plot(pydates1,xx1_t[flag],linewidth=mlw,color='r') 
#units='[mW m$^{-2}$nm$^{-1}$sr$^{-1}$]'
units='[Z-score]'
ax6.set_ylabel('CSIF %s'%units,  fontsize= fontsize*0.4,color='k')
 
#ax1.tick_params(axis='y', which='major', labelsize=fontsize*0.3)  #'both'
#ax1.locator_params(tight=False, axis='y',nticks=8)
#ax1.locator_params(tight=False, axis='x',nticks=15)

plt.setp(ax6.get_xticklabels(), visible=True)
 
DEV=max( abs(np.nanmin(np.array(xx2_csif))), abs(np.nanmax(np.array(xx2_csif)))  )
ax6.set_ylim([-DEV*1.1,DEV*1.1])
  


#ax1.set_ylim([-2,2])
#ax2.set_ylim([-2.5,2.5])
#ax3.set_ylim([-2.5,2.5])
#ax4.set_ylim([-2.5,2.5])


yy = np.zeros(len(pydates1))

if flag == 1:
  ax1.plot(ind+1*width,yy, color='gray',linestyle='-',linewidth=1)
  ax2.plot(ind+1*width,yy, color='gray',linestyle='-',linewidth=1)
  ax3.plot(ind+1*width,yy, color='gray',linestyle='-',linewidth=1)
  ax4.plot(ind+1*width,yy, color='gray',linestyle='-',linewidth=1)
  ax5.plot(ind+1*width,yy, color='gray',linestyle='-',linewidth=1)
  ax6.plot(ind+1*width,yy, color='gray',linestyle='-',linewidth=1)
    
#frame1 = plt.gca()
#frame1.axes.get_xaxis().set_visible(True)

 
ax1.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth=1)
ax2.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth=1)
ax3.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth=1)
ax4.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth=1)
ax5.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth=1)
ax6.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth=1)
    
#dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax1.get_xticklabels()]
dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax1.get_yticklabels()]
dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax2.get_yticklabels()]
dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax3.get_yticklabels()]
dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax4.get_yticklabels()]
dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax5.get_yticklabels()]
dummy = [lab.set_fontsize(0.45*fontsize) for lab in ax6.get_yticklabels()]

dummy = [lab.set_fontsize(0.5*fontsize) for lab in ax6.get_xticklabels()]


#ax1.yaxis.get_major_formatter().set_powerlimits((0,2))
#ax2.yaxis.get_major_formatter().set_powerlimits((0,2))
#ax3.yaxis.get_major_formatter().set_powerlimits((0,2))
#ax4.yaxis.get_major_formatter().set_powerlimits((0,2))
  
        
#del xx1_t 
#del xx1_p
#
#del xx1_smroot    
#del xx1_ewt
 
#ax1.xaxis.set_major_locator(pltdt.MonthLocator([1]))
#ax1.xaxis.set_major_formatter(pltdt.DateFormatter('%Y'))
#
#ax2.xaxis.set_major_locator(pltdt.MonthLocator([1]))
#ax2.xaxis.set_major_formatter(pltdt.DateFormatter('%Y'))
#
#ax3.xaxis.set_major_locator(pltdt.MonthLocator([1]))
#ax3.xaxis.set_major_formatter(pltdt.DateFormatter('%Y'))
#
#ax4.xaxis.set_major_locator(pltdt.MonthLocator([1]))
#ax4.xaxis.set_major_formatter(pltdt.DateFormatter('%Y'))


#ax4.tick_params(axis='x',rotation=30)

#ax = plt.gca()
ax1.xaxis.grid(True,color ="lightgray",linestyle='--')
ax2.xaxis.grid(True,color ="lightgray",linestyle='--')
ax3.xaxis.grid(True,color ="lightgray",linestyle='--')
ax4.xaxis.grid(True,color ="lightgray",linestyle='--')
ax5.xaxis.grid(True,color ="lightgray",linestyle='--')
ax6.xaxis.grid(True,color ="lightgray",linestyle='--')

#ax1.yaxis.grid(True,color ="lightgray",linestyle='--')
#ax2.yaxis.grid(True,color ="lightgray",linestyle='--')
#ax3.yaxis.grid(True,color ="lightgray",linestyle='--')
#ax4.yaxis.grid(True,color ="lightgray",linestyle='--')
  

plt.subplots_adjust(hspace = 0.1)
  
fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(12,7+3) 

#fig.tight_layout() 

#fn='climate_anomanly_yearly_2001-2015_Europe_20211101'
outp=os.path.join(dir1,fn+'_SPEI01_GS_supp.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'_SPEI01_GS_supp.eps')
print(outp)
plt.savefig(outp,dpi=600,format='eps')    
    
   







############################################################################################################

# EUROCOM rectangular boundary extent  (yearly)


#xx1_tair=np.array([-0.02116617,  0.1137169 ,  0.01164828, -0.04627171, -0.06358433,
#        0.05410195,  0.18560422,  0.11902397, -0.03263106, -0.49034042,
#        0.07406014, -0.2362357 , -0.08294058,  0.24397639,  0.17103811])

#xx1_tair=np.array([-0.08198599,  0.37807811,  0.05422055, -0.07647025, -0.38576541,
#       -0.42633992,  0.52939742,  0.05791124,  0.06561335, -0.21545829,
#        0.13401174,  0.13919017, -0.57416393,  0.35579175,  0.04596946])
    

# Month 1-5 
xx2_tair= array([ 0.26606531,  1.43270246, -0.73049016, -0.23705002, -0.56182747,
       -1.50239939,  1.36679707,  1.25189256, -0.03422065, -1.54851968,
       -0.35680329, -0.68751203, -0.79563188,  1.49456599,  0.64243119])
 
xx2_p = array([ 0.99856495, -0.06477403, -2.43016903,  0.58102336,  0.04271971,
       -0.88511625,  0.8601859 ,  0.45301065, -0.24076718,  1.24814209,
       -0.01333889,  1.13171556, -0.80165033,  0.58461405, -1.46416056])
 
xx2_vpd =  array([-0.79076799,  0.50041514,  2.66158898, -1.45635728, -0.91335755,
        0.66133606, -0.08791349, -0.81920053,  0.29984657, -1.20807957,
        0.31324798,  0.61436585, -0.01705635, -0.55517817,  0.79711038])
      
xx2_smroot = array([ 0.45871871, -0.58883113, -1.53883719,  1.79536737, -0.21165325,
       -0.49271523, -0.41628817, -0.42231387,  0.04327251,  2.19750266,
        0.89588057, -0.99959753,  0.57809862, -0.20085532, -1.09774876])
    

#xx2_lai = array([ 1.07179732,  1.35075921, -2.68887033,  0.17840437,  0.83670065,
#       -1.28678984,  0.51212119,  0.00849627,  0.03041568, -0.78912846,
#        0.60184101, -0.86532462,  0.24035101,  0.26767841,  0.53154813])
    
xx2_csif = array([ 0.03392641,  1.39651349, -2.06323655,  0.90184543, -0.24735494,
       -1.05123745,  0.26939754,  0.03667437,  0.63583857, -0.10659838,
        1.15681193, -1.23439646,  0.05291679,  1.41974791, -1.20084866])

xx2_nirv = array([-0.74688166,  1.2189839 , -0.81963908,  1.18057679, -0.25250691,
       -0.84648718,  0.15654494, -0.65656586,  0.54431649, -0.76632705,
        1.26923613, -1.03574632,  1.06603988,  1.40988304, -1.72142712]) 
    

fluxsat = np.array([ 0.26003791,  1.07627851, -2.12021775,  1.036637  , -0.63649882,
        -0.80817583,  0.87476406,  0.66063278,  0.20005102, -0.5149928 ,
         0.68478335, -1.28156104, -0.39160552,  1.68690214, -0.727035  ])

gosifgpp = np.array([ 0.08042334,  1.1836332 , -1.74246687,  0.59750821, -0.42429958,
        -0.29585108,  0.76455718, -0.02735759,  0.53893993, -0.95846424,
         1.19887145, -1.78867297, -0.41402997,  1.76182978, -0.47462078])
  
   
#xx1_nirv= np.array([ 0.26813743,  0.68845201, -1.06370091,  0.82925323,  0.41391929,
#       -0.71100682,  1.16795062, -0.57954687, -0.54417644, -1.90307954,
#        0.74570927, -1.46607934,  0.10508905,  1.89759611,  0.15148292])

        
#xx1_spei=np.array([ 0.53005096,  0.14021018, -2.14352641,  0.70701526, -0.36682762,
#       -0.69343371,  0.67713118,  0.95364045,  0.62526444,  1.5540558 ,
#       -1.75015104,  0.93896238, -0.45616444,  0.09253059, -0.80875801])
 
#xx1_ewt=np.array([ 1.6626236 ,  0.15178945, -1.31274565,  0.07247854,  0.02049039,
#       -1.47222696, -0.64158079, -0.47098851, -0.63057627,  2.34160467,
#        1.00684471, -0.48628365,  0.51370506, -0.34010853, -0.41502607])
 

#pearsonr(xx1_csif[flag],xx1_nirv[flag])
#Out[118]: (0.8683556162272904, 2.6883967902347843e-05)
#
# 
#
#pearsonr(xx1_smroot[flag],xx1_spei[flag])
#Out[115]: (0.7093986722468171, 0.0030573348423114256)
#
#pearsonr(xx1_ewt[flag],xx1_smroot[flag])
#Out[117]: (0.7754710716927633, 0.0006811522893336032)
# 
#pearsonr(xx1_ewt[flag],xx1_spei[flag])
#Out[116]: (0.3328879756540325, 0.22536788986702416)
#
#
#
#
#pearsonr(xx1_csif[flag],xx1_smroot[flag])
#Out[122]: (0.2446478379877342, 0.37952024633836623)
#
#pearsonr(xx1_csif[flag],xx1_ewt[flag])
#Out[123]: (0.10210260520450358, 0.7172893179572029)
#
#pearsonr(xx1_csif[flag],xx1_spei[flag])
#Out[124]: (0.1653703039377622, 0.5558668577246484)
#
#pearsonr(xx1_nirv[flag],xx1_smroot[flag])
#Out[125]: (-0.07006766602380465, 0.8040320136081733)
#
#pearsonr(xx1_nirv[flag],xx1_ewt[flag])
#Out[126]: (-0.047866090378081275, 0.8654840971498852)
#
#pearsonr(xx1_nirv[flag],xx1_spei[flag])
#Out[127]: (-0.15266084846667738, 0.5870194717700081)
#
#

   
    


cte2019 = np.array([-0.06809985,  0.02320596, -0.16463827,  0.00152958,  0.0341887 ,
        -0.05612546,  0.15560942,  0.16821656,  0.17526091,
        
        -0.16851694,
         0.18520867, -0.12777879, -0.06998283, -0.0751083 , -0.01296937])
    
    
    
    
ct2019 = np.array([-0.07235407,  0.03719057,  0.10041716,  0.09296429,  0.01469185,
        -0.23994688, -0.05469818,  0.15312724, -0.09964366, 
        
        -0.09966622,
         0.29979735, -0.26371135,  0.22060409, -0.07290569, -0.01586649])
    
    
    
    
cams = np.array([-0.01607115, -0.18672789, -0.08769807,  0.11503489,  0.25514457,
         0.13926298, -0.09985728,  0.07037811,  0.08279048, 
         
         -0.20412863,
         0.00084564, -0.13624966, -0.03601567, -0.05785119,  0.16114287])
    
    
    
    
jena99v43 = np.array([ 0.01171644, -0.11077561, -0.1858765 ,  0.00692794,  0.24291837,
         0.25740639,  0.03339367,  0.10846295, -0.16906522,
         
         -0.15238051,
        -0.02598949,  0.01332   , -0.25433869,  0.00710164,  0.21717863])
    
    
jena99v43NEET = np.array([-0.03105824,  0.0433105 , -0.11816678,  0.09202262, -0.02168894,
        -0.08859283,  0.12654513,  0.07345238,  0.02283297, 
        
        -0.04720971,
         0.00296774, -0.02960854, -0.05265019,  0.04653433, -0.01869045])
    
    
    
    
jena99v43NEETW=  np.array([ 0.03464263,  0.03270925, -0.22071977,  0.10491675,  0.01843641,
        -0.04174402,  0.08800169,  0.06989118,  0.01314202, 
-0.13987459,
         0.04303823,  0.03091268, -0.07604809,  0.09014218, -0.04744656])
    
Trendynep = np.array([ 0.03773827,  0.01529954, -0.13602604,  0.10100564, -0.04313734,
       -0.08008579,  0.08330819,  0.0048757,  0.01180316,  0.09710952,
       -0.01785144, -0.12889138,  0.00811853,  0.13104703, -0.08431359])    
    

#Fluxcomnep = np.array([-0.01933919, -0.01947307, -0.02269379,  0.08139276,  0.02478938,
#        0.00166406, -0.00817153,  0.02632747,  0.04304065, -0.03074335,
#        0.00060811, -0.02253787, -0.03096929, -0.01049734, -0.01339699])

    
#ANN: JRA+ERA5
#Fluxcomnep = np.array([-0.00565191, -0.01010053, -0.0531066 ,  0.07561099,  0.02004914,
#        0.01300299, -0.02500568,  0.02027931,  0.03410706,  0.02639708,
#       -0.01025408, -0.02922055, -0.01406607, -0.00512667, -0.03691449])

    # JRA
Fluxcomnep1 = np.array([-0.00433007, -0.02503668, -0.06081315,  0.06735869,  0.01857612,
        0.00874353, -0.01574782,  0.0259082 ,  0.03947302,  0.04650378,
       -0.00464033, -0.03002585, -0.01186338, -0.01648363, -0.03762242])
    
    # ERA5
Fluxcomnep2 = np.array([-0.00697374,  0.00483561, -0.04540004,  0.0838633 ,  0.02152216,
        0.01726245, -0.03426353,  0.01465043,  0.0287411 ,  0.00629039,
       -0.01586783, -0.02841526, -0.01626876,  0.00623029, -0.03620656])    
  
    
#RF: JRA+ERA5
    
    # ERA5
#Fluxcomnep1 = np.array([-0.01856318, -0.03405157, -0.02088665,  0.04627107,  0.00646526,
#       -0.01669227,  0.02089448,  0.02116812,  0.04116855, -0.01636904,
#        0.01139911,  0.01259959, -0.03030944, -0.02020221, -0.00289182])
#    # JRA
#Fluxcomnep2=np.array([-0.03252524, -0.01035579,  0.01421517,  0.06013414,  0.00304454,
#       -0.00358337,  0.00275233,  0.0140961 ,  0.04558413, -0.04676741,
#       -0.00095394,  0.0052675 , -0.04221794, -0.00710127, -0.00158895])    
#    
#
#
##MARS: JRA+ERA5
#    
#    # ERA5
#Fluxcomnep1 = np.array([-0.02264459, -0.04462756, -0.02411996,  0.11452027,  0.05065146,
#       -0.0092859 , -0.00566562,  0.04916373,  0.05064405, -0.07014726,
#        0.00703241, -0.03951628, -0.04050604, -0.02203457,  0.00653587])
#
#    # JRA
#Fluxcomnep2=np.array([-0.03099831, -0.00760245,  0.00084189,  0.11620907,  0.04847675,
#        0.01353991, -0.01699899,  0.03297824,  0.05263304, -0.10397057,
#        0.00667927, -0.05513693, -0.04465021, -0.00339267, -0.00860804]) 
    


eurocom1 = array([-0.07796869,  0.16029364, -0.14973553,  0.17016346, -0.14245146,
        0.03969857])
eurocom2 =  array([ 0.03635836,  0.09542384, -0.16956765, -0.0642361 ,  0.07368795,
        0.02833359])
eurocom3 =  array([ 0.05698552,  0.02557602, -0.08046738, -0.3413053 ,  0.53678107,
       -0.19756992])
eurocom4 =  array([-0.10844791,  0.18667693, -0.15973426,  0.06867671,  0.1373812 ,
       -0.12455267])
eurocom5 =  array([ 0.13075221, -0.10266092, -0.12046228, -0.00315861,  0.12458666,
       -0.02905707])
eurocom6 =  array([ 0.23153015, -0.11316518, -0.13655968, -0.14684059, -0.00162981,
        0.16666511])
eurocom7 =  array([ 0.1477639 , -0.12292061, -0.17783597,  0.13337816,  0.01961452])

Byrne2020= array([ 0.07188867, -0.0792685 , -0.13504031,  0.12854501,  0.10566152,
       -0.09178641])
cmsflux = array([-0.04233668,  0.04380452, -0.13279528,  0.23493794, -0.03502472,
       -0.06858578])
gcasv2 = array([-0.01296051,  0.04737569, -0.17032468,  0.11297656,  0.16032072,
       -0.13738777])
ccdas=array([-0.12920947,  0.17686924, -0.11287834, -0.0060121 ,  0.2892296 ,
       -0.21799893]) 
    
    
    
    
X_arr=[]
X_arr.append(xx2_tair)
X_arr.append(xx2_p)
X_arr.append(xx2_vpd)
X_arr.append(xx2_smroot)
X_arr.append(xx2_nirv)
X_arr.append(xx2_csif)
X_arr.append(fluxsat)
X_arr.append(gosifgpp)
 
Y_arr=[]
Y_arr.append(cte2019)
Y_arr.append(ct2019)
Y_arr.append(cams)
Y_arr.append(jena99v43)
Y_arr.append(jena99v43NEET)
Y_arr.append(jena99v43NEETW)
Y_arr.append(Trendynep)
Y_arr.append(Fluxcomnep1)
Y_arr.append(Fluxcomnep2)

#Y_arr = np.array(Y_arr)
 
     
    

corr = []
pval = []

for x in X_arr:
    for y in Y_arr:
         ps = pearsonr(x,y)
         corr.append(ps[0])
         pval.append(ps[1])

    ps = pearsonr(x[-6:],eurocom1)
    corr.append(ps[0])
    pval.append(ps[1])
         
    ps = pearsonr(x[-6:],eurocom2)
    corr.append(ps[0])
    pval.append(ps[1])    
         
    ps = pearsonr(x[-6:],eurocom3)
    corr.append(ps[0])
    pval.append(ps[1])   
    
    ps = pearsonr(x[-6:],eurocom4)
    corr.append(ps[0])
    pval.append(ps[1])   
    
    ps = pearsonr(x[-6:],eurocom5)
    corr.append(ps[0])
    pval.append(ps[1])   
    
    ps = pearsonr(x[-6:],eurocom6)
    corr.append(ps[0])
    pval.append(ps[1])   
    
    ps = pearsonr(x[-5:],eurocom7)
    corr.append(ps[0])
    pval.append(ps[1])   
    
    
    ps = pearsonr(x[-6:],Byrne2020)
    corr.append(ps[0])
    pval.append(ps[1])  
    
    ps = pearsonr(x[-6:],gcasv2)
    corr.append(ps[0])
    pval.append(ps[1])  
    
    ps = pearsonr(x[-6:],cmsflux)
    corr.append(ps[0])
    pval.append(ps[1])      

    ps = pearsonr(x[-6:],ccdas)
    corr.append(ps[0])
    pval.append(ps[1])    
    
    
corr = np.array(corr)
pval = np.array(pval)    



corr_grid = corr.reshape(7+1,9+11)
pval_grid = pval.reshape(7+1,9+11)


corr_grid1 = np.round(corr_grid, 2) 
 

##################################################################   
         
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm  #, shiftgrid,maskoceans
 
fontsize = 32

#ax1=fig.add_subplot(4,1,1)


#plt.subplots(figsize=(15*2,12))

#ax = sns.heatmap(corr_grid, cmap = cm.GMT_no_green, annot = True, \
#                 square = True,linewidths=1, linecolor = 'k')
ax = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k'  , cbar_kws={'shrink': 0.6} )  # cmap = cm.GMT_no_green
 
#plt.setp(ax.get_legend().get_texts(), fontsize=fontsize*0.4) 


ax.set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard','JenaNEET','JenaNEETW','TRENDY_v6','FLUXCOM_CRUJRA','FLUXCOM_ERA5','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','GCASv2','CMS-Flux2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
         
ax.set_yticklabels( ['Ta','Prep','VPD','SM','NIRv','CSIF','FluxSat','GOSIFGPP'] , rotation=0, fontsize= fontsize*0.4, color='k')
         
#plt.legend(fontsize='50')


fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,7.5) 

fig.tight_layout() 

fn='Carbonflux_Environmentalvariables_correlation_anomanly_yearly_2001-2015_Europe_20211101_FLUXCOM_ANNOnly+Ta3_20221004'
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print(outp)
plt.savefig(outp,dpi=300,format='eps')    
    



#import seaborn as sns
#import numpy as np
#  
#  
#np.random.seed(0)
#  
#data = np.random.rand(12, 12)
#ax = sns.heatmap(data, cmap="PiYG")
#




############

plt.subplots(figsize=(15,10))
#ax = sns.heatmap(pval_grid, cmap = cm.GMT_no_green, annot = True, \
#                 square = True,linewidths=1, linecolor = 'k')
ax = sns.heatmap(pval_grid, fmt='.2f',vmin = 0, vmax = 1, center= 0.5, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k')  # cmap = cm.GMT_no_green


ax.set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard','JenaNEET','JenaNEETW','TRENDY_v6','FLUXCOM_CRUJRA','FLUXCOM_ERA5','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','GCASv2','CMS-Flux2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
         
ax.set_yticklabels( ['Ta','Prep','VPD','SM','NIRv','CSIF','FluxSat','GOSIFGPP'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='Carbonflux_Environmentalvariables_correlation_anomanly_yearly_2001-2015_Europe_20211101_pval_20221004'
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print(outp)
plt.savefig(outp,dpi=300,format='eps')    
        
   
                
''' 




################################################################
#
#   Use TBMs as reference
#
################################################################ 

# Table 2 & 3

corr_grid =[[-0.005,0.11,-0.42,-0.32,0.75,0.56,0.72,0.18,0.50],
            [0.35,-0.09,0.23,-0.34,0.68,0.63,0.58,0.63,0.74],
            [0.22,0.12,0.01,-0.05,0.73,0.80,0.89,0.77,0.58],
            [0.12,0.21,0.14,0.19,0.50,0.23,0.19,-0.04,0.61],
            [0.04,0.13,0.01,-0.03,0.64,0.45,0.57,0.41,0.49]
            ]

 

plt.subplots(figsize=(15,10))
ax = sns.heatmap(corr_grid, cmap = cm.GMT_no_green, annot = True, \
                 square = True,linewidths=1, linecolor = 'k')


ax.set_xticklabels( ['CTE2018','CT2019','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','FLUXCOM'] , rotation=90, fontsize= fontsize*0.4, color='k')
         
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='Inversions+Fluxcom_TRENDY_correlation_anomanly_yearly_2001-2015_Europe_20211101'
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print(outp)
plt.savefig(outp,dpi=600,format='eps')    
        




corr_grid1 =[[-0.33,0.33,0.44,0.07,-0.26,-0.06,-0.04,0.37,-0.15,0.19],
            [0.51,0.11,0.59,0.63,0.56,0.65,0.08,0.28,0.43,0.67],
            [-0.11,0.59,0.71,0.38,0.23,-0.40,0.68,0.51,0.71,0.77],
            [-0.35,-0.34,0.13,0.72,0.37,0.07,0.63,0.46,0.81,0.71],
            [-0.18,0.50,0.54,0.41,0.57,0.11,0.56,0.94,0.85,0.62] ]

plt.subplots(figsize=(15,10))
ax = sns.heatmap(corr_grid1, cmap = cm.GMT_no_green, annot = True, \
                 square = True,linewidths=1, linecolor = 'k')


ax.set_xticklabels( ['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
     
     
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='EUROCOM_TRENDY_correlation_anomanly_yearly_2001-2015_Europe_20211101'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        




fontsize = 32*2

corr_grid2 =[[-0.30,-0.44,0.71,0.35,0.36,0.59,-0.20,-0.56, 0,0.31],           
            [0.10,-0.24,0.60,0.58,0.78,0.54,0.32,-0.40,0,0.21],
            [-0.11,0.17,0.68,0.64,0.78,0.71,0.95,0.08,0,0.22],
            [-0.04,-0.39,0.94,0.91,0.84,0.96,0.98,0.37,0,-0.29],
            [-0.07,0.21,0.67,0.57,0.74,0.64,0.90,-0.15,0,0.56] ]

plt.subplots(figsize=(15,10))
#ax = sns.heatmap(  np.array(corr_grid2), cmap = cm.GMT_no_green, annot = True, \
#                 square = True,linewidths=1, linecolor = 'k')
ax = sns.heatmap(corr_grid2, vmin = -0.6, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k')   #cmap = cm.GMT_no_green


ax.set_xticklabels( ['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
     
     
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='EUROCOM_TRENDY_correlation_anomanly_yearly_2001-2015_Europe_20211101_prior'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        





plt.subplots(figsize=(15,10))
#ax = sns.heatmap(  np.array(corr_grid1)-np.array(corr_grid2), cmap = cm.GMT_no_green, annot = True, \
#                 square = True,linewidths=1, linecolor = 'k')

ax = sns.heatmap(np.array(corr_grid1)-np.array(corr_grid2), vmin = -1, vmax = 1, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k')   #cmap = cm.GMT_no_green

ax.set_xticklabels( ['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
     
     
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='EUROCOM_TRENDY_correlation_anomanly_yearly_2001-2015_Europe_20211101_diff'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        





################################################################
#
#   Use CSIF as reference
#
################################################################ 

# Table 2 & 3

corr_grid =[[0.40,0.15,0.32,0.06,-0.13,-0.03,0.01,0.19],
            [-0.009,-0.12,0.22,-0.39,0.71,0.71,0.81,0.76],
            [0.26,-0.01,-0.27,-0.29,0.57,0.68,0.88,0.69],
            [0.27,0.20,0.06,0.05,0.62,0.27,0.29,0.17],
            [0.25,0.18,0.10,-0.22,0.70,0.57,0.72,0.62]]

 

plt.subplots(figsize=(15,10))
ax = sns.heatmap(corr_grid, cmap = cm.GMT_no_green, annot = True, \
                 square = True,linewidths=1, linecolor = 'k')


ax.set_xticklabels( ['CTE2018','CT2019','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)'] , rotation=90, fontsize= fontsize*0.4, color='k')
         
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='Inversions+Fluxcom_CSIF_correlation_anomanly_yearly_2001-2015_Europe_20211101'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        



         
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm  #, shiftgrid,maskoceans
 

#corr_grid1 =[[0.23,0.08,0.08,0.62,0.26,-0.15,-0.46,0.21,0.48,0.17],
#            [0.53,-0.32,0.39,0.86,0.22,0.72,-0.61,-0.13,0.49,0.57],
#            [0.29,0.02,0.47,0.61,0.22,-0.34,0.52,0.27,0.47,0.66],
#            [-0.24,-0.24,0.08,0.69,0.28,-0.005,0.64,0.48,0.86,0.79],
#            [0.28,0.02,0.30,0.70,0.18,0.01,0.67,0.45,0.64,0.71] ]

corr_grid1 =[[0.29118436,0.09857707,0.10323954,0.77940415,0.32415406,-0.19407321,-0.45924322 , 0.24260029 ,0.54484755, -0.19656352],
            [5.48581556e-01, -3.32560584e-01,4.00345416e-01,8.85251148e-01,2.28697245e-01,7.41152155e-01,-6.50862390e-01,-1.45731262e-01,5.38304181e-01,-6.24386269e-01],
            [0.33140664,0.01901285,0.53038978,0.68517042,0.25051028,-0.38182493,0.55611734,0.31148717,0.53684352,-0.75814106],
            [-0.27376292,-0.2759206,0.08812699,0.80403329, 0.32844311,-0.0055882, 0.69526247, 0.49819787, 0.89799924,-0.82080986],
            [0.33822127, 0.02519285, 0.35514664, 0.83353934, 0.21417922, 0.00740873, 0.73753504, 0.52250905, 0.75398144,-0.83399207]]


#corr_grid1 =[[0.23,0.08,0.08,0.62,0.26,-0.15,-0.46,0.21,0.48,0.17],
#            [0.53,-0.32,0.39,0.86,0.22,0.72,-0.61,-0.13,0.49,0.57],
#            [0.29,0.02,0.47,0.61,0.22,-0.34,0.52,0.27,0.47,0.66],
#            [-0.24,-0.24,0.08,0.69,0.28,-0.005,0.64,0.48,0.86,0.79],
#            [0.28,0.02,0.30,0.70,0.18,0.01,0.67,0.45,0.64,0.71] ]



plt.subplots(figsize=(15,10))
ax = sns.heatmap(corr_grid1, cmap = cm.GMT_no_green, annot = True, \
                 square = True,linewidths=1, linecolor = 'k')


ax.set_xticklabels( ['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
     
     
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='EUROCOM_CSIF_correlation_anomanly_yearly_2001-2015_Europe_20211101_20220131'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        





corr_grid2 =[[0.28,0.17,0.17,0.59,-0.05,0.28,0.15,0.21,0,0.27],           
            [-0.06,-0.50,0.70,0.67,0.52,0.65,0.69,-0.13,0,-0.49],
            [-0.04,0.26,0.61,0.50,0.47,0.59,0.80,0.27,0,-0.50],
            [0.04,-0.25,0.78,0.80,0.82,0.80,0.81,0.48,0,0.49],
            [0.15,0.17,0.60,0.67,0.70,0.59,0.80,0.45,0,-0.15] ]

plt.subplots(figsize=(15,10))
#ax = sns.heatmap(  np.array(corr_grid2), cmap = cm.GMT_no_green, annot = True, \
#                 square = True,linewidths=1, linecolor = 'k')
ax = sns.heatmap(corr_grid2, vmin = -0.6, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k')   #cmap = cm.GMT_no_green



ax.set_xticklabels( ['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
     
     
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='EUROCOM_CSIF_correlation_anomanly_yearly_2001-2015_Europe_20211101_prior'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        







 

plt.subplots(figsize=(15,10))
#ax = sns.heatmap(  np.array(corr_grid1)-np.array(corr_grid2), cmap = cm.GMT_no_green, annot = True, \
#                 square = True,linewidths=1, linecolor = 'k')

ax = sns.heatmap(np.array(corr_grid1)-np.array(corr_grid2), vmin = -1, vmax = 1, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k')   #cmap = cm.GMT_no_green

ax.set_xticklabels( ['CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
     
     
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15,10) 

fig.tight_layout() 

fn='EUROCOM_CSIF_correlation_anomanly_yearly_2001-2015_Europe_20211101_diff_20220131'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        


 



'''
 

#RR_arr = np.array([[[ 0.05, -0.03,  0.5 , -0.09, -0.34, -0.04, -0.21,  0.17,  0.25,
#         -0.24, -0.2 ,  0.61, -0.15, -0.15,  0.06, -0.21, -0.21,  0.16,
#          0.02],
#        [-0.45,  0.25, -0.12,  0.18,  0.19,  0.21,  0.28,  0.13,  0.27,
#         -0.56,  0.12,  0.66, -0.04,  0.36, -0.68, -0.17, -0.23,  0.56,
#         -0.42],
#        [ 0.26,  0.27, -0.1 , -0.24,  0.32,  0.32,  0.38,  0.4 ,  0.6 ,
#         -0.53,  0.09,  0.8 , -0.06, -0.3 ,  0.19, -0.1 ,  0.08,  0.62,
#         -0.38],
#        [ 0.33, -0.08, -0.37, -0.34,  0.08, -0.11, -0.03,  0.24, -0.  ,
#          0.31,  0.38,  0.07,  0.34,  0.35,  0.11,  0.1 , -0.33,  0.17,
#         -0.39],
#        [ 0.28,  0.33,  0.18, -0.13,  0.22,  0.2 ,  0.24,  0.32,  0.67,
#         -0.42,  0.06,  0.76, -0.08,  0.12,  0.86, -0.05,  0.03,  0.62,
#         -0.19]],
#
#       [[-0.53, -0.71, -0.35, -0.43,  0.3 ,  0.61,  0.35,  0.17, -0.59,
#         -0.69, -0.61,  0.16, -0.45, -0.8 ,  0.84, -0.8 , -0.82, -0.66,
#         -0.67],
#        [-0.21, -0.39, -0.04, -0.51,  0.26,  0.23,  0.47,  0.45, -0.15,
#         -0.15,  0.05,  0.18,  0.38, -0.  ,  0.28, -0.41,  0.08, -0.23,
#         -0.34],
#        [-0.15, -0.45, -0.34, -0.24,  0.35,  0.38,  0.43,  0.21, -0.73,
#          0.53,  0.25, -0.41,  0.41, -0.28,  0.52,  0.45, -0.27, -0.6 ,
#          0.04],
#        [-0.42, -0.37,  0.03,  0.35,  0.24,  0.01,  0.12, -0.09, -0.76,
#         -0.63,  0.31,  0.26,  0.38,  0.36, -0.11,  0.9 ,  0.13, -0.18,
#          0.16],
#        [-0.1 , -0.32, -0.07, -0.25,  0.57,  0.38,  0.61,  0.38, -0.6 ,
#          0.31,  0.42, -0.1 ,  0.74,  0.2 , -0.38,  0.6 , -0.14, -0.52,
#         -0.18]],
#
#       [[ 0.56,  0.37,  0.53,  0.2 , -0.55, -0.56, -0.48, -0.15,  0.65,
#          0.25,  0.23,  0.5 ,  0.25,  0.18, -0.46,  0.31,  0.22,  0.49,
#          0.3 ],
#        [-0.38,  0.39, -0.46,  0.61, -0.55, -0.49, -0.64, -0.66, -0.29,
#         -0.25, -0.19, -0.27, -0.25, -0.42, -0.63,  0.04, -0.24,  0.64,
#         -0.15],
#        [ 0.19,  0.22,  0.07, -0.02, -0.14, -0.21, -0.41, -0.1 ,  0.64,
#         -0.66, -0.32,  0.31, -0.09,  0.36, -0.89, -0.7 , -0.16,  0.3 ,
#         -0.17],
#        [-0.05,  0.12, -0.5 , -0.29, -0.6 , -0.51, -0.45, -0.13,  0.48,
#          0.48, -0.14, -0.66, -0.36, -0.13, -0.59, -0.8 , -0.26, -0.36,
#         -0.37],
#        [ 0.01,  0.19, -0.34,  0.1 , -0.36, -0.29, -0.45, -0.19,  0.48,
#         -0.58, -0.14,  0.37, -0.4 ,  0.02,  0.44, -0.48,  0.05,  0.41,
#         -0.14]],
#
#       [[-0.39, -0.33, -0.42, -0.33,  0.32,  0.44,  0.4 , -0.1 , -0.62,
#         -0.19, -0.29, -0.22, -0.11, -0.67,  0.49, -0.46, -0.49, -0.5 ,
#         -0.47],
#        [ 0.16, -0.37,  0.14, -0.6 ,  0.37,  0.36,  0.43,  0.5 ,  0.21,
#          0.11,  0.31,  0.31,  0.61,  0.33,  0.4 , -0.04,  0.15, -0.52,
#         -0.05],
#        [ 0.04, -0.27, -0.15, -0.15,  0.24,  0.33,  0.32,  0.13, -0.64,
#          0.77,  0.29, -0.26,  0.02, -0.34,  0.73,  0.58,  0.21, -0.19,
#          0.16],
#        [-0.06, -0.04,  0.07,  0.55,  0.08, -0.18, -0.25, -0.48, -0.37,
#         -0.4 ,  0.15,  0.62,  0.56,  0.32,  0.69,  0.6 ,  0.26,  0.36,
#          0.56],
#        [-0.08, -0.12, -0.04, -0.06,  0.41,  0.26,  0.36,  0.16, -0.37,
#          0.58,  0.23, -0.02,  0.38, -0.1 ,  0.43,  0.72,  0.13, -0.13,
#          0.2 ]],
#
#       [[ 0.19,  0.  ,  0.33,  0.21, -0.2 ,  0.01, -0.04,  0.34,  0.18,
#         -0.3 , -0.15,  0.69,  0.28, -0.2 ,  0.02, -0.31, -0.18,  0.3 ,
#          0.19],
#        [-0.08, -0.06,  0.23, -0.36,  0.6 ,  0.59,  0.73,  0.66,  0.57,
#         -0.54,  0.23,  0.86,  0.12,  0.69, -0.61, -0.28, -0.3 ,  0.45,
#         -0.53],
#        [ 0.27,  0.05, -0.26, -0.38,  0.4 ,  0.53,  0.74,  0.58,  0.27,
#         -0.27,  0.29,  0.74,  0.07, -0.4 ,  0.6 ,  0.27,  0.19,  0.35,
#         -0.26],
#        [ 0.17,  0.05, -0.01,  0.07,  0.57,  0.2 ,  0.27,  0.11, -0.28,
#         -0.39,  0.07,  0.85,  0.41,  0.02,  0.59,  0.69,  0.48,  0.42,
#          0.51],
#        [ 0.15,  0.31,  0.14, -0.17,  0.51,  0.44,  0.58,  0.5 ,  0.47,
#         -0.26,  0.16,  0.91,  0.09, -0.  ,  0.74,  0.32,  0.2 ,  0.54,
#         -0.09]],
#
#       [[ 0.4 ,  0.34,  0.32,  0.06, -0.13, -0.04,  0.01,  0.19,  0.29,
#          0.1 ,  0.1 ,  0.78,  0.32, -0.19, -0.46,  0.14,  0.11,  0.51,
#          0.3 ],
#        [-0.01, -0.02,  0.22, -0.39,  0.71,  0.71,  0.81,  0.76,  0.55,
#         -0.33,  0.4 ,  0.89,  0.23,  0.74, -0.65, -0.23, -0.32,  0.5 ,
#         -0.54],
#        [ 0.26,  0.03, -0.27, -0.29,  0.57,  0.68,  0.88,  0.69,  0.33,
#          0.02,  0.53,  0.69,  0.25, -0.38,  0.56,  0.23, -0.01,  0.57,
#         -0.57],
#        [ 0.27,  0.14,  0.06,  0.05,  0.62,  0.27,  0.3 ,  0.18, -0.27,
#         -0.28,  0.09,  0.8 ,  0.33, -0.01,  0.7 ,  0.65,  0.45,  0.53,
#          0.51],
#        [ 0.25,  0.25,  0.1 , -0.22,  0.7 ,  0.57,  0.72,  0.62,  0.34,
#          0.03,  0.36,  0.83,  0.21,  0.01,  0.74,  0.49,  0.11,  0.61,
#         -0.18]],
#
#       [[ 0.39,  0.13,  0.01, -0.23,  0.17,  0.26,  0.37,  0.33,  0.03,
#          0.09,  0.1 ,  0.76,  0.18, -0.46, -0.36,  0.13, -0.05,  0.33,
#          0.07],
#        [-0.09,  0.14,  0.22, -0.3 ,  0.68,  0.66,  0.76,  0.69,  0.55,
#         -0.36,  0.26,  0.91,  0.14,  0.71, -0.54, -0.12, -0.23,  0.46,
#         -0.44],
#        [ 0.32,  0.12, -0.23, -0.2 ,  0.69,  0.76,  0.9 ,  0.76,  0.3 ,
#          0.1 ,  0.56,  0.67,  0.25, -0.4 ,  0.56,  0.22, -0.02,  0.63,
#         -0.61],
#        [ 0.22,  0.28, -0.1 ,  0.06,  0.49,  0.07,  0.11,  0.04, -0.19,
#         -0.24,  0.05,  0.84,  0.32, -0.06,  0.64,  0.6 ,  0.56,  0.48,
#          0.56],
#        [ 0.24,  0.32,  0.05, -0.18,  0.71,  0.57,  0.69,  0.62,  0.31,
#          0.07,  0.35,  0.83,  0.21, -0.04,  0.75,  0.56,  0.15,  0.61,
#         -0.15]],
#
#       [[ 0.64,  0.44,  0.49,  0.06, -0.28, -0.34, -0.21,  0.01,  0.45,
#          0.13,  0.15,  0.67,  0.56,  0.07, -0.66,  0.35,  0.37,  0.73,
#          0.55],
#        [-0.05,  0.05,  0.24, -0.23,  0.78,  0.76,  0.83,  0.76,  0.53,
#         -0.2 ,  0.55,  0.87,  0.26,  0.73, -0.71, -0.06, -0.4 ,  0.6 ,
#         -0.52],
#        [ 0.3 , -0.01, -0.18, -0.08,  0.61,  0.7 ,  0.94,  0.73,  0.35,
#          0.19,  0.73,  0.69,  0.33, -0.28,  0.57,  0.23,  0.  ,  0.67,
#         -0.62],
#        [ 0.4 ,  0.21,  0.12,  0.1 ,  0.65,  0.35,  0.32,  0.22, -0.2 ,
#         -0.1 ,  0.19,  0.8 ,  0.27, -0.02,  0.83,  0.52,  0.41,  0.71,
#          0.45],
#        [ 0.32,  0.23,  0.22, -0.07,  0.68,  0.53,  0.7 ,  0.59,  0.32,
#          0.11,  0.48,  0.82,  0.19,  0.08,  0.84,  0.47,  0.09,  0.75,
#         -0.19]]])




"""
RR_arr = np.array([[[ 0.05, -0.03,  0.5 , -0.09, -0.34, -0.04, -0.21,  0.17,  0.25,
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
"""


RR_arr = np.array([[[ 0.05, -0.03,  0.5 , -0.09, -0.34, -0.04, -0.21,  0.17,  0.25,
         -0.24, -0.2 ,  0.61, -0.15, -0.15,  0.06, -0.21, -0.21,  0.16,
         -0.14],
        [-0.45,  0.25, -0.12,  0.18,  0.19,  0.21,  0.28,  0.13,  0.27,
         -0.56,  0.12,  0.66, -0.04,  0.36, -0.68, -0.17, -0.23,  0.66,
          0.62],
        [ 0.26,  0.27, -0.1 , -0.24,  0.32,  0.32,  0.38,  0.4 ,  0.6 ,
         -0.53,  0.09,  0.8 , -0.06, -0.3 ,  0.19, -0.1 ,  0.08,  0.03,
          0.51],
        [ 0.33, -0.08, -0.37, -0.34,  0.08, -0.11, -0.03,  0.24, -0.  ,
          0.31,  0.38,  0.07,  0.34,  0.35,  0.11,  0.1 , -0.33,  0.06,
          0.15],
        [ 0.28,  0.33,  0.18, -0.13,  0.22,  0.2 ,  0.24,  0.32,  0.67,
         -0.42,  0.06,  0.76, -0.08,  0.12,  0.86, -0.05,  0.03,  0.28,
          0.46]],

       [[-0.53, -0.71, -0.35, -0.43,  0.3 ,  0.61,  0.35,  0.17, -0.59,
         -0.69, -0.61,  0.16, -0.45, -0.8 ,  0.84, -0.8 , -0.82, -0.57,
         -0.61],
        [-0.21, -0.39, -0.04, -0.51,  0.26,  0.23,  0.47,  0.45, -0.15,
         -0.15,  0.05,  0.18,  0.38, -0.  ,  0.28, -0.41,  0.08,  0.06,
          0.39],
        [-0.15, -0.45, -0.34, -0.24,  0.35,  0.38,  0.43,  0.21, -0.73,
          0.53,  0.25, -0.41,  0.41, -0.28,  0.52,  0.45, -0.27,  0.61,
         -0.21],
        [-0.42, -0.37,  0.03,  0.35,  0.24,  0.01,  0.12, -0.09, -0.76,
         -0.63,  0.31,  0.26,  0.38,  0.36, -0.11,  0.9 ,  0.13,  0.51,
          0.53],
        [-0.1 , -0.32, -0.07, -0.25,  0.57,  0.38,  0.61,  0.38, -0.6 ,
          0.31,  0.42, -0.1 ,  0.74,  0.2 , -0.38,  0.6 , -0.14,  0.34,
          0.27]],

       [[ 0.56,  0.37,  0.53,  0.2 , -0.55, -0.56, -0.48, -0.15,  0.65,
          0.25,  0.23,  0.5 ,  0.25,  0.18, -0.46,  0.31,  0.22,  0.51,
          0.21],
        [-0.38,  0.39, -0.46,  0.61, -0.55, -0.49, -0.64, -0.66, -0.29,
         -0.25, -0.19, -0.27, -0.25, -0.42, -0.63,  0.04, -0.24,  0.43,
          0.2 ],
        [ 0.19,  0.22,  0.07, -0.02, -0.14, -0.21, -0.41, -0.1 ,  0.64,
         -0.66, -0.32,  0.31, -0.09,  0.36, -0.89, -0.7 , -0.16, -0.76,
         -0.05],
        [-0.05,  0.12, -0.5 , -0.29, -0.6 , -0.51, -0.45, -0.13,  0.48,
          0.48, -0.14, -0.66, -0.36, -0.13, -0.59, -0.8 , -0.26, -0.8 ,
         -0.77],
        [ 0.01,  0.19, -0.34,  0.1 , -0.36, -0.29, -0.45, -0.19,  0.48,
         -0.58, -0.14,  0.37, -0.4 ,  0.02,  0.44, -0.48,  0.05, -0.13,
          0.12]],

       [[-0.39, -0.33, -0.42, -0.33,  0.32,  0.44,  0.4 , -0.1 , -0.62,
         -0.19, -0.29, -0.22, -0.11, -0.67,  0.49, -0.46, -0.49, -0.39,
         -0.3 ],
        [ 0.16, -0.37,  0.14, -0.6 ,  0.37,  0.36,  0.43,  0.5 ,  0.21,
          0.11,  0.31,  0.31,  0.61,  0.33,  0.4 , -0.04,  0.15, -0.07,
          0.22],
        [ 0.04, -0.27, -0.15, -0.15,  0.24,  0.33,  0.32,  0.13, -0.64,
          0.77,  0.29, -0.26,  0.02, -0.34,  0.73,  0.58,  0.21,  0.61,
          0.13],
        [-0.06, -0.04,  0.07,  0.55,  0.08, -0.18, -0.25, -0.48, -0.37,
         -0.4 ,  0.15,  0.62,  0.56,  0.32,  0.69,  0.6 ,  0.26,  0.61,
          0.54],
        [-0.08, -0.12, -0.04, -0.06,  0.41,  0.26,  0.36,  0.16, -0.37,
          0.58,  0.23, -0.02,  0.38, -0.1 ,  0.43,  0.72,  0.13,  0.5 ,
          0.22]],

       [[ 0.19,  0.  ,  0.33,  0.21, -0.2 ,  0.01, -0.04,  0.34,  0.18,
         -0.3 , -0.15,  0.69,  0.28, -0.2 ,  0.02, -0.31, -0.18,  0.22,
         -0.15],
        [-0.08, -0.06,  0.23, -0.36,  0.6 ,  0.59,  0.73,  0.66,  0.57,
         -0.54,  0.23,  0.86,  0.12,  0.69, -0.61, -0.28, -0.3 ,  0.43,
          0.53],
        [ 0.27,  0.05, -0.26, -0.38,  0.4 ,  0.53,  0.74,  0.58,  0.27,
         -0.27,  0.29,  0.74,  0.07, -0.4 ,  0.6 ,  0.27,  0.19,  0.42,
          0.7 ],
        [ 0.17,  0.05, -0.01,  0.07,  0.57,  0.2 ,  0.27,  0.11, -0.28,
         -0.39,  0.07,  0.85,  0.41,  0.02,  0.59,  0.69,  0.48,  0.88,
          0.81],
        [ 0.15,  0.31,  0.14, -0.17,  0.51,  0.44,  0.58,  0.5 ,  0.47,
         -0.26,  0.16,  0.91,  0.09, -0.  ,  0.74,  0.32,  0.2 ,  0.6 ,
          0.68]],

       [[ 0.4 ,  0.34,  0.32,  0.06, -0.13, -0.04,  0.01,  0.19,  0.29,
          0.1 ,  0.1 ,  0.78,  0.32, -0.19, -0.46,  0.14,  0.11,  0.54,
          0.2 ],
        [-0.01, -0.02,  0.22, -0.39,  0.71,  0.71,  0.81,  0.76,  0.55,
         -0.33,  0.4 ,  0.89,  0.23,  0.74, -0.65, -0.23, -0.32,  0.54,
          0.62],
        [ 0.26,  0.03, -0.27, -0.29,  0.57,  0.68,  0.88,  0.69,  0.33,
          0.02,  0.53,  0.69,  0.25, -0.38,  0.56,  0.23, -0.01,  0.54,
          0.76],
        [ 0.27,  0.14,  0.06,  0.05,  0.62,  0.27,  0.3 ,  0.18, -0.27,
         -0.28,  0.09,  0.8 ,  0.33, -0.01,  0.7 ,  0.65,  0.45,  0.9 ,
          0.82],
        [ 0.25,  0.25,  0.1 , -0.22,  0.7 ,  0.57,  0.72,  0.62,  0.34,
          0.03,  0.36,  0.83,  0.21,  0.01,  0.74,  0.49,  0.11,  0.75,
          0.83]],

       [[ 0.39,  0.13,  0.01, -0.23,  0.17,  0.26,  0.37,  0.33,  0.03,
          0.09,  0.1 ,  0.76,  0.18, -0.46, -0.36,  0.13, -0.05,  0.51,
          0.19],
        [-0.09,  0.14,  0.22, -0.3 ,  0.68,  0.66,  0.76,  0.69,  0.55,
         -0.36,  0.26,  0.91,  0.14,  0.71, -0.54, -0.12, -0.23,  0.71,
          0.75],
        [ 0.32,  0.12, -0.23, -0.2 ,  0.69,  0.76,  0.9 ,  0.76,  0.3 ,
          0.1 ,  0.56,  0.67,  0.25, -0.4 ,  0.56,  0.22, -0.02,  0.55,
          0.78],
        [ 0.22,  0.28, -0.1 ,  0.06,  0.49,  0.07,  0.11,  0.04, -0.19,
         -0.24,  0.05,  0.84,  0.32, -0.06,  0.64,  0.6 ,  0.56,  0.89,
          0.8 ],
        [ 0.24,  0.32,  0.05, -0.18,  0.71,  0.57,  0.69,  0.62,  0.31,
          0.07,  0.35,  0.83,  0.21, -0.04,  0.75,  0.56,  0.15,  0.82,
          0.87]],

       [[ 0.64,  0.44,  0.49,  0.06, -0.28, -0.34, -0.21,  0.01,  0.45,
          0.13,  0.15,  0.67,  0.56,  0.07, -0.66,  0.35,  0.37,  0.73,
          0.42],
        [-0.05,  0.05,  0.24, -0.23,  0.78,  0.76,  0.83,  0.76,  0.53,
         -0.2 ,  0.55,  0.87,  0.26,  0.73, -0.71, -0.06, -0.4 ,  0.69,
          0.72],
        [ 0.3 , -0.01, -0.18, -0.08,  0.61,  0.7 ,  0.94,  0.73,  0.35,
          0.19,  0.73,  0.69,  0.33, -0.28,  0.57,  0.23,  0.  ,  0.56,
          0.81],
        [ 0.4 ,  0.21,  0.12,  0.1 ,  0.65,  0.35,  0.32,  0.22, -0.2 ,
         -0.1 ,  0.19,  0.8 ,  0.27, -0.02,  0.83,  0.52,  0.41,  0.92,
          0.85],
        [ 0.32,  0.23,  0.22, -0.07,  0.68,  0.53,  0.7 ,  0.59,  0.32,
          0.11,  0.48,  0.82,  0.19,  0.08,  0.84,  0.47,  0.09,  0.77,
          0.85]]])

##########################################################################################
##########################################################################################
## 
##         global + EUROCOM + satellite
##
##########################################################################################
###########################################################################################

 
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm  #, shiftgrid,maskoceans
 
#   Use CSIF as reference
 
fontsize= 16  #32*2

# Table 2 & 3

#corr_grid =[[0.40,0.15,0.32,0.06,-0.13,-0.03,0.01,0.19,0.23,0.08,0.08,0.62,0.26,-0.15,-0.46,0.21,0.48,0.17],
#            [-0.009,-0.12,0.22,-0.39,0.71,0.71,0.81,0.76,0.53,-0.32,0.39,0.86,0.22,0.72,-0.61,-0.13,0.49,0.57],
#            [0.26,-0.01,-0.27,-0.29,0.57,0.68,0.88,0.69,0.29,0.02,0.47,0.61,0.22,-0.34,0.52,0.27,0.47,0.66],
#            [0.27,0.20,0.06,0.05,0.62,0.27,0.29,0.17,-0.24,-0.24,0.08,0.69,0.28,-0.005,0.64,0.48,0.86,0.79],
#            [0.25,0.18,0.10,-0.22,0.70,0.57,0.72,0.62,0.28,0.02,0.30,0.70,0.18,0.01,0.67,0.45,0.64,0.71]]

names = np.array(['Ta','Prep','VPD','SM','NIRv','CSIF','FluxSat GPP','GOSIF GPP'])

fig, axs = plt.subplots(figsize = (6*2,2*3+0.5),nrows=4, ncols=2)


for i in np.arange(0,RR_arr[:4].shape[0]):
     
    corr_grid  = RR_arr[i]

    #plt.subplots(figsize=(15,15))
    #axs[i] = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
    #                 square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.6})   #cmap = cm.GMT_no_green
    
    if i == 2: 
        sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,0],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "bwr_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
      
    else:
        sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,0],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
 
     
    axs[i,0].set_xticklabels([])         
    axs[i,0].set_yticklabels( ['Northern','Western','Central','Southern','whole'] , rotation=0, fontsize= fontsize*0.5, color='k')    
    axs[i,0].set_title(names[i], fontsize=fontsize*0.55)
     
for i in np.arange(0,RR_arr[4:].shape[0]):
     
    corr_grid  = RR_arr[4+i]

    #plt.subplots(figsize=(15,15))
    #axs[i] = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
    #                 square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.6})   #cmap = cm.GMT_no_green
    
    #if i == 2: 
    sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,1],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
      
    
    #axs[i].set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','GCASv2','CMS-Flux2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
       
    axs[i,1].set_xticklabels([])         
    axs[i,1].set_yticklabels( ['Northern','Western','Central','Southern','whole'] , rotation=0, fontsize= fontsize*0.5, color='k')    
    axs[i,1].set_title(names[i+4], fontsize=fontsize*0.55)
     
    

axs[3,0].set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','CMS-Flux2020','GCASv2','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.6, color='k')      
axs[3,1].set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','CMS-Flux2020','GCASv2','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.6, color='k')
    
#fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
#fig.set_size_inches(15*2,10+5) 

fig.tight_layout() 

fn='Inversions_environmental_variable_correlation_anomanly_yearly_2001-2015_Europe_20221008_combined_updated'   #+str(i)
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print(outp)
plt.savefig(outp,dpi=300,format='eps')    
        



 

RR_arr_prior = np.array([[[ 0.04, -0.02,  0.05,  0.32,  0.08,  0.1 ,  0.24, 0, -0.3 ,
          0.16,  0.12],
        [-0.06, -0.37,  0.44,  0.4 , -0.06,  0.48,  0.68, 0,  0.2 ,
          0.56,  0.22],
        [ 0.28,  0.31,  0.26,  0.24,  0.03,  0.28,  0.18, 0,  0.38,
          0.62,  0.55],
        [-0.02,  0.23, -0.23, -0.16,  0.08, -0.26,  0.19,  0,  0.02,
          0.17,  0.08],
        [ 0.29,  0.17,  0.16,  0.39,  0.27,  0.16,  0.58, 0,  0.45,
          0.62,  0.61]],

       [[-0.69, -0.79,  0.33,  0.16, -0.07,  0.15, -0.32, 0, -0.89,
         -0.66, -0.41],
        [-0.54, -0.58,  0.41,  0.34,  0.37,  0.21,  0.26, 0, -0.76,
         -0.23,  0.18],
        [-0.64, -0.38,  0.21,  0.  ,  0.42,  0.2 ,  0.36,  0, -0.89,
         -0.6 , -0.73],
        [-0.63, -0.82,  0.81,  0.41,  0.49,  0.77,  0.75,  0, -0.91,
         -0.18, -0.36],
        [-0.62, -0.42,  0.6 ,  0.13,  0.62,  0.47,  0.41,  0 , -0.95,
         -0.52, -0.73]],

       [[ 0.41,  0.49, -0.28,  0.14, -0.07, -0.13,  0.37,  0,  0.27,
          0.49,  0.44],
        [ 0.11,  0.16, -0.41, -0.41, -0.72, -0.31,  0.13,  0,  0.71,
          0.64,  0.16],
        [ 0.54,  0.46, -0.4 , -0.18, -0.39, -0.39, -0.81, 0,  0.57,
          0.3 ,  0.5 ],
        [ 0.2 ,  0.57, -0.93, -0.79, -0.85, -0.93, -0.94, 0,  0.55,
         -0.36, -0.22],
        [ 0.43,  0.28, -0.4 , -0.  , -0.44, -0.35, -0.45, 0,  0.72,
          0.41,  0.76]],

       [[-0.55, -0.38,  0.44,  0.27, -0.12,  0.28, -0.48, 0, -0.53,
         -0.5 ,  0.24],
        [-0.21, -0.49,  0.42,  0.41,  0.54,  0.27,  0.13, 0, -0.9 ,
         -0.52,  0.2 ],
        [-0.33, -0.11,  0.43,  0.35,  0.57,  0.42,  0.83,  0 , -0.37,
         -0.19,  0.37],
        [-0.09, -0.52,  0.87,  0.67,  0.76,  0.88,  0.88,  0, -0.47,
          0.36, -0.27],
        [-0.22, -0.14,  0.65,  0.36,  0.61,  0.66,  0.98,  0, -0.65,
         -0.13,  0.58]],

       [[ 0.05, -0.12, -0.05,  0.55, -0.4 ,  0.1 , -0.26, 0, -0.43,
          0.3 , -0.08],
        [-0.19, -0.6 ,  0.67,  0.55,  0.35,  0.62,  0.48, 0, -0.05,
          0.45, -0.29],
        [ 0.1 ,  0.37,  0.54,  0.51,  0.24,  0.5 ,  0.57,  0 ,  0.11,
          0.35, -0.17],
        [-0.02, -0.29,  0.89,  0.89,  0.9 ,  0.91,  0.97,  0, -0.37,
          0.42, -0.34],
        [ 0.26,  0.3 ,  0.47,  0.68,  0.53,  0.5 ,  0.87,  0,  0.19,
          0.54, -0.43]],

       [[ 0.31,  0.27,  0.27,  0.76, -0.03,  0.4 ,  0.06,  0, -0.01,
          0.51,  0.49],
        [-0.06, -0.51,  0.72,  0.66,  0.49,  0.66,  0.61, 0, -0.02,
          0.5 ,  0.24],
        [ 0.06,  0.27,  0.67,  0.56,  0.51,  0.64,  0.67,  0,  0.  ,
          0.57,  0.35],
        [ 0.05, -0.27,  0.9 ,  0.92,  0.92,  0.91,  0.94,  0 , -0.29,
          0.53,  0.45],
        [ 0.18,  0.22,  0.7 ,  0.8 ,  0.78,  0.69,  0.9 ,  0,  0.02,
          0.61,  0.34]],

       [[ 0.25,  0.12,  0.48,  0.71,  0.13,  0.53,  0.03,  0, -0.07,
          0.33,  0.6 ],
        [-0.1 , -0.57,  0.81,  0.72,  0.45,  0.78,  0.78,  0, -0.04,
          0.46,  0.12],
        [ 0.06,  0.29,  0.69,  0.59,  0.6 ,  0.68,  0.72,  0,  0.01,
          0.63,  0.4 ],
        [ 0.09, -0.17,  0.87,  0.92,  0.91,  0.9 ,  0.97,  0, -0.24,
          0.48,  0.43],
        [ 0.17,  0.2 ,  0.76,  0.84,  0.83,  0.75,  0.92,  0, -0.02,
          0.61,  0.31]],

       [[ 0.59,  0.32, -0.16,  0.39, -0.16, -0.  , -0.08,  0,  0.2 ,
          0.73,  0.62],
        [ 0.17, -0.4 ,  0.71,  0.72,  0.52,  0.7 ,  0.68,  0,  0.07,
          0.6 ,  0.12],
        [ 0.16,  0.21,  0.75,  0.66,  0.62,  0.79,  0.77,  0,  0.04,
          0.67,  0.43],
        [ 0.27, -0.15,  0.8 ,  0.96,  0.88,  0.82,  0.85,  0, -0.07,
          0.71,  0.65],
        [ 0.36,  0.19,  0.71,  0.83,  0.76,  0.74,  0.87,  0,  0.1 ,
          0.75,  0.41]]])
    
 
###########################################################################################

 
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm  #, shiftgrid,maskoceans
 
#   Use CSIF as reference
 
fontsize= 16  #32*2

# Table 2 & 3

#corr_grid =[[0.40,0.15,0.32,0.06,-0.13,-0.03,0.01,0.19,0.23,0.08,0.08,0.62,0.26,-0.15,-0.46,0.21,0.48,0.17],
#            [-0.009,-0.12,0.22,-0.39,0.71,0.71,0.81,0.76,0.53,-0.32,0.39,0.86,0.22,0.72,-0.61,-0.13,0.49,0.57],
#            [0.26,-0.01,-0.27,-0.29,0.57,0.68,0.88,0.69,0.29,0.02,0.47,0.61,0.22,-0.34,0.52,0.27,0.47,0.66],
#            [0.27,0.20,0.06,0.05,0.62,0.27,0.29,0.17,-0.24,-0.24,0.08,0.69,0.28,-0.005,0.64,0.48,0.86,0.79],
#            [0.25,0.18,0.10,-0.22,0.70,0.57,0.72,0.62,0.28,0.02,0.30,0.70,0.18,0.01,0.67,0.45,0.64,0.71]]

names = np.array(['Ta','Prep','VPD','SM','NIRv','CSIF','FluxSat GPP','GOSIF GPP'])

fig, axs = plt.subplots(figsize = (6*2-4,2*3+0.5),nrows=4, ncols=2)
 

for i in np.arange(0,RR_arr_prior[:4,].shape[0]):
     
    corr_grid  = RR_arr_prior[i]

    #plt.subplots(figsize=(15,15))
    #axs[i] = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
    #                 square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.6})   #cmap = cm.GMT_no_green
    
    if i == 2: 
        sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,0],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "bwr_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
      
    else:
        sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,0],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
 
     
    axs[i,0].set_xticklabels([])         
    axs[i,0].set_yticklabels( ['Northern','Western','Central','Southern','whole'] , rotation=0, fontsize= fontsize*0.5, color='k')    
    axs[i,0].set_title(names[i], fontsize=fontsize*0.55)
     
for i in np.arange(0,RR_arr_prior[4:].shape[0]):
     
    corr_grid  = RR_arr_prior[4+i]

    #plt.subplots(figsize=(15,15))
    #axs[i] = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
    #                 square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.6})   #cmap = cm.GMT_no_green
    
    #if i == 2: 
    sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,1],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
      
    
    #axs[i].set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','GCASv2','CMS-Flux2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
       
    axs[i,1].set_xticklabels([])         
    axs[i,1].set_yticklabels( ['Northern','Western','Central','Southern','whole'] , rotation=0, fontsize= fontsize*0.5, color='k')    
    axs[i,1].set_title(names[i+4], fontsize=fontsize*0.55)
     
    

axs[3,0].set_xticklabels( [ 'CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','CMS-Flux2020','GCASv2','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.6, color='k')      
axs[3,1].set_xticklabels( [ 'CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','CMS-Flux2020','GCASv2','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.6, color='k')
    
#fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
#fig.set_size_inches(15*2,10+5) 

fig.tight_layout() 

fn='Inversions_environmental_variable_correlation_anomanly_yearly_2001-2015_Europe_20221008_combined_prior_updated'   #+str(i)
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print(outp)
plt.savefig(outp,dpi=300,format='eps')    
        




###########################################################################################

 
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm  #, shiftgrid,maskoceans
 
#   Use CSIF as reference
 
fontsize= 16  #32*2

# Table 2 & 3

#corr_grid =[[0.40,0.15,0.32,0.06,-0.13,-0.03,0.01,0.19,0.23,0.08,0.08,0.62,0.26,-0.15,-0.46,0.21,0.48,0.17],
#            [-0.009,-0.12,0.22,-0.39,0.71,0.71,0.81,0.76,0.53,-0.32,0.39,0.86,0.22,0.72,-0.61,-0.13,0.49,0.57],
#            [0.26,-0.01,-0.27,-0.29,0.57,0.68,0.88,0.69,0.29,0.02,0.47,0.61,0.22,-0.34,0.52,0.27,0.47,0.66],
#            [0.27,0.20,0.06,0.05,0.62,0.27,0.29,0.17,-0.24,-0.24,0.08,0.69,0.28,-0.005,0.64,0.48,0.86,0.79],
#            [0.25,0.18,0.10,-0.22,0.70,0.57,0.72,0.62,0.28,0.02,0.30,0.70,0.18,0.01,0.67,0.45,0.64,0.71]]

names = np.array(['Ta','Prep','VPD','SM','NIRv','CSIF','FluxSat GPP','GOSIF GPP'])

fig, axs = plt.subplots(figsize = (6*2-4,2*3+0.5),nrows=4, ncols=2)

RR_diff = RR_arr[:,:,8:] - RR_arr_prior

for i in np.arange(0,RR_diff[:4,].shape[0]):
     
    corr_grid  = RR_diff[i]

    #plt.subplots(figsize=(15,15))
    #axs[i] = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
    #                 square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.6})   #cmap = cm.GMT_no_green
    
    if i == 2: 
        sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,0],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "bwr_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
      
    else:
        sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,0],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
 
     
    axs[i,0].set_xticklabels([])         
    axs[i,0].set_yticklabels( ['Northern','Western','Central','Southern','whole'] , rotation=0, fontsize= fontsize*0.5, color='k')    
    axs[i,0].set_title(names[i], fontsize=fontsize*0.55)
     
for i in np.arange(0,RR_diff[4:].shape[0]):
     
    corr_grid  = RR_diff[4+i]

    #plt.subplots(figsize=(15,15))
    #axs[i] = sns.heatmap(corr_grid, fmt='.2f', vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
    #                 square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.6})   #cmap = cm.GMT_no_green
    
    #if i == 2: 
    sns.heatmap(corr_grid, fmt='.2f', ax = axs[i,1],  vmin = -0.9, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                         square = True,linewidths=1, linecolor = 'k',cbar_kws={'shrink': 0.7})   #cmap = cm.GMT_no_green
      
    
    #axs[i].set_xticklabels( ['CTE2018','CT2019B','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','GCASv2','CMS-Flux2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
       
    axs[i,1].set_xticklabels([])         
    axs[i,1].set_yticklabels( ['Northern','Western','Central','Southern','whole'] , rotation=0, fontsize= fontsize*0.5, color='k')    
    axs[i,1].set_title(names[i+4], fontsize=fontsize*0.55)
     
    

axs[3,0].set_xticklabels( [ 'CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','CMS-Flux2020','GCASv2','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.6, color='k')      
axs[3,1].set_xticklabels( [ 'CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','Byrne2020','CMS-Flux2020','GCASv2','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.6, color='k')
    
#fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
#fig.set_size_inches(15*2,10+5) 

fig.tight_layout() 

fn='Inversions_environmental_variable_correlation_anomanly_yearly_2001-2015_Europe_20221008_combined_post-prior_updated'   #+str(i)
outp=os.path.join(dir1,fn+'.png')
print(outp)
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print(outp)
plt.savefig(outp,dpi=300,format='eps')    
        



''' 
#   Use TBMs as reference
 
# Table 2 & 3

corr_grid =[[-0.005,0.11,-0.42,-0.32,0.75,0.56,0.72,0.18,-0.33,0.33,0.44,0.07,-0.26,-0.06,-0.04,0.37,-0.15,0.19],
            [0.35,-0.09,0.23,-0.34,0.68,0.63,0.58,0.63,0.51,0.11,0.59,0.63,0.56,0.65,0.08,0.28,0.43,0.67],
            [0.22,0.12,0.01,-0.05,0.73,0.80,0.89,0.77,-0.11,0.59,0.71,0.38,0.23,-0.40,0.68,0.51,0.71,0.77],
            [0.12,0.21,0.14,0.19,0.50,0.23,0.19,-0.04,-0.35,-0.34,0.13,0.72,0.37,0.07,0.63,0.46,0.81,0.71],
            [0.04,0.13,0.01,-0.03,0.64,0.45,0.57,0.41,-0.18,0.50,0.54,0.41,0.57,0.11,0.56,0.94,0.85,0.62]
            ]

 

plt.subplots(figsize=(15,10))
ax = sns.heatmap(corr_grid, vmin = -0.6, vmax = 0.9, center= 0, cmap= "RdBu_r", annot = True, annot_kws={"fontsize":fontsize*0.35}, \
                 square = True,linewidths=1, linecolor = 'k')  # cmap = cm.GMT_no_green


ax.set_xticklabels( ['CTE2018','CT2019','CAMS','JenaStandard(s99)','JenaNEET(sEXT)','JenaNEET(sEXT10)','JenaNEETW(s99)','JenaNEETW(sEXT10)','CTE(EUROCOM)','EnKF-RAMS','FLEXINVERT','CarboScope-Regional','LUMIA','PYVAR-CHIMERE','NAME-HB','GCASv2','Byrne2020','CCDAS-SM+VOD'] , rotation=90, fontsize= fontsize*0.4, color='k')
             
ax.set_yticklabels( ['R1','R2','R3','R4','All'] , rotation=0, fontsize= fontsize*0.4, color='k')
         

fig = plt.gcf()
#fig.set_size_inches(8,5)  #(8,10)
fig.set_size_inches(15*2,10) 

fig.tight_layout() 

fn='Inversions+EUROCOM+Satellite_TRENDY_correlation_anomanly_yearly_2001-2015_Europe_20211101'
outp=os.path.join(dir1,fn+'.png')
print outp
plt.savefig(outp,dpi=300)
   
    
outp=os.path.join(dir1,fn+'.eps')
print outp
plt.savefig(outp,dpi=600,format='eps')    
        






 


##########################################################################################
##########################################################################################

    
pearsonr(xx1_csif[flag],jena99v43NEETW)
pearsonr(xx1_csif[flag],jena99v43NEET) 
pearsonr(xx1_csif[flag],jena99v43) 
pearsonr(xx1_smroot[flag],jena99v43NEETW)
pearsonr(xx1_smroot[flag],jena99v43NEET) 
pearsonr(xx1_smroot[flag],jena99v43) 
pearsonr(xx1_spei[flag],jena99v43NEETW)
pearsonr(xx1_spei[flag],jena99v43NEET) 
pearsonr(xx1_spei[flag],jena99v43) 
pearsonr(xx1_csif[flag],cte2019)
pearsonr(xx1_csif[flag],ct2019)
pearsonr(xx1_csif[flag],cams)
pearsonr(xx1_smroot[flag],cte2019)
pearsonr(xx1_smroot[flag],ct2019)
pearsonr(xx1_smroot[flag],cams)



pearsonr(xx1_csif[flag],jena99v43NEETW)
Out[146]: (0.763337400074782, 0.0009289839032302101)

pearsonr(xx1_csif[flag],jena99v43NEET)
Out[151]: (0.7692178432022878, 0.0008010897087804667)  

pearsonr(xx1_csif[flag],jena99v43) 
Out[152]: (0.039315770673945614, 0.8893618264390387)

pearsonr(xx1_csif[flag],cte2019)
Out[161]: (0.5294173862244869, 0.04240483147645869)

pearsonr(xx1_csif[flag],ct2019)
Out[162]: (0.2617365749851211, 0.3460236914065041)

pearsonr(xx1_csif[flag],cams)
Out[163]: (-0.023496874624211128, 0.9337571946444031)

pearsonr(xx1_csif[flag],Trendynep)
Out[171]: (0.7068740540968441, 0.003212269363824864)

pearsonr(xx1_csif[flag],Fluxcomnep)
Out[181]: (0.24942639216873563, 0.36998232307702383)




pearsonr(gosifgpp[flag],jena99v43NEETW)
Out[193]: (0.6908325494002437, 0.004349135969238747)

pearsonr(gosifgpp[flag],jena99v43NEET)
Out[194]: (0.680691176981193, 0.0052185727408681866)

pearsonr(gosifgpp[flag],jena99v43)
Out[195]: (0.024397128364504087, 0.9312246173700226)

pearsonr(gosifgpp[flag],Trendynep)
Out[196]: (0.6767777909134626, 0.005588666935166362)



pearsonr(fluxsat[flag],jena99v43NEETW)
Out[197]: (0.7957166146985183, 0.0003884830355094404)

pearsonr(fluxsat[flag],jena99v43NEET)
Out[198]: (0.8278789464951721, 0.00013858998521802226)

pearsonr(fluxsat[flag],jena99v43)
Out[199]: (-0.0057676800806775, 0.9837241838112423)

pearsonr(fluxsat[flag],Trendynep)
Out[200]: (0.808234709934931, 0.00026601964899287554)

 



pearsonr(xx1_smroot[flag],jena99v43NEETW)
Out[155]: (0.18076142834132783, 0.5191245364079886)

pearsonr(xx1_smroot[flag],jena99v43NEET)
Out[156]: (0.27751749943955484, 0.31662016102456814)

pearsonr(xx1_smroot[flag],jena99v43)
Out[157]: (-0.3005102430708391, 0.27646039178739606)

pearsonr(xx1_smroot[flag],Trendynep)
Out[202]: (0.7760835370909652, 0.0006702372222639157)




pearsonr(xx1_spei[flag],fluxsat[flag])
Out[204]: (0.3618869744754752, 0.1850210679033608)

pearsonr(xx1_smroot[flag],fluxsat[flag])
Out[205]: (0.44690224497702696, 0.09490011151753375)

pearsonr(xx1_ewt[flag],fluxsat[flag])
Out[212]: (0.2177334634233776, 0.43565984279172654)
 

pearsonr(xx1_spei[flag],gosifgpp[flag])
Out[209]: (0.06826558284380463, 0.8089852343299259)

pearsonr(xx1_smroot[flag],gosifgpp[flag])
Out[210]: (0.2378859648242396, 0.3932410313310993)

pearsonr(xx1_ewt[flag],gosifgpp[flag])
Out[211]: (0.09807299738641333, 0.7280484960302882)






pearsonr(xx1_spei[flag],jena99v43NEETW)
Out[158]: (0.41392544372841306, 0.12507673019621002)

pearsonr(xx1_spei[flag],jena99v43NEET)
Out[159]: (0.49473284303795456, 0.06081304099413577)

pearsonr(xx1_spei[flag],jena99v43)
Out[160]: (-0.03975041132215345, 0.8881457457519648)

pearsonr(xx1_spei[flag],Trendynep)
Out[201]: (0.5509387414089815, 0.03329157812761162)

#  TRENDY 对水分的敏感性高估？
pearsonr(xx1_spei[flag],fluxsat[flag])
Out[204]: (0.3618869744754752, 0.1850210679033608)
pearsonr(xx1_smroot[flag],fluxsat[flag])
Out[205]: (0.44690224497702696, 0.09490011151753375)

pearsonr(xx1_spei[flag],xx1_csif[flag])
Out[207]: (0.1653703039377622, 0.5558668577246484)
pearsonr(xx1_smroot[flag],xx1_csif[flag])
Out[208]: (0.2446478379877342, 0.37952024633836623)

pearsonr(xx1_spei[flag],gosifgpp[flag])
Out[209]: (0.06826558284380463, 0.8089852343299259)
pearsonr(xx1_smroot[flag],gosifgpp[flag])
Out[210]: (0.2378859648242396, 0.3932410313310993)



pearsonr(xx1_smroot[flag],cte2019)
Out[164]: (-0.07147602318587608, 0.8001657924195402)

pearsonr(xx1_smroot[flag],ct2019)
Out[165]: (-0.05073396515622095, 0.8574983769446476)

pearsonr(xx1_smroot[flag],cams)
Out[166]: (-0.3717000681316926, 0.17251132690883192)

pearsonr(xx1_spei[flag],cte2019)
Out[167]: (0.029526051732752318, 0.9168083735137817)

pearsonr(xx1_spei[flag],ct2019)
Out[168]: (-0.43935737238808203, 0.10129839448506339)

pearsonr(xx1_spei[flag],cams)
Out[169]: (-0.24127419236608433, 0.3863331641400764)





pearsonr(xx1_nirv[flag],fluxsat[flag])
Out[247]: (0.7457385244558344, 0.0014136240367792405)

pearsonr(xx1_nirv[flag],gosifgpp[flag])
Out[248]: (0.8269246984945131, 0.00014331090298744023)




pearsonr(fluxsat[flag],gosifgpp[flag])
Out[246]: (0.930830354692178, 4.788219112529045e-07)

pearsonr(xx1_csif[flag],fluxsat[flag])
Out[244]: (0.9636044201267525, 7.983416928737914e-09)

pearsonr(xx1_csif[flag],gosifgpp[flag])
Out[245]: (0.9603713996224228, 1.377310002153593e-08)



pearsonr(xx1_spei[flag],xx1_smroot[flag])
pearsonr(xx1_spei[flag],xx1_ewt[flag])


pearsonr(Trendynep,jena99v43NEETW)


'''
