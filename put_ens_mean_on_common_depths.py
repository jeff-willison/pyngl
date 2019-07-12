#!/usr/bin/env python
# coding: utf-8

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset

import numpy as np
import Nio,Ngl
import os,sys
import datetime


# In[2]:


path_to_output = "/scratch/rhe/GoM/Data-Assimilation/jan2009c/"
common_grid_file = Nio.open_file("/home/rhe/common_gom_domain.nc","r")
depths = common_grid_file.variables['z_rho'][:]*-1
ndepths = depths.shape[0]
#-- open file
f = Nio.open_file(path_to_output+"instance_0001/roms_posterior_0001_54939.nc","r")
gridf = Nio.open_file(path_to_output+"instance_0001/useast_ini_0001.nc","r")

# Start date and end date
start_date=datetime.datetime(2009,5,24)
end_date=datetime.datetime(2009,9,23)

#############################################################
days_since=datetime.datetime(1858,11,17)
start_date_ordinal=datetime.date.toordinal(start_date)
end_date_ordinal=datetime.date.toordinal(end_date)
days_since_ordinal=datetime.date.toordinal(days_since)
  
start_date_index=start_date_ordinal-days_since_ordinal
end_date_index=end_date_ordinal-days_since_ordinal
  
# Optional sanity check
#start_date_string=datetime.date.fromordinal(days_since_ordinal+start_date_index)
nfiles=int((end_date_index-start_date_index+1)/2)
filesi=list(range(start_date_index,end_date_index,2))    


# In[3]:


#ffrac=2*((2*pi)/86400)*sin(26.010683 degrees)
ffrac=0.00006378287
h=gridf.variables['h'] # (eta_rho,xi_rho)
hneg = h[:,:]*(-1.)
 
Cs_r = gridf.variables['Cs_r']  # s_rho
hc_nio   = gridf.variables['hc'] # single value
hc = hc_nio.get_value()
vtransform_nio = gridf.variables['Vtransform']
vtransform = vtransform_nio.get_value()
Sc_r = gridf.variables['s_rho'] # s_rho
hinv=1/h[:,:]
 
pm = gridf.variables['pm'][:,:] # these only change with latitude (non-conformal grid)
pn = gridf.variables['pn'][:,:]
dx = 1/pm
dy = 1/pn
lats=gridf.variables['lat_rho'][:,0]
lons=gridf.variables['lon_rho'][0,:]
nlats=lats.shape[0]
nlons=lons.shape[0]


# In[7]:


## Keep this in the same cell for easier troubleshooting


# for i in range(0,nfiles):
for time in range(0,nfiles):
  file_date_string=datetime.date.fromordinal(filesi[time]+days_since_ordinal)
  
  os.system("rm -rf ../ncs/common_depths_"+str(file_date_string)+".nc")     #-- delete file if it exists
  dataset=Dataset("../ncs/common_depths_"+str(file_date_string)+".nc",'w',format='NETCDF4_CLASSIC')  #-- create dimensions

  level = dataset.createDimension('level', ndepths)
  lat = dataset.createDimension('lat', nlats)
  lon = dataset.createDimension('lon', nlons)    
  ocean_time = dataset.createDimension('ocean_time', 1) 
  
  nc_u = dataset.createVariable('u','f',('ocean_time','level','lat','lon'),fill_value=1e37)
  nc_v = dataset.createVariable('v','f',('ocean_time','level','lat','lon'),fill_value=1e37)
  nc_t = dataset.createVariable('temp','f',('ocean_time','level','lat','lon'),fill_value=1e37)
  nc_s = dataset.createVariable('salt','f',('ocean_time','level','lat','lon'),fill_value=1e37)
  nc_z = dataset.createVariable('zeta','f',('ocean_time','lat','lon'),fill_value=1e37)    
  
  nc_bottom_u = dataset.createVariable('bottom_u_velocity','f',('ocean_time','lat','lon'),fill_value=1e37)
  nc_bottom_v = dataset.createVariable('bottom_v_velocity','f',('ocean_time','lat','lon'),fill_value=1e37)
  nc_bottom_t = dataset.createVariable('bottom_temp','f',('ocean_time','lat','lon'),fill_value=1e37)
  nc_bottom_s = dataset.createVariable('bottom_salt','f',('ocean_time','lat','lon'),fill_value=1e37)
  
  nc_time = dataset.createVariable('ocean_time','d',('ocean_time'))
  nc_lat = dataset.createVariable('lat','f',('lat'))
  nc_lon = dataset.createVariable('lon','f',('lon'))
  nc_level = dataset.createVariable('level','f',('level'))
  
  nc_level[:]=depths[:]
  nc_lat[:] = lats[:]
  nc_lon[:] = lons[:]
  nc_time[0] = filesi[time]

  print(nc_time)
  setattr(dataset.variables['ocean_time'],'units',"days since 1858-11-17 00:00:00")
  setattr(dataset.variables['ocean_time'],'calendar',"gregorian")
  
  setattr(dataset.variables['zeta'],'long_name','free-surface')
  setattr(dataset.variables['zeta'],'units','meter')
  setattr(dataset.variables['zeta'],'time','ocean_time')
  setattr(dataset.variables['zeta'],'coordinates','ocean_time lon lat')
  setattr(dataset.variables['zeta'],'field','free-surface, scalar, series')    

  setattr(dataset.variables['temp'],'long_name','potential temperature')
  setattr(dataset.variables['temp'],'units','Celsius')
  setattr(dataset.variables['temp'],'time','ocean_time')
  setattr(dataset.variables['temp'],'coordinates','ocean_time level lon lat')
  setattr(dataset.variables['temp'],'field','temperature, scalar, series')

  setattr(dataset.variables['salt'],'long_name','salinity')
  setattr(dataset.variables['salt'],'time','ocean_time')
  setattr(dataset.variables['salt'],'coordinates','ocean_time level lon lat')
  setattr(dataset.variables['salt'],'field','salinity, scalar, series')
  
  setattr(dataset.variables['u'],'long_name','u-momentum component')
  setattr(dataset.variables['u'],'units','meter second-1')    
  setattr(dataset.variables['u'],'time','ocean_time')
  setattr(dataset.variables['u'],'coordinates','ocean_time level lon lat')
  setattr(dataset.variables['u'],'field','u-velocity, scalar, series')
  
  setattr(dataset.variables['v'],'long_name','v-momentum component')
  setattr(dataset.variables['v'],'units','meter second-1')    
  setattr(dataset.variables['v'],'time','ocean_time')
  setattr(dataset.variables['v'],'coordinates','ocean_time level lon lat')
  setattr(dataset.variables['v'],'field','v-velocity, scalar, series')    
  
  setattr(dataset.variables['bottom_u_velocity'],'long_name','u-momentum component at sea floor')
  setattr(dataset.variables['bottom_u_velocity'],'units','meter second-1')
  setattr(dataset.variables['bottom_u_velocity'],'time','ocean_time')
  setattr(dataset.variables['bottom_u_velocity'],'coordinates','lon lat ocean_time')
  setattr(dataset.variables['bottom_u_velocity'],'field','bottom_u, scalar, series')

  setattr(dataset.variables['bottom_v_velocity'],'long_name','v-momentum component at sea floor')
  setattr(dataset.variables['bottom_v_velocity'],'units','meter second-1')
  setattr(dataset.variables['bottom_v_velocity'],'time','ocean_time')
  setattr(dataset.variables['bottom_v_velocity'],'coordinates','lon lat ocean_time')
  setattr(dataset.variables['bottom_v_velocity'],'field','bottom_v, scalar, series')    

  setattr(dataset.variables['bottom_temp'],'long_name','temperature at sea floor')
  setattr(dataset.variables['bottom_temp'],'units','Celsius')
  setattr(dataset.variables['bottom_temp'],'time','ocean_time')
  setattr(dataset.variables['bottom_temp'],'coordinates','lon lat ocean_time')
  setattr(dataset.variables['bottom_temp'],'field','bottom_temp, scalar, series')
  
  setattr(dataset.variables['bottom_salt'],'long_name','temperature at sea floor')
  setattr(dataset.variables['bottom_salt'],'units','Celsius')
  setattr(dataset.variables['bottom_salt'],'time','ocean_time')
  setattr(dataset.variables['bottom_salt'],'coordinates','lon lat ocean_time')
  setattr(dataset.variables['bottom_salt'],'field','bottom_temp, scalar, series')
  
  print(time,nfiles)
  print(file_date_string)
  
  f = Nio.open_file(path_to_output+"output_mean."+str(filesi[time])+".nc","r")
  u=f.variables['u'][0,:,:,:]
  v=f.variables['v'][0,:,:,:]
  t=f.variables['temp'][0,:,:,:]
  s=f.variables['salt'][0,:,:,:]
  tb=f.variables['temp']
  zb = f.variables['zeta']
  vin = tb[0,:,:,:]
  zeta = zb[0,:,:]
  njni=zeta.shape
  nj=njni[0]
  ni=njni[1]
  
  
  # Put u and v points onto rho points 
  dims=vin.shape
  
  dimsu  = u.shape
  dimY   = dimsu[1]
  dimX   = dimsu[2]
  ur = np.empty(dims)
  ur[:,:,1:dimX-1] = 0.5*(u[:,:,:dimX-2] + u[:,:,1:dimX-1])
  ur[:,:,0]=ur[:,:,1]
  ur[:,:,dimX]=ur[:,:,dimX-1]
  #print(ur.shape)
  #print(dimX)
  dimsv  = v.shape
  dimY   = dimsv[1]
  dimX   = dimsv[2]

  vr = np.empty(dims)
  vr[:,1:dimY-1,:] = 0.5*(v[:,:dimY-2,:] + v[:,1:dimY-1,:])
  vr[:,0,:]=vr[:,1,:]
  vr[:,dimY,:]=vr[:,dimY-1,:]
  
  # have all I need for depth calculation
 
  depth=np.empty(dims)
  uinterpolated=np.empty([ndepths,nj,ni])
  vinterpolated=np.empty([ndepths,nj,ni])
  tinterpolated=np.empty([ndepths,nj,ni])
  sinterpolated=np.empty([ndepths,nj,ni])
  N=dims[0]

  if vtransform == 2:
    for k in range(0,N):
      cff = 1/(hc + h[:,:])
      cffr = hc*Sc_r[k] + h[:,:]*Cs_r[k]
      depth[k,:,:]=(zeta + ( zeta + h )*cffr*cff)
  if vtransform == 1: 
    for k in range(0,N):
      cffr = hc*(Sc_r[k] - Cs_r[k])
      depth[k,:,:]=cffr+Cs_r[k]*h[:,:] + zeta*(1+(cffr+Cs_r[k]*h)*hinv)
      
  
  for i in range(0,ni):
    for j in range(0,nj):
      x=depth[:,j,i]
      uplane=ur[:,j,i]
      vplane=vr[:,j,i]
      tplane=t[:,j,i]
      splane=s[:,j,i]
      ui=np.interp(depths,x,uplane)
      uinterpolated[:,j,i]=ui
      vi=np.interp(depths,x,vplane)
      vinterpolated[:,j,i]=vi
      tinterpolated[:,j,i]=np.interp(depths,x,tplane)
      sinterpolated[:,j,i]=np.interp(depths,x,splane)
      
  
  for k in range (0,ndepths):
    tinterpolated[k,:,:][hneg>depths[k]] = 1e37

  # Use temp vals to filter out extrapolated areas below bathy
  uinterpolated[tinterpolated==1e37] = 1e37
  vinterpolated[tinterpolated==1e37] = 1e37
  sinterpolated[tinterpolated==1e37] = 1e37
  
  # set land mask for bottom u and v
  ur[t>1000] = 1e37
  vr[t>1000] = 1e37

  nc_bottom_u[0,:,:] = ur[0,:,:]
  nc_bottom_v[0,:,:] = vr[0,:,:]
  nc_bottom_t[0,:,:] = t[0,:,:]    
  nc_bottom_s[0,:,:] = s[0,:,:]    
  nc_z[0,:,:] = zb[0,:,:]    
  
  nc_u[0,:,:,:] = uinterpolated[:,:,:] # write numpy structured array to netcdf compound var
  nc_v[0,:,:,:] = vinterpolated[:,:,:] # write numpy structured array to netcdf compound var    
  nc_t[0,:,:,:] = tinterpolated[:,:,:] # write numpy structured array to netcdf compound var
  nc_s[0,:,:,:] = sinterpolated[:,:,:] # write numpy structured array to netcdf compound var

  dataset.close()
  #os.system("ncatted -O -a _FillValue,salt,c,f,1e37 ../ncs/common_depths_"+str(file_date_string)+".nc")
  #os.system("ncatted -O -a _FillValue,temp,c,f,1e37 ../ncs/common_depths_"+str(file_date_string)+".nc")


# In[ ]:




