
# coding: utf-8

# In[5]:

from matplotlib import pyplot as plt
get_ipython().magic(u'matplotlib inline')
import sys
sys.path.append('../plotting/')
from plotting import *
import netCDF4
import numpy as np


# In[6]:

# plot IC
path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ocean_only/outflow_ice_shelf/quiet/'
icfile = 'MOM_IC.nc'
layer_o = netCDF4.Dataset(path+'layer/outflow/'+icfile)
layer_i = netCDF4.Dataset(path+'layer/ice_shelf/'+icfile)
rho_o = netCDF4.Dataset(path+'rho/outflow/'+icfile)
rho_i = netCDF4.Dataset(path+'rho/ice_shelf/'+icfile)
sigma_z_o = netCDF4.Dataset(path+'sigma_z/outflow/'+icfile)
sigma_z_i = netCDF4.Dataset(path+'sigma_z/ice_shelf/'+icfile)
z_o = netCDF4.Dataset(path+'z/outflow/'+icfile)
z_i = netCDF4.Dataset(path+'z/ice_shelf/'+icfile)
sigma_o = netCDF4.Dataset(path+'sigma/outflow/'+icfile)
sigma_i = netCDF4.Dataset(path+'sigma/ice_shelf/'+icfile)
# Read the horizontal coordinate which is the same for all configurations 
xq = layer_o.variables['latq'][:] # This is the coordinate of the cell corners (u-points in 1D)
xq = np.concatenate(([0],xq)) # Inserts left most edge of domain in to coordinate
record = 0
plt.figure(figsize=(18,9))
plt.subplot(6,2,1); plot_section(layer_i, record, xq, variable='Temp', eta='eta'); plt.title('Layer', fontsize=18);
plt.subplot(6,2,2); plot_section(rho_i, record, xq, variable='Temp', eta='eta'); plt.title(r'$\rho$', fontsize=18);
plt.subplot(6,2,3); plot_section(layer_o, record, xq, variable='Temp', eta='eta'); plt.title(r'$\sigma$', fontsize=18);
plt.subplot(6,2,4); plot_section(rho_o, record, xq, variable='Temp', eta='eta'); plt.title(r'$\rho$', fontsize=18);
plt.subplot(6,2,5); plot_section(sigma_z_i, record, xq, variable='Temp', eta='eta'); plt.title('$\sigma z*$', fontsize=18);
plt.subplot(6,2,6); plot_section(z_i, record, xq, variable='Temp', eta='eta'); plt.title(r'$z*$', fontsize=18);
plt.subplot(6,2,7); plot_section(sigma_z_o, record, xq, variable='Temp', eta='eta'); plt.title(r'$\sigma z*$', fontsize=18);
plt.subplot(6,2,8); plot_section(z_o, record, xq, variable='Temp', eta='eta'); plt.title(r'$z*$', fontsize=18);
plt.tight_layout()


# In[7]:

# Open the output from the four experiments
prog_file='prog__0001_002.nc'
layer_file = netCDF4.Dataset(path+'layer/'+prog_file)
rho_file = netCDF4.Dataset(path+'rho/'+prog_file)
sigma_file = netCDF4.Dataset(path+'sigma/'+prog_file)
z_file = netCDF4.Dataset(path+'z/'+prog_file)


# In[8]:

# Read the horizontal coordinate which is the same for all configurations 
xq = layer_file.variables['yq'][:] # This is the coordinate of the cell corners (u-points in 1D)
xq = np.concatenate(([0],xq)) # Inserts left most edge of domain in to coordinate


# In[9]:

record = 9
print 'Time is:',layer_file.variables['Time'][record], 'days'
plt.figure(figsize=(18,9))
plt.subplot(2,2,1); plot_section(layer_file, record, xq); plt.title('Layer', fontsize=18);
plt.subplot(2,2,2); plot_section(rho_file, record, xq); plt.title(r'$\rho$', fontsize=18);
plt.subplot(2,2,3); plot_section(sigma_file, record, xq); plt.title(r'$\sigma$', fontsize=18);
plt.subplot(2,2,4); plot_section(z_file, record, xq); plt.title(r'$z^*$', fontsize=18);
plt.tight_layout()


# In[10]:

variable = 'u'
plt.figure(figsize=(18,9))
plt.subplot(2,2,1); plot_section(layer_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title('Layer', fontsize=18);
plt.subplot(2,2,2); plot_section(rho_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title(r'$\rho$', fontsize=18);
plt.subplot(2,2,3); plot_section(sigma_file, record, xq, variable=variable, clim=(-2e-5,2e-5)); plt.title(r'$\sigma$', fontsize=18);
plt.subplot(2,2,4); plot_section(z_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title(r'$z^*$', fontsize=18);
plt.tight_layout()


# In[11]:

record = -1
print 'Time is:',layer_file.variables['Time'][record], 'days'
plt.figure(figsize=(18,9))
plt.subplot(2,2,1); plot_section(layer_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title('Layer', fontsize=18);
plt.subplot(2,2,2); plot_section(rho_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title(r'$\rho$', fontsize=18);
plt.subplot(2,2,3); plot_section(sigma_file, record, xq, variable=variable, clim=(-2e-5,2e-5)); plt.title(r'$\sigma$', fontsize=18);
plt.subplot(2,2,4); plot_section(z_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title(r'$z^*$', fontsize=18);
plt.tight_layout()


# In[14]:

variable = 'v'
plt.figure(figsize=(18,9))
plt.subplot(2,2,1); plot_section(layer_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title('Layer', fontsize=18);
plt.subplot(2,2,2); plot_section(rho_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title(r'$\rho$', fontsize=18);
plt.subplot(2,2,3); plot_section(sigma_file, record, xq, variable=variable, clim=(-2e-5,2e-5)); plt.title(r'$\sigma$', fontsize=18);
plt.subplot(2,2,4); plot_section(z_file, record, xq, variable=variable, clim=(-2e-10,2e-10)); plt.title(r'$z^*$', fontsize=18);
plt.tight_layout()


# In[13]:

plot_stats((layer_file,rho_file,sigma_file,z_file))

