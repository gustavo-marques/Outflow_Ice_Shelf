#!/usr/bin/env python

# generate files for used in the comparison between the flow behavior in 
# two similar test cases: flow downslope (outflow) and flow upslope (ice shelf)
# Author: Gustavo Marques

#import matplotlib
#matplotlib.use('Agg')
import argparse
from netCDF4 import MFDataset, Dataset
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Generate files for outflow/ice shelf study.
      ''',
  epilog='Written by Gustavo Marques, Feb. 2017.')

  parser.add_argument('-nx', type=int, default=8,
      help='''The total number of grid points in the x direction (default = 8).''')

  parser.add_argument('-ny', type=int, default=200,
      help='''The total number of grid points in the y direction (default = 200).''')

  parser.add_argument('-nz', type=int, default=63,
      help='''Number of model layers (default = 63).''')
  
  parser.add_argument('-W', type=float, default=8.,
      help='''Domain width in the x direction (km). Default is 8.''')

  parser.add_argument('-cshelf_lenght', type=float, default=20.,
      help='''Continental shelf lenght in the y direction (km). Default is 20.''')

  parser.add_argument('-slope_lenght', type=float, default=50.,
      help='''Continental shelf slope lenght in the y direction (km). Default is 50.''')

  parser.add_argument('-L', type=float, default=200.,
      help='''Domain lenght in the y direction (km). Default is 200''')

  parser.add_argument('-max_depth', type=float, default=1.0e3,
      help='''Maximum ocean depth (m). Default is 1000.''')

  parser.add_argument('-min_depth', type=float, default=100.0,
      help='''Minimum ocean depth (m). Default is 100.''')

  parser.add_argument('-ts', type=float, default=0.0,
      help='''Surface temperature (C). Default is 0.''')

  parser.add_argument('-tb', type=float, default=-1.0,
      help='''Bottom temperature (C). Default is -1''')

  parser.add_argument('-tout', type=float, default=9999.0,
      help='''Initial outflow temperature (C). Default is 9999.0 which means tout is not used.''')

  parser.add_argument('-tshelf', type=float, default=9999.0,
      help='''Initial sub-ice-shelf temperature (C). Default is 9999.0 which means tshelf is not used.''')

  parser.add_argument('-sponge_width', type=float, default=5.0,
      help='''Widht of sponge layer (km). Default is 5.''')

  parser.add_argument('-coupled_run', help='''Generate all the files needed to run an ocean_SIS2 simulation.''', action="store_true")
  
  parser.add_argument('-debug', help='''Adds prints and plots to help debug.''', action="store_true")
  
  parser.add_argument('-ice_shelf', help='''Make an ice shelf file.''', action="store_true")
  
  parser.add_argument('-surface_pressure_scale', type=float, default=8996.40,
      help='''A scaling factor to convert ice shelf thickness into a surface pressure 
              (kg/ (s^2 m2)). Default is 8996.40''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   nx = args.nx ; ny = args.ny; nz = args.nz
   L = args.L; W = args.W; D = args.max_depth

   dy = L/ny; dx = W/nx
   y = np.arange(dy/2.,L,dy)
   x = np.arange(dx/2.,W,dx)

   # create dir. to place figures
   os.system('mkdir PNG')

   #if args.ts_restart:
      # use T/S from specified run
   #   make_ts_restart(x,y,args)
   #else: # default
      # initial T/S from climatology
   rho, z = make_ts(x,y,args)
   
   # create topography
   Ocean_Depth = make_topo(x,y,z,rho,args)

   if args.coupled_run:

     make_mosaic(x,y,Ocean_Depth,args) 

   # create forcing
   #make_forcing(x,y,args) 
   
   return

def make_ice_shelf(x,y,z,rho,depth_target,args):
   '''
   Write ice shelf netcdf file. 
   '''
   x = x * 1.0e3 # im m
   y = y * 1.0e3 # im m
   dy = y[1]-y[0]
   dx = x[1]-x[0]
   
   # depth_target is the depth where we want to place
   # the ice shelf. However, we need to take into
   # account the ocean density (hydrostatic pressure) 
   # to find the corresponding ice shelf thickness.
   thick = np.zeros((args.ny,args.nx))
   dz = z[1] - z[0]
   znew = np.linspace(0,z[-1]+dz,len(z)+1)
   p_tmp = np.zeros(len(z)+1)
   for i in range(args.nx):
      for j in range(args.ny):
         p_tmp[:] = 0.0
         for k in range(1,len(p_tmp)): # compute hydrostatic pressure
             p_tmp[k] = p_tmp[k-1] + (rho[k-1,j,i]*9.8*dz)
         # use linear interp to find ice shelf thickness
         p_target = np.interp(depth_target[j,i], znew, p_tmp)
         thick[j,i] = p_target/args.surface_pressure_scale
         
   if args.debug:
      plt.plot(y,thick[:,0],'k')
      plt.savefig('PNG/ice_shelf_profile.png')

   area = np.ones(thick.shape) * dx * dy
   area[thick == 0.0] = 0.0

   # Create a netcdf horizontal ocean-grid file
   name = 'IC_IS'
   ncfile = Dataset(name+'.nc','w')
   ncfile.createDimension('nx',args.nx)
   ncfile.createDimension('ny',args.ny)
   nthick = ncfile.createVariable('thick','double',('ny','nx',))
   nthick.units = 'm'
   nthick.standard_name =  'ice shelf thickness'
   narea = ncfile.createVariable('area','double',('ny','nx',))
   narea.units = 'm2'
   narea.standard_name =  'ice shelf area'  

   # write into nc file
   nthick[:,:] = thick[:,:]
   narea[:,:] = area[:,:]

   ncfile.sync()
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!') 

def make_mosaic(x,y,Ocean_Depth,args):
   '''
   Create files used in coupled ice/ocean runs.
   '''
   ni = args.nx;   nj = args.ny
   snj, sni = 2*nj, 2*ni
   dxo = (x[1] - x[0]) * 1.0e3 # in m
   dyo = (y[1] - y[0]) * 1.0e3 # in m

   # Number of land points used nl GM: There has to be at least ONE land point
   nl=0
   #Define land points:
   for j in range(nj):
      for i in range(ni):
          if Ocean_Depth[j,i]==0:
              nl=nl+1

   # Using a supergrid refinement of 2, the grid lengths, area and coordinates are:
   dx = (0.5 * dxo) * np.ones((snj+1,sni))
   dy = (0.5 * dyo) * np.ones((snj,sni+1))
   area = (0.25 * (dxo * dyo)) * np.ones((snj,sni))
   x = np.zeros((snj+1,sni+1))
   x[:,1:] = np.cumsum( dx, axis=-1 )
   y = np.zeros((snj+1,sni+1))
   y[1:,:] = np.cumsum( dy, axis=-2 )

   # Create a netcdf horizontal ocean-grid file
   name = 'ocean_hgrid'
   ncfile = Dataset(name+'.nc','w')
   ncfile.createDimension('nx',sni)
   ncfile.createDimension('ny',snj)
   ncfile.createDimension('nxp',sni+1)
   ncfile.createDimension('nyp',snj+1)
   ncfile.createDimension('string',255)
   dx_h = ncfile.createVariable('dx','f8',('nyp','nx',))
   dx_h.units = 'm'
   dy_h = ncfile.createVariable('dy','f8',('ny','nxp',))
   dy_h.units = 'm'
   area_h = ncfile.createVariable('area','f8',('ny','nx',))
   area_h.units = 'm2'
   x_h = ncfile.createVariable('x','f8',('nyp','nxp',))
   x_h.units = 'm'
   y_h = ncfile.createVariable('y','f8',('nyp','nxp',))
   y_h.units = 'm'
   tile = ncfile.createVariable('tile','c',('string',)) 
   dx_h[:,:] = dx[:,:]
   dy_h[:,:] = dy[:,:]
   area_h[:,:] = area[:,:]
   x_h[:,:] = x[:,:]
   y_h[:,:] = y[:,:]
   tile[:] = '\000' * 255
   tile[:5] = 'tile1'
   set_string(tile,'tile1')
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # create ocean_mask file
   name = 'ocean_mask'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('nx',ni)
   rg.createDimension('ny',nj)
   mask = rg.createVariable('mask','f4',('ny','nx',))
   mask[:,:] = 1.
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # Create a mosaic description for the grid file
   name = 'ocean_mosaic'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('ntiles',1)
   rg.createDimension('string',255)
   mosaic = rg.createVariable('mosaic','c',('string',))
   mosaic.standard_name = 'grid_mosaic_spec'
   mosaic.children = 'contacts'
   mosaic.grid_descriptor = ''
   gridlocation = rg.createVariable('gridlocation','c',('string',))
   gridlocation.standard_name = 'grid_file_location'
   gridfiles = rg.createVariable('gridfiles','c',('ntiles','string',))
   gridtiles = rg.createVariable('gridtiles','c',('ntiles','string',))
   rg.grid_version = '0.2'
   # Fill in data
   mosaic[:] = '\000' * 255
   mosaic[:12] = 'ocean_mosaic'
   gridlocation[:] = '\000' * 255
   gridlocation[:2] = './'
   gridfiles[:] = '\000' * 255
   gridfiles[0,:14] = 'ocean_hgrid.nc'
   gridtiles[:] = '\000' * 255
   gridtiles[0,:5] = 'tile1'
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   name = 'atmos_mosaic_tile1Xland_mosaic_tile1'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('ncells',nl)  #It is unclear whether this works when nl=0. It does work for nl>0
   rg.createDimension('two',2)
   contact = rg.createVariable('contact','c',('string',))
   contact.standard_name = 'grid_contact_spec'
   contact.contact_type = 'exchange'
   contact.parent1_cell = 'tile1_cell'
   contact.parent2_cell = 'tile2_cell'
   contact.xgrid_area_field = 'xgrid_area'
   contact.distant_to_parent1_centroid = 'tile1_distance'
   contact.distant_to_parent2_centroid = 'tile2_distance'
   tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two'))
   tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
   tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two'))
   tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
   xgrid_area = rg.createVariable('xgrid_area','f8',('ncells'))
   xgrid_area.standard_name = 'exchange_grid_area'
   xgrid_area.units = 'm2'
   tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
   tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
   tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
   tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
   rg.grid_version = '0.2'
   # Fill in data
   contact[:] = '\000' * 255
   contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
   for i in range(nl):
      tile1_cell[i,:] = [ni,nj]
      tile2_cell[i,:] = [ni,nj]

   xgrid_area[:] = dxo * dyo
   tile1_distance[:] = 0.
   tile2_distance[:] = 0.
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   name = 'atmos_mosaic_tile1Xocean_mosaic_tile1'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
   rg.createDimension('two',2)
   contact = rg.createVariable('contact','c',('string',))
   contact.standard_name = 'grid_contact_spec'
   contact.contact_type = 'exchange'
   contact.parent1_cell = 'tile1_cell'
   contact.parent2_cell = 'tile2_cell'
   contact.xgrid_area_field = 'xgrid_area'
   contact.distant_to_parent1_centroid = 'tile1_distance'
   contact.distant_to_parent2_centroid = 'tile2_distance'
   tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
   tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
   tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
   tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
   xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
   xgrid_area.standard_name = 'exchange_grid_area'
   xgrid_area.units = 'm2'
   tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
   tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
   tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
   tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
   rg.grid_version = '0.2'
   # Fill in data
   contact[:] = '\000' * 255
   contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
   for i in range(nl):
      tile1_cell[i,:] = [ni,nj]
      tile2_cell[i,:] = [ni,nj]
   
   xgrid_area[:] = dxo * dyo
   count=-1
   for j in range(nj):
       for i in range(ni):
          if Ocean_Depth[j,i]!=0:
            count=count+1
            tile1_cell[count] = [i,j]
            tile2_cell[count] = [i,j]
            tile1_distance[count] = [0,0]
            tile2_distance[count] = [0,0]
            xgrid_area[count] = dxo * dyo

   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!') 

   name = 'land_mosaic_tile1Xocean_mosaic_tile1'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
   rg.createDimension('two',2)
   contact = rg.createVariable('contact','c',('string',))
   contact.standard_name = 'grid_contact_spec'
   contact.contact_type = 'exchange'
   contact.parent1_cell = 'tile1_cell'
   contact.parent2_cell = 'tile2_cell'
   contact.xgrid_area_field = 'xgrid_area'
   contact.distant_to_parent1_centroid = 'tile1_distance'
   contact.distant_to_parent2_centroid = 'tile2_distance'
   tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
   tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
   tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
   tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
   xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
   xgrid_area.standard_name = 'exchange_grid_area'
   xgrid_area.units = 'm2'
   tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
   tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
   tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
   tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
   rg.grid_version = '0.2'
   # Fill in data
   contact[:] = '\000' * 255
   contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
   for i in range(nl):
      tile1_cell[i,:] = [ni,nj]
      tile2_cell[i,:] = [ni,nj]
   
   xgrid_area[:] = dxo * dyo
   count=-1
   for j in range(nj):
      for i in range(ni):
        if Ocean_Depth[j,i]!=0:
            count=count+1
            tile1_cell[count] = [i,j]
            tile2_cell[count] = [i,j]
            tile1_distance[count] = [0,0]
            tile2_distance[count] = [0,0]
            xgrid_area[count] = dxo * dyo
  
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')


   name = 'mosaic' # sometimes grid_spec is called mosaic.nc
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('nfile_aXo',1) # -1 is for a single land point
   rg.createDimension('nfile_aXl',1) # -1 is for a single land point
   rg.createDimension('nfile_lXo',1) # -1 is for a single land point
   atm_mosaic_dir = rg.createVariable('atm_mosaic_dir','c',('string',))
   atm_mosaic_dir.standard_name = 'directory_storing_atmosphere_mosaic'
   atm_mosaic_file = rg.createVariable('atm_mosaic_file','c',('string',))
   atm_mosaic_file.standard_name = 'atmosphere_mosaic_file_name'
   atm_mosaic = rg.createVariable('atm_mosaic','c',('string',))
   atm_mosaic.standard_name = 'atmosphere_mosaic_name'
   lnd_mosaic_dir = rg.createVariable('lnd_mosaic_dir','c',('string',))
   lnd_mosaic_dir.standard_name = 'directory_storing_land_mosaic'
   lnd_mosaic_file = rg.createVariable('lnd_mosaic_file','c',('string',))
   lnd_mosaic_file.standard_name = 'land_mosaic_file_name'
   lnd_mosaic = rg.createVariable('lnd_mosaic','c',('string',))
   lnd_mosaic.standard_name = 'land_mosaic_name'
   ocn_mosaic_dir = rg.createVariable('ocn_mosaic_dir','c',('string',))
   ocn_mosaic_dir.standard_name = 'directory_storing_ocean_mosaic'
   ocn_mosaic_file = rg.createVariable('ocn_mosaic_file','c',('string',))
   ocn_mosaic_file.standard_name = 'ocean_mosaic_file_name'
   ocn_mosaic = rg.createVariable('ocn_mosaic','c',('string',))
   ocn_mosaic.standard_name = 'ocean_mosaic_name'
   ocn_topog_dir = rg.createVariable('ocn_topog_dir','c',('string',))
   ocn_mosaic_dir.standard_name = 'directory_storing_ocean_topog'
   ocn_topog_file = rg.createVariable('ocn_topog_file','c',('string',))
   ocn_topog_file.standard_name = 'ocean_topog_file_name'
   aXo_file = rg.createVariable('aXo_file','c',('nfile_aXo','string',))
   aXo_file.standard_name = 'atmXocn_exchange_grid_file'
   aXl_file = rg.createVariable('aXl_file','c',('nfile_aXl','string',))
   aXl_file.standard_name = 'atmXlnd_exchange_grid_file'
   lXo_file = rg.createVariable('lXo_file','c',('nfile_lXo','string',))
   lXo_file.standard_name = 'lndXocn_exchange_grid_file'
   #Global attributes
   rg.grid_version = '0.2'
   rg.code_version = "$Name:  $"
   rg.history = "/net2/aja/workspace/MOM6-examples/ice_ocean_SIS2/OM4_025/preprocessing/fre_nctools/tools/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ocean_topog.nc"
   # Fill in data
   atm_mosaic_dir[:] = '\000' * 255
   atm_mosaic_dir[:2] = './'
   atm_mosaic_file[:] = '\000' * 255
   atm_mosaic_file[:15] = 'ocean_mosaic.nc'
   atm_mosaic[:] = '\000' * 255
   atm_mosaic[:12] = 'atmos_mosaic'
   lnd_mosaic_dir[:] = '\000' * 255
   lnd_mosaic_dir[:2] = './'
   lnd_mosaic_file[:] = '\000' * 255
   lnd_mosaic_file[:15] = 'ocean_mosaic.nc'
   lnd_mosaic[:] = '\000' * 255
   lnd_mosaic[:11] = 'land_mosaic'
   ocn_mosaic_dir[:] = '\000' * 255
   ocn_mosaic_dir[:2] = './'
   ocn_mosaic_file[:] = '\000' * 255
   ocn_mosaic_file[:15] = 'ocean_mosaic.nc'
   ocn_mosaic[:] = '\000' * 255
   ocn_mosaic[:12] = 'ocean_mosaic'
   ocn_topog_dir[:] = '\000' * 255
   ocn_topog_dir[:2] = './'
   ocn_topog_file[:] = '\000' * 255
   ocn_topog_file[:14] = 'ocean_topog.nc'
   aXo_file[:,:] = '\000' * 255
   aXo_file[:,:40] = 'atmos_mosaic_tile1Xocean_mosaic_tile1.nc'
   aXl_file[:,:] = '\000' * 255
   aXl_file[:,:39] = 'atmos_mosaic_tile1Xland_mosaic_tile1.nc'
   lXo_file[:,:] = '\000' * 255
   lXo_file[:,:39] = 'land_mosaic_tile1Xocean_mosaic_tile1.nc'
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')
   os.system('cp mosaic.nc grid_spec.nc')
   print ('*** Run make_quick_mosaic ***')
   os.system('module load fre')
   os.system('make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ice_topog.nc')

def get_profile(t,i,j,var,depth,z,args,vname):
   '''
   Get data (var), remove mask, smooth profile then interpolate to z grid.
   '''
   if args.mean_profile:
      data = np.mean(Dataset('WOA05_pottemp_salt.nc').variables[var][:,:,j,i], axis=0)
      
   else:
      data = Dataset('WOA05_pottemp_salt.nc').variables[var][t,:,j,i]

   # replace mask value with last good value
   tmp = np.nonzero(data.mask==False)[0][-1]
   data[tmp+1::] = data[tmp]
   # smooth  
   data = gaussian_filter(data,2) # filter std 2
   # interpo
   f1 = interp1d(depth, data)
   if args.debug:
      plt.figure()
      plt.plot(data,-depth,'b',f1(z),-z,'rx')
      plt.title(vname)
      plt.savefig('PNG/'+vname+'.png')

   return f1(z)

def set_freezing_temp(T,S,z,y,args):
    '''
    Set the temp. in the upper 50 m to the surface freezing point, then smooth data
    '''
    d = np.nonzero(z<=50.)[0]
    for k in range(len(d)):
        for j in range(args.ny):
           # set T to freezing point just near the coast
           if y[j] <= (args.ISL + args.cshelf_lenght):
              T[0,k,j,:] = eos.tfreeze(S[0,k,j,:],1.0e5) 

    #for j in range(args.ny):
    #  for i in range(args.nx):
    #     T[0,d,j,i] = gaussian_filter(T[0,d,j,i],1) # filter std 1

    # smooth T and S in the cross slope direction
    #for k in range(args.nz):
    #   for j in range(args.ny):
    #       T[0,k,:,i] = gaussian_filter(T[0,k,:,i],2)
    #       S[0,k,:,i] = gaussian_filter(T[0,k,:,i],2)

    return T,S

def make_ts_restart(x,y,args):
   '''
   Extract last T/S from ts_file(ocean_month_z) and save them in a ncfile. 
   '''
   if args.ts_file == '': 
      print
      raise MyError( 'Parameter -ts_file must be specified.')

   print 'Processing ', args.ts_file + ', this might take a while...'
   xh = Dataset(args.ts_file).variables['xh'][:]
   yh = Dataset(args.ts_file).variables['yh'][:]
   zl = Dataset(args.ts_file).variables['zt'][:]
   t = args.restart_time_indice
   temp = Dataset(args.ts_file).variables['temp'][t,:]
   salt = Dataset(args.ts_file).variables['salt'][t,:]

   XH,YH = np.meshgrid(xh,yh)
   # 3D fields where data will be interpolated 
   # x,y, are high-res and zl is from month_z
   temp3D = np.zeros((len(zl),len(y),len(x)))
   salt3D = np.zeros((len(zl),len(y),len(x)))
   temp3D_new = np.zeros((1,args.nz,len(y),len(x)))
   salt3D_new = np.zeros((1,args.nz,len(y),len(x)))

   dz = args.max_depth / args.nz
   z = np.linspace(0.5*dz,args.max_depth-0.5*dz,args.nz)

   # replace mask value with last good value
   # first in the vertical, then in the horizontal (done in the second loop)
   for i in range(len(xh)):
      for j in range(len(yh)):
         tmp = np.nonzero(temp.mask[:,j,i]==False)[0]
         if len(tmp)>0:
            tmp = tmp[-1]
            if tmp < (len(zl)-1):
               temp[tmp+1::,j,i] = temp[tmp,j,i]
               salt[tmp+1::,j,i] = salt[tmp,j,i]

   for k in range(len(zl)):
       print 'level', str(k) + ' out of ' + str(len(zl))
       # replace mask value with last good value
       # in the x dir
       for i in range(len(xh)):
          tmp = np.nonzero(temp.mask[k,:,i]==False)[0][0]
          temp[k,0:tmp,i] = temp[k,tmp,i]
          salt[k,0:tmp,i] = salt[k,tmp,i]
          
       #ftemp = interpolate.RectBivariateSpline(yh, xh, temp[k,:,:], bbox=[0,args.L,0,args.W])
       ftemp = interpolate.RectBivariateSpline(yh, xh, temp[k,:,:])
       temp3D[k,:] = ftemp(y,x)
       fsalt = interpolate.RectBivariateSpline(yh, xh, salt[k,:,:], bbox=[0,args.L,0,args.W])
       salt3D[k,:] = fsalt(y,x)

   # now interpolate on final grid
   for i in range(len(x)):
      for j in range(len(y)):
          ftemp = interp1d(zl,temp3D[:,j,i])
          fsalt = interp1d(zl,salt3D[:,j,i])
          temp3D_new[0,:,j,i] = ftemp(z)
          salt3D_new[0,:,j,i] = fsalt(z)

   # write ncfile
   write_ic_ncfile('ic_ts',x,y,z,temp3D_new,salt3D_new)   

def make_ts(x,y,args):
   '''
   Make initial conditions. 
   '''
   # model depth
   z = np.linspace(0,args.max_depth,args.nz) # positive downward
   # 3d fields
   temp3D = np.zeros((1,args.nz,len(y),len(x)))
   salt3D = np.ones((1,args.nz,len(y),len(x))) * 0.0

   dtdz = (args.ts - args.tb)/args.max_depth #
   dsdz = 1.0/args.max_depth #
   
   for j in range(args.ny):
      if (y[j] < args.cshelf_lenght) and (args.tout != 9999. or args.tshelf != 9999.):
         if args.tshelf != 9999.:
            temp3D[0,:,j,:] = args.tshelf 
            salt3D[0,:,j,:] = 1.0 
         else:
            temp3D[0,:,j,:] = args.tout
            salt3D[0,:,j,:] = 1.0
      else: 
         for k in range(args.nz):
            temp3D[0,k,j,:] = - dtdz * z[k]
            #salt3D[0,k,:,i] = dsdz * z[k]

   # compute linear eos
   rho = get_rho(salt3D,temp3D)
   print 'rho min/max', rho.min(), rho.max()
 
   if args.debug:
      plt.figure()
      plt.contourf(y,-z,rho[0,:,:,0])
      plt.colorbar()
      #plt.contour(y,-z,rho[0,:,:,0]-1000.,layers-1000,colors='k',linewidths=2)
      plt.title('Density - 1000.')
      plt.savefig('PNG/rho_section.png')

   layers = np.linspace(rho.min(),rho.max(),args.nz) 
   # create ncfiles
   # vertical coord. (layers)
   name = 'layers'
   ncfile = Dataset(name+'.nc','w')
   ncfile.createDimension('Layer',args.nz)
   Layer = ncfile.createVariable('Layer',np.dtype('double').char,('Layer'))
   Layer.units = 'kg/m^3'
   Layer[:] = layers[:]

 ##  ncfile.close()
 ##  print ('*** SUCCESS creating '+name+'.nc!')
   write_ic_ncfile('ic_ts',x,y,z,temp3D,salt3D)
   dz_half = (z[1] - z[0]) * 0.5
   layer_depth = np.arange(dz_half,z[-1]-dz_half,len(z))
   return rho[0,:], layer_depth

def get_rho(salt,temp):
   '''
   Compute linear EoS
   '''  
   alpha = 5.3e-5; rho0 = 1025.0; beta = 0.0
   return (rho0*(1+beta*salt-alpha*temp))

def write_ic_ncfile(name,x,y,z,T,S):
   '''
   Write the initial T/S condition into a netcdf file called name.

   '''
   # 1) The name of the z-space input file used to initialize
   #  the layer thicknesses, temperatures and salinities.
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('DEPTH',len(z))
   ncfile.createDimension('LAT',len(y))
   ncfile.createDimension('LON',len(x))
   ncfile.createDimension('TIME',1)

   # create variables
   LAT = ncfile.createVariable('LAT',np.dtype('double').char,('LAT'))
   LAT.units = 'degrees_north'
   LAT.long_name = 'h point nominal latitude'
   LAT.cartesian_axis = 'Y'
   LAT[:] = y[:]

   LON = ncfile.createVariable('LON',np.dtype('double').char,('LON'))
   LON.units = 'degrees_east'
   LON.long_name = 'h point nominal longitude'
   LON.cartesian_axis = 'X'
   LON[:] = x[:]

   TIME = ncfile.createVariable('TIME',np.dtype('double').char,('TIME'))
   TIME.units = 'days since 0001-01-01 00:00:00'
   TIME.calendar = 'noleap'
   TIME.cartesian_axis = 'T'
   TIME[0] = 0.0

   DEPTH = ncfile.createVariable('DEPTH',np.dtype('double').char,('DEPTH'))
   DEPTH.units = 'm'
   DEPTH.direction = -1
   DEPTH.cartesian_axis = 'Z'
   DEPTH[:] = z[:]

   PTEMP = ncfile.createVariable('PTEMP',np.dtype('float32').char,('TIME','DEPTH','LAT','LON'), fill_value = -1.e+34)
   PTEMP.missing_value = -1.e+34
   PTEMP[:] = T[:] 

   SALT = ncfile.createVariable('SALT',np.dtype('float32').char,('TIME','DEPTH','LAT','LON'), fill_value = -1.e+34)  
   SALT.missing_value = -1.e+34
   SALT[:] = S[:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # 2) The file from which the coordinate densities are read
   #name = 'ts_ic_profile'
   #ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   #ncfile.createDimension('Layer',args.nz)

   # create variables
   #PTEMP = ncfile.createVariable('PTEMP',np.dtype('double').char,('Layer'))
   #PTEMP.units = 'Celcius'
   #PTEMP.long_name = 'Potential Temperature'
   #PTEMP[:] = temp_int[:]

   #SALT = ncfile.createVariable('SALT',np.dtype('double').char,('Layer'))
   #SALT.units = 'PSU'
   #SALT.long_name = 'Salinity'
   #SALT[:] = salt_int[:]

   #ncfile.close()
   #print ('*** SUCCESS creating '+name+'.nc!')

def celcius_to_kelvin(tc):
    '''
    Convert temp in C to temp in K
    '''
    return tc + 273.15

def make_forcing(x,y,args):
   # forcing parameters
   Ly = args.L # domain size km
   W = args.W # domain width km
   CSL = args.cshelf_lenght  # km
   tau_asf = args.taux # default is 0.075 # N/m^2
   tauy_max = args.tauy_max # def is 0.05 N/m2
   tauy_min = args.tauy_min  # min northward wind
   sponge = args.sponge_width # 100 km
   ISL = args.ISL
   nx = len(x); ny = len(y)
   if args.add_seasonal_cycle:
     nt = 365*4 # one year every 6 hours
     print('About to create forcing file with',nt, ' records...')
   else:
     nt = 1
   
   # atmos/ice forcing
   time_days = np.zeros(nt)
   t_bot = np.zeros((nt,ny,nx)) 
   wind_x = np.zeros((nt,ny,nx)) 
   wind_y = np.zeros((nt,ny,nx)) 
   tau_x = np.zeros((nt,ny,nx)) 
   tau_y = np.zeros((nt,ny,nx))
   liq = np.zeros((nt,ny,nx))
   snow = np.zeros((nt,ny,nx))
   salt = np.zeros((nt,ny,nx))

   # atm params
   t1min = celcius_to_kelvin(args.t1min) # -40
   dt1 = args.dt1 # 30
   t2min = celcius_to_kelvin(args.t2min)  # -15
   dt2 = args.dt2 # 20
   t3min = celcius_to_kelvin(args.t3min) # 0
   dt3 = args.dt3 # 10

   w1min = args.w1min
   dw1 = args.dw1
   w2min = args.w2min
   dw2 = args.dw2
   w2_y_lim = args.w2_y_lim
   w3min = args.w3min
   dw3 = args.dw3
   w3_y_lim = args.w3_y_lim
   w4min = args.w4min
   dw4 = args.dw4

   # loop in time
   for t in range(nt):
     # wind and heat
     time_days[t] = t * 0.25
     temp_period = args.temp_forcing_period # in days
     wind_period = args.wind_forcing_period # in days
     #print'Time:',time_days[t]
     temp_cos = 1 #np.cos(np.pi*time_days[t]/temp_period)**2
     wind_cos = 1 #np.cos(np.pi*time_days[t]/wind_period)**2
     temp_season_cos = np.cos(np.pi*time_days[t]/temp_period)**2
     temp_season_sin = np.sin(np.pi*time_days[t]/temp_period)**2
     wind_season_cos = np.cos(np.pi*time_days[t]/wind_period)**2
     wind_season_sin = np.sin(np.pi*time_days[t]/wind_period)**2
   
     # temperatures, start with summer 
     t1 = t1min + temp_season_cos*dt1
     t2 = t2min + temp_season_cos*dt2
     t3 = t3min + temp_season_cos*dt3

     # winds
     w1 = w1min + wind_season_cos*dw1 # katabatic
     w2 = w2min + wind_season_cos*dw2 # along-slope south of CS break
     w3 = w3min + wind_season_cos*dw3 # ASF
     w4 = w4min + wind_season_cos*dw4 # ACC

     # wind x-dir
     #if w1 == 0.0: # mode 2 and mode 3
     #   wind_x_pos = args.ISL #- (season_cos * 200.) # x wind moves with season
     #   wind_x_pos = args.wind_x_pos
     #else: # mode1
     #   wind_x_pos = args.wind_x_pos

     Lw2 = w2_y_lim - ISL
     Lw3 = w3_y_lim - w2_y_lim 
     Lw4 = Ly - w3_y_lim 
     for j in range(ny):
       if y[j] <= ISL:
          tau_x[t,j,:] = 0.0; wind_x[t,j,:] = 0.0
       elif y[j] > ISL and y[j] <= w2_y_lim: # shelf wind
          tmp = 1* np.pi*(y[j]-ISL)/(Lw2)
	  tau_x[t,j,:] = (tau_asf * np.sin(tmp)**2) # not time dependent
	  wind_x[t,j,:] = (w2 * np.sin(tmp)**2)
       elif y[j] > w2_y_lim and y[j] <= w3_y_lim: # shelf break wind (ASF)
          tmp = 1* np.pi*(y[j]-w2_y_lim)/(Lw3)
	  tau_x[t,j,:] = (tau_asf * np.sin(tmp)**2) # not time dependent
	  wind_x[t,j,:] = (w3 * np.sin(tmp)**2) 
       else: # ACC
          tmp = 1* np.pi*(y[j]-w3_y_lim)/(Lw4)
	  tau_x[t,j,:] = 0.0
	  wind_x[t,j,:] = (w4 * np.sin(tmp)**2)

     # heat
     # Follow http://onlinelibrary.wiley.com/doi/10.1002/wea.436/pdf
     # Q is ~ linear
     
     # katabatic
     # follow ~ http://journals.ametsoc.org/doi/pdf/10.1175/1520-0493(1994)122%3C0671%3ADOATDM%3E2.0.CO%3B2
     # # has a gaussian shape: tauy = tauy_max * np.exp(((x-W/2.)**2)/(2*W_v10))
     W_v10 = 20000. # km
     # tau_y and heat are zero at y = 300 and 400 km, respect.
     tmp = 400.0 - ISL; tmp1 = 400.0 - ISL
     tmp_inv = 1.0/tmp; tmp1_inv = 1.0/tmp1

     #delta_tauy = tauy_max - tauy_min
     #delta_wind_y = wind_y_max - wind_y_min
     # uncomment below to add variations in tauy
     #delta_wind_y = delta_wind_y - (wind_season_cos * delta_wind_y*0.75) # gets weaker in summer by 75%
     efold = args.tauy_efold  # km

     # exp decay forcing
     for j in range(ny):
	 if y[j] < ISL+efold:
	    t_bot[t,j,:] = t1
            if args.tauy_confined:
                for i in range(nx):
                   tmp1 =  np.exp((-(x[i]-W*0.5)**2)/(2*W_v10))
                   wind_y[t,j,i] = (w1 * tmp1) 
            else:
	        wind_y[t,j,:] = w1 #(delta_wind_y * wind_cos) + wind_y_min
	 #elif y[j] >= ISL+efold and y[j] <= wind_x_pos:
	 else:
	    t_bot[t,j,:] = (t1 - t2) * np.exp(-(y[j]-ISL-efold)/(efold)) + t2 
            if args.tauy_confined:
               for i in range(nx):
                   tmp1 =  np.exp((-(x[i]-W*0.5)**2)/(2*W_v10))
                   #tau_y[t,j,i] = (((delta_tauy * np.exp(-(y[j]-ISL-efold)/(2*efold))) \
                   #         * tmp1) * wind_cos) + tauy_min
                   #wind_y[t,j,i] = (((delta_wind_y * np.exp(-(y[j]-ISL-efold)/(2*efold))) \
                   #         * tmp1) * wind_cos) + wind_y_min
                   wind_y[t,j,i] = w1 * np.exp(-(y[j]-ISL-efold)/(2*efold)) * tmp1
            else:
               #tau_y[t,j,:] = ((delta_tauy * np.exp(-(y[j]-ISL-efold)/(2*efold))) \
               #               * wind_cos) + tauy_min
               #wind_y[t,j,:] = ((delta_wind_y * np.exp(-(y[j]-ISL-efold)/(2*efold))) \
               #               * wind_cos) + wind_y_min
               wind_y[t,j,:] = w1 * np.exp(-(y[j]-ISL-efold)/(2*efold))
         #else: # linear
         #   t2_tmp = (t1 - t2) * np.exp(-(500.-ISL-efold)/(efold)) + t2
         #   t_bot[t,j,:] =  ((t3-t2_tmp)/(Ly-500))*(y[j]-500.)+ t2_tmp
         #   wind_y[t,j,:] = 0.0

     # linear temp
     if args.linear_forcing:
       for j in range(ny):
         if y[j] < ISL+efold:
            t_bot[t,j,:] = t1
         else:
            t_bot[t,j,:] = (t3 - t1)*((y[j]-(ISL+efold))/(Ly-(ISL+efold))) + t1


     # lprec, fprec
     lprec = args.liq_prec # lprec 
     fprec = args.frozen_prec # lprec 
     #tmp = args.cshelf_lenght
     for j in range(ny):
	if y[j] < w2_y_lim:
	   liq[t,j,:] = 0.0; snow[t,j,:] = 0.0
	elif y[j]>= w2_y_lim and y[j]< (Ly-sponge):
           #tmp = (Ly-sponge) - wind_x_pos
	   #liq[t,j,:] = season_cos * lprec *2./3. + lprec * 1./3. #* np.sin((np.pi * (y[j]-wind_x_pos))/ tmp)
	   liq[t,j,:] = lprec 
           #snow[t,j,:] = season_sin * fprec  #* np.sin((np.pi * (y[j]-wind_x_pos))/ tmp)
           snow[t,j,:] = fprec  
        else:
           liq[t,j,:] = 0.0
           snow[t,j,:] = 0.0

     for j in range(ny):
        if y[j] > ISL and y[j] <= ISL + 50.:
           salt[t,j,:] = 5.0e-6
   #
   # End of time loop
   #

   # power lost due to heat loss
   #grid_area = (y[1]-y[0]) * (x[1]-x[0])
   #power = np.sum(heat*grid_area)
   #print 'Power due to sensible heat (J/s, negative means loss):', power

   # plots
   if args.debug:

      plt.figure()
      u=wind_x[0,::5,::5]; v=wind_y[0,::5,::5]; mag = np.sqrt(u**2 + v**2)
      plt.quiver(x[::5],y[::5],u, v, mag, cmap = plt.cm.seismic)
      plt.colorbar()
      plt.plot(x,np.ones(nx)*ISL,'k')
      plt.plot(x,np.ones(nx)*(ISL+CSL),'k')
      plt.title('Wind (m/s)')
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/wind_vector.png')

      #plt.figure()
      #plt.subplot(211)
      #plt.plot(y,u10[0,:,1])
      #plt.title('u10'); plt.xlabel('y [km]'); plt.ylabel('m/s')
      #plt.subplot(212)
      #plt.contourf(x,y,u10[0,:,:]); plt.colorbar()
      #plt.xlabel('x [km]'); plt.ylabel('y [km]')
      #plt.savefig('PNG/u10.png')

      #plt.figure()
      #plt.subplot(211)
      #plt.plot(y,v10[0,:,1])
      #plt.title('v10'); plt.xlabel('y [km]'); plt.ylabel('m/s')
      #plt.subplot(212)
      #plt.contourf(x,y,v10[0,:,:]); plt.colorbar()
      #plt.xlabel('x [km]'); plt.ylabel('y [km]')
      #plt.savefig('PNG/v10.png')

      plt.figure()
      plt.subplot(211)
      plt.plot(y,tau_x[0,:,1])
      plt.title('taux'); plt.xlabel('y [km]'); plt.ylabel('Pa')
      plt.subplot(212)
      plt.contourf(x,y,tau_x[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/taux.png')

      plt.figure()
      plt.subplot(211)
      plt.plot(y,tau_y[0,:,args.nx/2.])
      plt.title('tauy'); plt.xlabel('y [km]'); plt.ylabel('Pa')
      plt.subplot(212)
      plt.contourf(x,y,tau_y[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/tauy.png') 

   # create ncfile
   # open a new netCDF file for writing.
   # # used when forcing is applied in the atm
   name = 'forcing_10'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('LON',nx)
   ncfile.createDimension('LAT',ny)
   ncfile.createDimension('TIME',None)
   # create variables
   LON = ncfile.createVariable('LON',np.dtype('double').char,('LON'))
   LON.units = 'km'
   LON.long_name = 'h point nominal longitude'
   LON.cartesian_axis = 'X'
   LON[:] = x[:]

   LAT = ncfile.createVariable('LAT',np.dtype('double').char,('LAT'))
   LAT.units = 'km'
   LAT.long_name = 'h point nominal latitude'
   LAT.cartesian_axis = 'Y'
   LAT[:] = y[:]

   time = ncfile.createVariable('TIME',np.dtype('double').char,('TIME'))
   time.long_name = 'time'
   time.units = 'days since 0001-01-01 00:00:00'
   time.cartesian_axis = 'T'
   time.calendar_type = 'NOLEAP'
   time.calendar = 'NOLEAP'
   time.bounds = 'time_bounds'
   if args.add_seasonal_cycle:
      time.modulo = '' # make it cyclic
      time[:] = time_days[:] + (time_days[1] - time_days[0])
   else:
      time[0] = 0.0

   t_10 = ncfile.createVariable('T_10',np.dtype('float32').char,('TIME','LAT','LON')) 
   t_10.long_name = 'Air Temperature'
   t_10.units = 'Kelvin'
   t_10[:] = t_bot[:]

   u_10 = ncfile.createVariable('U_10',np.dtype('float32').char,('TIME','LAT','LON'))
   u_10.long_name = 'U wind'
   u_10.units = 'm/s'
   u_10[:] = wind_x[:]

   v_10 = ncfile.createVariable('V_10',np.dtype('float32').char,('TIME','LAT','LON'))
   v_10.long_name = 'V wind'
   v_10.units = 'm/s'
   v_10[:] = wind_y[:]
  
   salt_flux = ncfile.createVariable('salt_flux',np.dtype('float32').char,('TIME', 'LAT', 'LON'))
   salt_flux.units = 'kg/(m^2 s)'
   salt_flux.long_name = 'salt flux'
   salt_flux[:] = salt # + adds salt from ocean
 
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # used when forcing is applied in the sea ice
   name = 'forcing'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('xh',nx)
   ncfile.createDimension('yh',ny)
   ncfile.createDimension('xq',nx)
   ncfile.createDimension('yq',ny)
   ncfile.createDimension('time',None)
   ncfile.createDimension('nv',2)

   # create variables
   xh = ncfile.createVariable('xh',np.dtype('double').char,('xh'))
   xh.units = 'km'
   xh.long_name = 'h point nominal longitude'
   xh.cartesian_axis = 'X'
   xh[:] = x[:]
   
   xq = ncfile.createVariable('xq',np.dtype('double').char,('xq'))
   xq.units = 'km'
   xq.long_name = 'q point nominal longitude'
   xq.cartesian_axis = 'X'
   xq[:] = x[:]

   yh = ncfile.createVariable('yh',np.dtype('double').char,('yh'))
   yh.units = 'km'
   yh.long_name = 'h point nominal latitude'
   yh.cartesian_axis = 'Y'
   yh[:] = y[:]
   
   yq = ncfile.createVariable('yq',np.dtype('double').char,('yq'))
   yq.units = 'km'
   yq.long_name = 'q point nominal latitude'
   yq.cartesian_axis = 'Y'
   yq[:] = y[:]

   time = ncfile.createVariable('time',np.dtype('double').char,('time'))
   time.long_name = 'time'
   time.units = 'days since 0001-01-01 00:00:00'
   time.cartesian_axis = 'T'
   time.calendar_type = 'NOLEAP'
   time.calendar = 'NOLEAP'
   time.bounds = 'time_bounds'
   if args.add_seasonal_cycle:
      time.modulo = '' # make it cyclic
      time[:] = time_days[:] + (time_days[1] - time_days[0])
   else:
      time[0] = 0

   nv = ncfile.createVariable('nv',np.dtype('double').char,('nv'))   
   nv.long_name = 'vertex number'
   nv.units = 'none'
   nv.cartesian_axis = 'N'
   nv[:] = [1,2]

   if args.coupled_run:
     u_flux = ncfile.createVariable('u_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     u_flux.units = 'Pa'
     u_flux.missing_value = 1.e+20
     u_flux.long_name = 'i-direction wind stress'
     u_flux[:] = -tau_x[:] # change sign in ice

     v_flux = ncfile.createVariable('v_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     v_flux.units = 'Pa'
     v_flux.missing_value = 1.e+20
     v_flux.long_name = 'j-direction wind stress'
     v_flux[:] = -tau_y[:] # change sign in ice

     t_flux = ncfile.createVariable('t_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     t_flux.units = 'Watt meter-2'
     t_flux.missing_value = 1.e+20
     t_flux.long_name = 'Sensible heat flux'
     t_flux[:] = 0.0 # -heat[:] # change sign in ice

     # latent heat
     lt = ncfile.createVariable('latent',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     lt.units = 'Watt meter-2'
     lt.missing_value = 1.e+20
     lt.long_name = 'Latent heat flux'
     lt[:] = 0.0 # change sign in ice

     salt_flux = ncfile.createVariable('salt_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     salt_flux.units = 'kg/(m^2 s)'
     salt_flux.missing_value = 1.e+20
     salt_flux.long_name = 'salt flux'
     salt_flux[:] = 1.0e-5 # + adds salt from ocean

     lprec = ncfile.createVariable('lprec',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     lprec.units = 'kg/(m^2 s)'
     lprec.missing_value = 1.e+20
     lprec.long_name = 'liquid precipitation'
     lprec[:] = liq[:] # positive is adding water into the ocean

     fprec = ncfile.createVariable('fprec',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     fprec.units = 'kg/(m^2 s)'
     fprec.missing_value = 1.e+20
     fprec.long_name = 'froze precipitation'
     fprec[:] = snow[:] # positive is adding water into the ocean

   else:
     SW = ncfile.createVariable('SW',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     SW.units = 'Watt meter-2'
     SW.missing_value = 1.e+20
     SW.long_name = 'surface_net_downward_shortwave_flux'
     SW[:] = 0.0

     LW = ncfile.createVariable('LW',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     LW.units = 'Watt meter-2'
     LW.missing_value = 1.e+20
     LW.long_name = 'surface_net_downward_longwave_flux'
     LW[:] = 0.0

     latent = ncfile.createVariable('latent',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     latent.units = 'Watt meter-2'
     latent.missing_value = 1.e+20
     latent.long_name = 'Latent heat flux into ocean due to fusion and evaporation'
     latent[:] = 0.0

     sensible = ncfile.createVariable('sensible',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     sensible.units = 'Watt meter-2'
     sensible.missing_value = 1.e+20
     sensible.long_name = 'surface_downward_sensible_heat_flux'
     sensible[:] = heat[:]

     evap = ncfile.createVariable('evap',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     evap.units = 'kilogram meter-2 second-1'
     evap.missing_value = 1.e+20
     evap.long_name = 'Evaporation at ocean surface (usually negative)'
     evap[:] = 0.0

     taux = ncfile.createVariable('taux',np.dtype('float32').char,('time', 'yh', 'xq'), fill_value = 1.e+20)
     taux.units = 'Pascal'
     taux.missing_value = 1.e+20
     taux.long_name = 'Zonal Wind Stress'
     taux[:] = tau_x[:]   

     tauy = ncfile.createVariable('tauy',np.dtype('float32').char,('time', 'yq', 'xh'), fill_value = 1.e+20)
     tauy.units = 'Pascal'
     tauy.missing_value = 1.e+20
     tauy.long_name = 'Meridional Wind Stress'
     tauy[:] = 0.0

     ustar = ncfile.createVariable('ustar',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     ustar.units = 'meter second-1'
     ustar.missing_value = 1.e+20
     ustar.long_name = 'Surface friction velocity'
     ustar[:] = 0.0

     SST = ncfile.createVariable('SST',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = -1.e+34)
     SST.units = 'Celsius'
     SST.missing_value = -1.e+34
     SST.long_name = 'Sea Surface Temperature'
     SST[:] = 0.0

     SSS = ncfile.createVariable('SSS',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = -1.e+34)
     SSS.units = 'PSU'
     SSS.missing_value = -1.e+34
     SSS.long_name = 'Sea Surface Salinity'
     SSS[:] = 0.0

   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
  
   print ('*** Run make_quick_mosaic ***')
   os.system('module load fre')
   os.system('make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ice_topog.nc')
   return

def make_topo(x,y,z,rho,args):
   # parameters
   name = 'ocean_topog'
   Hmin = args.min_depth  # shelf's depth
   Ys = args.cshelf_lenght # len of shelf in km
   Ws = args.slope_lenght # len of slope in km
   Hmax = args.max_depth # max depth
   W = args.W  # domain width in km
   [X,Y] = np.meshgrid(x,y) # x,y in km
   nx = len(x); ny = len(y)
   D = np.zeros((ny,nx))

   a = (Hmax - 2*Hmin)/(Ws*1.0e3) # slope
   for j in range(ny):
       if y[j] <= Ys:
          D[j,:] = Hmin
       elif y[j] > Ys and y[j] <= (Ys + Ws):  
          D[j,:] = a * (y[j] - Ys)*1.0e3 + Hmin
       elif y[j] > (Ys + Ws) and y[j] <= (2*Ys + Ws):
          D[j,:] = (Hmax - Hmin)
       else:
          D[j,:] = Hmax
 
   if args.ice_shelf:
      thick = -(D-Hmax)
      D[:,:] = Hmax
      make_ice_shelf(x,y,z,rho,thick,args)
      # lenght of ice shelf
      args.ISL = 2*Ys + Ws
 
   # to avoid sea ice formation under ice shelves,
   # two topography files need to be constructed.
   # The coupler topo is where the cavity is masked.

   #if not args.ice_shelf:
   #  D[0,:] = 0.
   #else:
   #  D[0,0] = 0.
   D[0,:] = 0.0

   Dice = D.copy()
   #for j in range(ny):
   #  if y[j]<= ymask:
   #  #if y[j]<= ISL:
   #    Dice[j,:] = 0.0

   # 1) topography used in the coupler
   # open a new netCDF file for writing.
   name = 'ice_topog'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('ntiles',1)
 
   # create variables
   nxx = ncfile.createVariable('nx',np.dtype('double').char,('nx'))
   nxx.units = 'km'
   nxx.description = 'x location of cell centers'
   nxx.long_name = 'x location of cell centers'
   nxx[:] = x[:]

   nyy = ncfile.createVariable('ny',np.dtype('double').char,('ny'))
   nyy.units = 'km'
   nyy.description = 'y location of cell centers'
   nyy.long_name = 'y location of cell centers'
   nyy[:] = y[:]

   depth = ncfile.createVariable('depth',np.dtype('float32').char,('ny','nx'))
   depth.units = 'm'
   depth.description = 'depth at h points'
   depth.long_name = 'depth at h points'
   depth[:,:] = Dice[:,:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # 2) topography used ny the ocean
   # open a new netCDF file for writing.
   name = 'ocean_topog'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('ntiles',1)

   # create variables
   nxx = ncfile.createVariable('nx',np.dtype('double').char,('nx'))
   nxx.units = 'km'
   nxx.description = 'x location of cell centers'
   nxx.long_name = 'x location of cell centers'
   nxx[:] = x[:]

   nyy = ncfile.createVariable('ny',np.dtype('double').char,('ny'))
   nyy.units = 'km'
   nyy.description = 'y location of cell centers'
   nyy.long_name = 'y location of cell centers'
   nyy[:] = y[:]

   depth = ncfile.createVariable('depth',np.dtype('float32').char,('ny','nx'))
   depth.units = 'm'
   depth.description = 'depth at h points'
   depth.long_name = 'depth at h points'
   depth[:,:] = D[:,:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
 
   return D

# A helper function to use when writing strings in a netcdf file
def set_string(variable, value):
   """Sets "variable" to "value" padded with blanks where
   "variable" is a netcdf variable object and "value" is a string."""
   variable[:] = '\000' * variable.shape[0]
   variable[:len(value)] = value
   return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
