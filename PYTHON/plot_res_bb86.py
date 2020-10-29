#!python3
# plot results from HYCOM-BB86
import sys
import os
from os.path import join
### might not need next line in older or newer version of proplot
os.environ['PROJ_LIB'] = '/Users/abozec/opt/anaconda3/share/proj'
import proplot as plot
import numpy as np

iodir='/Users/abozec/Documents/GitHub/BB86_PACKAGE/PYTHON/'
sys.path.append(iodir)
from hycom.info import read_field_names,read_field_grid_names
from hycom.io import read_hycom_fields, read_hycom_grid

## PATH
io = iodir+'../'
ps_dir=iodir+'PNG/'


## size of the domain
idm = 101 ; jdm = 101  ## size of the domain
kdm = 2                ##  number of layer in bb86
dp0 = 500.             ## thickness of the first layer (m)
tdm = 1800             ## 12 month climatology
tplot1 = 1800 ; tplot2 = 1800 ## time-stamp to plot (starts from 1)

normeref_hyc=0.05      ## vector norm of HYCOM-BB86 (in m/s)
onevectout = 2         ## plot one vector out of 'onevectorout'
min_dp = -450. ;  max_dp = 450. ## layer thickness anomaly in m (min and max)

## Save figure
PS=0

## constants
rho = 1000.    ## reference density
g = 9.806      ## gravity

## Read grid
filet='regional.grid.BB86.python.a'
grid_data=join(io,'topo/'+filet)
#print(F"The fields available are: {read_field_grid_names(grid_data)}")
fieldg=['plon','plat']
grid_field= read_hycom_grid(grid_data, fieldg)
plon=grid_field['plon'][:,:]
plat=grid_field['plat'][:,:]

##  READ the HYCOM files

io_hycom=io+'expt_01.0/data/'
for t in np.arange(tplot2-tplot1+1)+tplot1:
   print('t=',t)
   if (t <= 344):
      year = '0001'
      day = '{:03d}'.format(t+16) ## 16 u and v == 0. so we start at day 17.
   elif np.logical_and((t > 344), (t <= 704)):
      year = '0002'
      day = '{:03d}'.format(t-344)
   elif np.logical_and((t > 704),(t <= 1064)):
      year = '0003'
      day = '{:03d}'.format(t-704)
   elif np.logical_and((t > 1064), (t <=  1424)):
      year = '0004'
      day = '{:03d}'.format(t-1064)
   elif np.logical_and((t > 1424), (t <= 1784)):
      year = '0005'
      day = '{:03d}'.format(t-1424)
   elif np.logical_and((t > 1784), (t <= 2144)):
      year = '0006'
      day = '{:03d}'.format(t-1784)
   else:
      print('Day over 2144, no case available')



   file1='archv.'+year+'_'+day+'_00.a'
   file_data = join(io_hycom,file1)
   print(file_data)
   # get metadata of .[ab]
   # Printing the fields available in the file
   print(F"The fields available are: {read_field_names(file_data)}")
   # get barotropic velocities
   fields=['u_btrop','v_btrop']
   layers=[0]
   hycom_field= read_hycom_fields(file_data, fields, layers)
   ubaro=hycom_field[fields[0]][0,0:jdm,0:idm]
   vbaro=hycom_field[fields[1]][0,0:jdm,0:idm]
   # get baroclinic velocities
   ## Note that u-vel and v-vel are total velocities in time-average file (archm)
   ## Note that u-vel and v-vel are baroclinic velocities in instantaneous file (archv)
   fields=['u-vel.','v-vel.']
   layers=np.arange(kdm+1)
   hycom_field= read_hycom_fields(file_data, fields, layers)
   ubac=np.zeros([jdm,idm,kdm+1])
   vbac=np.zeros([jdm,idm,kdm+1])
   for k in np.arange(kdm+1):
      ubac[:,:,k]=hycom_field[fields[0]][k]
      vbac[:,:,k]=hycom_field[fields[1]][k]

   ## Get utot & vtot
   uhyc=np.zeros([jdm,idm,kdm+1,tdm])
   vhyc=np.zeros([jdm,idm,kdm+1,tdm])
   for k in np.arange(kdm+1):
      uhyc[:, :, k, t-1] = ubaro[:, :]+ubac[:, :,k]
      vhyc[:, :, k, t-1] = vbaro[:, :]+vbac[:, :,k]

   ## get layer thickness
   fields=['thknss'] ; layers=np.arange(kdm+1)
   hycom_field= read_hycom_fields(file_data, fields, layers)
   dp=np.zeros([jdm,idm,kdm+1])
   for k in np.arange(kdm+1):
      dp[:,:,k]= hycom_field[fields[0]][k,0:jdm,0:idm]

   dphyc= np.zeros([jdm,idm,kdm+1,tdm])
   dphyc[:,:,:,t-1]=dp[:,:,:]/(rho*g) ## convert from pressure to m

   ## put the u/v vel on the p-grid
   maskt = np.ones([jdm, idm, kdm+1, tdm])
   maskt[uhyc == 0.] = 0.
   maxval = maskt+np.roll(maskt, -1, 1)
   maxval[maxval == 0.] = np.nan

   ## get the average value of u at the p-point
   uthyc = (uhyc+np.roll(uhyc, -1, 1))/maxval
   uthyc[:,idm-1, :, :] = uthyc[:,idm-2, :, :]

   ## put the v vel on the p-grid
   ## get the vmask
   maskt = np.ones([jdm, idm, kdm+1, tdm])
   maskt[vhyc == 0.] = 0.
   maxval = maskt+np.roll(maskt, -1, 0)
   maxval[maxval == 0.] = np.nan

   ## get the average value of v at the p-point
   vthyc = (vhyc+np.roll(vhyc, -1, 0))/maxval
   vthyc[jdm-1, :, :, :] = vthyc[jdm-2,:, :, :]

   ## get the vectors to plot
   uu = uthyc[::onevectout, ::onevectout, 1, t-1]
   vv = vthyc[::onevectout, ::onevectout, 1, t-1]
   lon = plon[::onevectout, ::onevectout]
   lat = plat[::onevectout, ::onevectout]

   ## we normalize the vectors
   normeref = normeref_hyc ## in m
   print(F'min: {np.nanmin(uu)}, max:{np.nanmax(uu)}')
   print(F'min: {np.nanmin(vv)}, max:{np.nanmax(vv)}')

   ## plot
   levels2=np.linspace(min_dp,max_dp,41)
   day_plot='{:04d}'.format(t)
   ## plot the field
   fig, axs = plot.subplots(ncols=2,nrows=1,width='12in',aspect=1.5,\
                            bottom='5em',span=False)
   axs[0].format(title='HYCOM-BB86 velocity field (m/s) day='+day_plot,xlim=(0,20),\
                 xlabel='Longitude (E)',ylabel='Latitude (N)')
   q = axs[0].quiver(lon,lat, uu,vv,scale=normeref_hyc,scale_units='x')
   axs[0].quiverkey(q, X=0.1, Y=-0.1, U=normeref_hyc,\
             label='0.05 m/s', labelpos='E')

   axs[1].format(title='HYCOM-BB86 thickness (m) ano day='+day_plot,\
                 xlabel='Longitude (E)',ylabel='Latitude (N)')
   axs[1].contour(plon,plat,dphyc[:,:,1,t-1]-500., levels=levels2,color='k',lw=0.5)
   m = axs[1].contourf(plon,plat,dphyc[:,:,1,t-1]-dp0, cmap='jet', extend='both',levels=levels2)
   fig.colorbar(m, loc='r')
   if (PS == 1):
      fig.savefig(join(ps_dir,'uv_dp_d'+day_plot+'_bb86-hycom.png'),dpi=150,\
               facecolor='w', edgecolor='w',transparent=False) ## .pdf,.eps,.png
      plot.close()
   else:
      plot.show()
