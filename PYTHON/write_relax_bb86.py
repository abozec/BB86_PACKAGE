#!python3

# create relax files
import sys
import os
from os.path import join
### might not need next line in older or newer version of proplot
os.environ['PROJ_LIB'] = '/Users/abozec/opt/anaconda3/share/proj'
import proplot as plot
import numpy as np

iodir='/Users/abozec/Documents/GitHub/BB86_PACKAGE/PYTHON/'
from hycom.io import read_hycom_depth, write_hycom_relax
from UTILITIES.density import tofsig

file_topo = 'depth_BB86_02_python.a'
file_int_new = 'relax_int_BB86_python'  ## !! without .a or .b !!
file_tem_new = 'relax_tem_BB86_python'
file_sal_new = 'relax_sal_BB86_python'
io=iodir+'../relax/011/'
pl = 0 ##  1 or 0 for plot or not

## size of the domain
idm = 101 ; jdm = 101
kdm = 3  ## number of layers (first layer very thin)
tdm = 12 ## 12 month climatology

## interface depth (!!!first interface depth always 0. !!!)
id = np.zeros([kdm])
id = np.array([0., 1., 500.])

## target density
## from bb86: rho=27.01037, rho=27.22136
d = np.zeros([kdm])
d = np.array([27.0100, 27.01037,27.22136])

## salinity profile
sa = np.zeros(kdm)
sa = np.array([37., 37.,  37.])

##### END of the USER inputs ######
## constants
rho = 1000.
g = 9.806
vmiss = 2.**100

## temperature profile for sigma0
sigma = 0
te = np.zeros([kdm])
for k in np.arange(kdm):
   te[k] = tofsig(sa[k], d[k], sigma)

print(te)

## read bathy for mask
depth_data=join(iodir+'../topo/',file_topo)
bathy=read_hycom_depth(depth_data,idm,jdm)
ind=np.where(np.isnan(bathy) == True)
bathy[ind]=0.

## create mask
mask2d = np.ones([jdm, idm])
mask2d[ind] = 0.
mask = np.zeros([jdm, idm, kdm, tdm])
for t in np.arange(tdm):
   for k in np.arange(kdm):
      mask[:, :, k, t] = mask2d[:,:]


## interface depth files (first layer always the surface i.e. 0.)
int = np.zeros([jdm, idm, kdm, tdm])
for k in np.arange(kdm):
   int[:, :, k,:] = rho*g * id[k] ## 1st layer tiny to fake a two layer config

## make sure that the interface depths are not lower than the bathy
for t in np.arange(tdm):
   for k in np.arange(kdm):
      for i in np.arange(idm):
         for j in np.arange(jdm):
            if (int[j, i, k, t]/9806. > bathy[j, i]):
               int[j, i, k, t] = bathy[j,i]*9806.



## mask
ind3d = np.where(mask == 0.)
int[ind3d] = vmiss

## write the field in relax file
write_hycom_relax(join(io, file_int_new), 'intf', d, int)
print('Interface depth OK')

## Temperature
tem = np.zeros([jdm, idm, kdm, tdm])
for k in np.arange(kdm):
   tem[:, :, k,:] = te[k]

## mask
tem[ind3d] = vmiss

## write the field in relax file
write_hycom_relax(join(io, file_tem_new), 'temp', d, tem) ;
print('Temperature OK')

## Salinity
sal = np.zeros([jdm, idm, kdm, tdm])
for k in np.arange(kdm):
   sal[:,:, k,:] = sa[k]

## mask
sal[ind3d] = vmiss

## write the field in relax file
write_hycom_relax(join(io, file_sal_new), 'saln', d, sal);
print('Salinity OK')


## Plot
if (pl == 1):
   levels2=np.linspace(16,18,12)
   ## plot the field
   fig, axs = plot.subplots(nrows=1,figsize=(0.5*11,0.5*8.5))
   axs[0].format(title='Temperature BB86')
   m = axs[0].contourf(tem[:,:,0,0], cmap='jet', extend='neither',levels=levels2)
   fig.colorbar(m, loc='b')
   plot.show()
