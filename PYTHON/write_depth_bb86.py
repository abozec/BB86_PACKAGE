#!python3

## create and write depth file for HYCOM
## include utilities function
import sys
import os
from os.path import join
### might not need next line in older or newer version of proplot
os.environ['PROJ_LIB'] = '/Users/abozec/opt/anaconda3/share/proj'
import proplot as plot
import numpy as np

iodir='/Users/abozec/Documents/GitHub/BB86_PACKAGE/PYTHON/'
sys.path.append(iodir)
from hycom.io import write_hycom_depth

# define path
io=join(iodir+'../topo/')
file_bat_new='depth_BB86_02_python'
pl=0

# size of the domain
idm=101 ; jdm=101

# definition of the bathymetry
depth = 5000. ## constant bathy everywhere

# closed boundary (latbdy=0) or open boundaries (latbdy=1) or cyclic (latbdy =2)
latbdy = 0

######### END of the USER inputs ##########
vmiss = 2.**100  ## HYCOM missing values

## get the Depth
bathy= np.zeros([jdm, idm])
bathy[:,:] = depth
print('Depth Ok')

## Plot
if pl == 1:
   levels2=np.linspace(4000,6000,12)
   ## plot the field
   fig, axs = plot.subplots(nrows=1,figsize=(0.5*11,0.5*8.5))
   axs[0].format(title='Bathymetry BB86')
   m = axs[0].contourf(bathy, cmap='Mako', extend='max',levels=levels2)
   fig.colorbar(m, loc='b')
   plot.show()

## Mask by missing values if any
bathy[bathy == 0.] = vmiss

## Get the edge of the domain right
if latbdy == 0:
   bathy[jdm-1,:] = vmiss
   bathy[:,idm-1] = vmiss
   bathy[0, :] = vmiss
   bathy[:, 0] = vmiss
elif latbdy == 1:
   bathy[jdm-1,:] = vmiss
   bathy[:,idm-1] = vmiss
elif latbdy == 2:
   bathy[jdm-1,:] = vmiss
   bathy[0,:] = vmiss


## write file
write_hycom_depth(bathy,join(io,file_bat_new+'.a'),'Bathymetry BB86')
