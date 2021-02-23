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
sys.path.append(iodir)
from hycom.io import read_hycom_grid

io=iodir+'../topo/'
pl = 0 ##  1 or 0 for plot or not

## size of the domain
idm = 101 ; jdm = 101
tdm = 12 ## 12 month climatology

## name the new files
file_E = 'forcing.tauewd.BB86.python'
file_N = 'forcing.taunwd.BB86.python'

## Read grid
grid_data=join(io,'regional.grid.BB86.python.a')
#print(F"The fields available are: {read_field_grid_names(grid_data)}")
fieldg=['plon','plat']
grid_field= read_hycom_grid(grid_data, fieldg)
plon=grid_field['plon'][:,:]
plat=grid_field['plat'][:,:]

## Calculation of the analytical wind-stress (N/m2)
ustress = np.zeros([jdm])
stressa = -1.
sconv = 1.e-1   ## scale factor form dyn/cm2 to N/m2

## BB86 formulation
for j in np.arange(jdm):
   ustress[j] = stressa*np.cos(int(j-1)/int(jdm-1)*6.28318530718)*sconv

## figure
fig, ax = plot.subplots()
ax.plot(ustress, plat[:,0], lw=2)
ax.format(suptitle='Wind-stress', xlim=(-0.2,0.2),ylim=(0,20))
plot.show(block=False)

## Write the wind-stress files
tte = np.zeros([jdm, idm, tdm])
ttn = np.zeros([jdm, idm, tdm])

## get size of pad for writing the file
layer_size=idm*jdm
npad = 4096-np.mod(layer_size, 4096)

## we apply  a taux , no tauy
for j in np.arange(jdm):
   tte[j, :, :] = ustress[j]

## open .a file eastward wind
fout=open(iodir+'../force/'+file_E+'.a','wb')

## write .a file
for t in np.arange(tdm):
   field32=np.array(tte[:,:,t],dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

fout.close()

# writing .b file
fbout=open(iodir+'../force/'+file_E+'.b','w')
fbout.write("Analytical Eastward Wind-stress \n")
fbout.write("  \n")
fbout.write("  \n")
fbout.write("  \n")
fbout.write(F"i/jdm = {idm},{jdm}  \n")
for t in np.arange(tdm):
   fbout.write("tauewd:  month,range = {0:02d} {1:10.5E} {2:10.5E} \n".\
             format(t+1,np.nanmin(tte[:,:,t]),np.nanmax(tte[:,:,t])))

fbout.close()


## open .a file northward wind
fout=open(iodir+'../force/'+file_N+'.a','wb')

## write .a file
for t in np.arange(tdm):
   field32=np.array(ttn[:,:,t],dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

fout.close()

# writing .b file
fbout=open(iodir+'../force/'+file_N+'.b','w')
fbout.write("Analytical Northward Wind-stress \n")
fbout.write("  \n")
fbout.write("  \n")
fbout.write("  \n")
fbout.write(F"i/jdm = {idm},{jdm}  \n")
for t in np.arange(tdm):
   fbout.write("taunwd:  month,range = {0:02d} {1:10.5E} {2:10.5E} \n".\
             format(t+1,np.nanmin(ttn[:,:,t]),np.nanmax(ttn[:,:,t])))

fbout.close()
