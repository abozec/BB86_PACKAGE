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
from hycom.io import write_hycom_grid

# define path
io=join(iodir+'../topo/')
file_grid_new='regional.grid.BB86.python'

## size of the domain
idm = 101 ; jdm = 101

## longitude/latitude starting point + resolution  (in degrees) (NB: those
## variables are  not used in HYCOM)
ini_lon = 0. ; ini_lat = 0.
res = 0.20

## scale dx and dy (in m) (used in HYCOM)
dx = 20e3 ; dy = dx

## grid type
mapflg=0  ##  uniform or mercator


######### END of the USER inputs #############
## missing value in HYCOM
vmiss = 2.**100

## Get the p-point the grid
plon = np.zeros([jdm, idm])
plat = np.zeros([jdm, idm])

## Longitude
plon[:, 0] = ini_lon
for i in np.arange(idm-1)+1:
   plon[:,i] = plon[:, i-1] + res

 ## Latitude
plat[0, :] = ini_lat
for j in np.arange(jdm-1)+1:
   plat[j, :] = plat[j-1, :] + res

print('p-points grid OK')

## Declaration of the grid tabs
qlon = np.zeros([jdm,idm])
qlat = np.zeros([jdm,idm])
ulon = np.zeros([jdm,idm])
ulat = np.zeros([jdm,idm])
vlon = np.zeros([jdm,idm])
vlat = np.zeros([jdm,idm])

pang = np.zeros([jdm,idm])

pscx = np.zeros([jdm,idm])
pscy = np.zeros([jdm,idm])

qscx = np.zeros([jdm,idm])
qscy = np.zeros([jdm,idm])

uscx = np.zeros([jdm,idm])
uscy = np.zeros([jdm,idm])

vscx = np.zeros([jdm,idm])
vscy = np.zeros([jdm,idm])

cori = np.zeros([jdm,idm])
pasp = np.zeros([jdm,idm])

## Longitude/Latitude for each point (Q-points stay at the same location)
## longitude
vlon = plon
for i in np.arange(idm-1)+1:
   qlon[:, i] = 0.5*(plon[:,i]+plon[:, i-1])

diff = plon[0, 3]-plon[0, 2]
qlon[:, 0] = plon[:, 0]-0.5*diff
ulon = qlon

## latitude
ulat = plat
for j in np.arange(jdm-1)+1:
   qlat[j, :] = 0.5*(plat[j, :]+plat[j-1, :])

diff = plat[3, 0]-plat[2, 0]
qlat[0, :] = plat[0, :]-0.5*diff
vlat = qlat

## simplified grid with prescribed dx and dy
pscx[:,:] = dx ; pscy[:,:] = dy
uscx[:,:] = dx ; uscy[:,:] = dy
vscx[:,:] = dx ; vscy[:,:] = dy
qscx[:,:] = dx ; qscy[:,:] = dy

## Coriolis
beta = 2.e-11
for  j in np.arange(jdm):
   for  i in np.arange(idm):
      cori[j, i]=.93e-4+(j-int((jdm-2)/2))*dx*beta


## pang ( p-angle  for rotated grid)
## here uniform or mercator so pang = 0.

## pasp= paspect: pscx/pscy
pasp=pscx/pscy

##  Writing new grid file
write_hycom_grid(join(io+file_grid_new+'.a'), plon, plat, qlon,  qlat,  ulon, ulat, \
 vlon, vlat, pang,pscx, pscy, qscx,  qscy, uscx, uscy,  vscx, vscy, cori, pasp,mapflg)

print('Writing grid file done ')
