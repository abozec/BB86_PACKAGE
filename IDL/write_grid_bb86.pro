PRO write_grid_bb86

   ;; include the Utilities functions
   iodir = '/Net/yucatan/abozec/BB86_PACKAGE/IDL/' ;; where you are
   !path=expand_path('+/'+iodir+'UTILITIES')+':'+expand_path('+'+!dir)

   close, /all
   ;; PATH
   io = iodir+'/../topo/'
   file_grid_new = 'regional.grid.BB86' ;; !! without .a or .b !!

   ;; size of the domain
   idm = 101
   jdm = 101 

   ;; longitude/latitude starting point + resolution  (in degrees) (NB: those
   ;; variables are  not used in HYCOM)
   ini_lon = 0.
   ini_lat = 0.
   res = 0.20

   ;; scale dx and dy (in m) (used in HYCOM)
   dx = 20e3
   dy = dx


   ;;;;;;; END of the USER inputs ;;;;;;;;;;;;;;;

   ;; missing value in HYCOM
   vmiss = 2.^100

   ;; Get the p-point the grid
   plon = fltarr(idm, jdm)
   plat = fltarr(idm, jdm)

   ;; Longitude
   plon(0, *) = ini_lon
   FOR i = 1, idm-1 DO plon(i, *) = plon(i-1, *) + res
   ;; Latitude
   plat(*, 0) = ini_lat
   FOR j = 1, jdm-1 DO plat(*, j) = plat(*, j-1) + res

   print, 'p-points grid OK'


   ;; Declaration of the grid tabs
   idm_hd = idm
   jdm_hd = jdm
   
   plon_hd = plon
   plat_hd = plat

   qlon_hd = fltarr(idm_hd,  jdm_hd)
   qlat_hd = fltarr(idm_hd,  jdm_hd)
   ulon_hd = fltarr(idm_hd,  jdm_hd)
   ulat_hd = fltarr(idm_hd,  jdm_hd)
   vlon_hd = fltarr(idm_hd,  jdm_hd)
   vlat_hd = fltarr(idm_hd,  jdm_hd)

   pang_hd = fltarr(idm_hd,  jdm_hd)

   pscx_hd = fltarr(idm_hd,  jdm_hd)
   pscy_hd = fltarr(idm_hd,  jdm_hd)

   qscx_hd = fltarr(idm_hd,  jdm_hd)
   qscy_hd = fltarr(idm_hd,  jdm_hd)

   uscx_hd = fltarr(idm_hd,  jdm_hd)
   uscy_hd = fltarr(idm_hd,  jdm_hd)

   vscx_hd = fltarr(idm_hd,  jdm_hd)
   vscy_hd = fltarr(idm_hd,  jdm_hd)

   cori_hd = fltarr(idm_hd,  jdm_hd)
   pasp_hd = fltarr(idm_hd,  jdm_hd)

   
   ;; Longitude/Latitude for each point 
   ;; longitude
   vlon_hd = plon_hd
   FOR i = 1, idm-1 DO BEGIN 
      qlon_hd(i, *) = 0.5*(plon_hd(i, *)+plon_hd(i-1, *))
   ENDFOR 
   diff = plon_hd(2, 0)-plon_hd(1, 0)
   qlon_hd(0, *) = plon_hd(0, *)-0.5*diff
   ulon_hd = qlon_hd
   
   ;; latitude
   ulat_hd = plat_hd  
   FOR j = 1, jdm-1 DO BEGIN 
      qlat_hd(*, j) = 0.5*(plat_hd(*, j)+plat_hd(*, j-1))
   ENDFOR 
   diff = plat_hd(0, 2)-plat_hd(0, 1)
   qlat_hd(*, 0) = plat_hd(*, 0)-0.5*diff
   vlat_hd = qlat_hd

   ;; simplified grid with prescribed dx and dy
   pscx_hd(*, *) = dx &  pscy_hd(*, *) = dy
   uscx_hd(*, *) = dx &  uscy_hd(*, *) = dy
   vscx_hd(*, *) = dx &  vscy_hd(*, *) = dy
   qscx_hd(*, *) = dx &  qscy_hd(*, *) = dy

   ;; Coriolis
   beta = 2.E-11
   FOR  j= 0,jdm_hd-1 DO BEGIN 
      FOR  i= 0,idm_hd-1 DO BEGIN 
         cori_hd(i, j)=.93e-4+float(j-(jdm_hd-2)/2)*dx*beta
      ENDFOR 
   ENDFOR 
   

   ;; pang ( p-angle  for rotated grid)
   ;; here uniform or mercator so pang = 0.


   ;; pasp= paspect: pscx/pscy
   
   FOR  j= 0,jdm_hd-1 DO BEGIN 
      FOR  i= 0,idm_hd-1 DO BEGIN 
         pasp_hd(i,j) = pscx_hd(i,j)/pscy_hd(i,j)
      ENDFOR 
   ENDFOR 
   
   
   ;; Writing new grid file
   write_grid_hycom, idm_hd, jdm_hd, io, file_grid_new, plon_hd, plat_hd, ulon_hd,  ulat_hd,  vlon_hd, vlat_hd, qlon_hd, qlat_hd, pang_hd,pscx_hd, pscy_hd, qscx_hd,  qscy_hd, uscx_hd, uscy_hd,  vscx_hd, vscy_hd, cori_hd, pasp_hd

   print,  'Writing grid file done '



   stop

END 
