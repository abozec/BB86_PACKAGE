PRO read_depth_hycom, im, jm, file, bathy

   ;; Script to read the HYCOM bathymetry file
   ;; A. Bozec Aug, 2011

   close,/all 
   ;; Dimensions of the domain
   idm = im
   jdm = jm
   idm1 = float(idm)
   ijdm = idm1*jdm


   ;; NPAD size and Tabs definition
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   toto = fltarr(npad)
   bathy = fltarr(idm, jdm)


   ;; Grid Directory and file 
   file1 = file

   ;; READING the file
   openu, 1, file1, /swap_endian
   readu, 1,  rr2
   FOR j = 0, jdm-1 DO BEGIN 
      FOR i = 0, idm-1 DO bathy(i, j) = rr2(j*idm1+i)
   ENDFOR 
   close, 1

   ;; Mask the bathymetry
   bathy(where(bathy GT 1e20)) = !values.f_nan

END
