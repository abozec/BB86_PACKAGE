PRO write_depth_hycom, im, jm, file, bathy


   close,/all 
   ;; Dimensions of the domain
   idm = im
   jdm = jm
   idm1 = float(idm)
   ijdm = idm1*jdm

   ;; NPAD size
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   toto = fltarr(npad)

   ;; mask for min max
   tt = bathy
   ind = where(tt GT 1e20)
   IF (ind NE [-1] ) THEN tt(ind) = !values.f_nan
   
   file1 = file

   ;; READING the file
   openw, 1, file+'.a', /swap_endian
   FOR j = 0, jdm-1 DO BEGIN
      FOR i = 0, idm-1 DO rr2(j*idm1+i) = bathy(i, j) 
   ENDFOR
   writeu, 1, rr2,  toto
   close, 1

   ;; REWRITE .B
   
    openw,1,file+'.b'
    printf, 1, format = '(A10)', 'Bathymetry'
    printf, 1, format = '(A5,1x,i4,1x,i3)', 'i/jdm ', idm, jdm
    printf, 1, format = '(A1)', ' '
    printf, 1, format = '(A1)', ''
    printf, 1, format = '(A1)', ''
    printf, 1, format = '(A18,1x,f8.3,2x,f8.3)','min,max depth =   ', min(tt, /nan), max(tt, /nan)
         
    close,1
   




END
