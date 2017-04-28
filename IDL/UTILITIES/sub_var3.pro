PRO sub_var3, im, jm, km, fName, num2, num3, ivar, var


   ;; Script to read a 3-D variable from hycom archive file
   ;;----------------------------------------------------------------------
   ;;   Input:
   ;;        fName     : file name '.a'
   ;;    im, jm, km    : dimension for full domain  
   ;;        num2/3    : number of 2d/3d variables -- .b 
   ;;        ivar      : index of the 3d variable  (starting from 1st num3)
   ;;----------------------------------------------------------------------
   ;; 
   ;; A. Bozec May, 2013

   close,/all 

   ;; Dimensions of the domain
   idm = im
   jdm = jm
   idm1 = float(im)
   ijdm = idm1*jdm
   kdm = km

   tt = fltarr(idm, jdm, num2)
   tt1 = fltarr(idm, jdm, num3)

   si_2d = num2
   si_3d = num3

   ;; NPAD size
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   toto = fltarr(npad)

   ;; Directory and file
   iodir = fName

   ;; Tab Declaration
   var = fltarr(idm, jdm, kdm)

   ;; Open file
   openr, 1, iodir, /swap_endian

   ;; Read 2D fields
   FOR jk = 0, si_2d-1 DO BEGIN 
      readu, 1, rr2,  toto
      FOR j = 0, jdm-1 DO BEGIN 
         FOR i = 0, idm-1 DO tt(i, j, jk) = rr2(j*idm1+i)
      ENDFOR 
   ENDFOR 
     
   ;; Read 3D fields 
   FOR ll = 0, kdm-1 DO BEGIN 
      FOR jk = 0, si_3d-1 DO BEGIN 
         readu, 1, rr2,  toto
         FOR j = 0, jdm-1 DO BEGIN 
            FOR i = 0, idm-1 DO tt1(i, j, jk) = rr2(j*idm1+i)
         ENDFOR 
      ENDFOR 
      var(*, *, ll) = tt1(*, *, ivar-1) 
   ENDFOR 
   close, 1
   ;; Put nan on missing values
   var(where(var GT 1e20)) = !values.f_nan

   
END 

