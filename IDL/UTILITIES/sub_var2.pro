PRO sub_var2, im, jm, fName, num2, ivar, var


   ;; Script to read a 2-D variable from hycom archive file
   ;;----------------------------------------------------------------------
   ;;   Input:
   ;;        fName     : file name '.a'
   ;;       im, jm     : dimension for full domain  
   ;;        num2      : number of 2d variables -- .b 
   ;;        ivar      : index of the 2d variable  
   ;;----------------------------------------------------------------------
   ;; 
   ;; A. Bozec May, 2013

   close,/all 

   ;; Dimensions of the domain
   idm = im
   jdm = jm
   idm1 = float(im)
   ijdm = idm1*jdm
   
   tt = fltarr(idm, jdm, num2)

   si_2d = num2

   ;; NPAD size
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   toto = fltarr(npad)

   ;; Directory and file
   iodir = fName

   ;; Tab Declaration
   var = fltarr(idm, jdm)

   ;; Open file
   openr, 1, iodir, /swap_endian

   ;; Read 2D fields
   FOR jk = 0, si_2d-1 DO BEGIN 
      readu, 1, rr2,  toto
      FOR j = 0, jdm-1 DO BEGIN 
         FOR i = 0, idm-1 DO tt(i, j, jk) = rr2(j*idm1+i)
      ENDFOR 
   ENDFOR 
      var(*, *) = tt(*, *, ivar-1) 
     
   close, 1
   
   ;; Put nan on missing values
   var(where(var GT 1e20)) = !values.f_nan
   
END 

