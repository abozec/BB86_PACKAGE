PRO sub_var3, im, jm, km, fName, num2, num3, ivar, var, ii, jj, kk


   ;; Script to read a 3-D variable from hycom archive file
   ;;----------------------------------------------------------------------
   ;;   Input:
   ;;        fName     : file name '.a'
   ;;    im, jm, km    : dimension for full domain  
   ;;        num2/3    : number of 2d/3d variables -- .b 
   ;;        ivar      : index of the 3d variable  (starting from 1st num3)
   ;;       [ii jj kk] : dimension for sub domain [optional] 
   ;;----------------------------------------------------------------------
   ;; 
   ;; A. Bozec Feb, 2014

   compile_opt idl2, strictarrsubs

   ; im = 320 &  jm = 384
   ; km = 41
   ; fName = '/Users/abozec/IDL/IDL_DATA/POP1v6/expt_36.9/RAW_adjusted/369_archm.0003_m01.a'
   ; num2 = 29
   ; num3 = 7
   ; ivar = 6


   close,/all 

   ;; Dimensions of the domain
   idm = im
   jdm = jm
   ijdm = idm*jdm
   kdm = km
   

   ;; NPAD size
   npad = ceil(float(idm)*jdm/4096)*4096
;   print, n_elements(ii)
   IF (n_elements(ii) NE 0.) THEN BEGIN     
      ni = n_elements(ii)
      nj = n_elements(jj)
      nk = n_elements(kk)
      var = dblarr(ni, nj, nk)
      FOR k = 0, nk-1 DO BEGIN
         offs = (num2+num3*kk[k]+ivar-1.D)
         for j = 0, nj-1 do begin 
            dstart = (offs*npad + jj[j]*idm + ii[0])*4.
            dstart = long64(dstart)
            result = read_binary(fName, data_start = dstart, data_type = 4, data_dims = [ni], endian = 'big')
            var[*, j, k] = result[*]
         endfor 
      endfor 
      ;; Put nan on missing values
      var[where(var GT 1e20)] = !values.f_nan
   ENDIF ELSE BEGIN
      
      ;; Tab Declaration
      var = dblarr(idm, jdm, kdm)
      FOR k = 0, kdm-1 DO BEGIN 
         
         dstart = (num2+num3*k+ivar-1.D)*npad*4.D
         dstart = long64(dstart)
         result = read_binary(fName, data_start = dstart, data_type = 4, data_dims = [ijdm], endian = 'big')
         FOR j = 0, jdm-1 DO BEGIN 
            FOR i = 0, idm-1 DO var[i, j, k] = result[j*idm+i]
         ENDFOR 
         
      ENDFOR 
      ;; Put nan on missing values
      var[where(var GT 1e20)] = !values.f_nan
   ENDELSE 
   
END 

