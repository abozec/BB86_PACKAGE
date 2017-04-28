PRO sub_var2, im, jm, fName, num2, ivar, var, ii, jj


   ;; Script to read a 2-D variable from hycom archive file
   ;;----------------------------------------------------------------------
   ;;   Input:
   ;;        fName     : file name '.a'
   ;;       im, jm     : dimension for full domain  
   ;;        num2      : number of 2d variables -- .b 
   ;;        ivar      : index of the 2d variable  
   ;;        [ii jj]   : optional subset domain
   ;;----------------------------------------------------------------------
   ;; 
   ;; A. Bozec Feb, 2014

   compile_opt idl2, strictarrsubs

   close,/all 

   ;; Dimensions of the domain
   idm = im
   jdm = jm
   ijdm = idm*jdm
   si_2d = num2

   ;; NPAD size
;   npad=4096 - ijdm MOD 4096
   npad = ceil(double(idm)*jdm/4096)*4096
   npad=double(npad)
   IF (n_elements(ii) NE 0.) THEN BEGIN     
      ni = n_elements(ii)
      nj = n_elements(jj)
      var = dblarr(ni, nj)
      offs = ivar-1.D
      for j = 0, nj-1 do begin 
         dstart = (offs*npad + jj[j]*idm + ii[0])*4.
         dstart = long64(dstart)
         result = read_binary(fName, data_start = dstart, data_type = 4, data_dims = [ni], endian = 'big')
         var[*, j] = result[*]
      endfor 
      
      ;; Put nan on missing values
      var[where(var GT 1e20)] = !values.f_nan
   ENDIF ELSE BEGIN

      ;; find data offset
      dstart = (ivar-1.D)*npad*4.D
      dstart=long64(dstart)
      
      ;; Tab Declaration
      var = fltarr(idm, jdm) + 100.
      result = read_binary(fName, data_start = dstart, data_type = 4, data_dims = [ijdm], endian = 'big')
      FOR j = 0L, jdm-1L DO BEGIN 
         FOR i = 0L, idm-1L DO var[i, j] = result[j*idm+i]
      ENDFOR 
      
      ;; Put nan on missing values
      index=where(var GT 1e20)
      if (index NE [-1] ) then var[index] = !values.f_nan
   ENDELSE 
END 

