PRO write_relax_hycom, im, jm, km, io, file, field, d, res


   ;; field='intf'
   ;; field='temp'
   ;; field='saln'
   ;; d(kdm); target density
 
   CASE (field) OF 
      'intf': BEGIN 
         long_name = 'Interface Depths'
         short_name = 'int'
      END 
      'temp': BEGIN 
         long_name = 'Potential Temperature'
         short_name = 'tem'
      END 
      'saln': BEGIN 
         long_name = 'Salinity'
         short_name = 'sal'
      END 

   ENDCASE 

   close,/all 

   ;; Dimensions of the domain
   kdm = km
   idm = im
   jdm = jm
   idm1 = float(idm)
   ijdm = idm1*jdm

   ;; NPAD size
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   if (npad NE 4096) then toto = fltarr(npad)

   ;; Grid Directory and file 
   file1 = io+file

   ;; WRITING .a binary file
   openw, 1, file1+'.a', /swap_endian
   For t=0,11 do begin 
      FOR k = 0, kdm-1 DO BEGIN 
         FOR j = 0, jdm-1 DO BEGIN 
            FOR i = 0, idm-1 DO rr2(j*idm1+i) = res(i, j, k, t)
         ENDFOR 
         if (npad NE 4096) then begin 
            writeu, 1, rr2,  toto
         endif else begin 
            writeu, 1, rr2
         endelse 
      ENDFOR
   ENDFOR   
   close, 1


   ;; WRITING .b text file
   mo = findgen(12)+1
   months = string(mo, format ='(i2.2)')
   l = findgen(kdm)+1
   layer = string(l, format ='(i2.2)')
   density = d

   ;; mask the field by nan
   index = where(res GT 1e20)
   IF (index NE [-1]) THEN res(index) = !values.f_nan
   
   openw,1,file1+'.b'
   printf, 1, format = '(A)', ''+long_name 
   printf, 1, format = '(A)', ''
   printf, 1, format = '(A)', ''
   printf, 1, format = '(A)', ''
   printf, 1, format = '(A5,1x,i3,1x,i3)', 'i/jdm ', idm, jdm

   
   FOR m = 0, 11 DO BEGIN 
      FOR k = 0, kdm-1 DO  BEGIN 
         printf, 1, format = '(A,1x,i2.2,1x,i2.2,2x,f6.3,2x,E14.7,2x,E14.7)',' '+short_name+': month,layer,dens,range =', m+1, k+1, density(k), min(res(*, *, k), /nan), max(res(*, *, k), /nan)
         
      ENDFOR 
   ENDFOR 
   close,1



   

END 



