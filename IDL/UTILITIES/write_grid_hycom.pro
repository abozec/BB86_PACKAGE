PRO  write_grid_hycom, im, jm,  io, file, plon, plat, ulon,  ulat,  vlon, vlat, qlon, qlat, pang, pscx, pscy, qscx,  qscy, uscx, uscy,   vscx, vscy, cori, pasp


   close,/all 

;   im = 271 &  jm = 193
;   im = 541 &  jm = 385
;   io = '/Users/abozec/IDL/IDL_DATA/GOMl0.08/topo/'
;   file = 'regional.grid.a'
;   file = 'regional.grid.GOMl0.04.a'

   ;; Dimensions of the domain
   idm = im
   jdm = jm
   idm1 = float(im)
   ijdm = idm1*jdm

   tt2 = fltarr(idm, jdm, 19)

   ;; NPAD size
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   if (npad NE 4096) then toto = fltarr(npad)

   ;; Grid Directory and file 
   iodir = io
;   file = 'regional.grid.a'

   tt2(*, *, 0) = plon(*, *)
   tt2(*, *, 1) = plat(*, *)

   tt2(*, *, 2) = qlon(*, *)
   tt2(*, *, 3) = qlat(*, *)

   tt2(*, *, 4) = ulon(*, *)
   tt2(*, *, 5) = ulat(*, *)

   tt2(*, *, 6) = vlon(*, *)
   tt2(*, *, 7) = vlat(*, *)

   tt2(*, *, 8) = pang(*, *)

   tt2(*, *, 9)  = pscx(*, *)
   tt2(*, *, 10) = pscy(*, *)

   tt2(*, *, 11) = qscx(*, *)
   tt2(*, *, 12) = qscy(*, *)

   tt2(*, *, 13) = uscx(*, *)
   tt2(*, *, 14) = uscy(*, *)

   tt2(*, *, 15) = vscx(*, *)
   tt2(*, *, 16) = vscy(*, *)

   tt2(*, *, 17) = cori(*, *)
   tt2(*, *, 18) = pasp(*, *)

   openw, 1, io+file+'.a', /swap_endian
   FOR jk = 0, 18 DO BEGIN 
      FOR j = 0, jm-1 DO BEGIN 
         FOR i = 0, im-1 DO  rr2(j*idm1+i) =tt2(i, j, jk)
      ENDFOR 
         if (npad NE 4096) then begin 
            writeu, 1, rr2,  toto
         endif else begin 
            writeu, 1, rr2
         endelse 
   ENDFOR 
   close, 1
  

   ;; Writing .b file
   openw,1,io+file+'.b'
   printf, 1, format = '(i8,1x,A33)', im,    '''idm   '' = longitudinal array size'
   printf, 1, format = '(i8,1x,A32)', jm,    '''jdm   '' = latitudinal array size'
   printf, 1, format = '(i8,1x,A53)', 0,     '''mapflg'' = map flag (0=mercator,10=panam,12=ulon-panam)'
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','plon:  min,max = ', min(tt2(*,*,0)), max(tt2(*,*,0))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','plat:  min,max = ', min(tt2(*,*,1)), max(tt2(*,*,1))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','qlon:  min,max = ', min(tt2(*,*,2)), max(tt2(*,*,2))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','qlat:  min,max = ', min(tt2(*,*,3)), max(tt2(*,*,3))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','ulon:  min,max = ', min(tt2(*,*,4)), max(tt2(*,*,4))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','ulat:  min,max = ', min(tt2(*,*,5)), max(tt2(*,*,5))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','vlon:  min,max = ', min(tt2(*,*,6)), max(tt2(*,*,6))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','vlat:  min,max = ', min(tt2(*,*,7)), max(tt2(*,*,7))
   printf, 1, format = '(A18,1x,E14.5,1x,E14.5)','pang:  min,max = ', min(tt2(*,*,8)), max(tt2(*,*,8))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','pscx:  min,max = ', min(tt2(*,*,9)), max(tt2(*,*,9))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','pscy:  min,max = ', min(tt2(*,*,10)), max(tt2(*,*,10))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','qscx:  min,max = ', min(tt2(*,*,11)), max(tt2(*,*,11))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','qscy:  min,max = ', min(tt2(*,*,12)), max(tt2(*,*,12))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','uscx:  min,max = ', min(tt2(*,*,13)), max(tt2(*,*,13))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','uscy:  min,max = ', min(tt2(*,*,14)), max(tt2(*,*,14))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','vscx:  min,max = ', min(tt2(*,*,15)), max(tt2(*,*,15))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','vscy:  min,max = ', min(tt2(*,*,16)), max(tt2(*,*,16))
   printf, 1, format = '(A18,1x,E14.5,1x,E14.5)','cori:  min,max = ', min(tt2(*,*,17)), max(tt2(*,*,17))
   printf, 1, format = '(A18,1x,f14.5,1x,f14.5)','pasp:  min,max = ', min(tt2(*,*,18)), max(tt2(*,*,18))
   close,1


END 

