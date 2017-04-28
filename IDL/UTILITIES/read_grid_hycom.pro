PRO  read_grid_hycom, im, jm,  io, file, plon, plat, ulon,  ulat,  vlon, vlat, qlon, qlat, pang, pscx, pscy, qscx,  qscy, uscx, uscy,   vscx, vscy, cori, pasp

   ;; Script to read the HYCOM grid file
   ;; A. Bozec Aug, 2011

   close,/all 

   ;; Dimensions of the domain
   idm = im
   jdm = jm
   idm1 = float(im)
   ijdm = idm1*jdm

   tt = fltarr(idm, jdm, 19)

   ;; NPAD size
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   if (npad NE 4096) then toto = fltarr(npad)

   ;; Grid Directory and file 
   iodir = io
   plon = fltarr(idm,  jdm)
   plat = fltarr(idm,  jdm)

   qlon = fltarr(idm,  jdm)
   qlat = fltarr(idm,  jdm)

   ulon = fltarr(idm,  jdm)
   ulat = fltarr(idm,  jdm)

   vlon = fltarr(idm,  jdm)
   vlat = fltarr(idm,  jdm)

   pang= fltarr(idm,  jdm)
   pscx = fltarr(idm,  jdm)
   pscy = fltarr(idm,  jdm)

   qscx = fltarr(idm,  jdm)
   qscy = fltarr(idm,  jdm)

   uscx = fltarr(idm,  jdm)
   uscy = fltarr(idm,  jdm)

   vscx = fltarr(idm,  jdm)
   vscy = fltarr(idm,  jdm)

   cori = fltarr(idm,  jdm)
   pasp = fltarr(idm,  jdm)

   openu, 1, iodir+file, /swap_endian
   FOR jk = 0, 18 DO BEGIN 
      if (npad NE 4096) then begin 
         readu, 1, rr2,  toto
      endif else begin
         readu, 1, rr2
      endelse
      FOR j = 0, jdm-1 DO BEGIN 
         FOR i = 0, idm-1 DO tt(i, j, jk) = rr2(j*idm1+i)
      ENDFOR 
   ENDFOR 
   
   plon(*, *) = tt(*, *, 0)
   plat(*, *) = tt(*, *, 1)

   qlon(*, *) = tt(*, *, 2)
   qlat(*, *) = tt(*, *, 3)

   ulon(*, *) = tt(*, *, 4)
   ulat(*, *) = tt(*, *, 5)

   vlon(*, *) = tt(*, *, 6)
   vlat(*, *) = tt(*, *, 7)

   pang(*, *) = tt(*, *, 8)

   pscx(*, *) = tt(*, *, 9)
   pscy(*, *) = tt(*, *, 10)

   qscx(*, *) = tt(*, *, 11)
   qscy(*, *) = tt(*, *, 12)

   uscx(*, *) = tt(*, *, 13)
   uscy(*, *) = tt(*, *, 14)

   vscx(*, *) = tt(*, *, 15)
   vscy(*, *) = tt(*, *, 16)

   cori(*, *) = tt(*, *, 17)
   pasp(*, *) = tt(*, *, 18)

   close, 1 


END 
