Forward_FUNCTION tofsig

PRO write_relax_bb86

   ;; include the Utilities functions
   iodir = '/Net/yucatan/abozec/BB86_PACKAGE/IDL/' ;; where you are
   !path=expand_path('+'+!dir)+':'+expand_path('+/'+iodir+'UTILITIES')

   close, /all
   ;; PATH
   io = iodir+'../relax/010/'
   file_topo = 'depth_BB86_01.a'
   file_int_new = 'relax_int_BB86' ;; !! without .a or .b !!
   file_tem_new = 'relax_tem_BB86'
   file_sal_new = 'relax_sal_BB86'
   pl = 0 ;; 1 or 0 for plot or not

   ;; size of the domain
   idm = 101
   jdm = 101 
   kdm = 3  ;; number of layers (first layer very thin)
   tdm = 12 ;; 12 month climatology

   ;; interface depth 
   id = fltarr(kdm)
   id = [0., 1.0, 500.] ;; first interface depth always 0.


   ;; target density
   ;; from bb86: rho=27.01037, rho=27.22136
   d = fltarr(kdm)
   d = [27.0100, 27.01037,27.22136]
   
   ;; salinity profile
   sa = fltarr(kdm)
   sa = [37., 37.,  37.]

   ;;;;;;; END of the USER inputs ;;;;;;;;;;;;;;;

   ;; constants
   rho = 1000.
   g = 9.806
   vmiss = 2.^100

   ;; temperature profile for sigma0
   sigma = 0
   te = fltarr(kdm)
   FOR k = 0, kdm-1 DO te(k) = tofsig(sa(k), d(k), sigma)
   print,  te

   
   ;; read bathy for mask
   read_depth_hycom, idm, jdm, io+'../../topo/'+file_topo,  bathy
   ind = where(bathy GT 1e20)
   IF (ind NE [-1]) THEN bathy(ind) = 0.
   ;; create mask
   mask2d = fltarr(idm, jdm)+1
   mask2d(ind) = 0.
   mask = fltarr(idm, jdm, kdm, tdm)
   FOR t = 0, tdm-1 DO BEGIN 
      FOR k = 0, kdm-1 DO mask(*, *, k, t) = mask2d
   ENDFOR 

   ;; interface depth files (first layer always the surface i.e. 0.)
   int = fltarr(idm, jdm, kdm, tdm)
   FOR k= 1, kdm-1 DO int(*, *, k, *) = rho*g * id(k) ;; 1st layer tiny to fake a two layer config

   ;; make sure that the interface depths are not lower than the bathy
   FOR t = 0, tdm-1 DO BEGIN 
      FOR k = 0, kdm-1 DO BEGIN 
         FOR i = 0, idm-1 DO BEGIN 
            FOR j = 0, jdm-1 DO BEGIN 
               IF (int(i, j, k, t)/9806. GT bathy(i, j)) THEN int(i, j, k, t) = bathy(i,j)*9806. 
            ENDFOR 
         ENDFOR 
      ENDFOR 
   ENDFOR 
   
   ;; mask 
   ind3d = where(mask EQ 0.)
   IF (ind3d NE [-1]) THEN int(ind3d) = vmiss

   ;; write the field in relax file
   write_relax_hycom, idm, jdm, kdm, io, file_int_new, 'intf', d, int
   print,  'Interface depth OK'

   ;; Temperature
   tem = fltarr(idm, jdm, kdm, tdm)
   FOR k = 0, kdm-1 DO tem(*, *, k, *) = te(k)
   
   ;; mask 
   IF (ind3d NE [-1]) THEN tem(ind3d) = vmiss

   ;; write the field in relax file
   write_relax_hycom, idm, jdm, kdm, io, file_tem_new, 'temp', d, tem
   print,  'Temperature OK'

   ;; Salinity
   sal = fltarr(idm, jdm, kdm, tdm)
   FOR k = 0, kdm-1 DO sal(*, *, k, *) = sa(k)

   ;; mask 
   IF (ind3d NE [-1]) THEN sal(ind3d) = vmiss

   ;; write the field in relax file
   write_relax_hycom, idm, jdm, kdm, io, file_sal_new, 'saln', d, sal
   print,  'Salinity OK'



   ;; Plot
   IF (pl EQ 1) THEN BEGIN 


      device, decomposed = 0
      
      ;; defined colors
      levels = 15
      Minss = 15. &  Maxss = 18.
      step = (Maxss - Minss) / levels
      num_level = IndGen(levels) * step + Minss 
      loadct, 33, ncolors = levels
      
      ;; let's make a plot
      window, 1, xsize = 800, ysize = 800
      contour, reform(tem(*, *, 1, 0)), xstyle = 1, ystyle = 1, levels = num_level, c_colors = indgen(levels), /fill, $
       Color = !P.Background, Background = !P.color, /follow, position = [0.15, 0.15, 0.85, 0.85], title = 'Temperature 2st layer'
      colorbar2, ncolors = levels, divisions = levels, color = 0, range = [Minss, Maxss], $
       position = [0.17, 0.08, 0.83, 0.10], format = '(f7.1)' 

      device, decomposed = 1
   ENDIF 

   stop

END 
