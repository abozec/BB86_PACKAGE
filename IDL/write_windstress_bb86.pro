PRO write_windstress_bb86

   ;; include the Utilities function
   iodir = '/Net/yucatan/abozec/BB86_PACKAGE/IDL/' ;; where you are
   !path=expand_path('+/'+iodir+'UTILITIES')+':'+expand_path('+'+!dir)

   close, /all
   ;; PATH
   io = iodir+'/../force/'
   file_grid = 'regional.grid.BB86.a'

   ;; domain
   idm = 101 &  jdm = 101
   tdm = 12 ;; monthly files

   ;; name the new files
   file_E = 'forcing.tauewd.BB86'   ;; !! without .a or .b !!
   file_N = 'forcing.taunwd.BB86'


   ;; Read grid
   read_grid_hycom, idm, jdm, io+'../topo/', file_grid, plon, plat

   ;; Calculation of the analytical wind-stress (N/m2)
   ustress = fltarr(jdm)
   ustressb = fltarr(jdm)
   stressa = -1. 
   sconv = 1.e-1   ;; scale factor form dyn/cm2 to N/m2

   ;; BB86 formulation
   FOR j = 0, jdm-1 DO ustress(j) = stressa*cos(float(j)/float(jdm-1)*6.28318530718)*sconv


   ;; plot the stress 
   device, decomposed = 0
   zeros = fltarr(jdm)
   window, 1, xsize = 800, ysize = 800
   plot, ustress, reform(plat(0, *)), xstyle = 1, ystyle = 1,  ytitle = 'Latitude (N) ',  xrange = [-0.2, 0.2], xtitle = 'Taux (N/m2)',  $
    Color = !P.Background, Background = !P.color
   oplot, zeros,  reform(plat(0, *)),  linestyle = 1, color = 0


   device, decomposed = 1
   

   ;; Write the wind-stress files
   tte = fltarr(idm, jdm, tdm)
   ttn = fltarr(idm, jdm, tdm)
  
   ijdm = float(idm)*jdm
   idm1 = float(idm)
   npad=4096. - ijdm MOD 4096
   rr2 = fltarr(ijdm)
   toto = fltarr(npad)
 
   ;; we apply  a taux , no tauy
   FOR j = 0, jdm-1 DO tte(*, j, *) = ustress(j)
   
   ;; Taux
   openw, 1, io+file_E+'.a',  /swap_endian
   
   FOR ll = 0, tdm-1 DO BEGIN
      FOR j = 0, jdm-1 DO BEGIN 
         FOR i = 0, idm-1 DO  rr2(j*idm1+i) = tte(i, j, ll)
      ENDFOR 
      writeu, 1, rr2, toto
      print,  'll= ',  ll
   ENDFOR 
   close, 1
   
   ;; Tauy
   openw, 1, io+file_N+'.a',  /swap_endian
   FOR ll = 0, tdm-1 DO BEGIN
      FOR j = 0, jdm-1 DO BEGIN 
         FOR i = 0, idm-1 DO rr2(j*idm1+i) = ttn(i, j, ll) 
      ENDFOR 
      writeu, 1,  rr2, toto
      print,  'll= ',  ll
   ENDFOR 
   close, 1


   ;; create .b file 
   mo = findgen(12)+1
   months = string(mo, format ='(i2.2)')
   time = findgen(12)+1
   span = 0.250
   

;   openw,1,io+'tauewd_agul05_42S.b'
   openw,1,io+file_E+'.b'
   printf, 1, format = '(A22)', 'Analytical Wind-stress'
   printf, 1, format = '(A1)', ''
   printf, 1, format = '(A1)', ''
   printf, 1, format = '(A1)', ''
   printf, 1, format = '(A5,1x,i3,1x,i3)', 'i/jdm ', idm, jdm
   FOR m = 0, tdm-1 DO BEGIN 
      printf,1,format = '(A23,1x,i2.2,1x,E14.7,2x,E14.7)',' tauewd: month,range =', time(m), min(tte(*, *, m)), max(tte(*, *, m))
      
   ENDFOR 
;; tau_nwd: month,range = 01  -1.5917161E-01   1.5640029E-01
   
   close,1

   openw,1,io+file_N+'.b'
   printf, 1, format = '(A22)', 'Analytical Wind-stress'
   printf, 1, format = '(A1)', ''
   printf, 1, format = '(A1)', ''
   printf, 1, format = '(A1)', ''
   printf, 1, format = '(A5,1x,i3,1x,i3)', 'i/jdm ', idm, jdm
   FOR m = 0, tdm-1 DO BEGIN 
      printf,1,format = '(A23,1x,i2.2,1x,E14.7,2x,E14.7)',' taunwd: month,range =', time(m), min(ttn(*, *, m)), max(ttn(*, *, m))
      
   ENDFOR 
   
   close,1


stop

END 
