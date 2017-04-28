PRO write_depth_bb86

   ;; include the Utilities functions
   iodir = '/Net/yucatan/abozec/BB86_PACKAGE/IDL/' ;; where you are
   !path=expand_path('+/'+iodir+'UTILITIES')+':'+expand_path('+'+!dir)

   close, /all
   ;; PATH
   io = iodir+'/../topo/'
   file_bat_new = 'depth_BB86_01' ;; !! without .a or .b !!
   pl = 0  ;; 1 or 0 for plot or not

   ;; size of the domain
   idm = 101
   jdm = 101 

   ;; definition of the bathymetry
   depth = 5000. ;; constant bathy everywhere

   ;; closed boundary (latbdy=0) or open boundaries (latbdy=1) or cyclic (latbdy =2)
   latbdy = 0

   ;;;;;;; END of the USER inputs ;;;;;;;;;;;;;;;

   idm_hd = idm
   jdm_hd = jdm
   vmiss = 2.^100  ;; HYCOM missing values

   ;; get the Depth 
   bathy_hd = fltarr(idm, jdm)
   bathy_hd(*, *) = depth
   print, 'Depth Ok'

   ;; Plot
   IF (pl EQ 1) THEN BEGIN 
      device, decomposed = 0
      
      ;; defined colors
      levels = 12
      Minss = 4000. &  Maxss = 6000.
      step = (Maxss - Minss) / levels
      num_level = IndGen(levels) * step + Minss 
      loadct, 33, ncolors = levels
      
      ;; let's make a plot
      window, 1, xsize = 800, ysize = 800
      contour, bathy_hd, xstyle = 1, ystyle = 1, levels = num_level, c_colors = indgen(levels), /fill, $
        Color = !P.Background, Background = !P.color, /follow, position = [0.15, 0.15, 0.85, 0.85]
      colorbar2, ncolors = levels, divisions = levels, color = 0, range = [Minss, Maxss], $
       position = [0.17, 0.08, 0.83, 0.10], format = '(f7.1)' 

      device, decomposed = 1
   ENDIF 
      
   ;; Mask by missing values if any
   ind = where(bathy_hd EQ 0.)
   IF (ind NE [-1]) THEN bathy_hd (ind) = vmiss


   ;; Get the edge of the domain right
   CASE (latbdy) OF 
      0: BEGIN          
         bathy_hd(*, jdm_hd-1) = vmiss 
         bathy_hd(idm_hd-1, *) = vmiss 
         bathy_hd(*, 0) = vmiss 
         bathy_hd(0, *) = vmiss
      END 
      1: BEGIN 
         bathy_hd(*, jdm_hd-1) = vmiss 
         bathy_hd(idm_hd-1, *) = vmiss 
      END 
      2: BEGIN 
         bathy_hd(*, jdm_hd-1) = vmiss 
         bathy_hd(*, 0) = vmiss 
      END 
   ENDCASE


   ;; write the file
   write_depth_hycom,  idm_hd, jdm_hd, io+file_bat_new, bathy_hd
   print,  'Writing depth file done '

   stop
END 
