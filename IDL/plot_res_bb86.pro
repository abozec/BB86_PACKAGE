PRO plot_res_bb86

   close, /all
   ;; include the Utilities function
   iodir = '/Net/yucatan/abozec/BB86_PACKAGE/IDL/' ;; where you are
   !path=expand_path('+/'+iodir+'UTILITIES')+':'+expand_path('+'+!dir)

   ;; PATH
   io = '/Net/yucatan/abozec/BB86_PACKAGE/'
   ;; domain
   idm = 101 &  jdm = 101       ;; size of the domain
   kdm = 2                      ;; number of vertical layer in BB86
   tdm = 1800                   ;; number of time-stamp in bb86
   dp0 = 500.                   ;; thickness of the 1st layer  (m) 
   eps = 0.0001                 ;; epsilon to avoid dividing by 0.
   tplot1 = 1800 &  tplot2 = 1800 ;; time-stamp to plot (starts from 1.)

   ;; postscript parameters
   PS = 1                       ;; save as postscript
   key_portrait = 0
   page_size = [21.5900, 27.9400] ;; letter size
   file_ps = '../PS/uv_dp_d'+string(tplot2, format = '(i4.4)')+'_bb86-hycom.ps'

   ;; constants
   rho = 1000.   ;; reference density
   g = 9.806     ;; gravity
   

   ;; Read grid
   file_grid = 'regional.grid.BB86.a'
   read_grid_hycom, idm, jdm, io+'topo/', file_grid, plon, plat


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;; READ the HYCOM files 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   io_hycom = io+'expt_01.0/data/output/'  

   uhyc = fltarr(idm, jdm, kdm+1, tdm) ;; 1 more layer for the hybrid!!
   vhyc = fltarr(idm, jdm, kdm+1, tdm) 
   dphyc = fltarr(idm, jdm, kdm+1, tdm) 

   FOR t = tplot1, tplot2 DO BEGIN 
      CASE 1  OF 
         (t LE 344): BEGIN 
            year = '0001'
            day = string(t+16, format = '(i3.3)') ;; 16 u and v == 0. so we start at day 17.
         END 
         (t GT 344) AND (t LE 704) : BEGIN 
            year = '0002'
            day = string(t-344, format = '(i3.3)')
         END 
         (t GT 704) AND (t LE  1064) : BEGIN 
            year = '0003'
            day = string(t-704, format = '(i3.3)')
         END 
         (t GT 1064) AND (t LE  1424) : BEGIN 
            year = '0004'
            day = string(t-1064, format = '(i3.3)')
         END 
         (t GT 1424) AND (t LE  1784) : BEGIN 
            year = '0005'
            day = string(t-1424, format = '(i3.3)')
         END 
         (t GT 1784) AND (t LE  2144) : BEGIN 
            year = '0006'
            day = string(t-1784, format = '(i3.3)')
         END 
      ENDCASE 
      file = 'archv.'+year+'_'+day+'_00.a'
      print,  file, t

      ;; archive files 2D variable units (see archv.*.b for the list of available var.)
      ;;      montg1  : *1./g to get in m
      ;;      srfhgt  : *1./g to get in m
      ;;      surflx  : W/m2
      ;;      salflx  : kg/m2/s
      ;;      bl_dpth : *1./(rho*g) to get in m
      ;;      mix_dpth: *1./(rho*g) to get in m
      ;;      u_btrop : m/s
      ;;      v_btrop : m/s

      ;; extract barotropic velocities (ubarot & vbarot)
      num2 = 8     ;; number of 2D variables in the archive file
      ivar_ub = 7  ;; index of ubaro (7th 2D variable)
      ivar_vb = 8  ;; index of vbaro (8th 2D variable)
      sub_var2, idm, jdm, io_hycom+file, num2, ivar_ub, ubaro    
      sub_var2, idm, jdm, io_hycom+file, num2, ivar_vb, vbaro    

      
      ;; archive files 3D variable units (see archv.*.b for the list of available var.)
      ;;      u-vel.  : m/s
      ;;      v-vel.  : m/s
      ;;      thknss  : *1./(rho*g) to get in m
      ;;      temp    : C
      ;;      salin   : psu
      
      ;; extract velocities (ubaroc & vbaroc)
      num3 = 5   ;; number of 3D variables in the archive file
      ivar_u = 1 ;; index of u baroclinic
      ivar_v = 2 ;; index of v baroclinic
      sub_var3, idm, jdm, kdm+1, io_hycom+file, num2, num3, ivar_u, ubac    
      sub_var3, idm, jdm, kdm+1, io_hycom+file, num2, num3, ivar_v, vbac    

      ;; Get utot & vtot
      FOR k = 0, kdm DO uhyc(*, *, k, t-1) = ubaro(*, *)+ubac(*, *, k)
      FOR k = 0, kdm DO vhyc(*, *, k, t-1) = vbaro(*, *)+vbac(*, *, k)


      ;; extract the layer thickness
      ivar_dp = 3 
      sub_var3, idm, jdm, kdm+1, io_hycom+file, num2, num3, ivar_dp, dp
      dphyc(*, *, *, t-1) = dp/(rho*g) 


      ;; put the u on the p-grid

      ;; get the umask
      index = where(finite(uhyc) EQ 0.)
      maskt = fltarr(idm, jdm, kdm+1, tdm)+1
      if (index NE [-1]) then maskt(index) = 0.
      maxval = maskt+shift(maskt, -1, 0, 0, 0)
      maxval(where(maxval EQ 0.)) = !values.f_nan

      ;; get the average value of u at the p-point
      uthyc = (uhyc+shift(uhyc, -1, 0, 0, 0))/maxval
      uthyc(idm-1, *, *, *) = uthyc(idm-2, *, *, *)

      ;; put the v on the p-grid

      ;; get the vmask
      index = where(finite(vhyc) EQ 0.)
      maskt = fltarr(idm, jdm, kdm+1, tdm)+1
      if (index NE [-1]) then maskt(index) = 0.
      maxval = maskt+shift(maskt, 0, -1, 0, 0)
      maxval(where(maxval EQ 0.)) = !values.f_nan

      ;; get the average value of v at the p-point
      vthyc = (vhyc+shift(vhyc, 0, -1, 0, 0))/maxval
      vthyc(*, jdm-1, *, *) = vthyc(*, jdm-2, *, *)
     

   ENDFOR  

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;; plot 
   for t = tplot1,  tplot2 do begin

      ;; norme, min and max of the plots
      normeref_hyc  = 0.05    ;; vector norm of HYCOM-BB86 (in m)
      onevectout = 2          ;; plot one vector out of 'onevectorout'
      min_dp = -450. &  max_dp = 450. ;; layer thickness anomaly in m (min and max)


      device, true_color = 24, decomposed = 0   
      IF (PS EQ 0 ) THEN BEGIN 
         dimensions = get_screen_size(RESOLUTION=resolution)
         coef = floor(1./resolution[0])
         windowsize_scale = 1.
         coef = windowsize_scale * coef
         
         mipgsz = min(page_size, max = mapgsz)

         xsize = coef * (mipgsz*key_portrait + mapgsz*(1-key_portrait))
         ysize = coef * (mipgsz*(1-key_portrait) + mapgsz*key_portrait)
         
         window, 1, xsize = xsize, ysize = ysize 
         foreground = !P.Background &  background = !P.Color
         char = 1.
      ENDIF 
      IF (PS EQ 1 ) THEN BEGIN 
         IF !d.name EQ 'PS' then device,/close
         set_plot,'ps'
         device, /color, /helvetica, filename = file_ps $
          , LANDSCAPE = 1 - key_portrait, PORTRAIT = key_portrait $
          , xsize = max(page_size), ysize = min(page_size), xoffset = 0., yoffset = max(page_size) $
          , bits_per_pixel = 8 
         foreground = !P.Color &  background = !P.Background
         char = 0.75
      ENDIF 


      

      loadct, 39 ;; color palette


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; HYCOM-BB86 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; get a blank plot
      contour, reform(uthyc(*, *, 1, 0)), plon,  plat, xstyle = 1, ystyle = 1, $
       Color = foreground, Background = background,  position = [0.1, 0.6, 0.45, 0.90], $
       title = 'HYCOM-BB86 Velocity field  day = '+string(t, format = '(i4.4)'), /nodata, $
       xtitle = 'Longitude (E)',  ytitle = 'Latitude (N)',  charsize = char

      ;; add the velocities vectors
      ;; definition of ref vector (legend)
      uref = fltarr(2, 2) &  vref =fltarr(2, 2)
      uref(0, 0) = 1 

      ;; get the vectors to plot
      uu = reform(uthyc(0:idm-1:onevectout, 0:jdm-1:onevectout, 1, t-1)) 
      vv = reform(vthyc(0:idm-1:onevectout, 0:jdm-1:onevectout, 1, t-1)) 
      lon = reform(plon(0:idm-1:onevectout, 0))
      lat = reform(plat(0, 0:jdm-1:onevectout))

      ;; we normalize the vectors
      normeref = normeref_hyc ;; in m
      norme =sqrt(uu^2.+vv^2.)
      normemax = max(norme, /nan)/normeref

      ;; plotting
      velovect, uu, vv, lon, lat, /overplot, color = 0., length = 3*normemax ;; ref vector over 3*dx

      ;; add the ref vector in the legend take 3*dx (= 1.2º) as reference
      velovect, uref, vref,  [lon(3, 0), lon(4, 0)], [-1.2, 0.], /overplot, color = 0, length = 3
      xyouts, 1.2, -2.5, string(normeref, format = '(f4.2)')+' m', color = 0



      ;; plot the dp's 
      ;; defined colors
      levels = 40
      Minss = min_dp &  Maxss = max_dp
      step = (Maxss - Minss) / levels
      num_level = IndGen(levels) * step + Minss 
      loadct, 33, ncolors = levels, bottom = 1
      
      ;; get anomaly dp
      diff_dp = reform(dphyc(*, *, 1, t-1)) -(dp0-1.)
      diff_dp(where(finite(diff_dp) EQ 0)) = 0.
      contour, diff_dp, plon,  plat, xstyle = 1, ystyle = 1, levels = num_level, c_colors = indgen(levels), /fill, $
       Color = foreground, Background =  background, /follow, position = [0.5, 0.6, 0.85, 0.90], $
       title = 'HYCOM-BB86 thickness ano. day = '+string(t, format = '(i4.4)'), $
       xtitle = 'Longitude (E)',  ytitle = 'Latitude (N)', /noerase,  charsize = char
      contour, diff_dp, plon,  plat,levels = num_level, /overplot, color = 0

      colorbar2, ncolors = levels, divisions = levels/4. , color = 0, range = [Minss, Maxss], $
       position = [0.89, 0.60, 0.92, 0.9], format = '(f7.1)', /vertical,  charsize = char


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (PS EQ 1 ) THEN BEGIN 
         device, /close
         set_plot, 'x'
      ENDIF 

      device, decomposed = 1
   endfor 
   stop

END  
