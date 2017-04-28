FUNCTION spherdist, lon1, lat1, lon2, lat2

;;
;; --- ------------------------------------------------
;; --- Computes the distance between geo. pos.
;; --- lon1,lat1 and lon2,lat2. 
;; --- input is in degrees.
;;
;; --- output is real*4 for better global consistancy,
;; --- by truncating double precision roundoff errors.
;; --- real*4 is not in f90, but is widely supported.
;;
;; --- Based on m_spherdist.F90 from Geir Evanson.
;; --- ------------------------------------------------
;;
   invradian=double(0.017453292)
   rearth=double(6371001.0)       ; Radius of earth
   deg360 = double(360.)
   deg0 = double(0.)
   deg90 = double(90.)
   one = double(1.)
;;
;;
;;     ensure that spherdist(ax,ay,bx,by) == spherdist(bx,by,ax,ay)
;;
      dlon1 = double(lon1)
      dlon1 = dlon1 MOD deg360
      IF (dlon1 LT deg0) THEN dlon1 = dlon1 + deg360
      lat1 = double(lat1)
       
      dlon2 = double(lon2)
      dlon2 = dlon2 MOD deg360
      IF(dlon2 LT deg0) THEN dlon2 = dlon2 + deg360
      lat2 = double(lat2)
      
       CASE 1 OF 
       (lat1 LT lat2) : BEGIN 
         rlon1=dlon1*invradian            ;lon1 in rad
         rlat1=(deg90-lat1)*invradian     ;90-lat1 in rad 
         rlon2=dlon2*invradian            ;lon2 in rad
         rlat2=(deg90-lat2)*invradian     ;90-lat2 in rad 
       END
       (lat1 EQ lat2) AND (dlon1 LE dlon2) : BEGIN 
         rlon1=dlon1*invradian            ;lon1 in rad
         rlat1=(deg90-lat1)*invradian     ;90-lat1 in rad 
         rlon2=dlon2*invradian            ;lon2 in rad
         rlat2=(deg90-lat2)*invradian     ;90-lat2 in rad 
       END
       ELSE : BEGIN 
         rlon2=dlon1*invradian            ;lon1 in rad
         rlat2=(deg90-lat1)*invradian     ;90-lat1 in rad 
         rlon1=dlon2*invradian            ;lon2 in rad
         rlat1=(deg90-lat2)*invradian     ;90-lat2 in rad 
       END
       ENDCASE 

      
       
;;;
      x1= sin(rlat1)*cos(rlon1)        ;x,y,z of pos 1.
      y1= sin(rlat1)*sin(rlon1)
      z1= cos(rlat1) 
;;
      x2= sin(rlat2)*cos(rlon2)        ;x,y,z of pos 2.
      y2= sin(rlat2)*sin(rlon2)
      z2= cos(rlat2) 
;;
      dr=acos(min([one,x1*x2+y1*y2+z1*z2]))  ; Arc length
;;
      spher=dr*rearth
      spher = float(spher)
;;
return, spher



; lon1 = -60.08
; lon2 = -60.00
; lat1 = 35.08
; lat2 = 35.08
END 
