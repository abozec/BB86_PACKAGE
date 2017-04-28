FUNCTION sofsig, T, R, sig

   IF (sig EQ 2) THEN BEGIN 
;; --- coefficients for sigma-2 (based on Brydon & Sun fit)
      C1= 9.77093E+00 & C2=-2.26493E-02
      C3= 7.89879E-01 & C4=-6.43205E-03
      C5=-2.62983E-03 & C6= 2.75835E-05
      C7= 3.15235E-05
   ENDIF 

   IF (sig EQ 0) THEN BEGIN 
;; --- coefficients for sigma-0 (based on Brydon & Sun fit)
      C1=-1.36471E-01 &  C2= 4.68181E-02
      C3= 8.07004E-01 &  C4=-7.45353E-03
      C5=-2.94418E-03 &  C6= 3.43570E-05
      C7= 3.48658E-05
   ENDIF 
;
; --- salinity (mil) as a function of sigma and temperature (deg c)
      sof=(R-C1-T*(C2+T*(C4+C6*T)))/(C3+T*(C5+C7*T))

return,sof
end
