FUNCTION  tofsig, S, R, sig

   IF (sig EQ 2.) THEN BEGIN 
;; --- coefficients for sigma-2 (based on Brydon & Sun fit)
      C1= 9.77093E+00 & C2=-2.26493E-02
      C3= 7.89879E-01 & C4=-6.43205E-03
      C5=-2.62983E-03 & C6= 2.75835E-05
      C7= 3.15235E-05
   ENDIF 
   IF (sig EQ 0.) THEN BEGIN 
;; --- coefficients for sigma-0 (based on Brydon & Sun fit)
      C1=-1.36471E-01 &  C2= 4.68181E-02
      C3= 8.07004E-01 &  C4=-7.45353E-03
      C5=-2.94418E-03 &  C6= 3.43570E-05
      C7= 3.48658E-05
   ENDIF 


      AZERO  =0.d
      AHALF  =1/2.d
      ATHIRD =1/3.d
      A1P5   =3.0/2d
;; --- auxiliary statements for finding root of 3rd degree polynomial
      A0=(C1+C3*S -R)/C6
      A1=(C2+C5*S)/C6
      A2=(C4+C7*S)/C6
      CUBQ=ATHIRD*A1-(ATHIRD*A2)^2
      CUBR=ATHIRD*(AHALF*A1*A2-A1P5*A0) -(ATHIRD*A2)^3
;; --- if q**3+r**2>0, water is too dense to yield real root at given
;; --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
;; --- lowering sigma until a double real root is obtained.
      CUBAN=ATHIRD*ATAN(SQRT(MAX([AZERO,-(CUBQ^3+CUBR^2)])),CUBR)
      CUBRL=SQRT(-CUBQ)*COS(CUBAN)
      CUBIM=SQRT(-CUBQ)*SIN(CUBAN)

; --- temp (deg c) as a function of sigma and salinity (mil)
      TOF=-CUBRL+SQRT(3.)*CUBIM-ATHIRD*A2



return,TOF
end
