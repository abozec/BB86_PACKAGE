FUNCTION sigma2_hycom, T, S


;;
;; --- coefficients for sigma-2 (based on Brydon & Sun fit)
      C1= 9.77093E+00 & C2=-2.26493E-02
      C3= 7.89879E-01 & C4=-6.43205E-03
      C5=-2.62983E-03 & C6= 2.75835E-05
      C7= 3.15235E-05

;; --- sigma-theta as a function of temp (deg c) and salinity (mil)
      SIG=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))

return,SIG
end
