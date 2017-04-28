function sig=sigma0_hycom(T,S)
%%
%% --- coefficients for sigma-0 (based on Brydon & Sun fit)
      C1= -1.36471E-01 ;
      C2=  4.68181E-02 ;
      C3=  8.07004E-01 ;  
      C4= -7.45353E-03 ;
      C5= -2.94418E-03 ;
      C6=  3.43570E-05 ;
      C7=  3.48658E-05 ;

%% --- sigma-theta as a function of temp (deg c) and salinity (mil)
      sig=(C1+C3*S+T.*(C2+C5*S+T.*(C4+C7*S+C6*T)))

