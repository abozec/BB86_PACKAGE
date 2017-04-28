                             PROGRAM RIGIDSF
c
C       ****************************************************************
C      *  THIS PROGRAM IS AN ADAPTATION OF THE BLECK AND BOUDRA (1986) *
C      *       ISOPYCNIC COORDINATE GENERAL CIRCULATION MODEL.         *
C      *  RECTANGULAR BASIN VERSION  WITH OR WITHOUT NO-SLIP BOUNDARY  *
C      *                      CONDITION                                *
C      *****************************************************************
C
C       II  =  NUMBER OF GRID POINTS IN x-DIRECTION                 C
C       JJ  =  NUMBER OF GRID POINTS IN y-DIRECTION                 C
C       KK  =  NUMBER OF LAYERS                                     C
C       II1 = II - 1                                                C
C       JJ1 = JJ - 1                                                C
C	DP          --> LAYER THICKNESS in pressure units
C	PBOT	    --> PRESSURE AT THE BOTTOM  [dyn/(cm**2)]
C	UFLUX       --> LAYER VELOCITY TIMES LAYER THICKNESS
C	UTROP       --> BAROTROPIC VELOCITY
C
C======================================================================C
C======================================================================C
C	RECALL THAT THE # OF GRID POINTS HAS TO BE A PRODUCT OF SMALL
C       PRIMES IN ORDER TO MAKE USE OF THE FFT TO SOLVE THE POISSON 
C	EQUATION.
C

      PARAMETER(II=101,JJ=II,KK=2,II2=II-2,JJ2=JJ-2,MAX=II2)
      PARAMETER(NW=INT(2.5*MAX+II))
C
C======================================================================C
C
      REAL U(II,JJ,2*KK) ,V(II,JJ,2*KK) ,DP(II,JJ,2*KK),MONTG(II,JJ,KK)
     .    ,P(II,JJ,KK+1) ,CORIO(II,JJ)  ,ABSVOR(II,JJ)   
     .    ,UTROP(II,JJ)  ,VTROP(II,JJ)  ,VORT(II,JJ)   ,PS(II,JJ)
     .    ,DEL2U(II,JJ)  ,STRESY(II,JJ) ,STRMF(II,JJ)  ,WSAVE(NW)
     .    ,EIG(II,JJ)    ,UFLUX(II,JJ)   
     .    ,VFLUX(II,JJ)  ,UOLD(II,JJ,KK), DEL2V(II,JJ)
     .    ,VOLD(II,JJ,KK),PBOT(JJ)      ,THETA(KK)     ,DP0(KK)
     .    ,DISVIS(KK)    ,STRESS(KK)    ,XCONT(KK)     ,ZETA(II,JJ)
     .    ,UTIL1(II,JJ)  ,DPM(II,JJ,KK)    ,UM(II,JJ,KK)
     .    ,UTIL2(II,JJ)  , VM(II,JJ,KK)
     .    ,STRESX(II,JJ), PU(II,JJ,KK+1), PV(II,JJ,KK+1), ULD(II,JJ,KK)
     .    ,VLD(II,JJ,KK)
c
c
C
      EQUIVALENCE   (VORT,STRMF)
C
C     ----------------------------------------------------------------------
C
c	-------------------------------------------------------------
c
C ---           FUNCTIONS TO MIMIC CRAY VECTOR MASK FUNCTIONS
c
              CVMGP(A,B,C)=A*(.5+SIGN(.5,C))+B*(.5-SIGN(.5,C))
c
c	-------------------------------------------------------------     
C
C
C======================================================================C
C                   DEFINE PHYSICAL PARAMETERS                         C
C======================================================================C
C
C --- SUBTRACT FOLLOWING VALUES FROM 1. TO OBTAIN LAYER SPECIFIC VOLUME
c
              DATA THETA/0.0263,0.0265/   
c
C
c	RHO=1./(1-THETA) GIVES THE RELATION BETWEEN RHO AND THETA
c
C
C --- HORIZONTAL GRID SCALE

               DATA SCALE/20.E5/             ! in cm
C
c
C --- PRESSURE AND DP's  ARE GIVEN IN UNITS OF dyn/(cm**2)
C --- EAST-WEST DIRECTION ARRAY OF BOTTOM PRESSURES

               DATA PBOT/II*5000.E5/
	       DATA PBOTTM/5000.E5/         

C
C --- GRAVITATIONAL ACCELERATION
               DATA G/980.6/
C
C
C --- VELOCITY TIME SMOOTHING WEIGHTS
               DATA WGT1,WGT2/.75,.125/
C
C
C --- DEPTH OF WIND STRESS FORCING AND LINEAR BOTTOM DRAG COEFFICIENT
               DATA PSTRES/50.E5/,DRAG/1.E-7/
C
C
C --- TWO LEAP FROG TIME LEVELS
               DATA M,N/2,1/
C
               DATA VISCOS/1.E6/
C
               DATA EPS/1.E5/,EPS1/1.E5/,EPSIL/1.E-9/
C
C --- 'EPS' IS USED WHEREVER LAYER THICKNESS MUST 
C ---       BE BOUNDED AWAY FROM ZERO.
C

c
C --- MODEL IS TO BE INTEGRATED FROM TIME STEP 'NSTEP1' TO 'NSTEP2'
C
        READ(5,*) NDAY1,NDAY2
          print*,'nday1= ', nday1, ' nday2 = ', nday2
C
	READ(5,*)DELT               ! TIME STEP OF INTEGRATION
C
        READ(5,*)(DP0(K),K=1,KK)
          print*,'dp0(1) = ',dp0(1),' dp0(2) = ',dp0(2)
C
        READ(5,*)STRESSA,IFREE
          print*,'stressa = ', stressa
	if(ifree.eq.1) print*,'Model run with free-slip boundary cond.'
	if(ifree.eq.2) print*,'Model run with no-slip boundary cond.'
C
C --------------------------------------------------------------------
C
      NSTEP1=NDAY1*86400./DELT + .0001
      NSTEP2=NDAY2*86400./DELT + .0001
C
C -------------------------------------------------------------------
C
	       SCALE2=SCALE**2	
               X=1./SCALE
               II1 = II-1
               JJ1 = JJ-1
               BETA=2.E-13
C
C
C
C======================================================================C
C                        INITIALIZATION                                C
C======================================================================C
C
C --- CALL INITIALIZATION SUBROUTINE FOR POISSON SOLVER:
C
          CALL POINIT(II-2,JJ-2,II,PBOT,EIG,WSAVE)                   
C
C --- SPECIFY INITIAL AND BOUNDARY CONDITIONS
C
C
      DO 10 I=1,II
          DO 10 J=1,JJ
C
C --- CORIOLIS PARAMETER
          CORIO(I,J)=.93E-4+FLOAT(J-JJ1/2)*SCALE*BETA
C
C======================================================================C
C                    WIND STRESS PROFILE                               C
C======================================================================C
C
          STRESX(I,J)=STRESSA*COS(FLOAT(J-1)/FLOAT(JJ1)
     .                *6.28318530718)*G/PSTRES
          STRESY(I,J)=0.
C
C======================================================================C
c	The wind stress is prescribed as a body force and is assumed 
c	to decrease linearly to the depth corresponding to PSTRES.
c	If the upper layer thickness (actually corresponding pressure)
c	becomes less than PSTRES, then the forcing is distributed 
c	to the second layer and so on according to the linear law.
C
C
          ABSVOR(I,J)=CORIO(I,J)
          UFLUX(I,J)=0.
          VFLUX(I,J)=0.
          STRMF(I,J)=0.
          MONTG(I,J,KK)=0.
          P(I,J,1)=0.
C
             DO 10 K=1,KK
             U(I,J,K   )=0.
             U(I,J,K+KK)=0.
             V(J,J,K   )=0.
             V(I,J,K+KK)=0.
             DP(I,J,K   )=DP0(K)
10           DP(I,J,K+KK)=DP0(K)
C======================================================================C
C======================================================================C
C
                   DELT1=DELT
C
        WRITE(6,551)II
551	FORMAT(2X'NUMBER OF GRID ELEMENTS IN EACH DIRECTION : 'I3)
       WRITE(6,552)SCALE
552	FORMAT(2X'GRID SIZE : 'E6.1)
       WRITE(6,553)DELT
553	FORMAT(2X'TIME STEP : 'F8.1)
       WRITE(6,554)VISCOS
554	FORMAT(2X'Viscosity : 'e10.4)
C
C======================================================================C
C --- RESTART FROM RESTART FILES IF COMPUTATION IS NOT TO START         C
C               FROM REST                                              C
C======================================================================C
C
       IF(NSTEP1.EQ.0) GO TO 75
c
C     ---------------------------------------------------------------
C
      DELT1=DELT+DELT
C
C     ---------------------------------------------------------------
C
      NO=10
      READ (NO) NSTEP0,TIME0,U
      NO=11
      READ (NO) NSTEP0,TIME0,V
      NO=12
      READ (NO) NSTEP0,TIME0,DP
C
C======================================================================C
C
C --- FOLLOWIND DO LOOP NEEDED ONLY IF HISTORY FILE USED TO RESTART
C
C           DO 2 K = 1,KK
C           DO 2 I = 1,II
C           DO 2 J = 1,JJ
C      V(I,J,K+KK) = V(I,J,K)
C      U(I,J,K+KK) = U(I,J,K)
C    2 DP(I,J,K+KK)= DP(I,J,K)
C      DELT1=DELT
C
C======================================================================C
C
       WRITE (6,1112) NSTEP0,NSTEP1
1112    FORMAT (15X'TIME STEP IN HISTORY FILE - 'I7,9X'WANTED - 'I7)
C
C
      IF(NSTEP0.NE.NSTEP1) THEN
      WRITE (6,103) NSTEP0,NSTEP1
103   FORMAT (/' STEP NO.'I5' SHOULD BE'I5)
      ENDIF
C
75       	CONTINUE
      IREC=0
      NSTEP=NSTEP1
      WRITE (6,102) NSTEP
102   FORMAT (/' MODEL INTEGRATION STARTS FROM TIME STEP ',I6/)
C     
C
C
C=======================================================================C
C                                                                       C
C                     MAIN LOOP STARTS HERE                             C
C                     =====================                             C
C                           =========                                   C
C=======================================================================C
C
15    MM=(M-1)*KK
      NN=(N-1)*KK
c
c
      NSTEP=NSTEP+1
      TIME=NSTEP*DELT/86400.
      INDEX=0
      IF (AMOD(TIME+.0001,1.).LT..0002) INDEX=1
C
C
C=======================================================================C
C           CONTINUITY EQUATION (LAYER THICKNESS PREDICTION)            C
C=======================================================================C
C
C     ****************************************************************
C --- ----------------------------------------------------------------
C
C ---    (USING 'FLUX CORRECTED TRANSPORT' ALGORITHM)
C
C --- ----------------------------------------------------------------
C     ****************************************************************
C
      DO 39 J=1,JJ1
      DO 39 I=1,II1
      UTROP(I,J)=0.
39    VTROP(I,J)=0.
      DO 1 K=1,KK
      KM=K+MM
      KN=K+NN
C
C --- COMPUTE LOW-ORDER (DIFFUSIVE) FLUXES
C
      DO 11 J=1,JJ1
      DO 11 I=2,II1
11    UFLUX(I,J)=U(I,J,KM)*CVMGP(DP(I-1,J,KN),DP(I,J,KN),U(I,J,KM))
      DO 12 J=2,JJ1
      DO 12 I=1,II1
12    VFLUX(I,J)=V(I,J,KM)*CVMGP(DP(I,J-1,KN),DP(I,J,KN),V(I,J,KM))
C
C --- ADVANCE  DP  FIELD USING LOW-ORDER (DIFFUSIVE) FLUX VALUES
C
      DO 19 J=1,JJ1
      DO 19 I=1,II1
      UOLD(I,J,K)=DP(I,J,KN)
19    DP(I,J,KN)=DP(I,J,KN)-(UFLUX(I+1,J)-UFLUX(I,J)
     .                      +VFLUX(I,J+1)-VFLUX(I,J))*X*DELT1
C
C --- COMPUTE 'ANTIDIFFUSIVE' (I.E., HIGH-ORDER MINUS LOW-ORDER) FLUXES
C
C
C --- SECOND-ORDER FLUXES FOR CONTINUITY EQUATION
C
      DO 16 J=1,JJ1
      UFLUX(  2,J)=U(  2,J,KM)*.5*(DP(  2,J,KM)+DP(  1,J,KM))-
     .UFLUX(  2,J)
16    UFLUX(II1,J)=U(II1,J,KM)*.5*(DP(II2,J,KM)+DP(II1,J,KM))-
     .UFLUX(II1,J)
      DO 17 I=1,II1
      VFLUX(I,  2)=V(I,  2,KM)*.5*(DP(I,  2,KM)+DP(I,  1,KM))-
     .VFLUX(I,  2)
17    VFLUX(I,JJ1)=V(I,JJ1,KM)*.5*(DP(I,JJ2,KM)+DP(I,JJ1,KM))-
     .VFLUX(I,JJ1)
      DO 20 J=1,JJ1
      DO 20 I=3,II2
   20 UFLUX(I,J)=U(I,J,KM)*.5*(DP(I,J,KM)+DP(I-1,J,KM))-
     .UFLUX(I,J)
      DO 21 J=3,JJ2
      DO 21 I=1,II1
   21 VFLUX(I,J)=V(I,J,KM)*.5*(DP(I,J,KM)+DP(I,J-1,KM))-
     .VFLUX(I,J)
C
C --- AT EACH GRID POINT, DETERMINE THE RATIO OF THE LARGEST PERMISSIBLE
C --- POS. (NEG.) CHANGE IN  DP  TO THE SUM OF ALL INCOMING (OUTGOING) FLUXES
C
C --- FIRST, THE FOUR CORNER POINTS...
      UTIL1(  1,  1)=-(DP(  1,  1,KN)-AMAX1(
     .DP(  1,  1,KN),DP(  2,  1,KN),DP(  1,  2,KN)))
     ./((-AMIN1(0.,UFLUX(  2,  1))
     .   -AMIN1(0.,VFLUX(  1,  2))+EPSIL)*X*DELT1)
      UTIL2(  1,  1)=-(DP(  1,  1,KN)-AMAX1(0.,AMIN1(
     .DP(  1,  1,KN),DP(  2,  1,KN),DP(  1,  2,KN))))
     ./((-AMAX1(0.,UFLUX(  2,  1))
     .   -AMAX1(0.,VFLUX(  1,  2))-EPSIL)*X*DELT1)
C
      UTIL1(  1,JJ1)=-(DP(  1,JJ1,KN)-AMAX1(
     .DP(  1,JJ1,KN),DP(  2,JJ1,KN),DP(  1,JJ2,KN)))
     ./((-AMIN1(0.,UFLUX(  2,JJ1))
     .   +AMAX1(0.,VFLUX(  1,JJ1))+EPSIL)*X*DELT1)
      UTIL2(  1,JJ1)=-(DP(  1,JJ1,KN)-AMAX1(0.,AMIN1(
     .DP(  1,JJ1,KN),DP(  2,JJ1,KN),DP(  1,JJ2,KN))))
     ./((-AMAX1(0.,UFLUX(  2,JJ1))
     .   +AMIN1(0.,VFLUX(  1,JJ1))-EPSIL)*X*DELT1)
C
      UTIL1(II1,  1)=-(DP(II1,  1,KN)-AMAX1(
     .DP(II1,  1,KN),DP(II2,  1,KN),DP(II1,  2,KN)))
     ./(( AMAX1(0.,UFLUX(II1,  1))
     .   -AMIN1(0.,VFLUX(II1,  2))+EPSIL)*X*DELT1)
      UTIL2(II1,  1)=-(DP(II1,  1,KN)-AMAX1(0.,AMIN1(
     .DP(II1,  1,KN),DP(II2,  1,KN),DP(II1,  2,KN))))
     ./(( AMIN1(0.,UFLUX(II1,  1))
     .   -AMAX1(0.,VFLUX(II1,  2))-EPSIL)*X*DELT1)
C
      UTIL1(II1,JJ1)=-(DP(II1,JJ1,KN)-AMAX1(
     .DP(II1,JJ1,KN),DP(II2,JJ1,KN),DP(II1,JJ2,KN)))
     ./(( AMAX1(0.,UFLUX(II1,JJ1))
     .   +AMAX1(0.,VFLUX(II1,JJ1))+EPSIL)*X*DELT1)
      UTIL2(II1,JJ1)=-(DP(II1,JJ1,KN)-AMAX1(0.,AMIN1(
     .DP(II1,JJ1,KN),DP(II2,JJ1,KN),DP(II1,JJ2,KN))))
     ./(( AMIN1(0.,UFLUX(II1,JJ1))
     .   +AMIN1(0.,VFLUX(II1,JJ1))-EPSIL)*X*DELT1)
C
C --- NOW THE REMAINING LATERAL BOUNDARY POINTS...
      DO 25 I=2,II2
      UTIL1(I,  1)=-(DP(I,  1,KN)-AMAX1(
     .DP(I,  1,KN),DP(I-1,  1,KN),DP(I+1,  1,KN),DP(I,  2,KN)))
     ./((AMAX1(0.,UFLUX(I,  1))-AMIN1(0.,UFLUX(I+1,  1))
     .                         -AMIN1(0.,VFLUX(I  ,  2))+EPSIL)*X*DELT1)
      UTIL2(I,  1)=-(DP(I,  1,KN)-AMAX1(0.,AMIN1(
     .DP(I,  1,KN),DP(I-1,  1,KN),DP(I+1,  1,KN),DP(I,  2,KN))))
     ./((AMIN1(0.,UFLUX(I,  1))-AMAX1(0.,UFLUX(I+1,  1))
     .                         -AMAX1(0.,VFLUX(I  ,  2))-EPSIL)*X*DELT1)
      UTIL1(I,JJ1)=-(DP(I,JJ1,KN)-AMAX1(
     .DP(I,JJ1,KN),DP(I-1,JJ1,KN),DP(I+1,JJ1,KN),DP(I,JJ2,KN)))
     ./((AMAX1(0.,UFLUX(I,JJ1))-AMIN1(0.,UFLUX(I+1,JJ1))
     .  +AMAX1(0.,VFLUX(I,JJ1))                         +EPSIL)*X*DELT1)
 25   UTIL2(I,JJ1)=-(DP(I,JJ1,KN)-AMAX1(0.,AMIN1(
     .DP(I,JJ1,KN),DP(I-1,JJ1,KN),DP(I+1,JJ1,KN),DP(I,JJ2,KN))))
     ./((AMIN1(0.,UFLUX(I,JJ1))-AMAX1(0.,UFLUX(I+1,JJ1))
     .  +AMIN1(0.,VFLUX(I,JJ1))                         -EPSIL)*X*DELT1)
C
      DO 26 J=2,JJ2
      UTIL1(  1,J)=-(DP(  1,J,KN)-AMAX1(
     .DP(  1,J,KN),DP(  2,J,KN),DP(  1,J-1,KN),DP(  1,J+1,KN)))
     ./((                      -AMIN1(0.,UFLUX(  2,J  ))
     .  +AMAX1(0.,VFLUX(  1,J))-AMIN1(0.,VFLUX(  1,J+1))+EPSIL)*X*DELT1)
      UTIL2(  1,J)=-(DP(  1,J,KN)-AMAX1(0.,AMIN1(
     .DP(  1,J,KN),DP(  2,J,KN),DP(  1,J-1,KN),DP(  1,J+1,KN))))
     ./((                      -AMAX1(0.,UFLUX(  2,J  ))
     .  +AMIN1(0.,VFLUX(  1,J))-AMAX1(0.,VFLUX(  1,J+1))-EPSIL)*X*DELT1)
      UTIL1(II1,J)=-(DP(II1,J,KN)-AMAX1(
     .DP(II1,J,KN),DP(II2,J,KN),DP(II1,J-1,KN),DP(II1,J+1,KN)))
     ./((AMAX1(0.,UFLUX(II1,J))
     .  +AMAX1(0.,VFLUX(II1,J))-AMIN1(0.,VFLUX(II1,J+1))+EPSIL)*X*DELT1)
 26   UTIL2(II1,J)=-(DP(II1,J,KN)-AMAX1(0.,AMIN1(
     .DP(II1,J,KN),DP(II2,J,KN),DP(II1,J-1,KN),DP(II1,J+1,KN))))
     ./((AMIN1(0.,UFLUX(II1,J))
     .  +AMIN1(0.,VFLUX(II1,J))-AMAX1(0.,VFLUX(II1,J+1))-EPSIL)*X*DELT1)
C
C --- FINALLY, THE INTERIOR GRID POINTS...
      DO 27 J=2,JJ2
      DO 27 I=2,II2
      UTIL1(I,J)=-(DP(I,J,KN)-AMAX1(
     .DP(I,J,KN),DP(I-1,J,KN),DP(I+1,J,KN),DP(I,J-1,KN),DP(I,J+1,KN)))
     ./((AMAX1(0.,UFLUX(I,J))-AMIN1(0.,UFLUX(I+1,J))
     .  +AMAX1(0.,VFLUX(I,J))-AMIN1(0.,VFLUX(I,J+1))+EPSIL)*X*DELT1)
 27   UTIL2(I,J)=-(DP(I,J,KN)-AMAX1(0.,AMIN1(
     .DP(I,J,KN),DP(I-1,J,KN),DP(I+1,J,KN),DP(I,J-1,KN),DP(I,J+1,KN))))
     ./((AMIN1(0.,UFLUX(I,J))-AMAX1(0.,UFLUX(I+1,J))
     .  +AMIN1(0.,VFLUX(I,J))-AMAX1(0.,VFLUX(I,J+1))-EPSIL)*X*DELT1)
C
C --- LIMIT ANTIDIFFUSIVE FLUXES
C --- (RETAIN INFORMATION ABOUT FLUX CLIPPING IN -UTROP,VTROP-. THIS
C --- WILL ALLOW US LATER TO RESTORE NONDIVERGENCE OF BAROTROPIC FLOW.)
C
      DO 28 J=1,JJ1
      DO 28 I=2,II1
      UTROP(I,J)=UTROP(I,J)+UFLUX(I,J)*
     .(1.-CVMGP(AMIN1(1.,UTIL1(I,J),UTIL2(I-1,J)),
     .          AMIN1(1.,UTIL2(I,J),UTIL1(I-1,J)),UFLUX(I,J)))
28    UFLUX(I,J)=UFLUX(I,J)*CVMGP(AMIN1(1.,UTIL1(I,J),UTIL2(I-1,J)),
     .                            AMIN1(1.,UTIL2(I,J),UTIL1(I-1,J)),
     .UFLUX(I,J))
      DO 29 J=2,JJ1
      DO 29 I=1,II1
      VTROP(I,J)=VTROP(I,J)+VFLUX(I,J)*
     .(1.-CVMGP(AMIN1(1.,UTIL1(I,J),UTIL2(I,J-1)),
     .          AMIN1(1.,UTIL2(I,J),UTIL1(I,J-1)),VFLUX(I,J)))
29    VFLUX(I,J)=VFLUX(I,J)*CVMGP(AMIN1(1.,UTIL1(I,J),UTIL2(I,J-1)),
     .                            AMIN1(1.,UTIL2(I,J),UTIL1(I,J-1)),
     .VFLUX(I,J))
C
C --- ADD ANTIDIFFUSIVE FLUXES TO  DP  FIELD
C
      DO 1 J=1,JJ1
      DO 1 I=1,II1
1     DP(I,J,KN)=DP(I,J,KN)-(UFLUX(I+1,J)-UFLUX(I,J)
     .                      +VFLUX(I,J+1)-VFLUX(I,J))*X*DELT1
C
C --- RESTORE NONDIVERGENCE OF VERTICALLY INTEGRATED FLOW
C
      DO 5 J=1,JJ1
      DO 5 K=1,KK
      DO 5 I=1,II1
5     P(I,J,K+1)=P(I,J,K)+DP(I,J,K+NN)
      DO 14 K=1,KK
      KN=K+NN
      DO 44 J=1,JJ1
      DO 44 I=2,II1
44    UFLUX(I,J)=UTROP(I,J)*CVMGP(DP(I-1,J,KN)/P(I-1,J,KK+1),
     .                            DP(I  ,J,KN)/P(I  ,J,KK+1),UTROP(I,J))
      DO 45 J=2,JJ1
      DO 45 I=1,II1
45    VFLUX(I,J)=VTROP(I,J)*CVMGP(DP(I,J-1,KN)/P(I,J-1,KK+1),
     .                            DP(I,J  ,KN)/P(I,J  ,KK+1),VTROP(I,J))
      DO 14 J=1,JJ1
      DO 14 I=1,II1
14    DP(I,J,KN)=DP(I,J,KN)-(UFLUX(I+1,J)-UFLUX(I,J)
     .                      +VFLUX(I,J+1)-VFLUX(I,J))*X*DELT1
C
C     ****************************************************************
C    
C ---                   END OF  F C T  CALCULATIONS
C
C     ****************************************************************
C
C=======================================================================C
C                       HYDROSTATIC EQUATION                            C
C=======================================================================C
C
      DO 83 J=1,JJ1
      DO 8 K=1,KK
      DO 8 I=1,II1
    8      P(I,J,K+1)=P(I,J,K)+DP(I,J,K+MM)
      DO 13 KI=2,KK
      K=KK+1-KI
      DO 13 I=1,II1
   13    MONTG(I,J,K)=MONTG(I,J,K+1)+P(I,J,K+1)*(THETA(K+1)-THETA(K))
C
C --- TIME SMOOTHING OF THICKNESS FIELD
C
      DO 83 K=1,KK
      DO 83 I=1,II1
   83      DP(I,J,K+MM)=DP(I,J,K+MM)*.98+(UOLD(I,J,K)+DP(I,J,K+NN))*.01
C
C
C=======================================================================C
C                       MOMENTUM EQUATIONS                              C
C=======================================================================C
C
      DO 9 K=1,KK

      KM=K+MM
      KN=K+NN
      XCONT(K)=0.
      DISVIS(K)=0.
      STRESS(K)=0.
C
c
C
c
C=======================================================================C
C                        VORTICITY COMPUTATION                          C
C=======================================================================C
C
      DO 80 I=2,II1
      DO 80 J=2,JJ1
   80 ABSVOR(I,J)=(CORIO(I,J)+(V(I,J,KM)-V(I-1,J,KM)-U(I,J,KM)
     .+U(I,J-1,KM))*X)
c
c     Impose the boundary conditions on the calculation of ABSVOR
c     near the boundaries
c
C
	IF(IFREE.EQ.2) THEN      ! for no-slip

      DO 81 J=2,JJ1
      ABSVOR(1,J)=(CORIO(1,J)+2.*V(1,J,KM)*X)
   81 ABSVOR(II,J)=(CORIO(II,J)-2.*V(II1,J,KM)*X)
C
      DO 82 I=2,II1
      ABSVOR(I,1)=(CORIO(I,1)-2.*U(I,1,KM)*X)
 82   ABSVOR(I,JJ)=(CORIO(I,JJ)+2.*U(I,JJ1,KM)*X)
      
	ENDIF                                  
C
C
C======================================================================C
C                       OOOOOOOOOOOOOOOOOOOOO                           C
C=======================================================================C
C                            U EQUATION                                 C
C=======================================================================C
C
C=======================================================================C
C                THICKNESS FIELD FOR LATERAL DIFFUSION                  C
C=======================================================================C

      DO 412 J=1,JJ1
      PS(1,J)=DP(1,J,KM)
      PS(II,J)=DP(II1,J,KM)
      UOLD(1,J,K)=0.
      UOLD(II,J,K)=0.
C
      DO 412 I=2,II1
      PS(I,J)=.5*(DP(I,J,KM)+DP(I-1,J,KM))
      UOLD(I,J,K)=U(I,J,KN)
  412 CONTINUE
C
      DO 37 J=1,JJ1
      JA=MAX0(  1,J-1)
      JB=MIN0(JJ1,J+1)
c
c	Set up the free-slip or no-slip boundary conditions
c	(if.ifree.eq.1 ---> free slip boundary condition)
C
             SIG1 = 1.
             SIG2 = 1.
	if(ifree.eq.2) then        ! no-slip boundary condition
             IF(J.EQ.  1) SIG1 = -1.
             IF(J.EQ.JJ1) SIG2 = -1.
	endif
C
C --- RECTANGULAR BASIN
c
C	...this DEL2U term is neccesary for the laplacian viscosity term...
C
      DO 37 I=2,II1
   37 DEL2U(I,J)=
     .((PS(I+1,J)+PS(I,J))*(UOLD(I+1,J,K)-UOLD(I,J,K))-
     . (PS(I,J)+PS(I-1,J))*(UOLD(I,J,K)-UOLD(I-1,J,K))+
     . (PS(I,JB)+PS(I,J))*(SIG2*UOLD(I,JB,K)-UOLD(I,J,K))-
     . (PS(I,J )+PS(I,JA))*(UOLD(I,J,K)-SIG1*UOLD(I,JA,K)))*.5

C=======================================================================C
C            BEGIN INTEGRATION OF U MOMENTUM EQUATION                   C
C=======================================================================C
      DO 61 J=1,JJ1
CDIR$ IVDEP
      DO 6 I=2,II1
6     U(I,J,KN)=
C
C --- HORIZONTAL PRESSURE FORCE (X-DIRECTION)
     .U(I,J,KN)-DELT1*((MONTG(I,J,K)-MONTG(I-1,J,K)
C
C --- GRADIENT OF KINETIC ENERGY
     .+.25*(U(I+1,J,KM)**2+V(I  ,J,KM)**2+V(I  ,J+1,KM)**2
     .     -U(I-1,J,KM)**2-V(I-1,J,KM)**2-V(I-1,J+1,KM)**2))*X
C
C --- ABSOLUTE VORTICITY FLUX 
     .-.125*(V(I,J,KM)+V(I,J+1,KM)+V(I-1,J,KM)+V(I-1,J+1,KM))*
     .(ABSVOR(I,J)+ABSVOR(I,J+1))
C
C --- WIND STRESS FORCING
     .-STRESX(I,J)*(AMIN1(P(I,J,K+1)+P(I-1,J,K+1)+1.E3,PSTRES+PSTRES)
     .             -AMIN1(P(I,J,K  )+P(I-1,J,K  )   ,PSTRES+PSTRES))/
     .             (      P(I,J,K+1)+P(I-1,J,K+1)+1.E3
     .                   -P(I,J,K  )-P(I-1,J,K  ))
C
C --- VISCOUS DIFFUSION
     .-DEL2U(I,J)*VISCOS*X*X
     ./AMAX1(.5*(DP(I,J,KM)+DP(I-1,J,KM)),EPS))
C
      IF (INDEX.EQ.0) GO TO 61
C --- EVALUATE TERMS IN KINETIC ENERGY EQUATION
      DO 59 I=2,II1
C
      XCONT(K)=XCONT(K)-U(I,J,KM)*.5*(DP(I,J,KM)+DP(I-1,J,KM))*
     .(MONTG(I,J,K)-MONTG(I-1,J,K))/(G*SCALE)
C
      STRESS(K)=STRESS(K)+U(I,J,KM)*.5*(DP(I,J,KM)+DP(I-1,J,KM))*
     .STRESX(I,J)*(AMIN1(P(I,J,K+1)+P(I-1,J,K+1)+1.E3,PSTRES+PSTRES)
     .    -AMIN1(P(I,J,K  )+P(I-1,J,K  )   ,PSTRES+PSTRES))/
     .    (      P(I,J,K+1)+P(I-1,J,K+1)+1.E3
     .    -P(I,J,K  )-P(I-1,J,K  ))/G
C
      IF (K.EQ.KK) STRESS(K)=STRESS(K)-
     .DRAG*.5*(DP(I,J,KM)+DP(I-1,J,KM))*U(I,J,KM)*U(I,J,KN)/G
C
59    DISVIS(K)=DISVIS(K)-U(I,J,KM)*.5*(DP(I,J,KM)+DP(I-1,J,KM))/G
     .*DEL2U(I,J)*VISCOS*X*X
     ./AMAX1(.5*(DP(I,J,KM)+DP(I-1,J,KM)),EPS)
61    CONTINUE
C
C=======================================================================C
C                       OOOOOOOOOOOOOOOOOOOOO                           C
C=======================================================================C
C                             V EQUATION                                C
C=======================================================================C
C
      DO 76 I=1,II1
      PS(I,1)=DP(I,1,KM)
      PS(I,JJ)=DP(I,JJ1,KM)
      VOLD(I,1,K)=0.
      VOLD(I,JJ,K)=0.
      DO 76 J=2,JJ1
      PS(I,J)=.5*(DP(I,J,KM)+DP(I,J-1,KM))
      VOLD(I,J,K)=V(I,J,KN)
   76 CONTINUE
C
C
      DO 420 I=1,II1
      IA=MAX0(1,I-1)
      IB=MIN0(II1,I+1)
c
c	Set up the free-slip or no-slip boundary conditions
c	(if.ifree.eq.1 ---> free slip boundary condition)
C
             SIG1 = 1.
             SIG2 = 1.
	if(ifree.eq.2) then       ! no-slip boundary condition
             IF(I.EQ.  1) SIG1 = -1.
             IF(I.EQ.II1) SIG2 = -1.
	endif
C
      DO 420 J=2,JJ1
  420 DEL2V(I,J)=
     .((PS(IB,J)+PS(I,J))*(SIG2*VOLD(IB,J,K)-VOLD(I,J,K))-
     . (PS(I,J)+PS(IA,J))*(VOLD(I,J,K)-SIG1*VOLD(IA,J,K))+
     . (PS(I,J+1)+PS(I,J))*(VOLD(I,J+1,K)-VOLD(I,J,K))-
     . (PS(I,J)+PS(I,J-1))*(VOLD(I,J,K)-VOLD(I,J-1,K)))*.5

C======================================================================C
C                 BEGIN INTEGRATION OF V-MOMENTUM-EQUATION             C
C======================================================================C
C
      DO 62 I=1,II1
CDIR$ IVDEP
      DO 7 J=2,JJ1
7     V(I,J,KN)=
C
C --- HORIZONTAL PRESSURE FORCE (Y-DIRECTION)
     .V(I,J,KN)-DELT1*((MONTG(I,J,K)-MONTG(I,J-1,K)
C
C --- GRADIENT OF KINETIC ENERGY
     .+.25*(V(I,J+1,KM)**2+U(I,J  ,KM)**2+U(I+1,J  ,KM)**2
     .     -V(I,J-1,KM)**2-U(I,J-1,KM)**2-U(I+1,J-1,KM)**2))*X
C
C --- ABSOLUTE VORTICITY FLUX
     .+.125*(U(I,J,KM)+U(I+1,J,KM)+U(I,J-1,KM)+U(I+1,J-1,KM))*
     .(ABSVOR(I,J)+ABSVOR(I+1,J))
C
C --- WIND STRESS FORCING
     .-STRESY(I,J)*(AMIN1(P(I,J,K+1)+P(I,J-1,K+1)+1.E3,PSTRES+PSTRES)
     .             -AMIN1(P(I,J,K  )+P(I,J-1,K  )   ,PSTRES+PSTRES))/
     .             (      P(I,J,K+1)+P(I,J-1,K+1)+1.E3
     .                   -P(I,J,K  )-P(I,J-1,K  ))
C
C --- VISCOUS DIFFUSION
     .-DEL2V(I,J)*X*X*VISCOS
     ./AMAX1(.5*(DP(I,J,KM)+DP(I,J-1,KM)),EPS))
      IF (INDEX.EQ.0) GO TO 62
C
C --- EVALUATE TERMS IN KINETIC ENERGY EQUATION
      DO 60 J=2,JJ1
      XCONT(K)=XCONT(K)-V(I,J,KM)*.5*(DP(I,J,KM)+DP(I,J-1,KM))*
     .(MONTG(I,J,K)-MONTG(I,J-1,K))/(G*SCALE)
C
      STRESS(K)=STRESS(K)+V(I,J,KM)*.5*(DP(I,J,KM)+DP(I,J-1,KM))*
     .STRESY(I,J)*(AMIN1(P(I,J,K+1)+P(I,J-1,K+1)+1.E3,PSTRES+PSTRES)
     .    -AMIN1(P(I,J,K  )+P(I,J-1,K  )   ,PSTRES+PSTRES))/
     .    (      P(I,J,K+1)+P(I,J-1,K+1)+1.E3
     .    -P(I,J,K  )-P(I,J-1,K  ))/G
C
      IF (K.EQ.KK) STRESS(K)=STRESS(K)-
     .DRAG*.5*(DP(I,J,KM)+DP(I,J-1,KM))*V(I,J,KM)*V(I,J,KN)/G
C
60    DISVIS(K)=DISVIS(K)-V(I,J,KM)*.5*(DP(I,J,KM)+DP(I,J-1,KM))/G
     .*DEL2V(I,J)*X*X*VISCOS
     ./AMAX1(.5*(DP(I,J,KM)+DP(I,J-1,KM)),EPS)
62    CONTINUE
9     CONTINUE
C
C======================================================================C
C                         BOTTOM DRAG                                  C
C======================================================================C
C
      KN=KK+NN
      DO 18 J=1,JJ1
      DO 18 I=1,II1
      U(I,J,KN)=U(I,J,KN)*(1.-DRAG*DELT1)
18    V(I,J,KN)=V(I,J,KN)*(1.-DRAG*DELT1)
C
C----------------------------------------------------------------------C
C
      DO 850 K=1,KK
      DO 850 I=1,II1
      DO 850 J=1,JJ1
  850 P(I,J,K+1)=P(I,J,K)+DP(I,J,K+NN)
C
C
C --- SUBSTITUTE DEPTH-WEIGHTED AVERAGES FOR (U,V) AT MASSLESS GRID POINTS
C
      DO 771 K=1,KK+1
      DO 771 J=1,JJ1
      DO 771 I=2,II1
  771 PU(I,J,K)=.5*(P(I,J,K)+P(I-1,J,K))
C
      DO 772 K=1,KK+1
      DO 772 J=2,JJ1
      DO 772 I=1,II1
  772 PV(I,J,K)=.5*(P(I,J,K)+P(I,J-1,K))
C
C
      DO 774 K=1,KK
      DO 668 J=2,JJ1
      DO 668 I=1,II1
  668 VLD(I,J,K)=V(I,J,K+NN)*(PV(I,J,K+1)-PV(I,J,K))
C
      DO 773 K1=1,K-1
      DO 773 J=2,JJ1
      DO 773 I=1,II1
  773 VLD(I,J,K)=VLD(I,J,K)+V(I,J,K1+NN)*
     .(AMIN1(PV(I,J,K  ),AMAX1(PV(I,J,K+1)-EPS1,PV(I,J,K1+1)))
     .-AMIN1(PV(I,J,K  ),AMAX1(PV(I,J,K+1)-EPS1,PV(I,J,K1  ))))
C
      DO 774 K1=K+1,KK
      DO 774 J=2,JJ1
      DO 774 I=1,II1
  774 VLD(I,J,K)=VLD(I,J,K)+V(I,J,K1+NN)*
     .(AMAX1(PV(I,J,K+1),AMIN1(PV(I,J,K  )+EPS1,PV(I,J,K1+1)))
     .-AMAX1(PV(I,J,K+1),AMIN1(PV(I,J,K  )+EPS1,PV(I,J,K1  ))))
C
      DO 663 K=1,KK
      DO 663 J=2,JJ1
      DO 663 I=1,II1
  663 V(I,J,K+NN)=VLD(I,J,K)/
     .(AMAX1(PV(I,J,K+1),AMIN1(PV(I,J,K  )+EPS1,PV(I,J,KK+1)))
     .-AMIN1(PV(I,J,K  ),AMAX1(PV(I,J,K+1)-EPS1,PV(I,J,   1))))
C
      DO 863 J=1,JJ1
C
      DO 874 K=1,KK
      DO 888 I=2,II1
  888 ULD(I,J,K)=U(I,J,K+NN)*(PU(I,J,K+1)-PU(I,J,K))
      DO 873 K1=1,K-1
      DO 873 I=2,II1
  873 ULD(I,J,K)=ULD(I,J,K)+U(I,J,K1+NN)*
     .(AMIN1(PU(I,J,K  ),AMAX1(PU(I,J,K+1)-EPS1,PU(I,J,K1+1)))
     .-AMIN1(PU(I,J,K  ),AMAX1(PU(I,J,K+1)-EPS1,PU(I,J,K1  ))))
C
      DO 874 K1=K+1,KK
      DO 874 I=2,II1
  874 ULD(I,J,K)=ULD(I,J,K)+U(I,J,K1+NN)*
     .(AMAX1(PU(I,J,K+1),AMIN1(PU(I,J,K  )+EPS1,PU(I,J,K1+1)))
     .-AMAX1(PU(I,J,K+1),AMIN1(PU(I,J,K  )+EPS1,PU(I,J,K1  ))))
C
            DO 863 K=1,KK
      DO 863 I=2,II1
  863 U(I,J,K+NN)=ULD(I,J,K)/
     .(AMAX1(PU(I,J,K+1),AMIN1(PU(I,J,K  )+EPS1,PU(I,J,KK+1)))
     .-AMIN1(PU(I,J,K  ),AMAX1(PU(I,J,K+1)-EPS1,PU(I,J,1   ))))

C======================================================================C
C REMOVE DIVERGENT COMPONENT FROM MEAN MOTION (RIGID LID APPROXIMAT.)  C
C======================================================================C
C
      DO 36 J=1,JJ1
      DO 36 I=1,II1
        P(I,J,KK+1)=0.
         UTROP(I,J)=0.
36        VTROP(I,J)=0.
      DO 32 K=1,KK
       KN=K+NN
        DO 30 J=1,JJ1
         DO 30 I=1,II1
30        P(I,J,KK+1)=P(I,J,KK+1)+DP(I,J,KN)
      DO 46 J=1,JJ1
      UTROP(  2,J)=UTROP(  2,J)+U(  2,J,KN)*(DP(  2,J,KN)+DP(  1,J,KN))
46    UTROP(II1,J)=UTROP(II1,J)+U(II1,J,KN)*(DP(II2,J,KN)+DP(II1,J,KN))
      DO 31 J=1,JJ1
      DO 31 I=3,II2
   31 UTROP(I,J)=UTROP(I,J)+U(I,J,KN)*(DP(I,J,KN)+DP(I-1,J,KN))
      DO 47 I=1,II1
      VTROP(I,  2)=VTROP(I,  2)+V(I,  2,KN)*(DP(I,  2,KN)+DP(I,  1,KN))
47    VTROP(I,JJ1)=VTROP(I,JJ1)+V(I,JJ1,KN)*(DP(I,JJ2,KN)+DP(I,JJ1,KN))
      DO 32 J=3,JJ2
      DO 32 I=1,II1
   32 VTROP(I,J)=VTROP(I,J)+V(I,J,KN)*(DP(I,J,KN)+DP(I,J-1,KN))
C
C
      DO 34 J=1,JJ1
      DO 34 I=2,II1
34    UTROP(I,J)=UTROP(I,J)/(P(I,J,KK+1)+P(I-1,J,KK+1))
      DO 35 J=2,JJ1
      DO 35 I=1,II1
35    VTROP(I,J)=VTROP(I,J)/(P(I,J,KK+1)+P(I,J-1,KK+1))
      DO 33 J=2,JJ1
      DO 33 I=2,II1
33    VORT(I,J)=(VTROP(I,J)-VTROP(I-1,J)-UTROP(I,J)+UTROP(I,J-1))
C
C=======================================================================C
C --- SOLVE POISSON EQUATION TO DETERMINE ROTATIONAL PART OF MEAN MOTION
C --- FOR TESTING THE FFT ROUTINE, ACTIVATE THE TWO STATEMENTS ..PRINT 110...
C --- AND THE LOOP DO 52.... FOR PBOTTM=CONST., THE RESULTS SHOULD AGREE
C=======================================================================C
C
      IF (INDEX.GT.0) PRINT 110,(VORT(I,I),I=2,13)
110   FORMAT (1X,1P12E10.3)

       CALL POISSON(VORT(2,2),II,II2,JJ2,PBOT,EIG,WSAVE)
C
      DO 52 I=2,13
52    UTIL1(I,I)=(STRMF(I-1,I)+STRMF(I+1,I)
     .         +STRMF(I,I-1)+STRMF(I,I+1)-4.*STRMF(I,I))/PBOTTM
      IF (INDEX.GT.0) PRINT 110,(UTIL1(I,I),I=2,13)
C
C
C
C=======================================================================C
C     COMPUTE THE NONDIVERGENT BAROTROPIC VELOCITY FIELD (UTROP) FROM   C
C     THE  STREMFUNCTION AND SUBTRACT THIS FROM THE TOTAL TO FIND THE   C
C     DIVERGENT PART OF THE BAROTROPIC FLOW.                            C
C=======================================================================C
C
      DO 43 J=1,JJ1
       DO 43 I=2,II1
43    UTROP(I,J)=UTROP(I,J)-(STRMF(I,J)-STRMF(I,J+1))
     ./(.5*(P(I,J,KK+1)+P(I-1,J,KK+1)))
C
      DO 40 J=2,JJ1
       DO 40 I=1,II1
40    VTROP(I,J)=VTROP(I,J)-(STRMF(I+1,J)-STRMF(I,J))
     ./(.5*(P(I,J,KK+1)+P(I,J-1,KK+1)))
C
      DO 41 K=1,KK
        KN=K+NN
         DO 42 J=1,JJ1
          DO 42 I=2,II1
42         U(I,J,KN)=U(I,J,KN)-UTROP(I,J)
            DO 41 J=2,JJ1
             DO 41 I=1,II1
41    V(I,J,KN)=V(I,J,KN)-VTROP(I,J)
C
C=======================================================================C
C                 SMOOTH  U  AND  V  FIELDS IN TIME                     C
C=======================================================================C
C
      DO 22 K=1,KK
      KM=K+MM
      KN=K+NN
      DO 324 J=1,JJ1
      DO 324 I=2,II1
      U(I,J,KM)=U(I,J,KM)*WGT1+(UOLD(I,J,K)+U(I,J,KN))*WGT2
324   CONTINUE

      DO 222 J=2,JJ1
      DO 222 I=1,II1
      V(I,J,KM)=V(I,J,KM)*WGT1+(VOLD(I,J,K)+V(I,J,KN))*WGT2
222	CONTINUE
C
22	CONTINUE
c
C=======================================================================C
C                                                                       C
C                  OUTPUT AND DIAGNOSTIC CALCULATIONS                   C
C                  ==================================                   C
C                                                                       C
C=======================================================================C
C
      IF (INDEX.EQ.0) GO TO 23
      PRINT 100,NSTEP,TIME
100   FORMAT (' T I M E   S T E P'I9,25X'D A Y'F8.1)
C
C --- ENERGY DIAGNOSTICS
C
      PBAR=0.
      SUMWGT=0.
      SUMPOT=0.
      SUMKIN=0.
      SUMXGR=0.
      SUMSTR=0.
      SUMDIS=0.
      SUMEKT=0.
      DO 53 K=1,KK
      KM=K+MM
      EPOT=0.
      EKIN=0.
      WEIGHT=PBOTTM
      DO 51 I=1,II1
      DO 51 J=1,JJ1
      IF (K.EQ.1)
     .EPOT=EPOT+.5*(STRMF(I,J)*CORIO(I,J)*SCALE/P(I,J,KK+1))**2/G
      IF (K.GT.1)
     .EPOT=EPOT+.5*(P(I,J,K)-PBAR)**2*(THETA(K)-THETA(K-1))/G
      EKIN=EKIN+.25*DP(I,J,KM)*(U(I,J,KM)**2+U(I+1,J,KM)**2
     .                         +V(I,J,KM)**2+V(I,J+1,KM)**2)/G
C --- CORRECT KINETIC ENERGY BUDGET FOR REMOVAL OF DIVERGENT PART OF MEAN MOTION
      XCONT(K)=XCONT(K)-.5*DP(I,J,KM)*
     .(U(I  ,J,KM)*UTROP(I  ,J)+V(I,J  ,KM)*VTROP(I,J  )
     .+U(I+1,J,KM)*UTROP(I+1,J)+V(I,J+1,KM)*VTROP(I,J+1))/(G*DELT1)
51    WEIGHT=AMIN1(WEIGHT,DP(I,J,KM))
      EPOT=EPOT/(II1*JJ1)
      EKIN=EKIN/(II1*JJ1)
      DISVIS(K)=-DISVIS(K)/(II1*JJ1)
      XCONT(K)=XCONT(K)/(II1*JJ1)
      STRESS(K)=STRESS(K)/(II1*JJ1)
      EKINT=STRESS(K)+XCONT(K)+DISVIS(K)
      PRINT 101,WEIGHT,EPOT,EKIN,EKINT,XCONT(K),STRESS(K),DISVIS(K),K
101   FORMAT (1P' THKN'E10.2'  EPOT'E10.2'  EKIN'E10.2'  EKINT'E10.2
     .'  XCONT'E10.2'  STRESS'E10.2'  DISSIP'E10.2,I4)
      SUMWGT=SUMWGT+WEIGHT
      SUMKIN=SUMKIN+EKIN
      SUMPOT=SUMPOT+EPOT
      SUMXGR=SUMXGR+XCONT(K)
      SUMSTR=SUMSTR+STRESS(K)
      SUMDIS=SUMDIS+DISVIS(K)
      SUMEKT=SUMEKT+EKINT
53    PBAR=PBAR+DP0(K)
      PRINT 101,SUMWGT,SUMPOT,SUMKIN,SUMEKT,SUMXGR,SUMSTR,SUMDIS
C
C --- OUTPUT TO RESTART FILE
C
      NO=20
      REWIND (NO)
      WRITE (NO) NSTEP,DAY,U
      NO=21
      REWIND (NO)
      WRITE (NO) NSTEP,DAY,V
      NO=22
      REWIND (NO)
      WRITE (NO) NSTEP,DAY,DP
C
C --- OUTPUT TO HISTORY FILE
C      
      NO=30
      WRITE (NO) (((U(I,J,K),I=1,II),J=1,JJ),K=1,KK),NSTEP
      NO=31
      WRITE (NO) (((V(I,J,K),I=1,II),J=1,JJ),K=1,KK),NSTEP
      NO=32
      WRITE (NO) (((DP(I,J,K),I=1,II),J=1,JJ),K=1,KK),NSTEP
C
C --- OUTPUT TO AVERAGED FILES
C
C
C
      IREC=IREC+1
      DO 180 K=1,KK
      DO 180 J=1,JJ
      DO 180 I=1,II
      DPM(I,J,K)=DP(I,J,K+NN)+DPM(I,J,K)
      UM(I,J,K)=U(I,J,K+NN)+UM(I,J,K)
 180  VM(I,J,K)=V(I,J,K+NN)+VM(I,J,K)
C
      NO=23
C      REWIND (NO)
C      WRITE (NO) UM,IREC
      NO=24
C      REWIND (NO)
C      WRITE (NO) VM,IREC
      NO=25
C      REWIND (NO)
C      WRITE (NO) DPM,IREC
C
C
23    L=M
      M=N
      N=L
      IF (NSTEP.GE.NSTEP2) STOP '(NORMAL)'
      DELT1=DELT+DELT
      GO TO 15
      END
C=======================================================================C
C                                                                       C
          SUBROUTINE POINIT(M,N,L,PBOT,EIG,WSAVE)                      
C         =======================================                       C
C    INITIALIZATION ROUTINE FOR FFTPSSN:                                C
C    THIS VERSION OF POINIT PERMITS PBOTTOM TO VARY IN Y DIRECTION      C
C=======================================================================C

          REAL EIG(L,N), WORK(1), PBOT(L)                              
           PI = 4.*ATAN(1.)
C--- INITIALIZE FFT ROUTINE
           CALL SINTI(M,WSAVE)

C--- EIGENVALUES OF TRIDIAGONAL MATRIX
           DO 1 J=1,M
           DO 1 K=1,N
    1      EIG(J,K) = -1./PBOT(K) - 1./PBOT(K+1)
     .                -8.*SIN(J*PI/2./(M+1))**2/(PBOT(K)+PBOT(K+1))
C--- GAUSSIAN ELIMINATION OF TRI-DIAGONAL SYSTEM - LEFT HAND SIDE
           DO 2 K = 2,N
           DO 2 J = 1,M 
    2      EIG(J,K) = EIG(J,K) - 1./(PBOT(K)**2*EIG(J,K-1))
c           
           RETURN
           END
c
C=================================================================C
C                                                                 C
       SUBROUTINE POISSON(ZETA,II,M,N,PBOT,EIG,WSAVE)
C      ==============================================             C  
C                                                                 C
C                                                                 C
C       This subroutine uses a Fast Sine Transform algorithm      C
C       to solve  the  POISSON  EQUATION:                         C
C                                                                 C
C                      PSI  + PSI  = ZETA                         C
C                         xx     yy                               C
C                                                                 C
C       using  the Dirichlet's boundary  conditions:              C
C                                                                 C
C         PSI(0,y) = PSI(Lx,y) = PSI(x,0) = PSI(x,Ly) = 0         C
C                                                                 C
C       in a rectangular domain  with dimensions Lx, Ly.          C
C                      ---------------------------                C
C       Note that the the array ZETA inputs the right hand        C
C       side of the equation and  returns the solution, PSI.      C
C                                                                 C
C       WSAVE is an working array of dimension int(2.5*MAX+15)    C
C           if II2=JJ2 then WSAVE does not need to be modified    C
C             by calling SINTI again.                             C
C                    ********************                         C
C       THIS SUBROUTINE CALLS SUBROUTINES SINT AND SINTI,  FROM   C
C       NCAR'S PACKAGE: FFTPACK (LINK WITH LIBRARY MYLIBRY.OLB)   C         
C=================================================================C
        
        DIMENSION ZETA(II,*), EIG(II,*), PBOT(*), WSAVE(*)

          PI = 4.*ATAN(1.)
          N1 = N+1
          M1 = M+1
C=================================================================C
C         PART I - COMPUTES TRANSFORM OF ZETA(I,J)                C
C=================================================================C 
          
          DO 60 K = 1,N
          DO 50 J = 1,M
   50     ZETA(J,K) = ZETA(J,K)/(2.*(M1))
C--- FORWARD TRANSFORM:
   60     CALL SINT(M,ZETA(1,K),WSAVE)  

C--- GAUSSIAN ELIMINATION - RIGHT HAND SIDE
          DO 2 K=2,N
          DO 2 J=1,M
    2     ZETA(J,K) = ZETA(J,K)-ZETA(J,K-1)/(PBOT(K)*EIG(J,K-1)) 

          DO 3 J=1,M
    3     ZETA(J,N) = ZETA(J,N)/EIG(J,N)
 
          DO 4 K=N-1,1,-1
          DO 4 J=1,M
    4     ZETA(J,K)=(ZETA(J,K)-ZETA(J,K+1)/PBOT(K+1))/EIG(J,K)  

C--- BACK TRANSFORM:
          DO 6 K=1,N
          CALL SINT(M,ZETA(1,K),WSAVE)
    6     ZETA(M1,K) = 0.

          RETURN
          END
C
c=========================================================================
c=========================================================================
c		INPUT FILE NEEDED WITH PARAMETERS
c=========================================================================
c
c 0,5                         ! nstep1, nstep2
c 1200                        ! delt
c 500.e5,4500.e5              ! dp(1),dp(2)
c 1.0,2                       !  wind stress, 1 --> free slip, 2 --> no-slip
c 5.e6                        !  viscosity in cgs units
c
c=========================================================================
