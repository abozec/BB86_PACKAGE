      subroutine icloan(m,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- 'energy loan' ice model. no advection, no dynamics. ice amount
c --- represents energy 'loaned' to water column to prevent wintertime
c --- cooling below freezing level. loan is paid back in summer.
c
c --- modified version for ice-ocean "coupling".
c --- freeze/melt energy from relaxation to the freezing temperature.
c --- the atmosphere/ice surface exchange is applied to the ocean
c --- (previously done in thermf).
c
      integer i,j,l
      real    tfrz,tsur,tmxl,smxl,hfrz,paybak,borrow,hice,thkimx,t2f
      real    radfl,tdif,wind,airt,rair,snsibl,emnp,dtrmui
      real    thkimxy(jdm)
c
c --- hice   = actual ice thickness (m), local variable
c
c --- thkice = average ice thickness, i.e. hice x covice (m)
c --- covice = ice coverage, i.e. cell fraction (0.0 to 1.0)
c --- temice = ice surface temperature          (degC)
c --- flxice = cell average heat flux under ice (W/m^2)
c --- fswice = cell average swv  flux under ice (W/m^2)
c --- sflice = cell average salt flux under ice
c
c --- icefrq = e-folding time scale back to tfrz (no. time steps)
c --- thkfrz = maximum thickness of near-surface freezing zone (m)
c --- tfrz_0 = ice melting point (degC) at S=0psu
c --- tfrz_s = gradient of ice melting point (degC/psu)
c --- ticegr = vertical temperature gradient inside ice (deg/m)
c ---            (0.0 to get ice surface temp. from atmos. surtmp)
c --- hicemn = minimum ice thickness (m)
c --- hicemx = maximum ice thickness (m)
c
      real       tfrz_n,ticemn,ticemx,salice,rhoice,fusion,meltmx
      parameter (tfrz_n= -1.79, ! nominal ice melting point (degC)
     &           ticemn=-50.0,  ! minimum ice surface temperature (degC)
     &           ticemx=  0.0,  ! maximum ice surface temperature (degC)
     &           salice=  4.0,  ! salinity of ice (psu) - same as CICE
     &           rhoice=917.0,  ! density  of ice (kg/m**3)
     &           fusion=334.e3, ! latent heat of fusion (J/kg)
     &           meltmx= 33.e-7)! max. ice melting rate (m/sec), 0.285 m/day
c
      real       fluxmx         !max. ice melting flux (W/m^2)
      parameter (fluxmx=meltmx*fusion*rhoice)    !~1000 W/m^2 - like CICE
c
      real       csice,csubp,pairc,rgas,tzero
      parameter (csice =0.0006,      !ice-air sensible exchange coefficient
     &           csubp =1005.7,      !specific heat of air (j/kg/deg)
     &           pairc=1013.0*100.0, !air pressure (mb) * 100
     &           rgas =287.1,        !gas constant (j/kg/k)
     &           tzero=273.16)       !celsius to kelvin temperature offset
c
      include 'stmt_fns.h'
c
      dtrmui = delt1/(1.0*86400.0)  !dt*1/1days
c
c --- energy loan: add extra energy to the ocean to keep SST from dropping
c --- below tfrz in winter. return this borrowed energy to the 'energy bank'
c --- in summer.
c
c --- salt loan: analogous to energy loan.
c
      margin = 0  !no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j,l,i,
!$OMP&                    t2f,hfrz,smxl,tmxl,tfrz,borrow,paybak)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        thkimxy(j)=0.0 !simplifies OpenMP parallelization
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c ---       relax to tfrz with e-folding time of icefrq time steps
c ---       assuming the effective surface layer thickness is hfrz
c ---       multiply by dpbl(i,j)/hfrz to get the actual e-folding time
            hfrz   = min( thkfrz*onem, dpbl(i,j) )
            t2f    = (spcifh*hfrz)/(baclin*icefrq*g)
            smxl   = saln(i,j,1,n)
            tmxl   = temp(i,j,1,n)
            tfrz   = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
            borrow = (tfrz-tmxl)*t2f       !W/m^2 into ocean
c
c ---       limit heat flux range (for both forming and melting ice)
            borrow=max( -fluxmx, min( fluxmx, borrow ) )
c
cdiag       if (i.eq.itest .and. j.eq.jtest) then
cdiag         write (lp,'(i9,2i5,a,5f9.3)')
cdiag&          nstep,i+i0,j+j0,'  t,tfrz,flx,hfrz,cov:',
cdiag&          tmxl,tfrz,borrow,hfrz*qonem,covice(i,j)
cdiag       endif
c
            if (tmxl.lt.tfrz) then
c
c ---         add energy to move tmxl towards tfrz (only if tmxl < tfrz)
c
              thkice(i,j)=thkice(i,j)+borrow*(delt1/(fusion*rhoice))
              flxice(i,j)=            borrow
              sflice(i,j)=           +borrow*(smxl-salice)*(1.0/fusion)
            elseif (thkice(i,j).gt.0.0) then  !tmxl > tfrz
c
c ---         ice, so return the borrowed amount whenever tmxl > tfrz
c
              paybak=min( -borrow, thkice(i,j)*(fusion*rhoice/delt1) )
              thkice(i,j)=thkice(i,j)-paybak*(delt1/(fusion*rhoice))
              flxice(i,j)=           -paybak
              sflice(i,j)=           -paybak*(smxl-salice)*(1.0/fusion)
            else !tmxl > tfrz & thkice(i,j) == 0.0
c
c ---         no ice.
c
              flxice(i,j)=0.0
              sflice(i,j)=0.0
c
              if (icmflg.eq.2) then
c
c ---           add extra cooling under the ice mask (tsur<=tfrz_n)
c ---           don't allow a new tsur maximum, to preserve sea ice
c
                if     (yrflag.lt.2) then
                  tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1),
     &                             surtmp(i,j,l2), surtmp(i,j,l3) ),
     &                        surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1+
     &                        surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3   )
                else
                  tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1) ),
     &                        surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1   )
                endif
                if     (tsur.le.tfrz_n) then
                  surflx(i,j)=surflx(i,j)+borrow
                endif
              endif !icmflg.eq.2
            endif
c
            util1(i,j)=max(thkice(i,j)-hicemx,0.0)  !icex = ice exceeding hicemx
            thkimxy(j)=max(thkimxy(j),thkice(i,j))
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
      thkimx=maxval(thkimxy(1:jj))
      call xcmaxr(thkimx)
c
c --- spread out portion of ice thicker than hicemx
      if (thkimx.gt.hicemx) then
        call psmooth(util1, margin)  !smooth icex
      endif
c
!$OMP PARALLEL DO PRIVATE(j,l,i,hice,smxl,tfrz,
!$OMP&                    radfl,tdif,wind,airt,rair,snsibl,emnp)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            thkice(i,j)=util1(i,j)+min(thkice(i,j),hicemx) !icex_sm+rest
c
c ---       compute fractional ice coverage for energy flux calculation
            if (thkice(i,j).lt.1.e-5*hicemn) then
              covice(i,j)=0.0
            else
              covice(i,j)=min(1.0,thkice(i,j)*(1.0/hicemn))
              hice=thkice(i,j)/covice(i,j)  !minimum of hicemn
            end if
c
            if     (icmflg.eq.3) then
c ---         relax to sea ice concentration from coupler
c ---         ice thickness is therefore always then 0 and hicemn
              covice(i,j)=covice(i,j)+dtrmui*(si_c(i,j)-covice(i,j))
              thkice(i,j)=covice(i,j)*hicemn
              hice=hicemn
            endif
c
c ---       compute ice surface temperature
            if     (covice(i,j).eq.0.0) then
              temice(i,j)=ticemx
            elseif (ticegr.eq.0.0) then  !use surtmp
              temice(i,j)=max( ticemn,
     &                         min( ticemx,
     &                              surtmp(i,j,l0)*w0+
     &                              surtmp(i,j,l1)*w1+
     &                              surtmp(i,j,l2)*w2+
     &                              surtmp(i,j,l3)*w3  ) )
            else
              temice(i,j)=max( ticemn, ticemx-ticegr*hice )
            endif
c
            if     (icmflg.eq.3) then
              if     (min(covice(i,j),si_c(i,j)).gt.0.0) then
                temice(i,j)=max( ticemn, min( ticemx, si_t(i,j) ) )
              endif
            endif
c
c ---       atmosphere to ice surface exchange, applied to the ocean.
c
            if     (covice(i,j).gt.0.0) then
c ---         net radiative thermal flux (w/m**2) +ve into ocean/ice
c ---         radflx's Qsw includes the atmos. model's surface albedo,
c ---         i.e. it already allows for ice&snow where it is observed.
              radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
     &             +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
              if     (lwflag.gt.1) then
c ---           longwave correction to radfl (Qsw+Qlw).
c ---           this will be ~zero for ticegr==0.0 (temice=surtmp)
                tdif = temice(i,j) -
     &               ( surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1
     &                +surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3)
                !correction is blackbody radiation from tdif at temice
                radfl = radfl - (4.506+0.0554*temice(i,j)) * tdif
              endif
              if     (flxflg.ne.3) then
c ---           wind speed (m/s)
                wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
     &              +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
c ---           air temperature (C)
                airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
     &              +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
                rair   = pairc/(rgas*(tzero+airt))
                snsibl = csubp*rair*wind*csice*(temice(i,j)-airt)
              else
                snsibl = 0.0 !already in total flux (i.e. in radfl)
              endif
              flxice(i,j) = flxice(i,j) +
     &                      covice(i,j)*(radfl - snsibl) !no evap
c
c ---         add a time-invarient net heat flux offset
              if     (flxoff) then
                flxice(i,j) = flxice(i,j) + covice(i,j)*offlux(i,j)
              endif
c
c ---         emnp = evaporation minus precipitation (m/sec) into atmos.
c ---         no evap (sublimation) over ice, all precip enters ocean
              if     (pcipf) then
                emnp = -( precip(i,j,l0)*w0+precip(i,j,l1)*w1
     &                   +precip(i,j,l2)*w2+precip(i,j,l3)*w3)
              else
                emnp =  0.0
              endif
              if     (priver) then
                emnp = emnp - ( rivers(i,j,lr0)*wr0+rivers(i,j,lr1)*wr1
     &                         +rivers(i,j,lr2)*wr2+rivers(i,j,lr3)*wr3)
              endif
c ---         sflice = salt flux (10**-3 kg/m**2/sec) into ocean under ice
              sflice(i,j) = sflice(i,j) +
     &                      covice(i,j)*emnp*(saln(i,j,1,n)*qthref)
            endif !covice
c
            fswice(i,j) = 0.0 !no penetrating Qsw under ice
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
      return
      end subroutine icloan
c
c
c> Revision history
c>
c> June 2000 - conversion to SI units
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> May  2003 - added option to impose an ice mask
c> June 2003 - added 8 time step e-folding time scale
c> June 2003 - limited rate of ice formation
c> June 2003 - replaced constant saldif with smxl-salice
c> Mar. 2005 - freezing point linearly dependent on salinity
c> Mar. 2005 - ice surface temperature optionally from surtmp
c> Jun. 2006 - modified version for ice-ocean "coupling"
