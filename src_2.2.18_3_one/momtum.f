      subroutine momtum_hs(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      use mod_tides ! HYCOM tides
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- -----------------------------------------
c --- hydrostatic equation (and surface stress)
c --- -----------------------------------------
c
      logical, parameter :: lpipe_momtum=.false.      !usually .false.
      real,    parameter :: dragw_rho=0.00536*1026.0  !ice-ocean drag from CICE
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &                 stress,stresx,stresy,dpmx,thkbop,oneta,
     &                 vis2u,vis4u,vis2v,vis4v,vort,
     &                 wgtia,wgtib,wgtja,wgtjb,
     &                 dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib,
     &                 tclmn,tclpn  !see momtum4
      common/momtumr4/ stress,stresx,stresy,dpmx,thkbop,oneta,
     &                 vis2u,vis4u,vis2v,vis4v,vort,
     &                 wgtia,wgtib,wgtja,wgtjb,
     &                 dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib,
     &                 tclmn,tclpn  !see momtum4
      save  /momtumr4/
c
      real    dpdn,dpup,q,simo,uimo,vimo,dpsur,psur,usur,vsur,
     &        sumdp,sumth
      integer i,j,k,l,mbdy,hlstep
c
*     real*8    wtime
*     external  wtime
*     real*8    wtime1(10),wtime2(20,kdm),wtimes
c
      include 'stmt_fns.h'
c
      mbdy = 6
c
      call xctilr(pu(    1-nbdy,1-nbdy,2  ),1,  kk, 6,6, halo_us)
      call xctilr(pv(    1-nbdy,1-nbdy,2  ),1,  kk, 6,6, halo_vs)
      call xctilr(dpmixl(1-nbdy,1-nbdy,  m),1,   1, 6,6, halo_ps)
      call xctilr(dp(    1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_ps)
      call xctilr(dpu(   1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_us)
      call xctilr(dpv(   1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_vs)
      call xctilr(pbavg( 1-nbdy,1-nbdy,1  ),1,   3, 6,6, halo_ps)
      call xctilr(dpoldm(1-nbdy,1-nbdy,1  ),1,  kk, 6,6, halo_ps)
      call xctilr(saln(  1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_ps)
      call xctilr(temp(  1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_ps)
      call xctilr(th3d(  1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_ps)
c
c --- tidal forcing
c
      if(tidflg.eq.2) then
         hlstep=0
         call tides_force(hlstep)
      endif
c
c --- hydrostatic equation (and surface stress)
c
*        wtime1( 1) = wtime()
c
c --- rhs: th3d.m, temp.m, saln.m, p, pbavg.m
c --- lhs: thstar, p, oneta, montg
c
      margin = mbdy
c
!$OMP PARALLEL DO PRIVATE(j,l,k,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            sumdp = 0.0
            sumth = 0.0
            do k=1,kk
              if (kapref.ne.0) then  !thermobaric
c
c ---   sigma-star is virtual potential density, as defined in 
c ---   Sun et.al. (1999), 'Inclusion of thermobaricity in 
c ---   isopycnic-coordinate ocean models', JPO 29 pp 2719-2729.
c
c ---   use upper interface pressure in converting sigma to sigma-star.
c ---   to avoid density variations in layers intersected by bottom
c
                if     (kapref.gt.0) then
                  thstar(i,j,k,1)=th3d(i,j,k,m)
     &                             +kappaf(temp(i,j,k,m),
     &                                     saln(i,j,k,m),
     &                                     th3d(i,j,k,m)+thbase,
     &                                        p(i,j,k),
     &                                     kapref)
                else
                  thstar(i,j,k,1)=th3d(i,j,k,m)
     &                             +kappaf(temp(i,j,k,m),
     &                                     saln(i,j,k,m),
     &                                     th3d(i,j,k,m)+thbase,
     &                                        p(i,j,k),
     &                                     2)
                  thstar(i,j,k,2)=th3d(i,j,k,m)
     &                             +kappaf(temp(i,j,k,m),
     &                                     saln(i,j,k,m),
     &                                     th3d(i,j,k,m)+thbase,
     &                                        p(i,j,k),
     &                                     kapi(i,j))
                endif !kapref
              else  !non-thermobaric
                thstar(i,j,k,1)=th3d(i,j,k,m)
              endif  !thermobaric:else
c
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,m)
c
              if     (sshflg.ne.0) then
                sumth = sumth + dp(i,j,k,m)*th3d(i,j,k,m)
                sumdp = sumdp + dp(i,j,k,m)
              endif !sshflg
            enddo !k
c
            if     (sshflg.ne.0) then
              sumth = sumth / max( sumdp, onemm )  !vertical mean of th3d
              sumdp = sumdp*qonem * g              !depth(m) * g
              steric(i,j) =  sshgmn(i,j) +
     &                      (sshgmn(i,j) + sumdp) *
     &                      (thmean(i,j) - sumth) /
     &                      (1000.0+thbase+sumth)
            endif !sshflg
c
c ---       store (1+eta) (= p_total/p_prime) in -oneta-
            oneta(i,j)=1.+pbavg(i,j,m)/p(i,j,kk+1)
c
c ---       m_prime in lowest layer:
            montg(i,j,kk,1)=psikk(i,j,1)+
     &        ( p(i,j,kk+1)*(thkk(i,j,1)-thstar(i,j,kk,1))
     &          -pbavg(i,j,m)*(thstar(i,j,kk,1)+thbase) )*thref**2
            if     (kapref.eq.-1) then
              montg(i,j,kk,2)=psikk(i,j,2)+
     &          ( p(i,j,kk+1)*(thkk(i,j,2)-thstar(i,j,kk,2))
     &            -pbavg(i,j,m)*(thstar(i,j,kk,2)+thbase) )*thref**2
            endif !kapref.eq.-1
c
c ---       m_prime in remaining layers:
            do k=kk-1,1,-1
              montg(i,j,k,1)=montg(i,j,k+1,1)+p(i,j,k+1)*oneta(i,j)
     &            *(thstar(i,j,k+1,1)-thstar(i,j,k,1))*thref**2
              if     (kapref.eq.-1) then
                montg(i,j,k,2)=montg(i,j,k+1,2)+p(i,j,k+1)*oneta(i,j)
     &              *(thstar(i,j,k+1,2)-thstar(i,j,k,2))*thref**2
              endif !kapref.eq.-1
            enddo !k
c
c ---       srfhgt (used diagnostically, in mxmyaij and for tidal SAL).
            if     (kapref.ne.-1) then
              montg1(i,j) =                 montg(i,j,1,1)
            else
              montg1(i,j) =      skap(i,j) *montg(i,j,1,1) +
     &                      (1.0-skap(i,j))*montg(i,j,1,2)
            endif !kapref
            srfhgt(i,j) = montg1(i,j) + thref*pbavg(i,j,m)
c
cdiag       if     (sshflg.ne.0) then
cdiag         if     (itest.gt.0 .and. jtest.gt.0) then
cdiag           write (lp,'(i9,2i5,3x,a,2f12.6,f12.2)')
cdiag&            nstep,itest+i0,jtest+j0,
cdiag&              'sssh =',
cdiag&              steric(i,j),sshgmn(i,j),sumdp
cdiag           write (lp,'(i9,2i5,3x,a,3f12.6)')
cdiag&            nstep,itest+i0,jtest+j0,
cdiag&              'thmn =',
cdiag&              sumth,thmean(i,j),1000.0+thbase+sumth
cdiag           write (lp,'(i9,2i5,3x,a,3f12.6)')
cdiag&            nstep,itest+i0,jtest+j0,
cdiag&              'ssh  =',
cdiag&              srfhgt(i,j),steric(i,j),srfhgt(i,j)-steric(i,j)
cdiag         endif !test
cdiag       endif !sshflg
c
c ---       tidal corrections, note that these are not part of montg
c ---       but are included here to simplify the pressure gradient
            if     (tidflg.gt.0 .and. sshflg.eq.0) then !tides
              do k=1,kk
                montg(i,j,k,1)=montg(i,j,k,1)
     &                            -g*etide(i,j)-tidsal*srfhgt(i,j)
                if     (kapref.eq.-1) then
                  montg(i,j,k,2)=montg(i,j,k,2)
     &                              -g*etide(i,j)-tidsal*srfhgt(i,j)
                endif !kapref.eq.-1
              enddo !k
            elseif (tidflg.gt.0 .and. sshflg.ne.0) then !tides
              do k=1,kk
                montg(i,j,k,1)=montg(i,j,k,1)
     &                            -g*etide(i,j)
     &                            -tidsal*(srfhgt(i,j)-steric(i,j))
                if     (kapref.eq.-1) then
                  montg(i,j,k,2)=montg(i,j,k,2)
     &                              -g*etide(i,j)
     &                              -tidsal*(srfhgt(i,j)-steric(i,j))
                endif !kapref.eq.-1
              enddo !k
            endif !tides (sshflg)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
*        wtime1( 2) = wtime()
      call dpudpv(dpu(1-nbdy,1-nbdy,1,m),
     &            dpv(1-nbdy,1-nbdy,1,m),
     &            p,depthu,depthv, max(0,margin-1))
*        wtime1( 3) = wtime()
c
c --- account for temporal smoothing of mid-time dpmixl. calculate the vertical
c --- excursions of the coordinates immediately above and below the mixed
c --- layer base, then vertically interpolate this motion to dpmixl(i,j,m)
c
      if(hybrid .and. mxlkta) then
c
c ---   rhs: dpoldm, dpmixl.m
c ---   lhs: util1, util2
c
        margin = mbdy
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,dpup,dpdn,q)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              util1(i,j)=0.
              util2(i,j)=0.
            enddo !i
            do k=1,kk
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                util1(i,j)=util2(i,j)
                util2(i,j)=util2(i,j)+dpoldm(i,j,k)
                if     (util2(i,j).ge.dpmixl(i,j,m).and.
     &                  util1(i,j).lt.dpmixl(i,j,m)     ) then
                  dpup=p(i,j,k  )-util1(i,j)
                  dpdn=p(i,j,k+1)-util2(i,j)
                  q=(util2(i,j)-dpmixl(i,j,m))/max(onemm,dpoldm(i,j,k))
                  dpmixl(i,j,m)=dpmixl(i,j,m)+(dpdn+q*(dpup-dpdn))
                endif
              enddo !i
            enddo !k
          enddo !l
        enddo !j
!$OMP   END PARALLEL DO
      endif
c
c --- --------------
c --- surface stress
c --- --------------
c
      if     (windf) then
        margin = 0
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,dpsur,psur,usur,vsur,uimo,vimo,simo)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              if     (wndflg.eq.2) then ! tau on p grid
                surtx(i,j)=( taux(i,j,l0)*w0+taux(i,j,l1)*w1
     &                      +taux(i,j,l2)*w2+taux(i,j,l3)*w3)
                surty(i,j)=( tauy(i,j,l0)*w0+tauy(i,j,l1)*w1
     &                      +tauy(i,j,l2)*w2+tauy(i,j,l3)*w3)
              else ! tau on u&v grids - NOT RECOMMEDED
                surtx(i,j)=( (taux(i,j,l0)+taux(i+1,j,l0))*w0
     &                      +(taux(i,j,l1)+taux(i+1,j,l1))*w1
     &                      +(taux(i,j,l2)+taux(i+1,j,l2))*w2
     &                      +(taux(i,j,l3)+taux(i+1,j,l3))*w3)*0.5
                surty(i,j)=( (tauy(i,j,l0)+tauy(i,j+1,l0))*w0
     &                      +(tauy(i,j,l1)+tauy(i,j+1,l1))*w1
     &                      +(tauy(i,j,l2)+tauy(i,j+1,l2))*w2
     &                      +(tauy(i,j,l3)+tauy(i,j+1,l3))*w3)*0.5
              endif
              if     (iceflg.eq.2 .and. si_c(i,j).gt.0.0) then
c ---           average currents over top 10m
                usur = 0.0
                vsur = 0.0
                psur = 0.0
                do k= 1,kk
                  dpsur = min( dp(i,j,k,n), max( 0.0, tenm - psur ) )
                  usur  = usur + dpsur*(u(i,j,k,n)+u(i+1,j,k,n))
                  vsur  = vsur + dpsur*(v(i,j,k,n)+v(i,j+1,k,n))
                  psur  = psur + dpsur
                  if     (dpsur.eq.0.0) then
                    exit
                  endif
                enddo !k
                usur  = 0.5*( usur/psur + ubavg(i,  j,n) +
     &                                    ubavg(i+1,j,n)  )
                vsur  = 0.5*( vsur/psur + vbavg(i,j,  n) +
     &                                    vbavg(i,j+1,n)  )
c               allow for ice-ocean stress
                uimo = si_u(i,j) - usur
                vimo = si_v(i,j) - vsur
                simo = sqrt( uimo**2 + vimo**2 )
                surtx(i,j)=(1.0-si_c(i,j))*surtx(i,j) +
     &                          si_c(i,j) *dragw_rho*simo*uimo
                surty(i,j)=(1.0-si_c(i,j))*surty(i,j) +
     &                          si_c(i,j) *dragw_rho*simo*vimo
              endif !ice-ocean stress
            enddo !i
          enddo !l
        enddo !j
!$OMP   END PARALLEL DO
c
        call xctilr(surtx,1,1, 6,6, halo_ps)
        call xctilr(surty,1,1, 6,6, halo_ps)
      endif !windf
c
      return
      end
c
      subroutine momtum(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      use mod_tides ! HYCOM tides
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ------------------------------------------------------
c --- momentum equations (2nd order version)
c
c --- Enstrophy conserving advection scheme (Sadourney, 1975)
c
c --- diffusion is Laplacian and/or biharmonic, both with
c --- "constant" and deformation dependent coefficients.
c
c --- hydrostatic equation and surface stress via momtum_hs
c --- ------------------------------------------------------
c
      logical, parameter :: lpipe_momtum=.false.  !usually .false.
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &                 stress,stresx,stresy,dpmx,thkbop,oneta,
     &                 vis2u,vis4u,vis2v,vis4v,vort,
     &                 wgtia,wgtib,wgtja,wgtjb,
     &                 dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib,
     &                 tclmn,tclpn  !see momtum4
      common/momtumr4/ stress,stresx,stresy,dpmx,thkbop,oneta,
     &                 vis2u,vis4u,vis2v,vis4v,vort,
     &                 wgtia,wgtib,wgtja,wgtjb,
     &                 dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib,
     &                 tclmn,tclpn  !see momtum4
      save  /momtumr4/
c
      integer ifirst
      save    ifirst
c
      real dpia,dpib,dpja,dpjb,vis2a,vis4a,vis2b,vis4b,
     &     scuya,scuyb,scvxa,scvxb,vmag,dall,adrlim,
     &     dpxy,ptopl,pbotl,cutoff,qcutoff,h1,q,deform,aspy2,aspx2,
     &     dt1inv,phi,plo,pbop,pthkbl,ubot,vbot,pstres,
     &     dmontg,dthstr,dragu,dragv,qdpu,qdpv,dpthin
      integer i,ia,ib,j,ja,jb,k,ka,l,mbdy
c
*     real*8    wtime
*     external  wtime
*     real*8    wtime1(10),wtime2(20,kdm),wtimes
c
      character text*12
      integer, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & mask
c
      real hfharm,a,b
      include 'stmt_fns.h'
c
c --- harmonic mean divided by 2
      hfharm(a,b)=a*b/(a+b)
c
      data ifirst / 0 /
c
      mbdy = 6
c
      call xctilr(u(    1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_uv)
      call xctilr(v(    1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_vv)
      call xctilr(ubavg(1-nbdy,1-nbdy,1  ),1,   3, 6,6, halo_uv)
      call xctilr(vbavg(1-nbdy,1-nbdy,1  ),1,   3, 6,6, halo_vv)
      call xctilr(uflx( 1-nbdy,1-nbdy,1  ),1,  kk, 6,6, halo_uv)
      call xctilr(vflx( 1-nbdy,1-nbdy,1  ),1,  kk, 6,6, halo_vv)
c
      if     (ifirst.eq.0) then
        ifirst=1
c ---   setup zero fill.
        margin = mbdy
c
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            vis2u(i,j)=0.0
            vis4u(i,j)=0.0
            vis2v(i,j)=0.0
            vis4v(i,j)=0.0
            dl2u( i,j)=0.0
            dl2v( i,j)=0.0
          enddo !i
        enddo !j
      endif
c
c --- ---------------------------------------
c --- hydrostatic equation and surface stress
c --- ---------------------------------------
c
      call momtum_hs(m,n)
c
c +++ ++++++++++++++++++
c +++ momentum equations
c +++ ++++++++++++++++++
c
*        wtime1( 4) = wtime()
c
c --- rhs: p, u.n+, v.n+, ubavg.n+, vbavg.n+, depthv+, pvtrop+
c --- rhs: dpmixl.m+, taux+, dpu, depthu+, dpv, tauy+
c --- lhs: util1, util2, drag, ubrhs, stresx, vbrhs, stresy
c
      if     (drglim.gt.0.0) then
        adrlim = drglim
      else
        adrlim = 0.125
      endif
      dt1inv = 1./delt1
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,
!$OMP&                    phi,plo,pbop,ubot,vbot,vmag,dall,pstres)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c ---       bottom drag (standard bulk formula)
c ---       bottom stress is applied over thickness dp00 for the kpp bottom
c ---       b.l. just as for the surface b.l. otherwise, bottom stress is
c ---       applied over thickness thkbot
c
            if (mxlkpp .and. bblkpp) then
              thkbop(i,j)=dp00         !bottom stress applied over this thknss
            else
              thkbop(i,j)=thkbot*onem  !bottom stress applied over this thknss
            endif
c
c --- the bottom stress term is estimated using velocity averaged over the
c --- bottom boundary layer. this thickness is dpbbl for the kpp boundary
c --- layer; otherwise, it is thkbop
            ubot=0.0
            vbot=0.0
            if (mxlkpp .and. bblkpp) then
              pthkbl=max(dpbbl(i,j),thkbop(i,j))  !thknss of bot. b.l.
            else
              pthkbl=thkbop(i,j)                  !thknss of bot. b.l.
            endif
            pbop=p(i,j,kk+1)-pthkbl               !top of bot. b.l.
            phi =max(p(i,j,1),pbop)
            do k=1,kk
              plo =phi  ! max(p(i,j,k),pbop)
              phi =max(p(i,j,k+1),pbop)
              ubot=ubot + (u(i,j,k,n)+u(i+1,j,k,n))*(phi-plo)
              vbot=vbot + (v(i,j,k,n)+v(i,j+1,k,n))*(phi-plo)
            enddo !k
            ubot=ubot/min(pthkbl,p(i,j,kk+1))
     &            + (ubavg(i,j,n)+ubavg(i+1,j,n))
            vbot=vbot/min(pthkbl,p(i,j,kk+1))
     &            + (vbavg(i,j,n)+vbavg(i,j+1,n))
            vmag=0.5*sqrt(ubot**2+vbot**2)
!!Alex add  linear bottom drag cbar2
!            dall=cb*(vmag+cbar)  !no tidal drag
            dall=cb*(vmag+cbar) + cbar2  

            drag(i,j)=dall/min(thkbop(i,j)*qonem,depths(i,j))
            if (mxlkpp .and. bblkpp) then
              ustarb(i,j)=sqrt(dall*vmag)
            endif
c
c ---       tidal bottom drag
c
            util6(i,j)=thkdrg*onem  !bottom stress applied over this thknss
            util5(i,j)=dragrh(i,j)/min(util6(i,j)*qonem,depths(i,j))
          enddo !i
        enddo !l
c
c ---   store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c ---   time-interpolate wind stress
c
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubrhs(i,j)=
     &        (vbavg(i  ,j,  m)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,m)*depthv(i  ,j+1)
     &        +vbavg(i-1,j,  m)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,m)*depthv(i-1,j+1))
     &        *(pvtrop(i,j)+pvtrop(i,j+1))*.125
c
            if     (windf) then
              if(hybrid .and. mxlkrt) then
                pstres=0.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
!!Alex add the use of pstrsi if BB86 config
              else if (pstrsi .gt. 0.) then
                pstres=pstrsi*onem
              else
                pstres=dpu(i,j,1,m)
              endif
c ---         units of surtx are N/m^2 (i.e. Pa)
              stresx(i,j)=(surtx(i,j)+surtx(i-1,j))*0.5*g
     &                    /(pstres*thref)
            else  ! no taux
              stresx(i,j)=0.
            endif !windf:else
          enddo !i
        enddo !l
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbrhs(i,j)=
     &      -(ubavg(i,  j  ,m)*depthu(i,j  )
     &       +ubavg(i+1,j  ,m)*depthu(i+1,j  )
     &       +ubavg(i,  j-1,m)*depthu(i,j-1)
     &       +ubavg(i+1,j-1,m)*depthu(i+1,j-1))
     &       *(pvtrop(i,j)+pvtrop(i+1,j))*.125
c
            if     (windf) then
              if(hybrid .and. mxlkrt) then
                pstres=0.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
!!Alex add the use of pstrsi if BB86 config
              else if (pstrsi .gt. 0.) then
                pstres=pstrsi*onem
              else
                pstres=dpv(i,j,1,m)
              endif
c ---         units of surty are N/m^2 (i.e. Pa)
              stresy(i,j)=(surty(i,j)+surty(i,j-1))*0.5*g
     &                    /(pstres*thref)
            else  ! no tauy
              stresy(i,j)=0.
            endif !windf:else
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'uba.n  n=',n
        call pipe_compare(ubavg(1-nbdy,1-nbdy,n),iu,text)
        write (text,'(a9,i3)') 'vba.n  n=',n
        call pipe_compare(vbavg(1-nbdy,1-nbdy,n),iv,text)
        write (text,'(a9,i3)') 'drag   k=',0
        call pipe_compare(drag,ip,text)
      endif
c
c --- the old  momeq2.f  starts here
c
      dpthin  = 0.001*onemm
      h1      =       tenm  !used in lateral weighting of hor.pres.grad.
      cutoff  =   0.5*onem
      qcutoff = 1.0/cutoff
c
*        wtime1( 5) = wtime()
c
c --- rhs: 0.0
c --- lhs: util1, util2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
c ---   spatial weighting function for pressure gradient calculation:
          util1(i,j)=0.
          util2(i,j)=0.
        enddo !i
      enddo !j
c
      do 9 k=1,kk
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. store minimum thickness
c --- values for use in pot.vort. calculation in -dpmx-.
c
*         wtime2( 1,k) = wtime()
c
c --- rhs: dpmx, dp.m+
c --- lhs: dpmx
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
        dpmx(i,1)=2.*cutoff
      enddo !i
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          dpmx(i,j+1)=2.*cutoff
        enddo !i
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpmx(i,j+1)=max(dpmx(i,j+1),dp(i,j,k,m)+dp(i-1,j,k,m))
          enddo !i
        enddo !l
      enddo !j
c
*         wtime2( 2,k) = wtime()
c
c --- rhs: ubavg.m, ubavg.n, dp.m+, dpu
c --- lhs: utotm, utotn, uflux, dpmx, pu
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
c
          i=ifu(j,l)-1
          if (i.ge.1-margin) then
            if (iuopn(i,j).ne.0) then
              utotm(i,j)=u(i+1,j,k,m)+ubavg(i,j,m)
              utotn(i,j)=u(i+1,j,k,n)+ubavg(i,j,n)
              uflux(i,j)=utotm(i,j)*max(dp(i,j,k,m),cutoff)
            endif
          endif
          i=ilu(j,l)+1
          if (i.le.ii+margin) then
            if (iuopn(i,j).ne.0) then
              utotm(i,j)=u(i-1,j,k,m)+ubavg(i,j,m)
              utotn(i,j)=u(i-1,j,k,n)+ubavg(i,j,n)
              uflux(i,j)=utotm(i,j)*max(dp(i-1,j,k,m),cutoff)
            endif
          endif  
c
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpmx(i,j)=max(dpmx(i,j),dp(i,j,k,m)+dp(i-1,j,k,m))
            utotm(i,j)=u(i,j,k,m)+ubavg(i,j,m)
            utotn(i,j)=u(i,j,k,n)+ubavg(i,j,n)
            uflux(i,j)=utotm(i,j)*max(dpu(i,j,k,m),cutoff)
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,m)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
*         wtime2( 3,k) = wtime()
c
c --- rhs: vbavg.m, vbavg.n, dp.m+, dpv
c --- lhs: vtotm, vtotn, vflux, dpmx, pv
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
        do l=1,jsv(i)
          j=jfv(i,l)-1
          if (j.ge.1-margin) then
            if (ivopn(i,j).ne.0) then
              vtotm(i,j)=v(i,j+1,k,m)+vbavg(i,j,m)
              vtotn(i,j)=v(i,j+1,k,n)+vbavg(i,j,n)
              vflux(i,j)=vtotm(i,j)*max(dp(i,j,k,m),cutoff)
            endif
          endif  
          j=jlv(i,l)+1
          if (j.le.jj+margin) then
            if (ivopn(i,j).ne.0) then
              vtotm(i,j)=v(i,j-1,k,m)+vbavg(i,j,m)
              vtotn(i,j)=v(i,j-1,k,n)+vbavg(i,j,n)
              vflux(i,j)=vtotm(i,j)*max(dp(i,j-1,k,m),cutoff)
            endif
          endif  
        enddo !l
      enddo !i
c
*         wtime2( 4,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            dpmx(i  ,j)=max(dpmx(i  ,j),dp(i,j,k,m)+dp(i,j-1,k,m))
            dpmx(i+1,j)=max(dpmx(i+1,j),dp(i,j,k,m)+dp(i,j-1,k,m))
            vtotm(i,j)=v(i,j,k,m)+vbavg(i,j,m)
            vtotn(i,j)=v(i,j,k,n)+vbavg(i,j,n)
            vflux(i,j)=vtotm(i,j)*max(dpv(i,j,k,m),cutoff)
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,m)
          enddo !i
        enddo !l
      enddo !j
c
c --- define auxiliary velocity fields (via,vib,uja,ujb) to implement
c --- sidewall friction along near-vertical bottom slopes. wgtja,wgtjb,wgtia,
c --- wgtib indicate the extent to which a sidewall is present.
c
*         wtime2( 5,k) = wtime()
c
c --- rhs: pu, depthu+, utotn+, wgtja
c --- lhs: wgtja, wgtjb, uja, ujb, dl2u
c --- rhs: pv, depthv+, vtotn+, wgtia
c --- lhs: wgtia, wgtib, via, vib, dl2v
c --- rhs: vtotm, vort+, corio+, dp.m+, dpmx+, vtotn
c --- lhs: vort, potvor, defor2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,ia,ib,aspy2,aspx2)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c ---   assume margin<nblk
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            wgtja(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,ja))
     &                /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
            wgtjb(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,jb))
     &                /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
            uja(i,j)=(1.-wgtja(i,j))*utotn(i,ja)+
     &               wgtja(i,j)*slip*utotn(i,j)
            ujb(i,j)=(1.-wgtjb(i,j))*utotn(i,jb)+
     &               wgtjb(i,j)*slip*utotn(i,j)
c ---       Laplacian of utotn scaled by -0.25*max(scux,scuy)**2
            aspx2 = aspux(i,j)**2
            aspy2 = aspuy(i,j)**2
            dl2u(i,j)=0.5*(aspx2+aspy2)*utotn(i,j)
     &                   -0.25*aspx2*(utotn(i+1,j)+utotn(i-1,j))
     &                   -0.25*aspy2*(  uja(i,  j)+  ujb(i,  j))
          enddo !i
        enddo !l
c
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c ---       assume margin<nblk
            ia=i-1
            ib=i+1
            wgtia(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(ia,j))
     &                /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
            wgtib(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(ib,j))
     &                /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
            via(i,j)=(1.-wgtia(i,j))*vtotn(ia,j)+
     &                wgtia(i,j)*slip*vtotn(i,j)
            vib(i,j)=(1.-wgtib(i,j))*vtotn(ib,j)+
     &                wgtib(i,j)*slip*vtotn(i,j)
c ---       Laplacian of vtotn scaled by -0.25*max(scvx,scvy)**2
            aspx2 = aspvx(i,j)**2
            aspy2 = aspvy(i,j)**2
            dl2v(i,j)=0.5*(aspx2+aspy2)*vtotn(i,j)
     &                 -0.25*aspy2*(vtotn(i,j+1)+vtotn(i,j-1))
     &                 -0.25*aspx2*(  via(i,j  )+  vib(i,j  ))
          enddo !i
        enddo !l
c
c --- vorticity, pot.vort., defor. at lateral boundary points
        do l=1,isv(j)
          i=ifv(j,l)
          if     (i.ge. 1-margin) then
            vort(i  ,j)= vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i  ,j)
            potvor(i  ,j)=(vort(i  ,j)+corio(i  ,j)) * 8.
     &                     /max(8.*cutoff,
     &                          4.*(dp(i,j,k,m)+dp(i,ja ,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i+1,j))
            defor2(i  ,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*
     &                      scq2i(i  ,j)
          endif
          i=ilv(j,l)
          if     (i.le.ii+margin) then
            vort(i+1,j)=-vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i+1,j)
            potvor(i+1,j)=(vort(i+1,j)+corio(i+1,j)) * 8.
     &                     /max(8.*cutoff,
     &                          4.*(dp(i,j,k,m)+dp(i,ja ,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i+1,j))
            defor2(i+1,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*
     &                      scq2i(i+1,j)
          endif
        enddo !i
      enddo !l
!$OMP END PARALLEL DO
c
*         wtime2( 6,k) = wtime()
c
c --- vorticity, pot.vort., defor. at lateral boundary points
c
c --- rhs: utotm, vort+, corio+, dp.m+, dpmx+, utotn
c --- lhs: vort, potvor, defor2
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
c ---   assume margin<nblk
        ia=i-1
        do l=1,jsu(i)
          j=jfu(i,l)
          if     (j.ge. 1-margin) then
            vort(i,j  )=-utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j  )
            potvor(i,j  )=(vort(i,j  )+corio(i,j  )) * 8.
     &                     /max(8.*cutoff,
     &                          4.*(dp(i,j,k,m)+dp(ia ,j,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i,j+1))
            defor2(i,j  )=(utotn(i,j)*(1.-slip)*scux(i,j))**2*
     &                      scq2i(i,j  )
          endif
          j=jlu(i,l)
          if     (j.le.jj+margin) then
            vort(i,j+1)= utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j+1)
            potvor(i,j+1)=(vort(i,j+1)+corio(i,j+1)) * 8.
     &                     /max(8.*cutoff,
     &                          4.*(dp(i,j,k,m)+dp(ia ,j,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i,j+1))
            defor2(i,j+1)=(utotn(i,j)*(1.-slip)*scux(i,j))**2*
     &                      scq2i(i,j+1)
          endif
        enddo !l
      enddo !i
c
*         wtime2( 7,k) = wtime()
c
c --- rhs: p, utotn+, vtotn+
c --- lhs: util3, defor1
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            util3(i,j)=0.5*(p(i,j,k+1)+p(i,j,k))*oneta(i,j)

            defor1(i,j)=((utotn(i+1,j)*scuy(i+1,j)
     &                   -utotn(i,  j)*scuy(i,  j))
     &                  -(vtotn(i,j+1)*scvx(i,j+1)
     &                   -vtotn(i,j  )*scvx(i,j  )))**2
     &                  *scp2i(i,j)
          enddo !i
        enddo !l
      enddo !j
c
c --- vorticity, pot.vort., defor. at interior points (incl. promontories)
*         wtime2( 8,k) = wtime()
c
c --- rhs: vtotm+, utotm+, vort, dp.m+, dpmx+, vib+, via, ujb+, uja
c --- lhs: vort, potvor, defor2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isq(j)
          do i=max(1-margin,ifq(j,l)),min(ii+margin,ilq(j,l))
            vort(i,j)=(vtotm(i,j)*scvy(i,j)-vtotm(i-1,j)*scvy(i-1,j)
     &                -utotm(i,j)*scux(i,j)+utotm(i,j-1)*scux(i,j-1))
     &                *scq2i(i,j)
            potvor(i,j)=(vort(i,j)+corio(i,j)) * 8.
     &         /max(8.*cutoff,
     &              2.*(dp(i,j  ,k,m)+dp(i-1,j  ,k,m)+
     &                  dp(i,j-1,k,m)+dp(i-1,j-1,k,m)),
     &              dpmx(i,j),dpmx(i-1,j),dpmx(i+1,j),
     &                        dpmx(i,j-1),dpmx(i,j+1))
            defor2(i,j)=(vib(i-1,j)*scvy(i,j)-via(i,j)*scvy(i-1,j)
     &                  +ujb(i,j-1)*scux(i,j)-uja(i,j)*scux(i,j-1))**2
     &                  *scq2i(i,j)
          enddo !i
        enddo !l
      enddo !j
c
c --- define auxiliary del2 fields (dl2via,dl2vib,dl2uja,dl2ujb) to imple-
c --- ment biharmonic sidewall friction along near-vertical bottom slopes.
c
*         wtime2( 9,k) = wtime()
c
c --- rhs: wgtja, wgtjb, dlu2+, wgtia, wgtib, dlv2+
c --- lhs: dl2uja, dl2ujb, dl2via, dl2vib
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dl2uja(i,j)=(1.-wgtja(i,j))*dl2u(i,ja)+
     &                  wgtja(i,j)*slip*dl2u(i,j)
            dl2ujb(i,j)=(1.-wgtjb(i,j))*dl2u(i,jb)+
     &                  wgtjb(i,j)*slip*dl2u(i,j)
          enddo !i
        enddo !l
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
            dl2via(i,j)=(1.-wgtia(i,j))*dl2v(ia,j)+
     &                  wgtia(i,j)*slip*dl2v(i, j)
            dl2vib(i,j)=(1.-wgtib(i,j))*dl2v(ib,j)+
     &                  wgtib(i,j)*slip*dl2v(i, j)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
*         wtime2(10,k) = wtime()
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        do j=1,jj
          do i=1,ii
            mask(i,j)=min(1,iq(i,j)+iu(i,j  )+iv(i,j  )
     &                             +iu(i,j-1)+iv(i-1,j))
          enddo !i
        enddo !j
        write (text,'(a9,i3)') 'potvor k=',k
        call pipe_compare(potvor,mask,text)
      endif
c
c --- ----------
c --- u equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
*         wtime2(11,k) = wtime()
c
c --- rhs: defor1+, defor2+
c --- lhs: vis2u,vis4u
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,deform)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            deform = sqrt(0.5*(defor1(i,j)+defor1(i-1,j)+
     &                         defor2(i,j)+defor2(i,j+1) ))
            vis2u(i,j)=max(veldf2u(i,j),visco2*deform)
            vis4u(i,j)=max(veldf4u(i,j),visco4*deform)
*           vis2u(i,j)=max(veldf2*aspux(i,j),visco2*deform)
*           vis4u(i,j)=max(veldf4*aspux(i,j),visco4*deform)
          enddo !i
        enddo !l
      enddo !j
c
*         wtime2(12,k) = wtime()
c
c --- rhs: vis2u+, vis4u+, dl2u+, dpu.m+, dl2uja, wgtja,, dl2ujb, wgtjb
c --- lhs: vis2u, vis4u, uflux1, uflux2, uflux3
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,
!$OMP&                    dpxy,dpja,dpjb,
!$OMP&                    vis2a,vis4a,vis2b,vis4b,scuya,scuyb)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          if     (ifu(j,l).gt. 1-nbdy) then
            vis2u(ifu(j,l)-1,j)=vis2u(ifu(j,l),j)
            vis4u(ifu(j,l)-1,j)=vis4u(ifu(j,l),j)
          endif
          if     (ilu(j,l).lt.ii+nbdy) then
            vis2u(ilu(j,l)+1,j)=vis2u(ilu(j,l),j)
            vis4u(ilu(j,l)+1,j)=vis4u(ilu(j,l),j)
          endif
c
c ---     longitudinal turb. momentum flux (at mass points)
c ---     note that hfharm is 0.5*harmonic mean
c
          do i=max(1-margin,ifu(j,l)-1),min(ii+margin,ilu(j,l))
            uflux1(i,j)=
     &        ((vis2u(i,j)+vis2u(i+1,j))*(utotn(i,j)-utotn(i+1,j))+
     &         (vis4u(i,j)+vis4u(i+1,j))*(dl2u( i,j)-dl2u( i+1,j)) )
     &                   *hfharm(max(dpu(i  ,j,k,m),onemm),
     &                           max(dpu(i+1,j,k,m),onemm))
     &                   *scp2(i,j)*2./(scux(i,j)+scux(i+1,j))
          enddo !i
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpxy=max(dpu(i,j ,k,m),onemm)
            dpja=max(dpu(i,ja,k,m),onemm)
            dpjb=max(dpu(i,jb,k,m),onemm)
c
c --- check whether variables along coast have been initialized correctly
cdiag      if (k.eq.kk) then
cdiag        if (iu(i,ja).eq.0 .and. dpu(i,ja,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag          write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &      '  error - nonzero dpu(ja):',dpu(i,ja,k,m)
cdiag!$OMP END CRITICAL
cdiag        endif
cdiag        if (iu(i,jb).eq.0 .and. dpu(i,jb,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag          write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &      '  error - nonzero dpu(jb):',dpu(i,jb,k,m)
cdiag!$OMP END CRITICAL
cdiag        endif
cdiag      endif
c
            if (iu(i,ja).eq.0) then
              vis2a=vis2u(i,j )
              vis4a=vis4u(i,j )
              scuya=scuy( i,j )
            else
              vis2a=vis2u(i,ja)
              vis4a=vis4u(i,ja)
              scuya=scuy( i,ja)
            endif
            if (iu(i,jb).eq.0) then
              vis2b=vis2u(i,j )
              vis4b=vis4u(i,j )
              scuyb=scuy( i,j )
            else
              vis2b=vis2u(i,jb)
              vis4b=vis4u(i,jb)
              scuyb=scuy( i,jb)
            endif
            uflux2(i,j)=((vis2u(i,j)+vis2a)*(uja(   i,j)-utotn(i,j))+
     &                   (vis4u(i,j)+vis4a)*(dl2uja(i,j)-dl2u( i,j)) )
     &                  *hfharm(dpja+wgtja(i,j)*(dpxy-dpja),dpxy)
     &                  *scq2(i,j )*2./(scuy(i,j)+scuya)
            uflux3(i,j)=((vis2u(i,j)+vis2b)*(utotn(i,j)-ujb(   i,j))+
     &                   (vis4u(i,j)+vis4b)*(dl2u( i,j)-dl2ujb(i,j)) )
     &                  *hfharm(dpjb+wgtjb(i,j)*(dpxy-dpjb),dpxy)
     &                  *scq2(i,jb)*2./(scuy(i,j)+scuyb)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
c --- pressure force in x direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
*         wtime2(13,k) = wtime()
c
c --- rhs: depthu, pu, montg+, thstar+, p+, dp.m+
c --- lhs: util1, pgfx
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,dmontg,dthstr,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            util1(i,j)=max(0.,min(depthu(i,j)-pu(i,j,k),h1))
            if     (kapref.ne.-1) then
              dmontg=montg( i,j,k,1)-montg( i-1,j,k,1)
              dthstr=thstar(i,j,k,1)-thstar(i-1,j,k,1)
            else !2 thermobaric reference states
              q=0.5*(skap(i,j)+skap(i-1,j))
              dmontg=     q *(montg( i,j,k,1)-montg( i-1,j,k,1))+
     &               (1.0-q)*(montg( i,j,k,2)-montg( i-1,j,k,2))
              dthstr=     q *(thstar(i,j,k,1)-thstar(i-1,j,k,1))+
     &               (1.0-q)*(thstar(i,j,k,2)-thstar(i-1,j,k,2))
            endif
            pgfx(i,j)=util1(i,j)*
     &          ( dmontg+
     &            dthstr*(p(i,j,k+1)*p(i-1,j,k+1)-
     &                    p(i,j,k  )*p(i-1,j,k  ))*thref**2
     &                      /(dp(i,j,k,m)+dp(i-1,j,k,m)+epsil) )


          enddo !i
        enddo !l
      enddo !j
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'pu     k=',k
        call pipe_compare(pu(1-nbdy,1-nbdy,k),iu,text)
        write (text,'(a9,i3)') 'depthu k=',k
        call pipe_compare(depthu,iu,text)
        write (text,'(a9,i3)') 'util1  k=',k
        call pipe_compare(util1, iu,text)
        util4(1:ii,1:jj) = montg( 1:ii,1:jj,k,1)-montg( 0:ii-1,1:jj,k,1)
        write (text,'(a9,i3)') 'montgX k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = thstar(1:ii,1:jj,k,1)-thstar(0:ii-1,1:jj,k,1)
        write (text,'(a9,i3)') 'thstaX k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k)*p(0:ii-1,1:jj,k)
        write (text,'(a9,i3)') 'pSQ    k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k+1)*p(0:ii-1,1:jj,k+1)
        write (text,'(a9,i3)') 'pSQ    k=',k+1
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = dp(1:ii,1:jj,k,m)+dp(0:ii-1,1:jj,k,m)+epsil
        write (text,'(a9,i3)') 'dp.m+  k=',k
        call pipe_compare(util4, iu,text)
        write (text,'(a9,i3)') 'pgfx   k=',k
        call pipe_compare(pgfx,  iu,text)
      endif
c
*         wtime2(14,k) = wtime()
c
c --- rhs: pgfx+, util1+, p+, dpmixl.m+, dpu.m, utotn, drag+
c --- rhs: stresx, u.m, u.n, dpu.n, gradx, utotm+, vtotm+,
c --- rhs: vflux+, potvor+, uflux1+, uflux2, uflux3
c --- lhs: gradx, stress, u.m, u.n
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,
!$OMP&                    ptopl,pbotl,pstres,pbop,qdpu,dragu)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
c --- check whether variables along coast have been initialized correctly
cdiag       if (k.eq.kk) then
cdiag         if (iu(i,ja).eq.0 .and. pgfx(i,ja).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &       '  error - nonzero pgfx(ja):',pgfx(i,ja)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag         if (iu(i,jb).eq.0 .and. pgfx(i,jb).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &       '  error - nonzero pgfx(jb):',pgfx(i,jb)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag       endif
c
            gradx(i,j)=(pgfx(i,j)+(h1-util1(i,j))*
     &        (pgfx (i-1,j)+pgfx (i+1,j)+pgfx (i,ja)+pgfx (i,jb))/
     &        (util1(i-1,j)+util1(i+1,j)+util1(i,ja)+util1(i,jb)+
     &          epsil))/h1
c
            ptopl=min(depthu(i,j),0.5*(p(i,j,k  )+p(i-1,j,k  )))
            pbotl=min(depthu(i,j),0.5*(p(i,j,k+1)+p(i-1,j,k+1)))
            if(hybrid .and. mxlkrt) then
              pstres=0.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
!!Alex add the use of pstrsi if BB86 config
            else if (pstrsi .gt. 0.) then
                pstres=pstrsi*onem
            else
              pstres=dpu(i,j,1,m)
            endif
c
c ---       top and bottom boundary layer stress
            pbop=depthu(i,j)-0.5*(thkbop(i,j)+thkbop(i-1,j)) !top of bot. b.l.
            qdpu=1.0/max(dpu(i,j,k,m),onemm)
            dragu=0.5*(drag(i,j)+drag(i-1,j))*
     &                ( max(pbop,    pbotl)
     &                 -max(pbop,min(ptopl,pbotl-onemm)))*qdpu
            if     (drglim.gt.0.0) then
c             explicit drag: u.t+1 - u.t-1 = - 2dt*dragu*u.t-1
c                   limiter: 2dt*dragu <= drglim
              stress(i,j)=-utotn(i,j)*min(dragu,drglim*dt1inv) +
     &             (stresx(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpu
            else
c             implicit drag: u.t+1 - u.t-1 = - 2dt*dragu*u.t+1
c                                          = - 2dt*dragu / (1+2dt*dragu) * u.t+1
              stress(i,j)=-utotn(i,j)*dragu/(1.0+delt1*dragu) +
     &             (stresx(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpu
            endif
c
            if     (tidflg.gt.0 .and. drgscl.ne.0.0) then
              pbop=depthu(i,j)-0.5*(util6(i,j)+util6(i-1,j)) !top of bot. b.l.
              qdpu=1.0/max(dpu(i,j,k,m),onemm)
              dragu=0.5*(util5(i,j)+util5(i-1,j))*
     &                  ( max(pbop,    pbotl)
     &                   -max(pbop,min(ptopl,pbotl-onemm)))*qdpu
              if     (drglim.gt.0.0) then
                stress(i,j)=stress(i,j)
     &                      -utotn(i,j)*min(dragu,drglim*dt1inv)
              else
                stress(i,j)=stress(i,j)
     &                      -utotn(i,j)*dragu/(1.0+delt1*dragu)
              endif
c ---         !anti-drag on detided velocity (always explicit)
              stress(i,j)=stress(i,j)
     &                      +untide(i,j)*min(dragu,adrlim*dt1inv)
            endif !tidflg
c
c ---       time smoothing of -u- field  (part 1)
            u(i,j,k,m)=u(i,j,k,m)*(wuv1*dpu(i,j,k,m)+onemm)
     &                +u(i,j,k,n)* wuv2*dpu(i,j,k,n)
c
cdiag       util4(i,j)=u(i,j,k,n)
            u(i,j,k,n)=u(i,j,k,n)+delt1*(-scuxi(i,j)*(gradx(i,j)
     &      +0.25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,j+1)**2
     &            -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,j+1)**2))
     &      +0.125*(vflux(i  ,j)+vflux(i  ,j+1)+
     &              vflux(i-1,j)+vflux(i-1,j+1))
     &            *(potvor(i,j)+potvor(i,j+1))
     &      - ubrhs(i,j) + stress(i,j)
     &      -(uflux1(i,j)-uflux1(i-1,j)
     &       +uflux3(i,j)-uflux2(i  ,j))*qdpu/scu2(i,j) )
          enddo !i
        enddo !l
      enddo !j
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'vtotm  k=',k
        call pipe_compare(vtotm, iv,text)
        write (text,'(a9,i3)') 'utotm  k=',k
        call pipe_compare(utotm, iu,text)
        write (text,'(a9,i3)') 'gradx  k=',k
        call pipe_compare(gradx, iu,text)
        write (text,'(a9,i3)') 'ubrhs  k=',k
        call pipe_compare(ubrhs, iu,text)
        write (text,'(a9,i3)') 'stress k=',k
        call pipe_compare(stress,iu,text)
        write (text,'(a9,i3)') 'u.n    k=',k
        call pipe_compare(u(1-nbdy,1-nbdy,k,n),iu,text)
      endif
c
cdiag if (k.eq.1) then
cdiag write (lp,100) nstep
cdiag do j=jtest-1,jtest+1
cdiag do i=itest-1,itest+1
cdiag write (lp,'(2i5,i3,8f8.3)') i,j,k,
cdiag&  util4(i,j),u(i,j,k,n),-delt1*gradx(i,j)*scuxi(i,j),
cdiag&  -delt1*
cdiag& 0.25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,j+1)**2
cdiag&      -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,j+1)**2)*
cdiag&       scuxi(i,j),
cdiag&   delt1*(vflux(i  ,j)+vflux(i  ,j+1)+vflux(i-1,j)+vflux(i-1,j+1))
cdiag&        *(potvor(i,j)+potvor(i,j+1))*.125,
cdiag&  -delt1*ubrhs(i,j),delt1*stress(i,j),
cdiag&  -delt1*(uflux1(i,j)-uflux1(i-1,j)
cdiag&         +uflux3(i,j)-uflux2(i,j))/
cdiag&         (scu2(i,j)*max(dpu(i,j,k,m),onemm))
cdiag enddo !i
cdiag enddo !j
cdiag endif
 100    format(i9,8x,'uold    unew   gradp  nonlin   corio',
     &3x,'ubrhs  stress    fric')
c
c --- ----------
c --- v equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
*         wtime2(15,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i,deform)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            deform = sqrt(0.5*(defor1(i,j)+defor1(i,j-1)+
     &                         defor2(i,j)+defor2(i+1,j) ))
            vis2v(i,j)=max(veldf2v(i,j),visco2*deform)
            vis4v(i,j)=max(veldf4v(i,j),visco4*deform)
*           vis2v(i,j)=max(veldf2*aspvy(i,j),visco2*deform)
*           vis4v(i,j)=max(veldf4*aspvy(i,j),visco4*deform)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
*         wtime2(16,k) = wtime()
      do i=1-margin,ii+margin
        do l=1,jsv(i)
          if     (jfv(i,l).gt. 1-nbdy) then
            vis2v(i,jfv(i,l)-1)=vis2v(i,jfv(i,l))
            vis4v(i,jfv(i,l)-1)=vis4v(i,jfv(i,l))
          endif
          if     (jlv(i,l).lt.jj+nbdy) then
            vis2v(i,jlv(i,l)+1)=vis2v(i,jlv(i,l))
            vis4v(i,jlv(i,l)+1)=vis4v(i,jlv(i,l))
          endif
        enddo !l
      enddo !i
c
!$OMP PARALLEL DO PRIVATE(j,l,i,ia,ib,
!$OMP&                    dpxy,dpia,dpib,
!$OMP&                    vis2a,vis4a,vis2b,vis4b,scvxa,scvxb)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c
c ---   longitudinal turb. momentum flux (at mass points)
c ---   note that hfharm is 0.5*harmonic mean
c
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            if     (iv(i,j)+iv(i,j+1).gt.0) then
            vflux1(i,j)=
     &        ((vis2v(i,j)+vis2v(i,j+1))*(vtotn(i,j)-vtotn(i,j+1))+
     &         (vis4v(i,j)+vis4v(i,j+1))*(dl2v( i,j)-dl2v( i,j+1)) )
     &                   *hfharm(max(dpv(i,j  ,k,m),onemm),
     &                           max(dpv(i,j+1,k,m),onemm))
     &                   *scp2(i,j)*2./(scvy(i,j)+scvy(i,j+1))
            endif
          enddo !i
        enddo !l
c
c ---   lateral turb. momentum flux (at vorticity points)
c ---   (left and right fluxes are evaluated separately because of sidewalls)
c
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
            dpxy=max(dpv(i ,j,k,m),onemm)
            dpia=max(dpv(ia,j,k,m),onemm)
            dpib=max(dpv(ib,j,k,m),onemm)
c
c --- check whether variables along coast have been initialized correctly
cdiag       if (k.eq.kk) then
cdiag         if (iv(ia,j).eq.0 .and. dpv(ia,j,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &       '  error - nonzero dpv(ia):',dpv(ia,j,k,m)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag         if (iv(ib,j).eq.0 .and. dpv(ib,j,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &       '  error - nonzero dpv(ib):',dpv(ib,j,k,m)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag       endif
            if (iv(ia,j).eq.0) then
              vis2a=vis2v(i ,j)
              vis4a=vis4v(i ,j)
              scvxa=scvx( i, j)
            else
              vis2a=vis2v(ia,j)
              vis4a=vis4v(ia,j)
              scvxa=scvx( ia,j)
            endif
            if (iv(ib,j).eq.0) then
              vis2b=vis2v(i ,j)
              vis4b=vis4v(i ,j)
              scvxb=scvx( i, j)
            else
              vis2b=vis2v(ib,j)
              vis4b=vis4v(ib,j)
              scvxb=scvx( ib,j)
            endif
            vflux2(i,j)=((vis2v(i,j)+vis2a)*(via(   i,j)-vtotn(i,j))+
     &                   (vis4v(i,j)+vis4a)*(dl2via(i,j)-dl2v( i,j)) )
     &                  *hfharm(dpia+wgtia(i,j)*(dpxy-dpia),dpxy)
     &                  *scq2(i ,j)*2./(scvx(i,j)+scvxa)
            vflux3(i,j)=((vis2v(i,j)+vis2b)*(vtotn(i,j)-vib(   i,j))+
     &                   (vis4v(i,j)+vis4b)*(dl2v( i,j)-dl2vib(i,j)) )
     &                  *hfharm(dpib+wgtib(i,j)*(dpxy-dpib),dpxy)
     &                  *scq2(ib,j)*2./(scvx(i,j)+scvxb)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
c --- pressure force in y direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
*         wtime2(17,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i,dmontg,dthstr,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            util2(i,j)=max(0.,min(depthv(i,j)-pv(i,j,k),h1))
            if     (kapref.ne.-1) then
              dmontg=montg( i,j,k,1)-montg( i,j-1,k,1)
              dthstr=thstar(i,j,k,1)-thstar(i,j-1,k,1)
            else !2 thermobaric reference states
              q=0.5*(skap(i,j)+skap(i,j-1))
              dmontg=     q *(montg( i,j,k,1)-montg( i,j-1,k,1))+
     &               (1.0-q)*(montg( i,j,k,2)-montg( i,j-1,k,2))
              dthstr=     q *(thstar(i,j,k,1)-thstar(i,j-1,k,1))+
     &               (1.0-q)*(thstar(i,j,k,2)-thstar(i,j-1,k,2))
            endif
            pgfy(i,j)=util2(i,j)*
     &          ( dmontg+
     &            dthstr*(p(i,j,k+1)*p(i,j-1,k+1)-
     &                    p(i,j,k  )*p(i,j-1,k  ))*thref**2
     &                      /(dp(i,j,k,m)+dp(i,j-1,k,m)+epsil) )

          enddo !i
        enddo !l
      enddo !j
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'pv     k=',k
        call pipe_compare(pv(1-nbdy,1-nbdy,k),iv,text)
        write (text,'(a9,i3)') 'depthv k=',k
        call pipe_compare(depthv,iv,text)
        write (text,'(a9,i3)') 'util2  k=',k
        call pipe_compare(util2, iv,text)
        util4(1:ii,1:jj) = montg( 1:ii,1:jj,k,1)-montg( 1:ii,0:jj-1,k,1)
        write (text,'(a9,i3)') 'montgY k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = thstar(1:ii,1:jj,k,1)-thstar(1:ii,0:jj-1,k,1)
        write (text,'(a9,i3)') 'thstaY k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k)  *p(1:ii,0:jj-1,k)
        write (text,'(a9,i3)') 'pSQ    k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k+1)*p(1:ii,0:jj-1,k+1)
        write (text,'(a9,i3)') 'pSQ    k=',k+1
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = dp(1:ii,1:jj,k,m)+dp(1:ii,0:jj-1,k,m)+epsil
        write (text,'(a9,i3)') 'dp.m+  k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'pgfy   k=',k
        call pipe_compare(pgfy,  iv,text)
      endif
c
*         wtime2(18,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
c
c --- check whether variables along coast have been initialized correctly
cdiag       if (k.eq.kk) then
cdiag         if (iv(ia,j).eq.0 .and. pgfy(ia,j).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &       '  error - nonzero pgfy(ia):',pgfy(ia,j)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag         if (iv(ib,j).eq.0 .and. pgfy(ib,j).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag     &       '  error - nonzero pgfy(ib):',pgfy(ib,j)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag       endif
c
            grady(i,j)=(pgfy(i,j)+(h1-util2(i,j))*
     &        (pgfy (ia ,j)+pgfy (ib ,j)
     &        +pgfy (i,j-1)+pgfy (i,j+1))/
     &        (util2(ia ,j)+util2(ib ,j)
     &        +util2(i,j-1)+util2(i,j+1)+epsil))/h1
          enddo !i
        enddo !l
      enddo !j
c
*         wtime2(19,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i,
!$OMP&                    ptopl,pbotl,pstres,pbop,qdpv,dragv)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
            ptopl=min(depthv(i,j),0.5*(p(i,j,k  )+p(i,j-1,k  )))
            pbotl=min(depthv(i,j),0.5*(p(i,j,k+1)+p(i,j-1,k+1)))
            if(hybrid .and. mxlkrt) then
              pstres=0.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
!!Alex add the use of pstrsi if BB86 config
            else if (pstrsi .gt. 0.) then
              pstres=pstrsi*onem
            else
              pstres=dpv(i,j,1,m)
            endif
c
c ---       top and bottom boundary layer stress
            pbop=depthv(i,j)-0.5*(thkbop(i,j)+thkbop(i,j-1))  !top of bot. b.l.
            qdpv=1.0/max(dpv(i,j,k,m),onemm)
            dragv=0.5*(drag(i,j)+drag(i,j-1))*
     &                ( max(pbop,    pbotl)
     &                 -max(pbop,min(ptopl,pbotl-onemm)))*qdpv
            if     (drglim.gt.0.0) then
              stress(i,j)=-vtotn(i,j)*min(dragv,drglim*dt1inv) +
     &             (stresy(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpv
            else
              stress(i,j)=-vtotn(i,j)*dragv/(1.0+delt1*dragv) +
     &             (stresy(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpv
            endif
c
            if     (tidflg.gt.0 .and. drgscl.ne.0.0) then
              pbop=depthv(i,j)-0.5*(util6(i,j)+util6(i,j-1))  !top of bot. b.l.
              qdpv=1.0/max(dpv(i,j,k,m),onemm)
              dragv=0.5*(util5(i,j)+util5(i,j-1))*
     &                  ( max(pbop,    pbotl)
     &                   -max(pbop,min(ptopl,pbotl-onemm)))*qdpv
              if     (drglim.gt.0.0) then
                stress(i,j)=stress(i,j)
     &                      -vtotn(i,j)*min(dragv,drglim*dt1inv)
              else
                stress(i,j)=stress(i,j)
     &                      -vtotn(i,j)*dragv/(1.0+delt1*dragv)
              endif
c ---         !anti-drag on detided velocity (always explicit)
              stress(i,j)=stress(i,j)
     &                      +vntide(i,j)*min(dragv,adrlim*dt1inv)
            endif !tidflg
c
c ---       time smoothing of -v- field  (part 1)
            v(i,j,k,m)=v(i,j,k,m)*(wuv1*dpv(i,j,k,m)+onemm)
     &                +v(i,j,k,n)* wuv2*dpv(i,j,k,n)
c
cdiag       util4(i,j)=v(i,j,k,n)
            v(i,j,k,n)=v(i,j,k,n)+delt1*(-scvyi(i,j)*(grady(i,j)
     &      +0.25*(vtotm(i,  j+1)**2+utotm(i,  j  )**2
     &            +utotm(i+1,j  )**2-vtotm(i,  j-1)**2
     &            -utotm(i,  j-1)**2-utotm(i+1,j-1)**2))
     &      -0.125*(uflux(i,  j  )+uflux(i+1,j  )
     &             +uflux(i,  j-1)+uflux(i+1,j-1))
     &            *(potvor(i,j)+potvor(i+1,j))
     &      - vbrhs(i,j) + stress(i,j)
     &      -(vflux1(i,j)-vflux1(i,j-1)
     &       +vflux3(i,j)-vflux2(i,j  ))*qdpv/scv2(i,j) )
          enddo !i
        enddo !l
      enddo !j
c
*     if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        util4(1:ii,1:jj) =  vtotm(1:ii,  2:jj+1)**2
     &                     +utotm(1:ii,  1:jj  )**2
     &                     +utotm(2:ii+1,1:jj  )**2
     &                     -vtotm(1:ii,  0:jj-1)**2
     &                     -utotm(1:ii,  0:jj-1)**2
     &                     -utotm(2:ii+1,0:jj-1)**2
        write (text,'(a9,i3)') 'totm   k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'grady  k=',k
        call pipe_compare(grady, iv,text)
        write (text,'(a9,i3)') 'vbrhs  k=',k
        call pipe_compare(vbrhs, iv,text)
        if     (k.eq.kk) then
          write (text,'(a9,i3)') 'p      k=',k-1
          call pipe_compare(p(1-nbdy,1-nbdy,k-1),ip,text)
          write (text,'(a9,i3)') 'p      k=',k
          call pipe_compare(p(1-nbdy,1-nbdy,k)  ,ip,text)
          write (text,'(a9,i3)') 'drag   k=',k
          call pipe_compare(drag,ip,text)
        endif
        write (text,'(a9,i3)') 'stress k=',k
        call pipe_compare(stress,iv,text)
        util4(1:ii,1:jj) =  uflux(1:ii,  1:jj  )
     &                     +uflux(2:ii+1,1:jj  )
     &                     +uflux(1:ii,  0:jj-1)
     &                     +uflux(2:ii+1,0:jj-1)
        write (text,'(a9,i3)') 'uflux  k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) =  potvor(1:ii,  1:jj  )
     &                     +potvor(2:ii+1,1:jj  )
        write (text,'(a9,i3)') 'potvor k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) =  vflux1(1:ii,  1:jj  )
     &                     -vflux1(1:ii,  0:jj-1)
     &                     +vflux3(1:ii,  1:jj  )
     &                     -vflux2(1:ii,  1:jj  )
        write (text,'(a9,i3)') 'vflux  k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'v.n    k=',k
        call pipe_compare(v(1-nbdy,1-nbdy,k,n),iv,text)
      endif
c
cdiag if (k.eq.1) then
cdiag write (lp,101) nstep
cdiag do j=jtest-1,jtest+1
cdiag do i=itest-1,itest+1
cdiag write (lp,'(2i5,i3,8f8.3)') i,j,k,
cdiag&  util4(i,j),v(i,j,k,n),-delt1*grady(i,j)*scvyi(i,j),
cdiag&  -delt1*
cdiag& 0.25*(vtotm(i,j+1)**2+utotm(i,j  )**2+utotm(i+1,j  )**2
cdiag&      -vtotm(i,j-1)**2-utotm(i,j-1)**2-utotm(i+1,j-1)**2)*
cdiag&       scvyi(i,j),
cdiag&  -delt1*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,j-1)+uflux(i+1,j-1))
cdiag&        *(potvor(i,j)+potvor(i+1,j))*.125,
cdiag&  -delt1*vbrhs(i,j),delt1*stress(i,j),
cdiag&  -delt1*(vflux1(i,j)-vflux1(i,j-1)
cdiag&         +vflux3(i,j)-vflux2(i,j))/
cdiag&         (scv2(i,j)*max(dpv(i,j,k,m),onemm))
cdiag enddo !i
cdiag enddo !j
cdiag endif
 101    format(i9,8x,'vold    vnew   gradp  nonlin   corio',
     &3x,'vbrhs  stress    fric')
c
*         wtime2(20,k) = wtime()
 9    continue  ! k=1,kk
c
*        wtime1( 6) = wtime()
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 14 j=1-margin,jj+margin
        do k=1,kk
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            enddo !i
          enddo !l
c
c ---     compute new -dpu,dpv- field. save old -dpu,dpv- values in -pu,pv-.
c
          do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              pu(i,j,k+1)=dpu(i,j,k,n)
            enddo !i
          enddo !l
c
          do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              pv(i,j,k+1)=dpv(i,j,k,n)
            enddo !i
          enddo !l
        enddo !k
 14   continue !j
!$OMP END PARALLEL DO
c
*        wtime1( 7) = wtime()
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n),
     &            dpv(1-nbdy,1-nbdy,1,n),
     &            p,depthu,depthv, max(0,margin-1))
*        wtime1( 8) = wtime()
c
c --- substitute zero for (u,v) at bottomed massless grid points, and
c --- substitute depth-weighted averages at other massless grid points.
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,ka,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 30 j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ka=1
            do k=kk,2,-1
              if     (dpu(i,j,k,n).ge.dpthin) then
                ka = k  !deepest layer with significant mass
                exit
              endif
              u(i,j,k,n)=0.0
            enddo !k
            utotn(i,j)=u(i,j,1,n)*dpu(i,j,1,n)
            do k=2,ka
              q=min(dpu(i,j,k,m),dpu(i,j,k,n),cutoff)
              u(i,j,k,n)=(u(i,j,k,n)*q+u(i,j,k-1,n)*(cutoff-q))*qcutoff
              utotn(i,j)=utotn(i,j)+u(i,j,k,n)*dpu(i,j,k,n)
            enddo !k
            utotn(i,j)=utotn(i,j)/depthu(i,j)
          enddo !i
        enddo !l
c         
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ka=1
            do k=kk,2,-1
              if     (dpv(i,j,k,n).ge.dpthin) then
                ka = k  !deepest layer with significant mass
                exit
              endif
              v(i,j,k,n)=0.0
            enddo !k
            vtotn(i,j)=v(i,j,1,n)*dpv(i,j,1,n)
            do k=2,ka
              q=min(dpv(i,j,k,m),dpv(i,j,k,n),cutoff)
              v(i,j,k,n)=(v(i,j,k,n)*q+v(i,j,k-1,n)*(cutoff-q))*qcutoff
              vtotn(i,j)=vtotn(i,j)+v(i,j,k,n)*dpv(i,j,k,n)
            enddo !k
            vtotn(i,j)=vtotn(i,j)/depthv(i,j)
          enddo !i
        enddo !l
 30   continue !j
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        do k= 1,kk
          write (text,'(a9,i3)') 'dpv.m  k=',k
          call pipe_compare(dpv(1-nbdy,1-nbdy,k,m),iv,text)
          write (text,'(a9,i3)') 'dpv.n  k=',k
          call pipe_compare(dpv(1-nbdy,1-nbdy,k,n),iv,text)
          write (text,'(a9,i3)') 'v.n(2) k=',k
          call pipe_compare(  v(1-nbdy,1-nbdy,k,n),iv,text)
        enddo !k
        write (text,'(a9,i3)') 'vtotn  k=',kk
        call pipe_compare(vtotn, iv,text)
      endif
c
c --- time smoothing of -u,v- fields  (part 2)
c
*        wtime1( 9) = wtime()
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 22 j=1-margin,jj+margin
        do k=1,kk
          do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              u(i,j,k,n)=u(i,j,k,n)-utotn(i,j)
              u(i,j,k,m)=(u(i,j,k,m)+u(i,j,k,n)*wuv2*dpu(i,j,k,n))/
     &         (wuv1*dpu(i,j,k,m)+onemm+wuv2*(pu(i,j,k+1)+dpu(i,j,k,n)))
            enddo !i
          enddo !l  
c
          do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              v(i,j,k,n)=v(i,j,k,n)-vtotn(i,j)
              v(i,j,k,m)=(v(i,j,k,m)+v(i,j,k,n)*wuv2*dpv(i,j,k,n))/
     &         (wuv1*dpv(i,j,k,m)+onemm+wuv2*(pv(i,j,k+1)+dpv(i,j,k,n)))
            enddo !i
          enddo !l  
        enddo !k
 22   continue !j
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 867 j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            utotn(i,j)=utotn(i,j)*dt1inv
          enddo !i
        enddo !l
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vtotn(i,j)=vtotn(i,j)*dt1inv
          enddo !i
        enddo !l
 867  continue !j
!$OMP END PARALLEL DO
c
*        wtime1(10) = wtime()
*        if     (mod(nstep,100).eq.0) then
* ---      timer printout.
*          do i= 2,10
*            write(lp,'(4x,a,i3,a,f9.5,a)')
*    &        'momtum   section',i,' time =',
*    &        wtime1(i) - wtime1(i-1),' seconds'
*          enddo !i
*          do i= 2,20
*            wtimes=0.0
*            do k=1,kk
*              wtimes = wtimes + (wtime2(i,k)-wtime2(i-1,k))
*            enddo !k
*            write(lp,'(4x,a,i3,a,f9.5,a)')
*    &        'momtum k section',i,' time =',
*    &        wtimes,' seconds'
*          enddo !i
*          call flush(lp)
*        endif
c
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'v.n(3) k=',1
        call pipe_compare(v(1-nbdy,1-nbdy,1,n),iv,text)
      endif
c
      return
      end
c
      subroutine momtum4(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      use mod_tides ! HYCOM tides
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ------------------------------------------------------
c --- momentum equations (4th order version)
c
c --- Wither, N.G., Y.G. Morel, G. Evensen (2006)
c --- Efficiency of high order numerical schemes
c --- for momentum advection.
c --- Ocean Modeling (submitted).
c
c --- diffusion is biharmonic with speed dependent and
c --- deformation dependent coefficients.
c
c --- hydrostatic equation and surface stress via momtum_hs
c --- ------------------------------------------------------
c
      logical, parameter :: lpipe_momtum=.false.  !usually .false.
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &                 stress,stresx,stresy,dpmx,thkbop,oneta,
     &                 visc,
     &                 uwkn,uwkm,advu,diffu,
     &                 vwkn,vwkm,advv,diffv,
     &                 tdun,tdum,tdvn,tdvm,
     &                 tclmm,tclpm,tclmn,tclpn
      common/momtumr4/ stress,stresx,stresy,dpmx,thkbop,oneta,
     &                 visc,
     &                 uwkn,uwkm,advu,diffu,
     &                 vwkn,vwkm,advv,diffv,
     &                 tdun,tdum,tdvn,tdvm,
     &                 tclmm,tclpm,tclmn,tclpn
      save  /momtumr4/
c
      integer, parameter :: ijmax=max(idm,jdm)
c
      real, save, dimension (1-nbdy:ijmax+nbdy) ::
     &       dun,dum,dvn,dvm
      real, save, dimension (1-nbdy:idm+nbdy,kk+1,2) ::
     &       p1
c
      integer i,ia,ib,im1,im2,ip1,ip2,
     &        j,ja,jb,jm1,jm2,jp1,jp2,k,ka,l,mbdy
c
      real    vmag,dall,adrlim,ptopl,pbotl,cutoff,qcutoff,h1,q,
     &        dt1inv,phi,plo,pbop,pthkbl,ubot,vbot,pstres,
     &        dmontg,dthstr,dragu,dragv,qdpu,qdpv,dpthin
c
      real    clpn,clmn,crpn,crmn,clpm,crmm
      real    velmodul,defortot
      real    ulimmax,hmindiff
      real    coefdiv2,coefdiv3,duadvjp1
      real    dvadvjp1
c
*     real*8    wtime
*     external  wtime
*     real*8    wtime1(10),wtime2(20,kdm),wtimes
c
      character text*12
      integer, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & mask
c
      integer ifirst
      save    ifirst
c
      real hfharm,a,b
      include 'stmt_fns.h'
c
c --- harmonic mean divided by 2
      hfharm(a,b)=a*b/(a+b)
c
      data ifirst / 0 /
c
      mbdy = 6
c
      call xctilr(u(      1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_vv)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1,   3, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1,   3, 6,6, halo_vv)
c
      if     (ifirst.eq.0) then
        ifirst=1
c ---   setup zero fill.
        margin = mbdy
c
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            pu(i,j,1) =0.0
            pv(i,j,1) =0.0
          enddo
        enddo
      endif !ifirst
c
c ---------------------------------
c     First initialize fields
c ---------------------------------
c
      coefdiv2=1.0/(32.0*baclin)
      coefdiv3=1.0/ 12.0
      dpthin  = 0.001*onemm
      hmindiff=       tenm
      h1      =   0.5*onem  !used in lateral weighting of hor.pres.grad.
      cutoff  =   0.5*onem
      qcutoff = 1.0/cutoff
c
c --- ---------------------------------------
c --- hydrostatic equation and surface stress
c --- ---------------------------------------
c
      call momtum_hs(m,n)
c
c +++ ++++++++++++++++++
c +++ momentum equations
c +++ ++++++++++++++++++
c
*        wtime1( 4) = wtime()
c
c --- rhs: p, u.n+, v.n+, ubavg.n+, vbavg.n+, depthv+, pvtrop+
c --- rhs: dpmixl.m+, taux+, dpu, depthu+, dpv, tauy+
c --- lhs: util1, util2, drag, ubrhs, stresx, vbrhs, stresy
c
      if     (drglim.gt.0.0) then
        adrlim = drglim
      else
        adrlim = 0.125
      endif
      dt1inv = 1./delt1
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,
!$OMP&                    phi,plo,pbop,ubot,vbot,vmag,dall,pstres)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c     
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c         
c ---       bottom drag (standard bulk formula)
c ---       bottom stress is applied over thickness dp00 for the kpp bottom
c ---       b.l. just as for the surface b.l. otherwise, bottom stress is
c ---       applied over thickness thkbot
c
            if (mxlkpp .and. bblkpp) then
              thkbop(i,j)=dp00         !bottom stress applied over this thknss
            else
              thkbop(i,j)=thkbot*onem  !bottom stress applied over this thknss
            endif
c             
c --- the bottom stress term is estimated using velocity averaged over the
c --- bottom boundary layer. this thickness is dpbbl for the kpp boundary
c --- layer; otherwise, it is thkbop
            ubot=0.0
            vbot=0.0
            if (mxlkpp .and. bblkpp) then
              pthkbl=max(dpbbl(i,j),thkbop(i,j))  !thknss of bot. b.l.
            else
              pthkbl=thkbop(i,j)                  !thknss of bot. b.l.
            endif
            pbop=p(i,j,kk+1)-pthkbl               !top of bot. b.l.
            phi =max(p(i,j,1),pbop)
            do k=1,kk
              plo =phi  ! max(p(i,j,k),pbop)
              phi =max(p(i,j,k+1),pbop)
              ubot=ubot + (u(i,j,k,n)+u(i+1,j,k,n))*(phi-plo)
              vbot=vbot + (v(i,j,k,n)+v(i,j+1,k,n))*(phi-plo)
            enddo !k
            ubot=ubot/min(pthkbl,p(i,j,kk+1))
     &            + (ubavg(i,j,n)+ubavg(i+1,j,n))
            vbot=vbot/min(pthkbl,p(i,j,kk+1))
     &            + (vbavg(i,j,n)+vbavg(i,j+1,n))
            vmag=0.5*sqrt(ubot**2+vbot**2)
            dall=cb*(vmag+cbar)  !no tidal drag
            drag(i,j)=dall/min(thkbop(i,j)*qonem,depths(i,j))
            if (mxlkpp .and. bblkpp) then
              ustarb(i,j)=sqrt(dall*vmag)
            endif
c
c ---       tidal bottom drag
c
            util6(i,j)=thkdrg*onem  !bottom stress applied over this thknss
            util5(i,j)=dragrh(i,j)/min(util6(i,j)*qonem,depths(i,j))
          enddo !i
        enddo !l
c         
c ---   store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c ---   time-interpolate wind stress
c
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubrhs(i,j)=
     &        (vbavg(i  ,j,  m)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,m)*depthv(i  ,j+1)
     &        +vbavg(i-1,j,  m)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,m)*depthv(i-1,j+1))
     &        *(pvtrop(i,j)+pvtrop(i,j+1))*.125
c    
            if     (windf) then
              if(hybrid .and. mxlkrt) then
                pstres=0.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
              else
                pstres=dpu(i,j,1,m)
              endif
c ---         units of surtx are N/m^2 (i.e. Pa)
              stresx(i,j)=(surtx(i,j)+surtx(i-1,j))*0.5*g
     &                    /(pstres*thref)
            else  ! no taux
              stresx(i,j)=0.
            endif !windf:else
          enddo !i
        enddo !l
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbrhs(i,j)=
     &      -(ubavg(i,  j  ,m)*depthu(i,j  )
     &       +ubavg(i+1,j  ,m)*depthu(i+1,j  )
     &       +ubavg(i,  j-1,m)*depthu(i,j-1)
     &       +ubavg(i+1,j-1,m)*depthu(i+1,j-1))
     &       *(pvtrop(i,j)+pvtrop(i+1,j))*.125
c    
            if     (windf) then
              if(hybrid .and. mxlkrt) then
                pstres=0.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
              else
                pstres=dpv(i,j,1,m)
              endif
c ---         units of surty are N/m^2 (i.e. Pa)
              stresy(i,j)=(surty(i,j)+surty(i,j-1))*0.5*g
     &                    /(pstres*thref)
            else  ! no tauy
              stresy(i,j)=0.
            endif !windf:else
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'uba.n  n=',n
        call pipe_compare(ubavg(1-nbdy,1-nbdy,n),iu,text)
        write (text,'(a9,i3)') 'vba.n  n=',n
        call pipe_compare(vbavg(1-nbdy,1-nbdy,n),iv,text)
        write (text,'(a9,i3)') 'drag   k=',0
        call pipe_compare(drag,ip,text)
      endif
c
c --- the old  momeq2.f  starts here
c
c --- rhs: 0.0
c --- lhs: util1, util2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
c ---   spatial weighting function for pressure gradient calculation:
          util1(i,j)=0.0
          util2(i,j)=0.0
        enddo !i
      enddo !j
c
      do 9 k=1,kk
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. store minimum thickness
c --- values for use in pot.vort. calculation in -dpmx-.
c
*         wtime2( 1,k) = wtime()
c
c --- rhs: dpmx, dp.m+
c --- lhs: dpmx
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
        dpmx(i,1)=2.*cutoff
      enddo !i
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          dpmx(i,j+1)=2.*cutoff
        enddo !i
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpmx(i,j+1)=max(dpmx(i,j+1),dp(i,j,k,m)+dp(i-1,j,k,m))
          enddo !i
        enddo !l
      enddo !j
c       
*         wtime2( 2,k) = wtime()
c
c --- rhs: ubavg.m, ubavg.n, dp.m+, dpu
c --- lhs: utotm, utotn, uflux, dpmx, pu
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
c       
          i=ifu(j,l)-1
          if (i.ge.1-margin) then
            if     (iuopn(i,j).ne.0) then
              utotm(i,j)=u(i+1,j,k,m)+ubavg(i,j,m)
              utotn(i,j)=u(i+1,j,k,n)+ubavg(i,j,n)
              uflux(i,j)=utotm(i,j)*max(dp(i,j,k,m),cutoff)
            endif
*           if     (iuopn(i,j).eq.2) then
*             utotm(i,j)=u(i+1,j,k,m)+ubavg(i,j,m)
*             utotn(i,j)=u(i+1,j,k,n)+ubavg(i,j,n)
*             uflux(i,j)=utotm(i,j)*max(dp(i,j,k,m),cutoff)
*           elseif (iuopn(i,j).eq.1) then
*             utotm(i,j)=unest(i,j,k,1)+ubavg(i,j,m)
*             utotn(i,j)=unest(i,j,k,1)+ubavg(i,j,n)
*             uflux(i,j)=utotm(i,j)*max(dp(i,j,k,m),cutoff)
*           endif
          endif
          i=ilu(j,l)+1
          if (i.le.ii+margin) then
            if     (iuopn(i,j).ne.0) then
              utotm(i,j)=u(i-1,j,k,m)+ubavg(i,j,m)
              utotn(i,j)=u(i-1,j,k,n)+ubavg(i,j,n)
              uflux(i,j)=utotm(i,j)*max(dp(i-1,j,k,m),cutoff)
            endif
*           if     (iuopn(i,j).eq.2) then
*             utotm(i,j)=u(i-1,j,k,m)+ubavg(i,j,m)
*             utotn(i,j)=u(i-1,j,k,n)+ubavg(i,j,n)
*             uflux(i,j)=utotm(i,j)*max(dp(i-1,j,k,m),cutoff)
*           elseif (iuopn(i,j).eq.1) then
*             utotm(i,j)=unest(i,j,k,1)+ubavg(i,j,m)
*             utotn(i,j)=unest(i,j,k,1)+ubavg(i,j,n)
*             uflux(i,j)=utotm(i,j)*max(dp(i-1,j,k,m),cutoff)
*           endif
          endif
c           
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpmx(i,j)=max(dpmx(i,j),dp(i,j,k,m)+dp(i-1,j,k,m))
            utotm(i,j)=u(i,j,k,m)+ubavg(i,j,m)
            utotn(i,j)=u(i,j,k,n)+ubavg(i,j,n)
            uflux(i,j)=utotm(i,j)*max(dpu(i,j,k,m),cutoff)
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,m)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c     
*         wtime2( 3,k) = wtime()
c
c --- rhs: vbavg.m, vbavg.n, dp.m+, dpv
c --- lhs: vtotm, vtotn, vflux, dpmx, pv
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
        do l=1,jsv(i)
          j=jfv(i,l)-1
          if (j.ge.1-margin) then
            if      (ivopn(i,j).ne.0) then
              vtotm(i,j)=v(i,j+1,k,m)+vbavg(i,j,m)
              vtotn(i,j)=v(i,j+1,k,n)+vbavg(i,j,n)
              vflux(i,j)=vtotm(i,j)*max(dp(i,j,k,m),cutoff)
            endif
*           if      (ivopn(i,j).eq.2) then
*             vtotm(i,j)=v(i,j+1,k,m)+vbavg(i,j,m)
*             vtotn(i,j)=v(i,j+1,k,n)+vbavg(i,j,n)
*             vflux(i,j)=vtotm(i,j)*max(dp(i,j,k,m),cutoff)
*           elseif  (ivopn(i,j).eq.1) then
*             vtotm(i,j)=vnest(i,j,k,1)+vbavg(i,j,m)
*             vtotn(i,j)=vnest(i,j,k,1)+vbavg(i,j,n)
*             vflux(i,j)=vtotm(i,j)*max(dp(i,j,k,m),cutoff)
*           endif
          endif
          j=jlv(i,l)+1
          if (j.le.jj+margin) then
            if      (ivopn(i,j).ne.0) then
              vtotm(i,j)=v(i,j-1,k,m)+vbavg(i,j,m)
              vtotn(i,j)=v(i,j-1,k,n)+vbavg(i,j,n)
              vflux(i,j)=vtotm(i,j)*max(dp(i,j-1,k,m),cutoff)
            endif
*           if      (ivopn(i,j).eq.2) then
*             vtotm(i,j)=v(i,j-1,k,m)+vbavg(i,j,m)
*             vtotn(i,j)=v(i,j-1,k,n)+vbavg(i,j,n)
*             vflux(i,j)=vtotm(i,j)*max(dp(i,j-1,k,m),cutoff)
*           elseif  (ivopn(i,j).eq.1) then
*             vtotm(i,j)=vnest(i,j,k,1)+vbavg(i,j,m)
*             vtotn(i,j)=vnest(i,j,k,1)+vbavg(i,j,n)
*             vflux(i,j)=vtotm(i,j)*max(dp(i,j,k,m),cutoff)
*           endif
          endif
        enddo !l
      enddo !i
c       
*         wtime2( 4,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            dpmx(i  ,j)=max(dpmx(i  ,j),dp(i,j,k,m)+dp(i,j-1,k,m))
            dpmx(i+1,j)=max(dpmx(i+1,j),dp(i,j,k,m)+dp(i,j-1,k,m))
            vtotm(i,j)=v(i,j,k,m)+vbavg(i,j,m)
            vtotn(i,j)=v(i,j,k,n)+vbavg(i,j,n)
            vflux(i,j)=vtotm(i,j)*max(dpv(i,j,k,m),cutoff)
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,m)
          enddo !i
        enddo !l
      enddo !j
c
c --- rhs: corio+, dp.m+, dpmx+
c --- lhs: potvor
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c ---   assume margin<nblk
        ja=j-1
        jb=j+1
c
c --- pot.vort., at lateral boundary points
        do l=1,isv(j)
          i=ifv(j,l)
          if     (i.ge. 1-margin) then
c
c input de vorticite aux limites ouvertes - a tester - frontiere ouest
c
            potvor(i  ,j)=corio(i  ,j) * 8.0
     &                     /max(8.0*cutoff,
     &                          4.0*(dp(i,j,k,m)+dp(i,ja ,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i+1,j))
          endif
          i=ilv(j,l)
          if     (i.le.ii+margin) then
c
c input de vorticite aux limites ouvertes - a tester - frontiere est
c
            potvor(i+1,j)=corio(i+1,j) * 8.0
     &                     /max(8.0*cutoff,
     &                          4.0*(dp(i,j,k,m)+dp(i,ja ,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i+1,j))
          endif
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
c
c --- pot.vort., at lateral boundary points
c
c --- rhs: corio+, dp.m+, dpmx+
c --- lhs: potvor
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
c ---   assume margin<nblk
        ia=i-1
        do l=1,jsu(i)
          j=jfu(i,l)
c
c input de vorticite aux limites ouvertes - a tester - frontiere sud
c
          if     (j.ge. 1-margin) then
            potvor(i,j  )=corio(i,j  ) * 8.0
     &                     /max(8.0*cutoff,
     &                          4.0*(dp(i,j,k,m)+dp(ia ,j,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i,j+1))
          endif
          j=jlu(i,l)
c
c input de vorticite aux limites ouvertes - a tester - frontiere nord
c
          if     (j.le.jj+margin) then
            potvor(i,j+1)=corio(i,j+1) * 8.0
     &                     /max(8.0*cutoff,
     &                          4.0*(dp(i,j,k,m)+dp(ia ,j,k,m)),
     &                          dpmx(i,j),
     &                          dpmx(i,j+1))
          endif
        enddo !l
      enddo !i
c
c --- pot.vort., at interior points (incl. promontories)
c
c --- rhs: corio, dp.m+, dpmx+
c --- lhs: potvor
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isq(j)
          do i=max(1-margin,ifq(j,l)),min(ii+margin,ilq(j,l))
            potvor(i,j)=corio(i,j) * 8.0
     &         /max(8.0*cutoff,
     &              2.0*(dp(i,j  ,k,m)+dp(i-1,j  ,k,m)+
     &                   dp(i,j-1,k,m)+dp(i-1,j-1,k,m) ),
     &              dpmx(i,j),dpmx(i-1,j),dpmx(i+1,j),
     &                        dpmx(i,j-1),dpmx(i,j+1)    )
          enddo !i
        enddo !l
      enddo !j
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        do j=1,jj
          do i=1,ii
            mask(i,j)=min(1,iq(i,j)+iu(i,j  )+iv(i,j  )
     &                             +iu(i,j-1)+iv(i-1,j))
          enddo
        enddo
        write (text,'(a9,i3)') 'potvor k=',k
        call pipe_compare(potvor,mask,text)
      endif
c
c ----------------------------------------
c --- initialize fields for QUICK SCHEME
c ----------------------------------------
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          uwkm( i,j)=0.0
          uwkn( i,j)=0.0
          advu( i,j)=0.0
          diffu(i,j)=0.0
          vwkm( i,j)=0.0
          vwkn( i,j)=0.0
          advv( i,j)=0.0
          diffv(i,j)=0.0
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            uwkm(i,j)=utotm(i,j)
            uwkn(i,j)=utotn(i,j)
          enddo
          i=ifu(j,l)-1
          if (i.ge.1-margin.and.iuopn(i,j).ne.0) then
              uwkm(i,j)=utotm(i,j)
              uwkn(i,j)=utotn(i,j)
          endif
          i=ilu(j,l)+1
          if (i.le.ii+margin .and. iuopn(i,j).ne.0) then
              uwkm(i,j)=utotm(i,j)
              uwkn(i,j)=utotn(i,j)
          endif
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vwkm(i,j)=vtotm(i,j)
            vwkn(i,j)=vtotn(i,j)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
      do i=1-margin,ii+margin
        do l=1,jsv(i)
          j=jfv(i,l)-1
          if (j.ge.1-margin  .and. ivopn(i,j).ne.0) then
            vwkm(i,j)=vtotm(i,j)
            vwkn(i,j)=vtotn(i,j)
          endif
          j=jlv(i,l)+1
          if (j.le.jj+margin .and. ivopn(i,j).ne.0) then
            vwkm(i,j)=vtotm(i,j)
            vwkn(i,j)=vtotn(i,j)
          endif
        enddo !l
      enddo !i
c
c ------------------------------------------------------
c --- deformation-dependent eddy viscosity coefficient
c --- defor. at interior points (incl. promontories)
c --- defor takes into account CLOSE BOUND COND
c ------------------------------------------------------
c
c
      margin = 3
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 809 j=1-margin,jj+margin
        do 809 l=1,isp(j)
          do 809 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            defor1(i,j)=((uwkn(i+1,j)-uwkn(i,j))
     &                  -(vwkn(i,j+1)-vwkn(i,j)))**2
 809  continue
!$OMP END PARALLEL DO
c
      margin = 3
c
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 810 j=1-margin,jj+margin
        do 810 l=1,isq(j)
          do 810 i=max(1-margin,ifq(j,l)),min(ii+margin,ilq(j,l))
c
            defor2(i,j)=(vwkn(i-1,j)-vwkn(i,j)
     &                  +uwkn(i,j-1)-uwkn(i,j))**2
 810  continue
!$OMP END PARALLEL DO
c
c --- ----------
c --- u equation
c --- ----------
c
c --------------------------------------------------------
c --- START QUICK ADVECTION SCHEME HERE
c --------------------------------------------------------
c
c -------------------------------------
c First compute advection along X
c -------------------------------------
c
c
      margin = 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,dum,dun,
!$OMP&               velmodul,defortot,clpm,crmm,
!$OMP&               clpn,clmn,crpn,crmn,im1,im2,ip1,ip2)
!$OMP&         SCHEDULE(STATIC,jblk)
c
      do 16 j=1-margin,jj+margin
        do 16 l=1,isu(j)
c
c ---     modify uwk for boundary conditions : begin
          uwkn(ifu(j,l)-1,j)=0.0
          uwkn(ifu(j,l)-2,j)=-uwkn(ifu(j,l),j)
          uwkn(ilu(j,l)+1,j)=0.0
          uwkn(ilu(j,l)+2,j)=-uwkn(ilu(j,l),j)
c
          uwkm(ifu(j,l)-1,j)=0.0
          uwkm(ifu(j,l)-2,j)=-uwkm(ifu(j,l),j)
          uwkm(ilu(j,l)+1,j)=0.0
          uwkm(ilu(j,l)+2,j)=-uwkm(ilu(j,l),j)
c
c ---     Open Boundary Conditions
c
          i = ifu(j,l)-1
          if (iuopn(i,j).eq.1) then
             uwkn(ifu(j,l)-1,j)=uwkn(ifu(j,l),j)
             uwkn(ifu(j,l)-2,j)=uwkn(ifu(j,l),j)
             uwkm(ifu(j,l)-1,j)=uwkm(ifu(j,l),j)
             uwkm(ifu(j,l)-2,j)=uwkm(ifu(j,l),j)
          endif
c
          i = ilu(j,l)+1
          if (iuopn(i,j).eq.1) then
             uwkn(ilu(j,l)+1,j)=uwkn(ilu(j,l),j)
             uwkn(ilu(j,l)+2,j)=uwkn(ilu(j,l),j)
             uwkm(ilu(j,l)+1,j)=uwkm(ilu(j,l),j)
             uwkm(ilu(j,l)+2,j)=uwkm(ilu(j,l),j)
          endif
c ---     modify uwk for boundary conditions : end
c
           do 161 i=max(1-margin,ifu(j,l)-1),min(ii+margin,ilu(j,l)+2)
             dum(i)=uwkm(i,j)-uwkm(i-1,j)
             dun(i)=uwkn(i,j)-uwkn(i-1,j)
 161       continue
c
c ---     eddy viscosity coefficient
c
           do 163 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
             velmodul=sqrt( (uwkn(i,j)**2+uwkn(i-1,j)**2)/2.0 +
     &                      (iv(i-1,j+1)*vwkn(i-1,j+1)**2 +
     &                       iv(i-1,j)  *vwkn(i-1,j  )**2  )
     &                         /max(1,iv(i-1,j+1)+iv(i-1,j))   )
c
             defortot=sqrt( defor1(i-1,j)+
     &                      (defor2(i,j)    *iq(i,j)+
     &                       defor2(i-1,j)  *iq(i-1,j)+
     &                       defor2(i,j+1)  *iq(i,j+1)+
     &                       defor2(i-1,j+1)*iq(i-1,j+1) )
     &                  /max(1,iq(i,j)+iq(i-1,j)+iq(i,j+1)+iq(i-1,j+1))
     &                    )
c
             visc(i,j)=max(veldf4,velmodul*facdf4,visco4*defortot)
             visc(i,j)=min(scux(i,j)*coefdiv2,visc(i,j))
 163       continue
c
           visc(ilu(j,l)+1,j)=visc(ilu(j,l),j)
c
c
           do 164 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
             clpm=dum(i)-dum(i-1)
             crmm=dum(i+2)-dum(i+1)
             clpn=dun(i)-dun(i-1)
             clmn=dun(i+1)-dun(i)
             crpn=clmn
             crmn=dun(i+2)-dun(i+1)
c
c ---        calcul de U dxU et de la diffusion biLap
             advu(i,j)=uwkm(i,j)*(0.5*(uwkm(i+1,j)-uwkm(i-1,j)) -
     &                            coefdiv3*(crmm-clpm)            )*
     &                 scuxi(i,j)
c
             diffu(i,j)=
     &        (visc(i+1,j)*(crpn-crmn)*max(hmindiff,dp(i,  j,k,m)) -
     &         visc(i  ,j)*(clpn-clmn)*max(hmindiff,dp(i-1,j,k,m))  )
     &           *scuxi(i,j)/max(hmindiff,dpu(i,j,k,m))
 164       continue
c
c ---     restitution de uwk pour conditions aux limites: begin
          im1=ifu(j,l)-1
          im2=ifu(j,l)-2
          ip1=ilu(j,l)+1
          ip2=ilu(j,l)+2

          uwkn(im1,j)=utotn(im1,j)
          uwkn(im2,j)=utotn(im2,j)
          uwkn(ip1,j)=utotn(ip1,j)
          uwkn(ip2,j)=utotn(ip2,j)

          uwkm(im1,j)=utotm(im1,j)
          uwkm(im2,j)=utotm(im2,j)
          uwkm(ip1,j)=utotm(ip1,j)
          uwkm(ip2,j)=utotm(ip2,j)
c ---     restitution de uwk pour conditions aux limites: end
   16 continue
!$OMP END PARALLEL DO
c
c
c -------------------------------------
c Now compute advection along Y
c -------------------------------------
c
      margin = 4
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 180 j=1-margin,jj+margin
        do 180 l=1,isu(j)
          do 180 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
            if(iu(i,j-1).ne.0) then
              tdum(i,j)=uwkm(i,j)-uwkm(i,j-1)
              tdun(i,j)=uwkn(i,j)-uwkn(i,j-1)
            endif
  180 continue
!$OMP END PARALLEL DO
c
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,l,i,ia,velmodul,defortot)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 183 j=1-margin,jj+margin
        ja=j-1
        do 183 l=1,isu(j)
          do 183 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
          ia=i-1
          velmodul=sqrt( (vwkn(i,j)**2+vwkn(i-1,j)**2)/2.0 +
     &                   (iu(i,j) *uwkn(i,j)  **2 +
     &                    iu(i,ja)*uwkn(i,j-1)**2  )
     &                     /max(1,iu(i,j)+iu(i,ja))          )
c
          defortot=sqrt( defor2(i,j)+
     &                   (defor1(i,j)  *ip(i,j)  +
     &                    defor1(ia,j) *ip(ia,j) +
     &                    defor1(i,ja) *ip(i,ja) +
     &                    defor1(ia,ja)*ip(ia,ja) )
     &                     /max(1,ip(i,j)+ip(ia,j)+ip(i,ja)+ip(ia,ja))
     &                 )
c
          visc(i,j)=max(veldf4,velmodul*facdf4,visco4*defortot)
          visc(i,j)=min(scuy(i,j)*coefdiv2,visc(i,j))
  183 continue
!$OMP END PARALLEL DO
c
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 1841 j=1-margin,jj+margin
        do 1841 l=1,isu(j)
          do 1841 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            if((iu(i,j-2).ne.0).and.(iu(i,j-1).ne.0)) then
              tclpm(i,j)=tdum(i,j)-tdum(i,j-1)
              tclpn(i,j)=tdun(i,j)-tdun(i,j-1)
            endif
 1841 continue
!$OMP END PARALLEL DO
c
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
c
      do 1842 j=1-margin,jj+margin
        do 1842 l=1,isu(j)
          do 1842 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            if((iu(i,j-1).ne.0).and.(iu(i,j+1).ne.0)) then
              tclmm(i,j)=tdum(i,j+1)-tdum(i,j)
              tclmn(i,j)=tdun(i,j+1)-tdun(i,j)
            endif
 1842 continue
!$OMP END PARALLEL DO
c
      margin = 1
c
!$OMP PARALLEL DO PRIVATE(j,jm1,jp1,l,i,duadvjp1,
!$OMP&                    clpm,crmm,clpn,clmn,crpn,crmn)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 185 j=1-margin,jj+margin
        jm1=j-1
        do 185 l=1,isu(j)
          do 185 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            duadvjp1=0.5*(uwkm(i,j+1)-uwkm(i,j-1))
            clpm=tclpm(i,j)
            crmm=tclmm(i,j+1)
            clpn=tclpn(i,j)
            clmn=tclmn(i,j)
            crpn=tclpn(i,j+1)
            crmn=tclmn(i,j+1)
            if     ((j.eq.1).or.(iu(i,j-1).eq.0)) then
              duadvjp1=uwkm(i,j+1)-uwkm(i,j) 
              clpm=0.0
              clpn=0.0
              clmn=0.0
              crpn=0.0
            elseif (iu(i,j-2).eq.0) then
              clpm=0.0
              clpn=0.0
            endif
c
            if     (iu(i,j+1).eq.0) then
              duadvjp1=uwkm(i,j)-uwkm(i,j-1)
              crmm=0.0
              clmn=0.0
              crpn=0.0
              crmn=0.0
            elseif(iu(i,j+2).eq.0) then
              crmm=0.0
              crmn=0.0
            endif
c
            if (iu(i,j+1).eq.0 .and. iu(i,j-1).eq.0) then
              duadvjp1=0.0
            endif
c 
            if     (ivopn(i,  j-1).ne.0 .or.
     &              ivopn(i-1,j-1).ne.0 .or. 
     &              ivopn(i,  j+2).ne.0 .or.
     &              ivopn(i-1,j+2).ne.0     ) then
              advu(i,j)=advu(i,j)+
     &         0.25*(vwkm(i-1,j)+vwkm(i-1,j+1)+vwkm(i,j)+vwkm(i,j+1))
     &             *scuxi(i,j)* duadvjp1
            elseif (ivopn(i,  j)  .eq.0 .and.
     &              ivopn(i-1,j)  .eq.0 .and.
     &              ivopn(i,  j+1).eq.0 .and.
     &              ivopn(i-1,j+1).eq.0      ) then
              advu(i,j)=advu(i,j)+
     &         0.25*(vwkm(i-1,j)+vwkm(i-1,j+1)+vwkm(i,j)+vwkm(i,j+1))
     &             *scuxi(i,j)*(duadvjp1-coefdiv3*(crmm-clpm))
            endif
c
            jp1=j+1
            if(iu(i,j+1).eq.0) then
              jp1=j
            else
              jp1=j+1
            endif
            diffu(i,j)=diffu(i,j)+
     &       ( visc(i,jp1)*(crpn-crmn)*
     &         max(hmindiff,0.5*(dpu(i,jp1,k,m)+dpu(i,j,k,m)))
     &        -visc(i,j)*(clpn-clmn)*
     &         max(hmindiff,0.5*(dpu(i,j,k,m)+dpu(i,jm1,k,m)))
     &       )*scuxi(i,j)/max(hmindiff,dpu(i,j,k,m))
c
  185 continue
!$OMP END PARALLEL DO
c
c ------------------------------------------------------
c --- ----------
c --- v equation
c --- ----------
c ------------------------------------------------------
c
c --------------------------------------------------------
c --- START QUICK ADVECTION SCHEME HERE
c --------------------------------------------------------
c
c -------------------------------------
c First compute advection along X
c -------------------------------------
c
c
      margin = 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,dvm,dvn,
!$OMP&             ia,ja,velmodul,defortot,clpm,crmm,
!$OMP&             clpn,clmn,crpn,crmn,im1,im2,ip1,ip2)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 17 j=1-margin,jj+margin
        do 17 l=1,isv(j)
c
c ---     modification de vwk pour conditions aux limites: begin
          if     (ifv(j,l)-2 .gt. 1-nbdy) then
            vwkn(ifv(j,l)-1,j)=2.0*vwkn(ifv(j,l),  j)-vwkn(ifv(j,l)+1,j)
            vwkn(ifv(j,l)-2,j)=2.0*vwkn(ifv(j,l)-1,j)-vwkn(ifv(j,l),  j)
            vwkm(ifv(j,l)-1,j)=2.0*vwkm(ifv(j,l),  j)-vwkm(ifv(j,l)+1,j)
            vwkm(ifv(j,l)-2,j)=2.0*vwkm(ifv(j,l)-1,j)-vwkm(ifv(j,l),  j)
c
c ---       Open Boundary Conditions
c
            if     ( iuopn(ifv(j,l),j  ).ne.0 .or.
     &               iuopn(ifv(j,l),j-1).ne.0     ) then
              vwkn(ifv(j,l)-1,j)=vwkn(ifv(j,l),j)
              vwkn(ifv(j,l)-2,j)=vwkn(ifv(j,l)-1,j)
              vwkm(ifv(j,l)-1,j)=vwkm(ifv(j,l),j)
              vwkm(ifv(j,l)-2,j)=vwkm(ifv(j,l)-1,j)
            endif 
          endif 
c
          if(ilv(j,l)+2 .lt. idm+nbdy) then
            vwkn(ilv(j,l)+1,j)=2.0*vwkn(ilv(j,l),  j)-vwkn(ilv(j,l)-1,j)
            vwkn(ilv(j,l)+2,j)=2.0*vwkn(ilv(j,l)+1,j)-vwkn(ilv(j,l),  j)
            vwkm(ilv(j,l)+1,j)=2.0*vwkm(ilv(j,l),  j)-vwkm(ilv(j,l)-1,j)
            vwkm(ilv(j,l)+2,j)=2.0*vwkm(ilv(j,l)+1,j)-vwkm(ilv(j,l),  j)
c
c ---       Open Boundary Conditions
c
            if     (iuopn(ilv(j,l)+1,j  ).ne.0 .or.
     &              iuopn(ilv(j,l)+1,j-1).ne.0     ) then
              vwkn(ilv(j,l)+1,j)=vwkn(ilv(j,l),  j)
              vwkn(ilv(j,l)+2,j)=vwkn(ilv(j,l)+1,j)
              vwkm(ilv(j,l)+1,j)=vwkm(ilv(j,l),  j)
              vwkm(ilv(j,l)+2,j)=vwkm(ilv(j,l)+1,j)
            endif 
          endif 
c ---     modification de vwk pour conditions aux limites: end
c
          do 171 i=max(1-margin,ifv(j,l)-1),min(ii+margin,ilv(j,l)+2)
            dvm(i)=vwkm(i,j)-vwkm(i-1,j)
            dvn(i)=vwkn(i,j)-vwkn(i-1,j)
 171      continue
c
c ---     eddy viscosity coefficient
c
          do 173 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ja=j-1
            velmodul=sqrt( (vwkn(i,j)**2+vwkn(i-1,j)**2)/2.0 +
     &                     (iu(i,j) *uwkn(i,j)  **2 +
     &                      iu(i,ja)*uwkn(i,j-1)**2  )
     &                       /max(1,iu(i,j)+iu(i,ja))          )
c
            defortot=sqrt( defor2(i,j)+
     &                     (defor1(i,j)  *ip(i,j)  +
     &                      defor1(ia,j) *ip(ia,j) +
     &                      defor1(i,ja) *ip(i,ja) +
     &                      defor1(ia,ja)*ip(ia,ja) )
     &                       /max(1,ip(i,j)+ip(ia,j)+ip(i,ja)+ip(ia,ja))
     &                   )
c
            visc(i,j)=max(velmodul*facdf4,visco4*defortot)
            visc(i,j)=min(scvx(i,j)*coefdiv2,visc(i,j))
 173      continue
          visc(ilv(j,l)+1,j)=visc(ilv(j,l),j)
c
          do 174 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            im1=i-1
            ip1=i+1
c
            clpm=dvm(i)-dvm(i-1)
            crmm=dvm(i+2)-dvm(i+1)
c
            clpn=dvn(i)-dvn(i-1)
            clmn=dvn(i+1)-dvn(i)
            crpn=clmn
            crmn=dvn(i+2)-dvn(i+1)
c
c ---       calcul de U dxV et de la diffusion biLap
            advv(i,j)=
     &      0.25*(uwkm(i,j)+uwkm(i+1,j)+uwkm(i,j-1)+uwkm(i+1,j-1))*
     &                (0.5*(vwkm(i+1,j)-vwkm(i-1,j))-
     &                 coefdiv3*(crmm-clpm))*scvyi(i,j)
c 
            if(iuopn(i,  j).ne.0 .or. iuopn(i,  j-1).ne.0 .or.
     &         iuopn(i+1,j).ne.0 .or. iuopn(i+1,j-1).ne.0     ) then
              advv(i,j)=0.0
            endif
c
            if( iuopn(i-1,j).ne.0 .or. iuopn(i-1,j-1).ne.0 .or. 
     &          iuopn(i+2,j).ne.0 .or. iuopn(i+2,j-1).ne.0     ) then
              advv(i,j)=
     &        0.25*(uwkm(i,j)+uwkm(i+1,j)+uwkm(i,j-1)+uwkm(i+1,j-1))*
     &                  (0.5*(vwkm(i+1,j)-vwkm(i-1,j)))*scvyi(i,j)
            endif
c
            diffv(i,j)=
     &      ( visc(i+1,j)*(crpn-crmn)*
     &        max(hmindiff,0.5*(dpv(ip1,j,k,m)+dpv(i,  j,k,m))) -
     &        visc(i,  j)*(clpn-clmn)*
     &        max(hmindiff,0.5*(dpv(i,  j,k,m)+dpv(im1,j,k,m)))  )
     &         *scvyi(i,j)/max(hmindiff,dpv(i,j,k,m))
  174     continue
c
c ---     restitution de vwk pour conditions aux limites: begin
          im1=ifv(j,l)-1
          im2=ifv(j,l)-2
          ip1=ilv(j,l)+1
          ip2=ilv(j,l)+2

          vwkn(im1,j)=vtotn(im1,j)
          vwkn(im2,j)=vtotn(im2,j)
          vwkn(ip1,j)=vtotn(ip1,j)
          vwkn(ip2,j)=vtotn(ip2,j)

          vwkm(im1,j)=vtotm(im1,j)
          vwkm(im2,j)=vtotm(im2,j)
          vwkm(ip1,j)=vtotm(ip1,j)
          vwkm(ip2,j)=vtotm(ip2,j)
c ---     restitution de vwk pour conditions aux limites: end
   17 continue
!$OMP END PARALLEL DO
c
c -------------------------------------
c Now compute advection along Y
c -------------------------------------
c
      margin = 3
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 190 j=1-margin,jj+margin
        do 190 l=1,isv(j)
          do 190 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            tdvm(i,j)=vwkm(i,j)-vwkm(i,j-1)
            tdvn(i,j)=vwkn(i,j)-vwkn(i,j-1)
 190   continue
!$OMP END PARALLEL DO
c
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,velmodul,defortot)
!$OMP&         SCHEDULE(STATIC,jblk)
c
      do 193 j=1-margin,jj+margin
        do 193 l=1,isv(j)
          do 193 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            velmodul=sqrt( (vwkn(i,j)**2+vwkn(i,j-1)**2)/2.0 +
     &                     (iu(i,  j-1)*uwkn(i,  j-1)**2 +
     &                      iu(i+1,j-1)*uwkn(i+1,j-1)**2  )
     &                       /max(1,iu(i,j-1)+iu(i+1,j-1))    )
c
            defortot=sqrt( defor1(i,j-1)+
     &                     (defor2(i,  j)  *iq(i,  j)   +
     &                      defor2(i+1,j)  *iq(i+1,j)   +
     &                      defor2(i,  j-1)*iq(i,  j-1) +
     &                      defor2(i+1,j-1)*iq(i+1,j-1)  )
     &                /max(1,iq(i,j)+iq(i+1,j)+iq(i,j-1)+iq(i+1,j-1))
     &                   )
c
            visc(i,j)=max(velmodul*facdf4,visco4*defortot)
            visc(i,j)=min(scvy(i,j)*coefdiv2,visc(i,j))
  193 continue
!$OMP END PARALLEL DO
c
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 1941 j=1-margin,jj+margin
        do 1941 l=1,isv(j)
          do 1941 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            if(iv(i,j-1).ne.0) then
              tclpm(i,j)=tdvm(i,j)-tdvm(i,j-1)
              tclpn(i,j)=tdvn(i,j)-tdvn(i,j-1)
            endif
 1941 continue
c
      margin = 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 1942 j=1-margin,jj+margin
        do 1942 l=1,isv(j)
          do 1942 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            if(iv(i,j+1).ne.0) then
              tclmm(i,j)=tdvm(i,j+1)-tdvm(i,j)
              tclmn(i,j)=tdvn(i,j+1)-tdvn(i,j)
            endif
 1942 continue
c
      margin = 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,jp1,dvadvjp1,
!$OMP&     clpm,crmm,clpn,clmn,crpn,crmn)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 195 j=1-margin,jj+margin
        do 195 l=1,isv(j)
          do 195 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            dvadvjp1=0.5*(vwkm(i,j+1)-vwkm(i,j-1))
c
            clpm=tclpm(i,j)
            crmm=tclmm(i,j+1)
c
            clpn=tclpn(i,j)
            clmn=tclmn(i,j)
            crpn=tclpn(i,j+1)
            crmn=tclmn(i,j+1)
c
c ---       Closed Boundary Conditions
c
            if(iv(i,j-1).eq.0) then
              clpm=0.0
              clpn=0.0
            endif
c
            if     (iv(i,j+1).eq.0) then
              crmm=0.0
              crmn=0.0
              clmn=-vwkn(i,j)-tdvn(i,j)
              crpn=clmn
            elseif (iv(i,j+2).eq.0) then
              crmm=-vwkm(i,j+1)-tdvm(i,j+1)
              crmn=-vwkn(i,j+1)-tdvn(i,j+1)
            endif
c
c ---       Open Boundary Conditions
c
            if(ivopn(i,j-1).ne.0) then
              clpm=0.0
              clpn=0.0
            endif
c
            if(ivopn(i,j+1).ne.0) then
              crmm=0.0
              crmn=0.0
              clmn=-tdvn(i,j)
              crpn=clmn
            elseif(ivopn(i,j+2).ne.0) then
              crmm=-tdvm(i,j+1)
              crmn=-tdvn(i,j+1)
            endif
c
            advv(i,j)= advv(i,j)+
     &                  vwkm(i,j)*
     &                  (dvadvjp1-coefdiv3*(crmm-clpm))*
     &                  scvyi(i,j)
c
            jp1=j+1
            if(iv(i,j+1).eq.0) then
              jp1=j
            else
              jp1=j+1
            endif
            diffv(i,j)=diffv(i,j)+
     &      ( visc(i,jp1)*(crpn-crmn)*max(hmindiff,dp(i,j,k,m))
     &       -visc(i,j)*(clpn-clmn)*max(hmindiff,dp(i,j-1,k,m)) )*
     &      scvyi(i,j)/max(hmindiff,dpv(i,j,k,m))
c
  195 continue
!$OMP END PARALLEL DO
c
      margin = mbdy - 1
c
c --- ---------------------------------------
c --- ---------------------------------------
c --- ----------
c --- u equation
c --- ----------
c --- ---------------------------------------
c --- ---------------------------------------
c
c
c --- pressure force in x direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
c --- rhs: depthu, pu, montg+, thstar+, p+, dp.m+
c --- lhs: util1, pgfx
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,dmontg,dthstr,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            util1(i,j)=max(0.0,min(depthu(i,j)-pu(i,j,k),h1))
c
            if     (kapref.ne.-1) then
              dmontg=montg( i,j,k,1)-montg( i-1,j,k,1)
              dthstr=thstar(i,j,k,1)-thstar(i-1,j,k,1)
            else !2 thermobaric reference states
              q=0.5*(skap(i,j)+skap(i-1,j))
              dmontg=     q *(montg( i,j,k,1)-montg( i-1,j,k,1))+
     &               (1.0-q)*(montg( i,j,k,2)-montg( i-1,j,k,2))
              dthstr=     q *(thstar(i,j,k,1)-thstar(i-1,j,k,1))+
     &               (1.0-q)*(thstar(i,j,k,2)-thstar(i-1,j,k,2))
            endif
            pgfx(i,j)=util1(i,j)*
     &       ( advu(i,j)  -
     &         diffu(i,j) -
     &         0.125*(vflux(i  ,j)+vflux(i  ,j+1)+
     &                vflux(i-1,j)+vflux(i-1,j+1) )
     &              *(potvor(i,j)+potvor(i,j+1))    +
     &         scuxi(i,j)*(dmontg+
     &                     dthstr*(p(i,j,k+1)*p(i-1,j,k+1)-
     &                             p(i,j,k  )*p(i-1,j,k  ))*thref**2
     &                        /(dp(i,j,k,m)+dp(i-1,j,k,m)+epsil)    )
     &       )
          enddo
        enddo
      enddo
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'pu     k=',k
        call pipe_compare(pu(1-nbdy,1-nbdy,k),iu,text)
        write (text,'(a9,i3)') 'depthu k=',k
        call pipe_compare(depthu,iu,text)
        write (text,'(a9,i3)') 'util1  k=',k
        call pipe_compare(util1, iu,text)
        util4(1:ii,1:jj) = montg( 1:ii,1:jj,k,1)-montg( 0:ii-1,1:jj,k,1)
        write (text,'(a9,i3)') 'montgX k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = thstar(1:ii,1:jj,k,1)-thstar(0:ii-1,1:jj,k,1)
        write (text,'(a9,i3)') 'thstaX k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k)*p(0:ii-1,1:jj,k)
        write (text,'(a9,i3)') 'pSQ    k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k+1)*p(0:ii-1,1:jj,k+1)
        write (text,'(a9,i3)') 'pSQ    k=',k+1
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = dp(1:ii,1:jj,k,m)+dp(0:ii-1,1:jj,k,m)+epsil
        write (text,'(a9,i3)') 'dp.m+  k=',k
        call pipe_compare(util4, iu,text)
        write (text,'(a9,i3)') 'pgfx   k=',k
        call pipe_compare(pgfx,  iu,text)
      endif
c
c
c --- rhs: pgfx+, util1+, p+, dpmixl.m+, dpu.m, utotn, drag+
c --- rhs: stresx, u.m, u.n, dpu.n, gradx, utotm+, vtotm+,
c --- rhs: vflux+, potvor+, uflux1+, uflux2, uflux3
c --- lhs: gradx, stress, u.m, u.n
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,
!$OMP&                    ptopl,pbotl,pstres,pbop,qdpu,dragu)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
            gradx(i,j)=(pgfx(i,j)+(h1-util1(i,j))*
     &        (pgfx (i-1,j)+pgfx (i+1,j)+pgfx (i,ja)+pgfx (i,jb))/
     &        (util1(i-1,j)+util1(i+1,j)+util1(i,ja)+util1(i,jb)+
     &          epsil))/h1
c
            ptopl=min(depthu(i,j),0.5*(p(i,j,k  )+p(i-1,j,k  )))
            pbotl=min(depthu(i,j),0.5*(p(i,j,k+1)+p(i-1,j,k+1)))
            if(hybrid .and. mxlkrt) then
              pstres=0.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
            else
              pstres=dpu(i,j,1,m)
            endif
c             
c ---       top and bottom boundary layer stress
            pbop=depthu(i,j)-0.5*(thkbop(i,j)+thkbop(i-1,j)) !top of bot. b.l.
            qdpu=1.0/max(dpu(i,j,k,m),onemm)
            dragu=0.5*(drag(i,j)+drag(i-1,j))*
     &                ( max(pbop,    pbotl)
     &                 -max(pbop,min(ptopl,pbotl-onemm)))*qdpu
            if     (drglim.gt.0.0) then
c             explicit drag: u.t+1 - u.t-1 = - 2dt*dragu*u.t-1
c                   limiter: 2dt*dragu <= drglim
              stress(i,j)=-utotn(i,j)*min(dragu,drglim*dt1inv) +
     &             (stresx(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpu
            else
c             implicit drag: u.t+1 - u.t-1 = - 2dt*dragu*u.t+1
c                                          = - 2dt*dragu / (1+2dt*dragu) * u.t+1
              stress(i,j)=-utotn(i,j)*dragu/(1.0+delt1*dragu) +
     &             (stresx(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpu
            endif
c
            if     (tidflg.gt.0 .and. drgscl.ne.0.0) then
              pbop=depthu(i,j)-0.5*(util6(i,j)+util6(i-1,j)) !top of bot. b.l.
              qdpu=1.0/max(dpu(i,j,k,m),onemm)
              dragu=0.5*(util5(i,j)+util5(i-1,j))*
     &                  ( max(pbop,    pbotl)
     &                   -max(pbop,min(ptopl,pbotl-onemm)))*qdpu
              if     (drglim.gt.0.0) then
                stress(i,j)=stress(i,j)
     &                      -utotn(i,j)*min(dragu,drglim*dt1inv)
              else
                stress(i,j)=stress(i,j)
     &                      -utotn(i,j)*dragu/(1.0+delt1*dragu)
              endif
c ---         !anti-drag on detided velocity (always explicit)
              stress(i,j)=stress(i,j)
     &                      +untide(i,j)*min(dragu,adrlim*dt1inv)
            endif !tidflg
c
c ---       time smoothing of -u- field  (part 1)
            u(i,j,k,m)=u(i,j,k,m)*(wuv1*dpu(i,j,k,m)+onemm)
     &                +u(i,j,k,n)* wuv2*dpu(i,j,k,n)
c
            u(i,j,k,n)=u(i,j,k,n)+delt1*(-gradx(i,j)
     &                                   -ubrhs(i,j)
     &                                   +stress(i,j) )
          enddo !i
        enddo !l
      enddo !j
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'vtotm  k=',k
        call pipe_compare(vtotm, iv,text)
        write (text,'(a9,i3)') 'utotm  k=',k
        call pipe_compare(utotm, iu,text)
        write (text,'(a9,i3)') 'gradx  k=',k
        call pipe_compare(gradx, iu,text)
        write (text,'(a9,i3)') 'ubrhs  k=',k
        call pipe_compare(ubrhs, iu,text)
        write (text,'(a9,i3)') 'stress k=',k
        call pipe_compare(stress,iu,text)
        write (text,'(a9,i3)') 'u.n    k=',k
        call pipe_compare(u(1-nbdy,1-nbdy,k,n),iu,text)
      endif
c
 100    format(i9,8x,'uold    unew   gradp  nonlin   corio',
     &3x,'ubrhs  stress    fric')
c
c --- ----------
c --- v equation
c --- ----------
c
c --- pressure force in y direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
!$OMP PARALLEL DO PRIVATE(j,l,i,dmontg,dthstr,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            util2(i,j)=max(0.,min(depthv(i,j)-pv(i,j,k),h1))
            if     (kapref.ne.-1) then
              dmontg=montg( i,j,k,1)-montg( i,j-1,k,1)
              dthstr=thstar(i,j,k,1)-thstar(i,j-1,k,1)
            else !2 thermobaric reference states
              q=0.5*(skap(i,j)+skap(i,j-1))
              dmontg=     q *(montg( i,j,k,1)-montg( i,j-1,k,1))+
     &               (1.0-q)*(montg( i,j,k,2)-montg( i,j-1,k,2))
              dthstr=     q *(thstar(i,j,k,1)-thstar(i,j-1,k,1))+
     &               (1.0-q)*(thstar(i,j,k,2)-thstar(i,j-1,k,2))
            endif
            pgfy(i,j)=util2(i,j)*
     &       ( advv(i,j)  -
     &         diffv(i,j) +
     &         0.125*(uflux(i,j  )+uflux(i+1,j  )+
     &                uflux(i,j-1)+uflux(i+1,j-1) )
     &              *(potvor(i,j)+potvor(i+1,j))    +
     &         scvyi(i,j)*(dmontg+
     &                     dthstr*(p(i,j,k+1)*p(i,j-1,k+1)-
     &                             p(i,j,k  )*p(i,j-1,k  ))*thref**2
     &                      /(dp(i,j,k,m)+dp(i,j-1,k,m)+epsil)      )
     &       )
          enddo
        enddo
      enddo
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'pv     k=',k
        call pipe_compare(pv(1-nbdy,1-nbdy,k),iv,text)
        write (text,'(a9,i3)') 'depthv k=',k
        call pipe_compare(depthv,iv,text)
        write (text,'(a9,i3)') 'util2  k=',k
        call pipe_compare(util2, iv,text)
        util4(1:ii,1:jj) = montg( 1:ii,1:jj,k,1)-montg( 1:ii,0:jj-1,k,1)
        write (text,'(a9,i3)') 'montgY k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = thstar(1:ii,1:jj,k,1)-thstar(1:ii,0:jj-1,k,1)
        write (text,'(a9,i3)') 'thstaY k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k)  *p(1:ii,0:jj-1,k)
        write (text,'(a9,i3)') 'pSQ    k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k+1)*p(1:ii,0:jj-1,k+1)
        write (text,'(a9,i3)') 'pSQ    k=',k+1
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = dp(1:ii,1:jj,k,m)+dp(1:ii,0:jj-1,k,m)+epsil
        write (text,'(a9,i3)') 'dp.m+  k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'pgfy   k=',k
        call pipe_compare(pgfy,  iv,text)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
c
            grady(i,j)=(pgfy(i,j)+(h1-util2(i,j))*
     &        (pgfy (ia ,j)+pgfy (ib ,j)
     &        +pgfy (i,j-1)+pgfy (i,j+1))/
     &        (util2(ia ,j)+util2(ib ,j)
     &        +util2(i,j-1)+util2(i,j+1)+epsil))/h1
          enddo
        enddo
      enddo
c
!$OMP PARALLEL DO PRIVATE(j,l,i,
!$OMP&                    ptopl,pbotl,pstres,pbop,qdpv,dragv)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
            ptopl=min(depthv(i,j),0.5*(p(i,j,k  )+p(i,j-1,k  )))
            pbotl=min(depthv(i,j),0.5*(p(i,j,k+1)+p(i,j-1,k+1)))
            if(hybrid .and. mxlkrt) then
              pstres=0.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
            else
              pstres=dpv(i,j,1,m)
            endif
c
c ---       top and bottom boundary layer stress
            pbop=depthv(i,j)-0.5*(thkbop(i,j)+thkbop(i,j-1))  !top of bot. b.l.
            qdpv=1.0/max(dpv(i,j,k,m),onemm)
            dragv=0.5*(drag(i,j)+drag(i,j-1))*
     &                ( max(pbop,    pbotl)
     &                 -max(pbop,min(ptopl,pbotl-onemm)))*qdpv
            if     (drglim.gt.0.0) then
              stress(i,j)=-vtotn(i,j)*min(dragv,drglim*dt1inv) +
     &             (stresy(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpv
            else
              stress(i,j)=-vtotn(i,j)*dragv/(1.0+delt1*dragv) +
     &             (stresy(i,j)*thref*(min(pstres,pbotl+onemm)
     &                                -min(pstres,ptopl      )))*qdpv
            endif
c
            if     (tidflg.gt.0 .and. drgscl.ne.0.0) then
              pbop=depthv(i,j)-0.5*(util6(i,j)+util6(i,j-1))  !top of bot. b.l.
              qdpv=1.0/max(dpv(i,j,k,m),onemm)
              dragv=0.5*(util5(i,j)+util5(i,j-1))*
     &                  ( max(pbop,    pbotl)
     &                   -max(pbop,min(ptopl,pbotl-onemm)))*qdpv
              if     (drglim.gt.0.0) then
                stress(i,j)=stress(i,j)
     &                      -vtotn(i,j)*min(dragv,drglim*dt1inv)
              else
                stress(i,j)=stress(i,j)
     &                      -vtotn(i,j)*dragv/(1.0+delt1*dragv)
              endif
c ---         !anti-drag on detided velocity (always explicit)
              stress(i,j)=stress(i,j)
     &                      +vntide(i,j)*min(dragv,adrlim*dt1inv)
            endif !tidflg
c
c ---       time smoothing of -v- field  (part 1)
            v(i,j,k,m)=v(i,j,k,m)*(wuv1*dpv(i,j,k,m)+onemm)
     &                +v(i,j,k,n)* wuv2*dpv(i,j,k,n)
c
            v(i,j,k,n)=v(i,j,k,n)+delt1*(-grady(i,j)
     &                                   -vbrhs(i,j)
     &                                   +stress(i,j) )
          enddo
        enddo
      enddo
c
*     if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        util4(1:ii,1:jj) =  vtotm(1:ii,  2:jj+1)**2
     &                     +utotm(1:ii,  1:jj  )**2
     &                     +utotm(2:ii+1,1:jj  )**2
     &                     -vtotm(1:ii,  0:jj-1)**2
     &                     -utotm(1:ii,  0:jj-1)**2
     &                     -utotm(2:ii+1,0:jj-1)**2
        write (text,'(a9,i3)') 'totm   k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'grady  k=',k
        call pipe_compare(grady, iv,text)
        write (text,'(a9,i3)') 'vbrhs  k=',k
        call pipe_compare(vbrhs, iv,text)
        if     (k.eq.kk) then
          write (text,'(a9,i3)') 'p      k=',k-1
          call pipe_compare(p(1-nbdy,1-nbdy,k-1),ip,text)
          write (text,'(a9,i3)') 'p      k=',k
          call pipe_compare(p(1-nbdy,1-nbdy,k)  ,ip,text)
          write (text,'(a9,i3)') 'drag   k=',k
          call pipe_compare(drag,ip,text)
        endif
        write (text,'(a9,i3)') 'stress k=',k
        call pipe_compare(stress,iv,text)
        util4(1:ii,1:jj) =  uflux(1:ii,  1:jj  )
     &                     +uflux(2:ii+1,1:jj  )
     &                     +uflux(1:ii,  0:jj-1)
     &                     +uflux(2:ii+1,0:jj-1)
        write (text,'(a9,i3)') 'uflux  k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) =  potvor(1:ii,  1:jj  )
     &                     +potvor(2:ii+1,1:jj  )
        write (text,'(a9,i3)') 'potvor k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) =  vflux1(1:ii,  1:jj  )
     &                     -vflux1(1:ii,  0:jj-1)
     &                     +vflux3(1:ii,  1:jj  )
     &                     -vflux2(1:ii,  1:jj  )
        write (text,'(a9,i3)') 'vflux  k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'v.n    k=',k
        call pipe_compare(v(1-nbdy,1-nbdy,k,n),iv,text)
      endif
c
cdiag if (k.eq.1) then
cdiag write (lp,101) nstep
cdiag do j=jtest-1,jtest+1
cdiag do i=itest-1,itest+1
cdiag write (lp,'(2i5,i3,8f8.3)') i,j,k,
cdiag&  util4(i,j),v(i,j,k,n),-delt1*grady(i,j)*scvyi(i,j),
cdiag&  -delt1*
cdiag& 0.25*(vtotm(i,j+1)**2+utotm(i,j  )**2+utotm(i+1,j  )**2
cdiag&      -vtotm(i,j-1)**2-utotm(i,j-1)**2-utotm(i+1,j-1)**2)*
cdiag&       scvyi(i,j),
cdiag&  -delt1*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,j-1)+uflux(i+1,j-1))
cdiag&        *(potvor(i,j)+potvor(i+1,j))*0.125,
cdiag&  -delt1*vbrhs(i,j),delt1*stress(i,j),
cdiag&  -delt1*(vflux1(i,j)-vflux1(i,j-1)
cdiag&         +vflux3(i,j)-vflux2(i,j))/
cdiag&         (scv2(i,j)*max(dpv(i,j,k,m),onemm))
cdiag enddo
cdiag enddo
cdiag endif
 101    format(i9,8x,'vold    vnew   gradp  nonlin   corio',
     &3x,'vbrhs  stress    fric')
c
 9    continue  ! k=1,kk
c
C --- -----------------------------------------------------
C --- DEAL WITH VANISHING LAYERS AND TIME SMOOTHING
C --- -----------------------------------------------------
C
      margin = 1
c
*        wtime1( 6) = wtime()
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 14 j=1-margin,jj+margin
        do k=1,kk
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            enddo !i
          enddo !l
c           
c ---     compute new -dpu,dpv- field. save old -dpu,dpv- values in -pu,pv-.
c
          do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              pu(i,j,k+1)=dpu(i,j,k,n)
            enddo !i
          enddo !l
c           
          do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              pv(i,j,k+1)=dpv(i,j,k,n)
            enddo !i
          enddo !l
        enddo !k
 14   continue !j
!$OMP END PARALLEL DO
c
*        wtime1( 7) = wtime()
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n),
     &            dpv(1-nbdy,1-nbdy,1,n),
     &            p,depthu,depthv, max(0,margin-1))
*        wtime1( 8) = wtime()
c 
      margin = 0
c    
c --- substitute zero for (u,v) at bottomed massless grid points, and
c --- substitute depth-weighted averages at other massless grid points.
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,ka,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 30 j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ka=1
            do k=kk,2,-1
              if     (dpu(i,j,k,n).ge.dpthin) then
                ka = k  !deepest layer with significant mass
                exit
              endif
              u(i,j,k,n)=0.0
            enddo !k
            utotn(i,j)=u(i,j,1,n)*dpu(i,j,1,n)
            do k=2,ka
              q=min(dpu(i,j,k,m),dpu(i,j,k,n),cutoff)
              u(i,j,k,n)=(u(i,j,k,n)*q+u(i,j,k-1,n)*(cutoff-q))*qcutoff
              utotn(i,j)=utotn(i,j)+u(i,j,k,n)*dpu(i,j,k,n)
            enddo !k
            utotn(i,j)=utotn(i,j)/depthu(i,j)
          enddo !i
        enddo !l
c         
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ka=1
            do k=kk,2,-1
              if     (dpv(i,j,k,n).ge.dpthin) then
                ka = k  !deepest layer with significant mass
                exit
              endif
              v(i,j,k,n)=0.0
            enddo !k
            vtotn(i,j)=v(i,j,1,n)*dpv(i,j,1,n)
            do k=2,ka
              q=min(dpv(i,j,k,m),dpv(i,j,k,n),cutoff)
              v(i,j,k,n)=(v(i,j,k,n)*q+v(i,j,k-1,n)*(cutoff-q))*qcutoff
              vtotn(i,j)=vtotn(i,j)+v(i,j,k,n)*dpv(i,j,k,n)
            enddo !k
            vtotn(i,j)=vtotn(i,j)/depthv(i,j)
          enddo !i
        enddo !l
 30   continue !j
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        do k= 1,kk
          write (text,'(a9,i3)') 'dpv.m  k=',k
          call pipe_compare(dpv(1-nbdy,1-nbdy,k,m),iv,text)
          write (text,'(a9,i3)') 'dpv.n  k=',k
          call pipe_compare(dpv(1-nbdy,1-nbdy,k,n),iv,text)
          write (text,'(a9,i3)') 'v.n(2) k=',k
          call pipe_compare(  v(1-nbdy,1-nbdy,k,n),iv,text)
        enddo !k
        write (text,'(a9,i3)') 'vtotn  k=',kk
        call pipe_compare(vtotn, iv,text)
      endif
c
c --- time smoothing of -u,v- fields  (part 2)
c
*        wtime1( 9) = wtime()
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 22 j=1-margin,jj+margin
        do k=1,kk
          do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              u(i,j,k,n)=u(i,j,k,n)-utotn(i,j)
              u(i,j,k,m)=(u(i,j,k,m)+u(i,j,k,n)*wuv2*dpu(i,j,k,n))/
     &         (wuv1*dpu(i,j,k,m)+onemm+wuv2*(pu(i,j,k+1)+dpu(i,j,k,n)))
            enddo !i
          enddo !l
c           
          do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              v(i,j,k,n)=v(i,j,k,n)-vtotn(i,j)
              v(i,j,k,m)=(v(i,j,k,m)+v(i,j,k,n)*wuv2*dpv(i,j,k,n))/
     &         (wuv1*dpv(i,j,k,m)+onemm+wuv2*(pv(i,j,k+1)+dpv(i,j,k,n)))
            enddo !i
          enddo !l
        enddo !k
 22   continue !j
!$OMP END PARALLEL DO
c
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 867 j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            utotn(i,j)=utotn(i,j)*dt1inv
          enddo !i
        enddo !l
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vtotn(i,j)=vtotn(i,j)*dt1inv
          enddo !i
        enddo !l
 867  continue !j
!$OMP END PARALLEL DO
c
*        wtime1(10) = wtime()
*        if     (mod(nstep,100).eq.0) then
* ---      timer printout.
*          do i= 2,10
*            write(lp,'(4x,a,i3,a,f9.5,a)')
*    &        'momtum   section',i,' time =',
*    &        wtime1(i) - wtime1(i-1),' seconds'
*          enddo
*          do i= 2,20
*            wtimes=0.0
*            do k=1,kk
*              wtimes = wtimes + (wtime2(i,k)-wtime2(i-1,k))
*            enddo
*            write(lp,'(4x,a,i3,a,f9.5,a)')
*    &        'momtum k section',i,' time =',
*    &        wtimes,' seconds'
*          enddo
*          call flush(lp)
*        endif
c
      if     (lpipe .and. lpipe_momtum) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'v.n(3) k=',1
        call pipe_compare(v(1-nbdy,1-nbdy,1,n),iv,text)
      endif
c
      return
      end
c
c
c> Revision history
c>
c> Mar. 1995 - changed min.depth in pot.vort. calculation from 1 mm to 1 m
c>             (loops 812,802,803)
c> July 1997 - transferred -visc- and -vort- from common to local
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - added some loops that used to be in barotp.f
c> Aug. 1997 - transferred -wgtia,wgtib,wgtja,wgtjb- from common to local
c> Oct. 1999 - ekman depth calculation added in loops 66-70 - based on ustar
c>             which is calculated from the stress
c> Oct. 1999 - surface stress terms handled differently depending on m.l.
c>             model
c> Dec. 1999 - switched to biharmonic friction
c> Jan. 2000 - added biharm to select between laplacian and biharmonic
c> Jan. 2000 - ekman depth calculation removed - surface momentum flux
c>             distributed within layer one for all mixed layer models
c> Jan. 2000 - added r. bleck's re-formulation of -pgfx- and -pgfy-
c>             to properly commute pressure torque from layer to layer
c> Jul. 2000 - loop reordering and logic changes for OpenMP
c> Aug. 2000 - loop 117 executed only when hybrid vertical coordinate is used
c> Oct. 2001 - replaced biharm with veldf[24] and visco[24]
c> Sep. 2004 - kapref selects one of three thermobaric reference states
c> Jun. 2006 - split out momtum_hs, for initial ssh calculation
c> Nov. 2006 - inserted volume-force for tide
c> Apr. 2007 - added drglim, implicit or CFL-limited explicit bottom drag
c> Apr. 2007 - added dragrh, linear tidal drag based on bottom roughness
c> Apr. 2007 - btrlfr: moved [uvp]bavg assignment to barotp
c> Jun. 2007 - added momtum4
c> Mar  2009 - more accurate kappaf, with potential density
