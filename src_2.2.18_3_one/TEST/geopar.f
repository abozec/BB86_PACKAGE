      subroutine geopar
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
c
c --- set up model parameters related to geography
c
c --- hycom version 2.1
c
      implicit none
c
      include 'common_blocks.h'
c
      real      dp0kf,dpm,dpms,ds0kf,dsm,dsms
      real      hmina,hminb,hmaxa,hmaxb
      integer   i,ios,j,k,ktr,l
      character preambl(5)*79,cline*80
c
      real       aspmax
      parameter (aspmax=2.0)  ! maximum grid aspect ratio for diffusion
*     parameter (aspmax=1.0)  ! ignore  grid aspect ratio in  diffusion
c
c --- read grid location,spacing,coriolis arrays
c
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        write (lp,'(3a)') ' reading grid file from ',
     &                    flnmgrd(1:len_trim(flnmgrd)),'.[ab]'
        open (unit=9,file=flnmgrd(1:len_trim(flnmgrd))//'.b',
     &        status='old')
      endif
      call xcsync(flush_lp)
      call zagetc(cline,ios, 9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'I/O error from zagetc, iunit,ios = ',9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read(cline,*) i
c
      call zagetc(cline,ios, 9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'I/O error from zagetc, iunit,ios = ',9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read (cline,*) j
c
      if     (i.ne.itdm .or. j.ne.jtdm) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)')
     &    'error - wrong array size in grid file'
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      call zagetc(cline,ios, 9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'I/O error from zagetc, iunit,ios = ',9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      if     (mnproc.eq.1) then
      write (lp,'(a)') cline(1:len_trim(cline))
      endif
      read (cline,*) mapflg
c
      call zaiopf(flnmgrd(1:len_trim(flnmgrd))//'.a','old', 9)
c
      do k= 1,11
        call zagetc(cline,ios, 9)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)')
     &        'I/O error from zagetc, iunit,ios = ',9,ios
          endif !1st tile
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*) hminb,hmaxb
        if     (mnproc.eq.1) then
        write (lp,'(a)') cline(1:len_trim(cline))
        endif
        call xcsync(flush_lp)
c
        if     (k.eq.1) then
          call zaiord(plon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.2) then
          call zaiord(plat, ip,.false., hmina,hmaxa, 9)
          do i= 1,7
            call zagetc(cline,ios, 9)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(/ a,i4,i9 /)')
     &            'I/O error from zagetc, iunit,ios = ',9,ios
              endif !1st tile
              call xcstop('(geopar)')
                     stop '(geopar)'
            endif
            call zaiosk(9)
          enddo
        elseif (k.eq.3) then
          call zaiord(scpx, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.4) then
          call zaiord(scpy, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.5) then
          call zaiord(scqx, iq,.false., hmina,hmaxa, 9)
        elseif (k.eq.6) then
          call zaiord(scqy, iq,.false., hmina,hmaxa, 9)
        elseif (k.eq.7) then
          call zaiord(scux, iu,.false., hmina,hmaxa, 9)
        elseif (k.eq.8) then
          call zaiord(scuy, iu,.false., hmina,hmaxa, 9)
        elseif (k.eq.9) then
          call zaiord(scvx, iv,.false., hmina,hmaxa, 9)
        elseif (k.eq.10) then
          call zaiord(scvy, iv,.false., hmina,hmaxa, 9)
        else
          call zaiord(corio,iq,.false., hmina,hmaxa, 9)
        endif
c
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          endif
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
      enddo
c
      call zaiocl(9)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close(unit=9)
      endif
c
      if (itest.gt.0 .and. jtest.gt.0) then
        i=itest
        j=jtest
        write (lp,'(/ a,2i5,a,f8.3,a,f12.9,2f10.2/)')
     &   ' i,j=',i+i0,j+j0,
     &   ' plat=',plat(i,j),
     &   ' corio,scux,vy=',corio(i,j),scux(i,j),scvy(i,j)
      endif
      call xcsync(flush_lp)
c
c --- read basin depth array
c
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        write (lp,'(3a)') ' reading bathymetry file from ',
     &                    flnmdep(1:len_trim(flnmdep)),'.[ab]'
        open (unit=9,file=flnmdep(1:len_trim(flnmdep))//'.b',
     &        status='old')
        read (     9,'(a79)')  preambl
      endif
      call xcsync(flush_lp)
      call zagetc(cline,ios, 9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'I/O error from zagetc, iunit,ios = ',9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close(unit=9)
        write (lp,'(/(1x,a))') preambl,cline
      endif
c
      call zaiopf(flnmdep(1:len_trim(flnmdep))//'.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.5*huge) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
c --- determine do-loop limits for u,v,p,q points, and update halo for depths
      call bigrid(depths, mapflg, util1,util2,util3)
ccc      call prtmsk(ip,depths,util1,idm,ii,jj,0.0,1.0,
ccc     &     'bottom depth (m)')
c
c     now safe to apply halo to arrays.
c
      vland = 1.0
      call xctilr(plon,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(plat,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(corio, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scpx,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(scpy,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(scqx,  1,1, nbdy,nbdy, halo_qs)
      call xctilr(scqy,  1,1, nbdy,nbdy, halo_qs)
      call xctilr(scux,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scuy,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scvx,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvy,  1,1, nbdy,nbdy, halo_vs)
      vland = 0.0
c
c --- area of grid cells (length x width) at u,v,p,q points resp.
c
******!$OMP PARALLEL DO PRIVATE(j,i)
******!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          scu2(i,j)=scux(i,j)*scuy(i,j)
          scv2(i,j)=scvx(i,j)*scvy(i,j)
          scp2(i,j)=scpx(i,j)*scpy(i,j)
          scq2(i,j)=scqx(i,j)*scqy(i,j)
c
          scuxi(i,j)=1.0/max(scux(i,j),epsil)
          scvyi(i,j)=1.0/max(scvy(i,j),epsil)
          scp2i(i,j)=1.0/max(scp2(i,j),epsil)
          scq2i(i,j)=1.0/max(scq2(i,j),epsil)
c
c ---     largest grid spacing (within limits) used in all diffusion
c ---     coefficients: min(max(sc?x,sc?y),sc?x*aspmax,sc?y*aspmax)
          aspux(i,j)=min(max(scux(i,j),scuy(i,j)),
     &                   min(scux(i,j),scuy(i,j))*aspmax)
     &               /max(scux(i,j),epsil)
          aspuy(i,j)=min(max(scux(i,j),scuy(i,j)),
     &                   min(scux(i,j),scuy(i,j))*aspmax)
     &               /max(scuy(i,j),epsil)
          aspvx(i,j)=min(max(scvx(i,j),scvy(i,j)),
     &                   min(scvx(i,j),scvy(i,j))*aspmax)
     &               /max(scvx(i,j),epsil)
          aspvy(i,j)=min(max(scvx(i,j),scvy(i,j)),
     &                   min(scvx(i,j),scvy(i,j))*aspmax)
     &               /max(scvy(i,j),epsil)
c
          util1(i,j)=depths(i,j)*scp2(i,j)
        enddo
      enddo
c
      call xcsum(avgbot, util1,ip)
      call xcsum(area,   scp2, ip)
      avgbot=avgbot/area
      if     (mnproc.eq.1) then
      write (lp,'(/a,f9.1,-12p,f10.2)')
     &       ' mean basin depth (m) and area (10^6 km^2):',
     &       avgbot,area
      endif
      call xcsync(flush_lp)
c
c --- logorithmic k-dependence of dp0 (deep z's)
      dp00 =onem*dp00
      dp00x=onem*dp00x
      if     (isopyc) then
        dp0k(1)=thkmin*onem
      else
        dp0k(1)=dp00
      endif
      dp0kp(1)=dp0k(1)+onem
      dpm  = dp0k(1)*qonem
      dpms = dpm
      if     (mnproc.eq.1) then
      write(lp,*)
      write(lp,135) 1,dp0k(1)*qonem,dpm,dpms
      endif
 135  format('dp0k(',i2,') =',f7.2,' m',
     &          '    thkns =',f7.2,' m',
     &          '    depth =',f8.2,' m')
      call xcsync(flush_lp)
c
      dp0kf=1.0
      do k=2,kk
        dp0kf=dp0kf*dp00f
        if     (k.le.nhybrd) then
          dp0k(k)=min(dp00*dp0kf,dp00x)
        else
          dp0k(k)=0.0
        endif
        dp0kp(k)=dp0k(k)+onem
        dpm  = dp0k(k)*qonem
        dpms = dpms + dpm
        if     (mnproc.eq.1) then
        write(lp,135) k,dp0k(k)*qonem,dpm,dpms
        endif
        if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
          write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
          write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
        endif
        call xcsync(flush_lp)
      enddo
c
c --- logorithmic k-dependence of ds0 (shallow z-s)
      ds00 =onem*ds00
      ds00x=onem*ds00x
      if     (isopyc) then
        ds0k(1)=thkmin*onem
      else
        ds0k(1)=ds00
      endif
      dsm  = ds0k(1)*qonem
      dsms = dsm
      if     (mnproc.eq.1) then
      write(lp,*)
      write(lp,130) 1,ds0k(1)*qonem,dsm,dsms
      endif
 130  format('ds0k(',i2,') =',f7.2,' m',
     &          '    thkns =',f7.2,' m',
     &          '    depth =',f8.2,' m')
      call xcsync(flush_lp)
c
      ds0kf=1.0
      do k=2,nsigma
        ds0kf=ds0kf*ds00f
        ds0k(k)=min(ds00*ds0kf,ds00x)
        dsm  = ds0k(k)*qonem
        dsms = dsms + dsm
        if     (mnproc.eq.1) then
        write(lp,130) k,ds0k(k)*qonem,dsm,dsms
        endif
        if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
          write(6,*) 'geopar: ds0kf  = ',ds0kf,    mnproc
          write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
        endif
        call xcsync(flush_lp)
      enddo
      if     (mnproc.eq.1) then
      write(lp,*)
      endif
c
c --- sigma-depth scale factors
      do k=1,nsigma
        dssk(k)=ds0k(k)/dsms  ! onem * fraction of depths in sigma layer k
      enddo
      do k= nsigma+1,kdm
        ds0k(k)=dp0k(k)
        dssk(k)=0.0           ! these layers are zero in sigma mode
      enddo
c
c --- initialize some arrays
c --- set depthu,dpu,utotn,pgfx,depthv,dpv,vtotn,pgfy to zero everywhere,
c --- so that they can be used at "lateral neighbors" of u and v points.
c --- similarly for pbot,dp at neighbors of q points.
c
!$OMP PARALLEL DO PRIVATE(j,i,k,ktr)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(     i,j,1)=0.0
          pu(    i,j,1)=0.0
          pv(    i,j,1)=0.0
          utotn( i,j)=0.0
          vtotn( i,j)=0.0
          pgfx(  i,j)=0.0
          pgfy(  i,j)=0.0
          depthu(i,j)=0.0
          depthv(i,j)=0.0
          pbot(  i,j)=0.0
c
          ubavg( i,j,1)=huge
          ubavg( i,j,2)=huge
          ubavg( i,j,3)=huge
          vbavg( i,j,1)=huge
          vbavg( i,j,2)=huge
          vbavg( i,j,3)=huge
          utotm( i,j)=huge
          vtotm( i,j)=huge
          uflux( i,j)=huge
          vflux( i,j)=huge
          uflux1(i,j)=huge
          vflux1(i,j)=huge
          uflux2(i,j)=huge
          vflux2(i,j)=huge
          uflux3(i,j)=huge
          vflux3(i,j)=huge
          uja(   i,j)=huge
          ujb(   i,j)=huge
          via(   i,j)=huge
          vib(   i,j)=huge
          do k=1,kk
            dp( i,j,k,1)=0.0
            dp( i,j,k,2)=0.0
            dpu(i,j,k,1)=0.0
            dpu(i,j,k,2)=0.0
            dpv(i,j,k,1)=0.0
            dpv(i,j,k,2)=0.0
c
            u(  i,j,k,1)=huge
            u(  i,j,k,2)=huge
            v(  i,j,k,1)=huge
            v(  i,j,k,2)=huge
c
            uflx(  i,j,k)=huge
            vflx(  i,j,k)=huge
c
            dpav(  i,j,k)=0.0
            uflxav(i,j,k)=0.0
            vflxav(i,j,k)=0.0
            diaflx(i,j,k)=0.0
c
            do ktr= 1,ntracr
              tracer(i,j,k,1,ktr)=0.0
              tracer(i,j,k,2,ktr)=0.0
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
            ubavg(i,j,1)=0.0
            ubavg(i,j,2)=0.0
            ubavg(i,j,3)=0.0
            utotm (i,j)=0.0
            uflux (i,j)=0.0
            uflux2(i,j)=0.0
            uflux3(i,j)=0.0
            uja(i,j)=0.0
            ujb(i,j)=0.0
c
            do k=1,kk
              uflx(i,j,k)=0.0
              u(i,j,k,1)=0.0
              u(i,j,k,2)=0.0
            enddo
          enddo
        enddo
      enddo
c
      call xctilr(ubavg,    1,   3, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(utotm,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux2,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux3,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uja,      1,   1, nbdy,nbdy, halo_us)
      call xctilr(ujb,      1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflx,     1,  kk, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(u,        1,2*kk, nbdy,nbdy, halo_us)  ! note scalar
c
!$OMP PARALLEL DO PRIVATE(i,l,j,k)
!$OMP&         SCHEDULE(STATIC)
      do i=1,ii
        do l=1,jsp(i)
          do j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
            vbavg(i,j,1)=0.0
            vbavg(i,j,2)=0.0
            vbavg(i,j,3)=0.0
            vtotm (i,j)=0.0
            vflux (i,j)=0.0
            vflux2(i,j)=0.0
            vflux3(i,j)=0.0
            via(i,j)=0.0
            vib(i,j)=0.0
c
            do k=1,kk
              vflx(i,j,k)=0.0
              v(i,j,k,1)=0.0
              v(i,j,k,2)=0.0
            enddo
          enddo
        enddo
      enddo
c
      call xctilr(vbavg,    1,   3, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vtotm,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux2,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux3,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(via,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vib,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflx,     1,  kk, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(v,        1,2*kk, nbdy,nbdy, halo_vs)  ! note scalar
c
      return
      end
c
c
c> Revision history:
c>
c> May  1997 - extended list of variables set to 'huge' on land
c> Oct. 1999 - added code that defines the vertical distribution of dp0
c>             used in hybgen
c> Jan. 2000 - added mapflg logic for different projections
c> Feb. 2000 - added dp00f for logorithmic z-level spacing
c> Mar. 2000 - added dp00s for sigma-spacing in shallow water
c> May  2000 - conversion to SI units (still wrong corio)
c> Feb. 2001 - removed rotated grid option
c> Jan. 2002 - more flexible Z-sigma-Z vertical configuration
c> Jan. 2002 - all grids now via array input
