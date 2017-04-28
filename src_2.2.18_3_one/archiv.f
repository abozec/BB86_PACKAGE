      subroutine archiv(n, kkout, iyear,iday,ihour, intvl)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer   n, kkout, iyear,iday,ihour
      real      sssc,sstc
      character intvl*3
c
      include 'stmt_fns.h'
c
c --- write an archive file.
c
      character*80 cformat
      integer      i,j,k,ktr,l,ldot,nop,nopa
      real         coord,xmin,xmax
c
      ldot = index(flnmarc,'.',back=.true.)
      if     (ldot.eq.0) then
        if     (mnproc.eq.1) then
        write (lp,*) 'need decimal point in flnmarc'
        write (lp,*) 'flnmarc = ',trim(flnmarc)
        endif
        call xcstop('(flnmarc)')
               stop '(flnmarc)'
      endif
      ldot = min(ldot,len(flnmarc)-11)  !need 11 characters for archive date
c
      if     ((kkout.eq.1 .and. dsurfq.ge.1.0/24.0) .or.
     &        (kkout.gt.1 .and. diagfq.ge.1.0/24.0)     ) then
c ---   indicate the archive date
        write(flnmarc(ldot+1:ldot+11),'(i4.4,a1,i3.3,a1,i2.2)') 
     &   iyear,'_',iday,'_',ihour
        ldot=ldot+11
      else
c ---   indicate the archive time step
        write(flnmarc(ldot+1:ldot+11),'(i11.11)') nstep
        ldot=ldot+11
      endif
      nopa=13
      nop =13+uoff
c
c --- no .[ab] files for 1-D cases (<=6x6) or for dsur1p surface cases.
c
      if     (max(itdm,jtdm).gt.6 .and.
     &        .not.(dsur1p .and. kkout.eq.1)) then  !not 1-D output
c
      call zaiopf(flnmarc(1:ldot)//'.a', 'new', nopa)
      if     (mnproc.eq.1) then
      open (unit=nop,file=flnmarc(1:ldot)//'.b',status='new') !uoff+13
      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm
      call flush(nop)
      endif !1st tile
 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
c
c --- surface fields
c
      coord=0.
c
      call zaiowr(montg1,ip,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
c --- identify the equation of state on the first record
      write (nop,117) 'montg1  ',nstep,time,sigver,thbase,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(srfhgt,ip,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'srfhgt  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      if     (sshflg.ne.0) then
c ---   write out steric SSH.
        call zaiowr(steric,ip,.true.,
     &              xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'steric  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      endif !sshflg
c
      call zaiowr(surflx,ip,.true., xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'surflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(salflx,ip,.true., xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'salflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
      call zaiowr(dpbl,ip,.true., xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'bl_dpth ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(dpmixl(1-nbdy,1-nbdy,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'mix_dpth',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      if     (iceflg.ne.0) then
        call zaiowr(covice,ip,.true., xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'covice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
        call zaiowr(thkice,ip,.true., xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'thkice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
        call zaiowr(temice,ip,.true., xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'temice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      endif  !write ice fields
c
c --- depth averaged fields
c
      call zaiowr(ubavg(1-nbdy,1-nbdy,n),iu,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'u_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(vbavg(1-nbdy,1-nbdy,n),iv,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'v_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=sigma(k)
      call zaiowr(u(1-nbdy,1-nbdy,k,n),iu,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'u-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(v(1-nbdy,1-nbdy,k,n),iv,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'v-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(dp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'thknss  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(temp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'temp    ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(saln(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'salin   ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
c --- no tracers or diffusion for single layer case
c
      if     (kkout.gt.1) then
        do ktr= 1,ntracr
          call zaiowr(tracer(1-nbdy,1-nbdy,k,n,ktr),ip,.true.,
     &                xmin,xmax, nopa, .false.)
          if     (mnproc.eq.1) then
          write (nop,117) 'tracer  ',nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
        enddo !ktr
        if     (difout) then
          call zaiowr(vcty(1-nbdy,1-nbdy,k+1),ip,.true.,
     &                xmin,xmax, nopa, .false.)
          if     (mnproc.eq.1) then
          write (nop,117) 'viscty  ',nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
          call zaiowr(dift(1-nbdy,1-nbdy,k+1),ip,.true.,
     &                xmin,xmax, nopa, .false.)
          if     (mnproc.eq.1) then
          write (nop,117) 't-diff  ',nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
          call zaiowr(difs(1-nbdy,1-nbdy,k+1),ip,.true.,
     &                xmin,xmax, nopa, .false.)
          if     (mnproc.eq.1) then
          write (nop,117) 's-diff  ',nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
        endif !difout
      endif !kkout>1
 75   continue
c
 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
c
c --- output time-averaged mass fluxes, if required
c
      if (.not. (mxlkpp .or. mxlmy .or. mxlgiss) .and. kkout.eq.kk) then
        do k=1,kk
          coord=sigma(k)
          call zaiowr(diaflx(1-nbdy,1-nbdy,k),ip,.true.,
     &                xmin,xmax, nopa, .false.)
          if     (mnproc.eq.1) then
          write (nop,118) 'diafx',intvl,nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
 118      format (a5,a3,' =',i11,f11.3,i3,f7.3,1p2e16.7)
        enddo
      endif !diaflx
c
      close (unit=nop)
      call zaiocl(nopa)
c
      call xcsync(no_flush)
c
      endif  !not 1-D
c
      if     (itest.gt.0 .and. jtest.gt.0) then
        if     (relaxf .and. sstflg.le.1) then
          sstc = twall(itest,jtest,1,lc0)*wc0+
     &           twall(itest,jtest,1,lc1)*wc1+
     &           twall(itest,jtest,1,lc2)*wc2+
     &           twall(itest,jtest,1,lc3)*wc3
        else !synoptic observed sst
          sstc = seatmp(itest,jtest,l0)*w0+
     &           seatmp(itest,jtest,l1)*w1+
     &           seatmp(itest,jtest,l2)*w2+
     &           seatmp(itest,jtest,l3)*w3
        endif
        sssc = swall(itest,jtest,1,lc0)*wc0+
     &         swall(itest,jtest,1,lc1)*wc1+
     &         swall(itest,jtest,1,lc2)*wc2+
     &         swall(itest,jtest,1,lc3)*wc3
        open (unit=nop,file=flnmarc(1:ldot)//'.txt',status='new') !uoff+13
        write (nop,'(3a / a,6i7,2f8.1,i7,i7.4,i7.3,i7.2)')
     &      '##   expt    idm    jdm    kdm',
     &        '  itest  jtest  lontst  lattst',
     &        ' yrflag   year    day     hr',
     &      '##',iexpt,  itdm,  jtdm,   kdm,
     &          ittest,jttest,
     &          mod(plon(itest,jtest),360.0),plat(itest,jtest),
     &          yrflag, iyear,  iday, ihour
        write (nop,'(7a / a,f10.3, f8.2,4f8.1, 2f9.2,2f8.4,
     &                    f9.5,4f9.3, 2f8.3, 3f8.3, 4f8.2)')
     &    '## model-day',
     &    '  srfhgt  sswflx  mixflx  surflx  sstflx',
     &    '      E-P   sssE-P  bhtflx  buoflx',
     &    '    ustar   hekman    dpbbl     dpbl   dpmixl',
     &    '   tclim   sclim',
     &    '    tmix    smix   thmix    umix    vmix',
     &    '   ubavg   vbavg',
     &    '#',time,                                              !model-day
     &    srfhgt(itest,jtest)*100.0/g,                           !cm
     &    sswflx(itest,jtest),                                   !W/m**2
     &    mixflx(itest,jtest),                                   !W/m**2
     &    surflx(itest,jtest),                                   !W/m**2
     &    sstflx(itest,jtest),                                   !W/m**2
     &    salflx(itest,jtest)*thref*8.64E7/saln(itest,jtest,1,n),!mm/day
     &    sssflx(itest,jtest)*thref*8.64E7/saln(itest,jtest,1,n),!mm/day
     &    bhtflx(itest,jtest)*1.e6,                         !1.e6*m**2/sec**3
     &    buoflx(itest,jtest)*1.e6,                         !1.e6*m**2/sec**3
     &     ustar(itest,jtest),                                   !m/s?
     &    min(hekman(itest,jtest),         9999.999),            !m
     &    min( dpbbl(itest,jtest)  *qonem, 9999.999),            !m
     &    min(  dpbl(itest,jtest)  *qonem, 9999.999),            !m
     &    min(dpmixl(itest,jtest,n)*qonem, 9999.999),            !m
     &      sstc,                                                !degC
     &      sssc,                                                !psu
     &      tmix(itest,jtest),                                   !degC
     &      smix(itest,jtest),                                   !psu
     &     thmix(itest,jtest)+thbase,                            !SigmaT
     &    max(-999.99,min(999.99,
     &        (umix(itest,jtest)+ubavg(itest,jtest,n))*100.0)),  !cm/s
     &    max(-999.99,min(999.99,
     &        (vmix(itest,jtest)+vbavg(itest,jtest,n))*100.0)),  !cm/s
     &    max(-999.99,min(999.99,
     &        ubavg(itest,jtest,n)*100.0)),                      !cm/s
     &    max(-999.99,min(999.99,
     &        vbavg(itest,jtest,n)*100.0))                       !cm/s
        if     (iceflg.ne.0) then
          write (nop,'(2a / a,f10.3, 3f8.2,2f8.1,f9.2)')
     &    '## model-day',
     &    '  covice  thkice  temice  flxice  fswice   iceE-P',
     &    '#',time,                                                !model-day
     &      covice(itest,jtest)*100.0,                             !%
     &      thkice(itest,jtest),                                   !m
     &      temice(itest,jtest),                                   !degC
     &      flxice(itest,jtest),                                   !W/m**2
     &      fswice(itest,jtest),                                   !W/m**2
     &      sflice(itest,jtest)*thref*8.64E7/saln(itest,jtest,1,n) !mm/day
        endif !iceflg
        if     (ntracr.eq.0) then
          write(cformat,'(a)')
     &      '(3a / (i4,2f8.2,3f8.3,f9.3,f10.3,2f8.2))'
        else
          write(cformat,'(a,i2,a,i2,a)')
     &      '(3a,', ntracr,
     &      'a / (i4,2f8.2,3f8.3,f9.3,f10.3,2f8.2,', ntracr,
     &      'f8.3))'
        endif
        write (nop,cformat)
     &      '#  k',
     &      '    utot    vtot    temp    saln    dens',
     &      '    thkns      dpth  viscty  t-diff',
     &      ('  tracer',ktr=1,ntracr),
     &      (k,
     &       max(-999.99,min(999.99,
     &           (u(itest,jtest,k,n)+ubavg(itest,jtest,n))*100.0)),  !cm/s
     &       max(-999.99,min(999.99,
     &           (v(itest,jtest,k,n)+vbavg(itest,jtest,n))*100.0)),  !cm/s
     &       temp(itest,jtest,k,n),                                  !degC
     &       saln(itest,jtest,k,n),                                  !psu
     &       th3d(itest,jtest,k,n)+thbase,                           !SigmaT
     &         dp(itest,jtest,k,n)*qonem,                            !m
     &         (p(itest,jtest,k+1)+p(itest,jtest,k))*0.5*qonem,      !m
     &       vcty(itest,jtest,k+1)*1.e4,                             !m**2/s*2
     &       dift(itest,jtest,k+1)*1.e4,                             !m**2/s*2
     &       (tracer(itest,jtest,k,n,ktr),ktr=1,ntracr),             !0-999?
     &       k=1,kk)
        close (unit=nop)
      endif !test point tile
c
      call xcsync(no_flush)
cccc
cccc --- output to line printer
cccc
ccc      call prtmsk(ip,srfhgt,util3,idm,ii,jj,0.,100.0/g,
ccc     .     'sea surface height (cm)')
ccc      if(mxlkpp) call prtmsk(ip,dpbl,util3,idm,ii,jj,0.,1.*qonem,
ccc     .     'turb. b.l. depth (m)')
ccc      call prtmsk(ip,dpmixl,util3,idm,ii,jj,0.,1.*qonem,
ccc     .     'mixed layer depth (m)')
ccc      call prtmsk(ip,tmix,util3,idm,ii,jj,0.,10.,
ccc     .     'mix.layer temp. (.1 deg)')
ccc      call prtmsk(ip,smix,util3,idm,ii,jj,35.,100.,
ccc     .     'mx.lay. salin. (.01 mil)')
ccc!$OMP PARALLEL DO PRIVATE(j,l,i)
ccc!$OMP&         SCHEDULE(STATIC,jblk)
ccc      do j=1-margin,jj+margin
ccc        do l=1,isu(j)
ccc          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
ccc            util1(i,j)=umix(i,j)+ubavg(i,j,n)
ccc          enddo
ccc        enddo
ccc        do l=1,isv(j)
ccc          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
ccc            util2(i,j)=vmix(i,j)+vbavg(i,j,n)
ccc          enddo
ccc        enddo
ccc      enddo
ccc!$OMP END PARALLEL DO
ccc      call prtmsk(iu(2,1),util1(2,1),util3,idm,ii-2,jj,0.,1000.,
ccc     .     'mix.layer u vel. (mm/s)')
ccc      call prtmsk(iv(1,2),util2(1,2),util3,idm,ii,jj-2,0.,1000.,
ccc     .     'mix.layer v vel. (mm/s)')
ccc      call prtmsk(iu(2,1),ubavg(2,1,n),util3,idm,ii-2,jj,0.,1000.,
ccc     .     'barotrop. u vel. (mm/s)')
ccc      call prtmsk(iv(2,1),vbavg(1,2,n),util3,idm,ii,jj-2,0.,1000.,
ccc     .     'barotrop. v vel. (mm/s)')
      return
      end subroutine archiv

      subroutine archiv_tile(n, kkout, iyear,iday,ihour, intvl)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer   n, kkout, iyear,iday,ihour
      real      sssc,sstc
      character intvl*3
c
      include 'stmt_fns.h'
c
c --- write a partial archive file on a tile by tile basis.
c
      character*12 cdir
      character*80 cformat
      logical      lexist
      integer      i,j,k,ktr,l,ldot,nop,nopa
      real         coord,xmin,xmax
c
c --- only write archive when the corresponing directory exists
c
      write(cdir,'(a6,i5.5,a1)') 'ARCHT/',mnproc,'/'
      inquire(file=cdir(1:11),exist=lexist)
      if     (.not.lexist) then
        call xcsync(no_flush)  !called on all tiles, see end of routine
        return
      endif
c
      ldot = index(flnmarct,'.',back=.true.)
      if     (ldot.eq.0) then
        if     (mnproc.eq.1) then
        write (lp,*) 'need decimal point in flnmarct'
        write (lp,*) 'flnmarct = ',trim(flnmarct)
        endif
        call xchalt('(flnmarct)')
               stop '(flnmarct)'
      endif
      ldot = min(ldot,len(flnmarct)-11)  !need 11 characters for archive date
c
      if     (tilefq.ge.1.0/24.0) then
c ---   indicate the archive date
        write(flnmarct(ldot+1:ldot+11),'(i4.4,a1,i3.3,a1,i2.2)') 
     &   iyear,'_',iday,'_',ihour
        ldot=ldot+11
      else
c ---   indicate the archive time step
        write(flnmarct(ldot+1:ldot+11),'(i11.11)') nstep
        ldot=ldot+11
      endif
      nopa=13
      nop =13+uoff
c
      call ztiopf(cdir//flnmarct(1:ldot)//'.A', 'new', nopa)
      open (unit=nop,file=cdir//flnmarct(1:ldot)//'.B',status='new') !uoff+13
      write(nop,116) ctitle,iversn,iexpt,yrflag,i0+1,j0+1,ii,jj
      call flush(nop)
 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''i1    '' = longitudinal array starting index'/
     & i5,4x,'''j1    '' = latitudinal  array starting index'/
     & i5,4x,'''ii    '' = longitudinal array size'/
     & i5,4x,'''jj    '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
c
c --- surface fields
c
      coord=0.
c
      call ztiowr(montg1,ip,.true.,
     &            xmin,xmax, nopa, .false.)
c --- identify the equation of state on the first record
      write (nop,117) 'montg1  ',nstep,time,sigver,thbase,xmin,xmax
      call flush(nop)
      call ztiowr(srfhgt,ip,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'srfhgt  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      if     (sshflg.ne.0) then
c ---   write out steric SSH.
        call ztiowr(steric,ip,.true.,
     &              xmin,xmax, nopa, .false.)
        write (nop,117) 'steric  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
      endif !sshflg
c
      call ztiowr(surflx,ip,.true., xmin,xmax, nopa, .false.)
      write (nop,117) 'surflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      call ztiowr(salflx,ip,.true., xmin,xmax, nopa, .false.)
      write (nop,117) 'salflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
c
      call ztiowr(dpbl,ip,.true., xmin,xmax, nopa, .false.)
      write (nop,117) 'bl_dpth ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      call ztiowr(dpmixl(1-nbdy,1-nbdy,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'mix_dpth',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      if     (iceflg.ne.0) then
        call ztiowr(covice,ip,.true., xmin,xmax, nopa, .false.)
        write (nop,117) 'covice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        call ztiowr(thkice,ip,.true., xmin,xmax, nopa, .false.)
        write (nop,117) 'thkice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        call ztiowr(temice,ip,.true., xmin,xmax, nopa, .false.)
        write (nop,117) 'temice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
      endif  !write ice fields
c
c --- depth averaged fields
c
      call ztiowr(ubavg(1-nbdy,1-nbdy,n),iu,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'u_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      call ztiowr(vbavg(1-nbdy,1-nbdy,n),iv,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'v_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=sigma(k)
      call ztiowr(u(1-nbdy,1-nbdy,k,n),iu,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'u-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      call ztiowr(v(1-nbdy,1-nbdy,k,n),iv,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'v-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      call ztiowr(dp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'thknss  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      call ztiowr(temp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'temp    ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      call ztiowr(saln(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nopa, .false.)
      write (nop,117) 'salin   ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      do ktr= 1,ntracr
        call ztiowr(tracer(1-nbdy,1-nbdy,k,n,ktr),ip,.true.,
     &              xmin,xmax, nopa, .false.)
        write (nop,117) 'tracer  ',nstep,time,k,coord,xmin,xmax
        call flush(nop)
      enddo !ktr
      if     (difout) then
        call ztiowr(vcty(1-nbdy,1-nbdy,k+1),ip,.true.,
     &              xmin,xmax, nopa, .false.)
        write (nop,117) 'viscty  ',nstep,time,k,coord,xmin,xmax
        call flush(nop)
        call ztiowr(dift(1-nbdy,1-nbdy,k+1),ip,.true.,
     &              xmin,xmax, nopa, .false.)
        write (nop,117) 't-diff  ',nstep,time,k,coord,xmin,xmax
        call flush(nop)
        call ztiowr(difs(1-nbdy,1-nbdy,k+1),ip,.true.,
     &              xmin,xmax, nopa, .false.)
        write (nop,117) 's-diff  ',nstep,time,k,coord,xmin,xmax
        call flush(nop)
      endif
 75   continue
c
 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
c
      close (unit=nop)
      call ztiocl(nopa)
c
      call xcsync(no_flush)  !called on all tiles, see lexist above
      return
      end subroutine archiv_tile
c>
c> Revision history
c>
c> Nov  2002 - additional surface data in .txt output
c> Jun  2006 - dsur1p for .txt only surface output
c> Jun  2006 - archi .txt output
c> May  2007 - no diaflx output for K-profile based mixed layer models
c> May  2007 - removed mixed layer fields and th3d from the archive file
c> Feb  2008 - optionally added steric SSH to the archive file
c> Jun  2008 - added archiv_tile for per-tile archive output
