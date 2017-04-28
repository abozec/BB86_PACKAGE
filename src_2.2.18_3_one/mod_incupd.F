      module mod_incupd
      use mod_xc  ! HYCOM communication interface
c
      implicit none
c
c --- HYCOM incremental updating (for data assimilation)
c
      integer, save, public  ::
     & incflg,   ! incremental update flag (0=no, 1=yes, 2=full-velocity)
     & incstp,   ! no. timesteps for full update (1=full insertion)
     & incupf    ! number of days of incremental updating input
c
      integer, save, private :: 
     & ncount,   ! increment time step counter
     & ncountd   ! increment day counter
c
      real*8, save, private :: 
     & dtimeu    ! next days increment field
c
      real,    allocatable, dimension(:,:),
     &         save, private ::
     & ubinc,    ! ubaro increment
     & vbinc     ! vbaro increment

c
      real,    allocatable, dimension(:,:,:),
     &         save, private ::
     &  tinc,    !     t increment
     & thinc,    !    th increment
     & dpinc,    !    dp increment
     &  uinc,    !     u increment
     &  vinc     !     v increment

      contains

      subroutine incupd_init(dtime0)
c
      real*8 dtime0
c
c --- subroutine used to calculate increment field for the incremental updating
c --- version: dec 2005
c
      integer i,j,l,k
      logical lopen
c
c --- allocate arrays
c
      allocate(  tinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &          thinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &          dpinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &           uinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &           vinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &          ubinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          vbinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     )
c
c --- set counter to zero
c
      ncount=0
      ncountd=0
      dtimeu=1.d-6
c
c --- read the target fields, and initialize the "inc" arrays.
c
      call incupd_read(dtime0)
c
      return
      end subroutine incupd_init

      subroutine incupd_rd(dtime0)
c
      real*8 dtime0
c
c --- subroutine used to calculate increment field for the incremental updating
c --- version: dec 2005
c
      integer i,j,l,k
      logical lopen
c
      if     (ncountd.gt.incupf) then
         if (mnproc.eq.1) then
         write(lp,*) '... ended updating fields with increments ...'
         write(lp,*) 'ncountd= ',ncountd
         write(lp,*)
         endif !1st tile
         call xcsync(flush_lp)
         return
      endif
c
c --- read the target fields, and initialize the "inc" arrays.
c
      call incupd_read(dtime0)
c
      return
      end subroutine incupd_rd

      subroutine incupd(n)
c
      include 'common_blocks.h'
C
      integer n
c
c**********
c*
c 1)  update hycom variables with increments.
c
c 2)  parameters:
c
c     output:
c      incremental updated model variables
c
c 4)  Ole Martin Smedstad (PSI), December 2005
c
c**********
c
      real       zero,one
      parameter (zero=0.0, one=1.0)
c
      integer i,j,k,l
      real    utotij,vtotij
c
      include 'stmt_fns.h'
c
c --- update counter
c
      if     (incstp.ne.1) then
        ncount=ncount+1
      endif
c
      margin=0
c
      if     (ncount.gt.incstp) then
        if     (ncount.eq.incstp+1) then
          if (mnproc.eq.1) then
          write(lp,*) '... ended updating fields with increments ...'
          write(lp,*) 'ncount= ',ncount
          write(lp,*)
          endif !1st tile
          call xcsync(flush_lp)
        endif !ncount==incstp+1
        return
      endif !ncount>incstp

      if     (ncountd.gt.incupf) then
        if (mnproc.eq.1) then
         write(lp,*) '... ended updating fields with increments ...'
         write(lp,*) 'ncountd= ',ncountd
         write(lp,*)
        endif !1st tile
         call xcsync(flush_lp)
         return
      endif
c
      if (mnproc.eq.1) then
      write(lp,*)
      if     (incflg.eq.1) then
        write(lp,'(2a)') 'update fields with increments, ',
     &                   'but not ubavg and vbavg'
      else   !incflg.eq.2
        write(lp,'(2a)') 'update fields with increments, ',
     &                   'including ubavg and vbavg'
      endif !incflg
      write(lp,*) 'ncount= ',ncount
      endif !1st tile
      call xcsync(flush_lp)
c
c --- incremental update of dp (dpu, dpv).
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do k=1,kk-1
              dp(i,j,k,n) = dp(i,j,k,n) + dpinc(i,j,k)
               p(i,j,k+1) = min( p(i,j,k) + dp(i,j,k,n),
     &                           p(i,j,kk+1) )  !min may be unnecessary
            enddo !k
c ---       dp must be non-negative.
            do k=1,kk-1
               p(i,j,k+1) = max( p(i,j,k+1),  p(i,j,k) )
              dp(i,j,k,n) =      p(i,j,k+1) - p(i,j,k)
            enddo ! k
c ---       layer kk always touches the bottom
            dp(i,j,kk,n) = p(i,j,kk+1) - p(i,j,kk)
          enddo !i
        enddo !l 
      enddo !j
!$OMP END PARALLEL DO
c
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n),
     &            dpv(1-nbdy,1-nbdy,1,n),
     &            p,depthu,depthv, 0)
c
c --- incremental update of the other fields.
c --- salinity from updated th&S.
c --- rebalance u and v via utotij and vtotij.
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,utotij,vtotij)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do k=1,kk
              if     ( tinc(i,j,k).ne.0.0 .or.
     &                thinc(i,j,k).ne.0.0     ) then
                temp(i,j,k,n) = temp(i,j,k,n) +  tinc(i,j,k)
                th3d(i,j,k,n) = th3d(i,j,k,n) + thinc(i,j,k)
                saln(i,j,k,n) = sofsig(th3d(i,j,k,n)+thbase,
     &                                 temp(i,j,k,n)       )
              endif !non-zero increment
            enddo ! k
          enddo !i
        enddo !l
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            utotij = 0.0
            do k=1,kk
              u(i,j,k,n) = u(i,j,k,n) + uinc(i,j,k)
              utotij = utotij + u(i,j,k,n)*dpu(i,j,k,n)
            enddo ! k
            utotij=utotij/depthu(i,j)
            do k=1,kk
              u(i,j,k,n) = u(i,j,k,n) - utotij
            enddo ! k
            if     (incflg.eq.2) then !update ubavg
              ubavg(i,j,n) = ubavg(i,j,n) + ubinc(i,j)
*             ubavg(i,j,n) = ubavg(i,j,n) + ubinc(i,j) + utotij
            endif !incflg==2
          enddo !i
        enddo !l
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vtotij = 0.0
            do k=1,kk
              v(i,j,k,n) = v(i,j,k,n) + vinc(i,j,k)
              vtotij = vtotij + v(i,j,k,n)*dpv(i,j,k,n)
            enddo ! k
            vtotij=vtotij/depthv(i,j)
            do k=1,kk
              v(i,j,k,n) = v(i,j,k,n) - vtotij
            enddo ! k
            if     (incflg.eq.2) then !update vbavg
              vbavg(i,j,n) = vbavg(i,j,n) + vbinc(i,j)
*             vbavg(i,j,n) = vbavg(i,j,n) + vbinc(i,j) + vtotij
            endif !incflg==2
          enddo !i
        enddo !l
      enddo ! j
!$OMP END PARALLEL DO
c
      if (mnproc.eq.1) then
       write(lp,*) 'finished incupdate',ncount
       write(lp,*)
      endif !1st tile
      call xcsync(flush_lp)
c
      return
      end subroutine incupd

      subroutine incupd_read(dtime)
      use mod_za  ! HYCOM I/O interface
c
      real*8    dtime
c
      include 'common_blocks.h'
c
c --- input 3-d HYCOM fields (from an archive file) on model day dtime.
c --- directly insert the input covice and thkice (if they exist).
c --- calculate the increment between the input and the initial state.
c
c --- filenames incup/incupd.iyear_iday_ihour.[ab].
c --- I/O and array I/O unit 925 used here, but not reserved.
c
      logical    ldebug_incupd_read
      parameter (ldebug_incupd_read=.false.)
c
      character flnm*24, cline*80, cvarin*6, cfield*8
      integer   i,idmtst,ios,j,jdmtst,k,l,layer,nskip
      integer   iyear,iday,ihour
      logical   nodens
      real      tincstp
*     real      sumdp1,sumdp2,sumdpinc
c
      integer   nstep0
      real*8    dtime0
c
      include 'stmt_fns.h'
c
      call forday(dtime, yrflag, iyear,iday,ihour)
c
      write(flnm,'("incup/incupd.",i4.4,"_",i3.3,"_",i2.2)')
     &                           iyear,iday,ihour
c
      if(dtime.ge.dtimeu) then
c
      ncountd=ncountd+1
      ncount=0
c
      if     (ncountd.gt.incupf) then
        if (mnproc.eq.1) then
         write(lp,*) '... ended updating fields with increments ...'
         write(lp,*) 'ncountd= ',ncountd
         write(lp,*)
        endif !1st tile
         call xcsync(flush_lp)
         return
      endif
c
      if (mnproc.eq.1) then
      write(lp,*) 'read incremental updating ...'
      write(lp,*) 'ncountd ...',ncountd
      write (lp,*) 'incupd_read: ',flnm
      write (lp,*) '       time: ',dtime
      write (lp,*) 'iyear,iday,ihour: ',iyear,iday,ihour
      endif !1st tile
      call xcsync(flush_lp)
c
      call zaiopf(flnm//'.a','old', 925)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+925,file=flnm//'.b',form='formatted',
     &        status='old',action='read')
c
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
c
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
      endif !1st tile
c
      call zagetc(cline,ios, uoff+925)
      read(cline,*) idmtst,cvarin
*     if     (mnproc.eq.1) then
*     write(lp,*) cvarin,' = ',idmtst
*     endif !1st tile
      if (cvarin.ne.'idm   ') then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read - input ',cvarin,
     &                        ' but should be idm   '
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read)')
               stop '(incupd_read)'
      endif
      call zagetc(cline,ios, uoff+925)
      read(cline,*) jdmtst,cvarin
*     if     (mnproc.eq.1) then
*     write(lp,*) cvarin,' = ',jdmtst
*     endif !1st tile
      if (cvarin.ne.'jdm   ') then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read - input ',cvarin,
     &                        ' but should be jdm   '
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read)')
               stop '(incupd_read)'
      endif
c
      if (idmtst.ne.itdm .or. jdmtst.ne.jtdm) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',itdm,  jtdm,  '  (dimensions.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read)')
               stop '(incupd_read)'
      endif
c
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        read (uoff+925,*)
      endif
c
c --- skip (most) surface fields.
c
      call zaiosk(925)
      call zagetc(cline,ios, uoff+925)
      i = index(cline,'=')
      read(cline(i+1:),*) nstep0,dtime0,layer
      if     (mnproc.eq.1) then
        write(lp,*) 'dtime0= ',dtime0
      endif
      if (dtime0.ne.dtime) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read - input ',dtime0,
     &                      ' but dtime should be ',dtime
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read)')
               stop '(incupd_read)'
      endif
      nodens = layer.ne.0  !new or original archive type
      if     (nodens .and. layer.ne.sigver) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read - input ',layer,
     &                     ' sigver but should be ',sigver
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read)')
               stop '(incupd_read)'
      endif
c
c assumes that there is a new incremental updating file once a day
c for "incupf" days, see blkdat.input
c
      dtimeu=dtime0+1.d0
c
      if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'dtime, dtime0, dtimeu = ',dtime,
     &               dtime0, dtimeu
        write(lp,*)
      endif !1st tile
      call xcsync(flush_lp)
c
      if     (nodens) then
        do i= 2,6
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo
      else
        do i= 2,11
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo
      endif
c
      call rd_archive(ubinc, cfield,layer, 925)  !u_btrop or covice or mix_dpth
      if     (cfield.eq.'mix_dpth') then
c ---   archive contains 'steric  '
        call rd_archive(ubinc, cfield,layer, 925)  !u_btrop or covice
      endif
      if     (mnproc.eq.1) then
      write(lp,'(2a)') "surface: ",cfield
      endif
      call xcsync(flush_lp)
      if     (cfield.eq.'covice  ') then
c
c ---   directly insert covice and thkice.
c
        call rd_archive(util5, cfield,layer, 925)  !thkice
        if     (mnproc.eq.1) then
        write(lp,'(2a)') "surface: ",cfield
        endif
        call xcsync(flush_lp)
!$OMP   PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              covice(i,j)=ubinc(i,j)
              thkice(i,j)=util5(i,j)
            enddo !i
          enddo !l 
        enddo !j
!$OMP   END PARALLEL DO
        call zaiosk(925)  !temice
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call rd_archive(ubinc, cfield,layer, 925)
        if     (mnproc.eq.1) then
        write(lp,'(2a)') "surface: ",cfield
        endif
        call xcsync(flush_lp)
      endif
      call rd_archive(vbinc, cfield,layer, 925)
      if     (mnproc.eq.1) then
      write(lp,'(2a)') "surface: ",cfield
      endif
      call xcsync(flush_lp)
c
           if     (mnproc.eq.1) then
           write (lp,*) 'start 3-D archive file read'
           endif
           call xcsync(flush_lp)
c
c --- 3-d fields.
c
      nskip = 0
      do k=1,kk
        call rd_archive(uinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'u-vel.  ' .and. k.ne.2) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield,
     &           'error in incupd_read - expected ','u-vel.  '
          endif !1st tile
          call xcstop('(incupd_read)')
                 stop '(incupd_read)'
        elseif (cfield.ne.'u-vel.  ') then !k==2
c
c ---     count "tracer" fields (to be skipped)
c
          if     (mnproc.eq.1) then
          write(lp,'(2a)') "counting tracers: ",cfield
          endif
          do nskip= 2,99
            call rd_archive(uinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
            if     (mnproc.eq.1) then
            write(lp,'(2a)') "counting tracers: ",cfield
            endif
            if     (cfield.eq.'u-vel.  ') then
              exit
            endif
          enddo !nskip
          nskip = nskip - 1
          write(lp,'(a,i3)') "nskip =",nskip
        endif
        call rd_archive(vinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'v-vel.  ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield,
     &           'error in incupd_read - expected ','v-vel.  '
          endif !1st tile
          call xcstop('(incupd_read)')
                 stop '(incupd_read)'
        endif
c          if     (mnproc.eq.1) then
c            write (lp,*) 'read v-vel archive file'
c          endif
        call rd_archive(dpinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'thknss  ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield,
     &             'error in incupd_read - expected ','thknss  '
          endif !1st tile
          call xcstop('(incupd_read)')
                 stop '(incupd_read)'
        endif
c          if     (mnproc.eq.1) then
c            write (lp,*) 'read dpinc archive file'
c          endif
        call rd_archive(tinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'temp    ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield,
     &           'error in incupd_read - expected ','temp    '
          endif !1st tile
          call xcstop('(incupd_read)')
                 stop '(incupd_read)'
        endif
        if     (nodens) then
c ---     read salinity into thinc, later convert it to density
          call rd_archive(thinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
          if     (cfield.ne.'salin   ') then
            if     (mnproc.eq.1) then
            write(lp,'(/ a / a,a /)') cfield,
     &             'error in incupd_read - expected ','salin   '
            endif !1st tile
            call xcstop('(incupd_read)')
                   stop '(incupd_read)'
          endif
        else
c ---     skip salinity, and read density
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
          call rd_archive(thinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
          if     (cfield.ne.'density ') then
            if     (mnproc.eq.1) then
            write(lp,'(/ a / a,a /)') cfield,
     &             'error in incupd_read - expected ','density '
            endif !1st tile
            call xcstop('(incupd_read)')
                   stop '(incupd_read)'
          endif
        endif !nodens:else
c
c ---   skip (nskip) tracers
c
        do l= 1,nskip
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo !l
      enddo !k
c
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      close( unit=uoff+925)
      endif
      call zaiocl(925)
c
c --- calculate increments
c --- the "inc" reads, above, are full HYCOM fields (not increments yet).
c
      if(incstp.eq.1) then
        tincstp=1.0
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'tincstp=1.0        ',tincstp,incstp
        endif
      else
        tincstp=2.0/real(incstp)
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'tincstp=2.0/incstp ',tincstp,incstp
        endif
      endif !incstp
c
c
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,*) 'calculate t,s,u,v and dp increments'
      endif !1st tile
      call xcsync(flush_lp)
c
*     if     (iutest.gt.0 .and. jutest.gt.0) then
*       write(lp,*) '*',' iutest= ',iutest+i0,' jutest= ',jutest+j0,' *'
*       write(lp,*) '*********** dpinc input ************'
*               sumdp1=0.0
*               sumdp2=0.0
*               sumdpinc=0.0
*             write(lp,'(a)')
*    &                'k,dp1,dp2,dpinc='
*               do k= 1,kk
*                sumdp1=sumdp1+dp(iutest,jutest,k,1)
*                sumdp2=sumdp2+dp(iutest,jutest,k,2)
*                sumdpinc=sumdpinc+dpinc(iutest,jutest,k)
*                   write(lp,'(a,i3,3f20.5)')
*    &                'k= ',
*    &                 k,dp(iutest,jutest,k,1)*qonem,
*    &                 dp(iutest,jutest,k,2)*qonem,
*    &                 dpinc(iutest,jutest,k)*qonem
*                   call flush(lp)
*               enddo !k
*            write(lp,*) 'sumdp1,sumdp2= ', sumdp1*qonem,sumdp2*qonem
*            write(lp,*) 'sumdpinc= ', sumdpinc*qonem
*            call flush(lp)
*     endif
c
      margin=0
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubinc(i,j)=(ubinc(i,j) - ubavg(i,j,1))*tincstp
            do k=1,kk
c             use an approximate 2*dpu
              if     (dpinc(i,j,k)+dpinc(i-1,j,k).gt.2.0*onem) then
                uinc(i,j,k)=(uinc(i,j,k) - u(i,j,k,1))*tincstp
              else
                uinc(i,j,k)=0.0  !thin target layer
              endif
            enddo !k
          enddo !i
        enddo !l 
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbinc(i,j)=(vbinc(i,j) - vbavg(i,j,1))*tincstp
            do k=1,kk
c             use an approximate 2*dpv
              if     (dpinc(i,j,k)+dpinc(i,j-1,k).gt.2.0*onem) then
                vinc(i,j,k)=(vinc(i,j,k) - v(i,j,k,1))*tincstp
              else  
                vinc(i,j,k)=0.0  !thin target layer
              endif
            enddo !k
          enddo !i
        enddo !l 
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            do k=1,kk
              if     (dpinc(i,j,k).gt.onem) then
                if     (nodens) then
                  thinc(i,j,k)=sig(tinc(i,j,k),thinc(i,j,k))-thbase
                endif
                thinc(i,j,k)=(thinc(i,j,k) - th3d(i,j,k,1))*tincstp
                 tinc(i,j,k)= (tinc(i,j,k) - temp(i,j,k,1))*tincstp
              else
                 tinc(i,j,k)=0.0  !thin target layer
                thinc(i,j,k)=0.0  !thin target layer
              endif
              dpinc(i,j,k)=(dpinc(i,j,k) - dp(i,j,k,1))*tincstp
            enddo !k
          enddo !i
        enddo !l 
      enddo !j
!$OMP END PARALLEL DO
c
*     if     (iutest.gt.0 .and. jutest.gt.0) then
*       write(lp,*) '*',' iutest= ',iutest+i0,' jutest= ',jutest+j0,' *'
*       write(lp,*) '*********** dpinc out ************'
*             write(lp,'(a)')
*    &                'k,dp1,dp2,dpinc='
*               sumdpinc=0.0
*               do k= 1,kk
*                sumdpinc=sumdpinc+dpinc(iutest,jutest,k)
*                   write(lp,'(a,i3,3f20.5)')
*    &                'k= ',
*    &                 k,dp(iutest,jutest,k,1)*qonem,
*    &                 dp(iutest,jutest,k,2)*qonem,
*    &                 dpinc(iutest,jutest,k)*qonem
*                   call flush(lp)
*               enddo !k
*            write(lp,*) 'inc sumdpinc= ', sumdpinc*qonem
*            call flush(lp)
*     endif
c
      if (mnproc.eq.1) then
       write(lp,*) '... finnished reading incupd',dtime,dtime0
      endif !1st tile
      call xcsync(flush_lp)
c
      endif ! dtime
c
      return
      end subroutine incupd_read
c
      end module mod_incupd
c
c
c> Revision history:
c>
c> Feb  2006 - 1st module version
c> May  2006 - changed to read multiple increment files
