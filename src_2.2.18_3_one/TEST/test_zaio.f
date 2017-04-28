      program testxc
      use mod_xc    ! HYCOM communication interface
      use mod_za    ! HYCOM I/O interface
      implicit none
c
      logical, parameter :: lregion = .true.   ! input regional.depth?
c
      real*4     spval4
      parameter (spval4=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
c     test array I/O.
c
      integer i,ierr,irec,j,jerr,l,ksea,nrecl
      integer mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    aorig(itdm,jtdm)
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    ao(itdm),at(itdm),atmax,atmin,vsave
      real*4  b4(n2drec),depth(itdm,jtdm)
c
c --- machine-specific initialization
c
      call machine
c
c --- initialize SPMD processsing
c
      call xcspmd
c
c     Initialize array I/O.
c
      call zaiost
c
c     read in land/sea map?
c
      if     (lregion) then
        inquire(iolength=nrecl) depth
c
        open( unit=11, file='regional.depth.a',
     &        form='unformatted', action='read',
     &        access='direct', recl=nrecl)
        read( unit=11, rec=1) depth
        close(unit=11)
      else
        depth = 1.0
      endif
c
c     initialize.
c
      do i= 1,n2drec
        b4(i) = spval4
      enddo
c
      l    = 0
      ksea = 0
      do j= 1,jtdm
        do i= 1,itdm
          l = l + 1
          if     (depth(i,j).gt.0.0 .and. depth(i,j).lt.2.0**99) then
            aorig(i,j) = i + (j-1)*100
            b4(l) = aorig(i,j)
            ksea = ksea + 1
          else
            aorig(i,j) = spval4
          endif
        enddo
      enddo
c
      do j= 1,jj
        do i= 1,ii
          atile(i,j) = aorig(i+i0,j+j0)
        enddo
      enddo
      if     (mnproc.eq.1) then
      write(lp,*) 'itdm,jtdm = ',itdm,jtdm
      write(lp,*) 'idm, jdm  = ',idm, jdm 
      write(lp,*) 'sea, land = ',ksea,itdm*jtdm-ksea
      write(lp,*) 'aorig = ',aorig(1,1),aorig(1,2),
     +                       aorig(2,1),aorig(2,2)
      write(lp,*) 'atile = ',atile(1,1),atile(1,2),
     +                       atile(2,1),atile(2,2)
      write(lp,*)
      call flush(lp)
      endif
c
c     create test file.
c
      if     (mnproc.eq.1) then
        inquire(iolength=nrecl) b4
        open( unit=29, file='fort.029',
     &                 form='unformatted',
     &                 status='new',
     &                 action='write',
     &                 access='direct',
     &                 recl=nrecl)
        write(unit=29, rec=1) b4
        write(unit=29, rec=2) b4
        close(unit=29, status='keep')
      endif
c
      call xcsync(flush_lp)
c
c     tiled i/o on test file.
c
      call zaiopf('fort.029',  'old', 929)
      call zaiopn(             'new',  29)
      do irec= 1,2
c
        atile(:,:) = 0.0
c
c       read.
c
        call xcsync(flush_lp)
        if     (mnproc.eq.1) then
          write(lp,*) 'reading record ',irec
        endif
        call xcsync(flush_lp)
        call zaiord(atile, mask,.false., atmin,atmax, 929)
        call xcsync(flush_lp)
        if     (mnproc.eq.1) then
          write(lp,*) '    read complete ',atmin,atmax
        endif
        call xcsync(flush_lp)
c
c       test.
c
        jerr = 0
        do j= 1,jtdm
          vsave = vland
          vland = spval4
          call xclget(at,itdm, atile,1,j,+1,0, 0)
          vland = vsave
          call xxlget(ao,itdm, aorig,1,j,+1,0)
          call yycomp(ao,at,itdm ,ierr)
          if     (ierr.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,*) '        yycomp - j,ierr = ',j,ierr
            endif
          endif
          call xcsync(flush_lp)
          jerr = jerr + ierr
        enddo
        if     (mnproc.eq.1) then
          if     (jerr.eq.0) then
            write(6,*) '    read successfull'
          else
            write(6,*) '    read generated ',jerr,' bad elements'
          endif
        endif
        call xcsync(flush_lp)
c
c       write.
c
        if     (mnproc.eq.1) then
          write(6,*) 'writing record ',irec
        endif
        call xcsync(flush_lp)
        call zaiowr(atile, mask,.false., atmin,atmax,  29, .false.)
        call xcsync(flush_lp)
        if     (mnproc.eq.1) then
          write(6,*) '    write complete',atmin,atmax
        endif
        call xcsync(flush_lp)
      enddo
c
c     close.
c
      call zaiocl( 29)
      call zaiocl(929)
c
      call xcstop('(normal')
             stop '(normal)'
      end
      subroutine xxlget(aline,nl, a, i1,j1,ic,jc)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer nl,i1,j1,ic,jc
      real    aline(nl)
      real    a(itdm,jtdm)
c
c**********
c*
c  1) extracts the line a(i1:i1+(nl-1)*ic:ic,j1:j1+(nl-1)*jc:jc),
c      w.r.t. the 2-d grid.
c
c  2) ic and jc can each be -1, 0, or +1.
c*
c**********
c
      integer i
c
*     write(lp,'(a,5i5)') 'xxlget - nl,i1,j1,ic,jc = ',nl,i1,j1,ic,jc
      if     (jc.eq.0) then
        do i= 1,nl
          aline(i) = a(i1+ic*(i-1),j1)
        enddo
      elseif (ic.eq.0) then
        do i= 1,nl
          aline(i) = a(i1,j1+jc*(i-1))
        enddo
      else
        do i= 1,nl
          aline(i) = a(i1+ic*(i-1),j1+jc*(i-1))
        enddo
      endif
      return
c     end of xxlget.
      end
      subroutine yycomp(a,b,n, ierr)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer n,ierr
      real    a(n),b(n)
c
c**********
c*
c  1) tests if a and b are identical.
c*
c**********
c
      integer i
c
      ierr = 0
      do 110 i= 1,n
        if     (a(i).ne.b(i)) then
          if     (mnproc.eq.1 .and. mod(ierr,20).eq.0) then
          write(lp,*) 'i,a,b = ',i,a(i),b(i)
          endif
          ierr = ierr + 1
        endif
  110 continue
*     if     (mnproc.eq.1) then
*       if     (ierr.eq.0) then
*         write(lp,*) 'arrays are identical'
*         write(lp,*) 
*       else
*         write(lp,*) 'arrays have ',ierr,' differing elements'
*         write(lp,*) 
*       endif
*     endif
      return
c     end of yycomp.
      end
