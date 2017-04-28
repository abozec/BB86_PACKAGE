      program testxc
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      logical, parameter :: lregion = .true.   ! input regional.depth?
c
c     test xclget.
c
      integer i,j,ksea,nrecl
      real*4  depth(itdm,jtdm)
      real    aorig(itdm,jtdm)
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    aoi(jtdm),aoj(itdm),ati(jtdm),atj(itdm)
c
c --- machine-specific initialization
c
      call machine
c
c --- initialize SPMD processsing
c
      call xcspmd
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
      do j= 1,jtdm
        do i= 1,itdm
          if     (depth(i,j).gt.0.0 .and. depth(i,j).lt.2.0**99) then
            aorig(i,j) = i + (j-1)*100
          else
            aorig(i,j) = 0.0
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
      ksea = count( aorig(:,:).ne.0.0 )
      write(lp,*) 'sea, land = ',ksea,itdm*jtdm-ksea
      write(lp,*) 'aorig = ',aorig(1,1),aorig(1,2),
     +                       aorig(2,1),aorig(2,2)
      write(lp,*) 'atile = ',atile(1,1),atile(1,2),
     +                       atile(2,1),atile(2,2)
      write(lp,*)
      call flush(lp)
      endif
c
c     test.
c
      call xclget(ati,jtdm, atile,itdm/2,1,0,+1, 0)
      call xxlget(aoi,jtdm, aorig,itdm/2,1,0,+1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(ati,jtdm, atile, itdm/2,1,0,+1)'
      endif
      call yycomp(aoi,ati,jtdm)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(ati,10, atile,itdm/2,jtdm/2-5,0,+1, 0)
      call xxlget(aoi,10, aorig,itdm/2,jtdm/2-5,0,+1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(ati,10, atile, itdm/2,jtdm/2-5,0,+1)'
      endif
      call yycomp(aoi,ati,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(ati,jtdm, atile,itdm/2,jtdm,0,-1, 0)
      call xxlget(aoi,jtdm, aorig,itdm/2,jtdm,0,-1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(ati,jtdm, atile, itdm/2,jtdm,0,-1)'
      endif
      call yycomp(aoi,ati,jtdm)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(ati,10, atile,itdm/2,jtdm/2+5,0,-1, 0)
      call xxlget(aoi,10, aorig,itdm/2,jtdm/2+5,0,-1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(ati,10, atile, itdm/2,jtdm/2+5,0,-1)'
      endif
      call yycomp(aoi,ati,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,itdm, atile,1,jtdm/2,+1,0, 0)
      call xxlget(aoj,itdm, aorig,1,jtdm/2,+1,0)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,itdm, atile, 1,jtdm/2,+1,0)'
      endif
      call yycomp(aoj,atj,itdm)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,10, atile,itdm/2-5,jtdm/2,+1,0, 0)
      call xxlget(aoj,10, aorig,itdm/2-5,jtdm/2,+1,0)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,10, atile, itdm/2-5,jtdm/2,+1,0)'
      endif
      call yycomp(aoj,atj,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,itdm, atile,itdm,jtdm/2,-1,0, 0)
      call xxlget(aoj,itdm, aorig,itdm,jtdm/2,-1,0)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,itdm, atile, itdm,jtdm/2,-1,0)'
      endif
      call yycomp(aoj,atj,itdm)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,10, atile,itdm/2+5,jtdm/2,-1,0, 0)
      call xxlget(aoj,10, aorig,itdm/2+5,jtdm/2,-1,0)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,10, atile, itdm/2+5,jtdm/2,-1,0)'
      endif
      call yycomp(aoj,atj,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,10, atile,itdm/2+5,jtdm/2+5,-1,-1, 0)
      call xxlget(aoj,10, aorig,itdm/2+5,jtdm/2+5,-1,-1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,10, atile, itdm/2+5,jtdm/2+5,-1,-1)'
      endif
      call yycomp(aoj,atj,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,10, atile,itdm/2+5,jtdm/2-5,-1,+1, 0)
      call xxlget(aoj,10, aorig,itdm/2+5,jtdm/2-5,-1,+1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,10, atile, itdm/2+5,jtdm/2-5,-1,+1)'
      endif
      call yycomp(aoj,atj,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,10, atile,itdm/2-5,jtdm/2+5,+1,-1, 0)
      call xxlget(aoj,10, aorig,itdm/2-5,jtdm/2+5,+1,-1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,10, atile, itdm/2-5,jtdm/2+5,+1,-1)'
      endif
      call yycomp(aoj,atj,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xclget(atj,10, atile,itdm/2-5,jtdm/2-5,+1,+1, 0)
      call xxlget(aoj,10, aorig,itdm/2-5,jtdm/2-5,+1,+1)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xclget(atj,10, atile, itdm/2-5,jtdm/2-5,+1,+1)'
      endif
      call yycomp(aoj,atj,10)
      if     (mnproc.eq.1) then
      call flush(lp)
      endif
c
      call xcstop('(normal)')
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
*         if     (i.le.5) then
*           write(lp,'(a,3i5,f10.2)') 'xxlget - l,i,j,aline = ',
*    +                               i,i1+ic*(i-1),j1,aline(i)
*         endif
        enddo
      elseif (ic.eq.0) then
        do i= 1,nl
          aline(i) = a(i1,j1+jc*(i-1))
*         if     (i.le.5) then
*           write(lp,'(a,3i5,f10.2)') 'xxlget - l,i,j,aline = ',
*    +                               i,i1,j1+jc*(i-1),aline(i)
*         endif
        enddo
      else
        do i= 1,nl
          aline(i) = a(i1+ic*(i-1),j1+jc*(i-1))
*         if     (i.le.5) then
*           write(lp,'(a,3i5,f10.2)') 'xxlget - l,i,j,aline = ',
*    +                               i,i1+ic*(i-1),j1+jc*(i-1),
*    +                               aline(i)
*         endif
        enddo
      endif
      return
c     end of xxlget.
      end
      subroutine yycomp(a,b,n)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer n
      real    a(n),b(n)
c
c**********
c*
c  1) tests if a and b are identical.
c*
c**********
c
      integer i,ierr
c
      ierr = 0
      do 110 i= 1,n
        if     (a(i).ne.b(i)) then
          if     (mnproc.eq.1) then
          write(lp,*) 'i,a,b = ',i,a(i),b(i)
          endif
          ierr = ierr + 1
        endif
  110 continue
      if     (mnproc.eq.1) then
        if     (ierr.eq.0) then
          write(lp,*) 'arrays are identical'
          write(lp,*) 
        else
          write(lp,*) 'arrays have ',ierr,' differing elements'
          write(lp,*) 
        endif
      endif
      return
c     end of yycomp.
      end
