      program testxc
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      logical, parameter :: lregion = .true.   ! input regional.depth?
c
c     test xcaget and xcaput.
c
      integer i,j,kbad,ksea,nrecl
      real*4  depth(itdm,jtdm)
      real    aorig(itdm,jtdm),aotgt(itdm,jtdm)
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &        atrgt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
      common/test/ aotgt  ! required by xcaput
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
          aotgt(i,j) = 0.0
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
c     xcaget onto one processor.
c
      call xcaget(aotgt,atile, 1)
c
      if     (mnproc.eq.1) then
        write(lp,*) 'call xcaget(aotgt,atile, 1)'
        kbad = count( aotgt(:,:).ne.aorig(:,:) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe1: arrays are identical'
        else
          write(lp,*) 'pe1: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
      if     (mnproc.eq.2) then
        kbad = count( aotgt(:,:).ne.0.0 )
        if     (kbad.eq.0) then
          write(lp,*) 'pe2: array is zero'
        else
          write(lp,*) 'pe1: array has ',kbad,' non-zeros'
        endif
      endif
      call xcsync(flush_lp)
c
c     xcaput from one processor.
c
      call xcaput(aotgt,atrgt, 1)
c
      if     (mnproc.eq.1) then
        write(lp,*) 'call xcaput(aotgt,atrgt, 1)'
        kbad = count( atrgt(1:ii,1:jj) .ne. atile(1:ii,1:jj) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe1: arrays are identical'
        else
          write(lp,*) 'pe1: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
      if     (mnproc.eq.2) then
        write(lp,*) 'call xcaput(aotgt,atrgt, 1)'
        kbad = count( atrgt(1:ii,1:jj) .ne. atile(1:ii,1:jj) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe2: arrays are identical'
        else
          write(lp,*) 'pe2: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
c
c     xcaget onto all processors.
c
      call xcaget(aotgt,atile, 0)
c
      if     (mnproc.eq.1) then
        write(lp,*) 'call xcaget(aotgt,atile, 0)'
        kbad = count( aotgt(:,:).ne.aorig(:,:) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe1: arrays are identical'
        else
          write(lp,*) 'pe1: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
      if     (mnproc.eq.2) then
        write(lp,*) 'call xcaget(aotgt,atile, 0)'
        kbad = count( aotgt(:,:).ne.aorig(:,:) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe2: arrays are identical'
        else
          write(lp,*) 'pe2: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
c
c     xcaput from all processors.
c
      atrgt = 0.0
      call xcaput(aotgt,atrgt, 0)
c
      if     (mnproc.eq.1) then
        write(lp,*) 'call xcaput(aotgt,atrgt, 0)'
        kbad = count( atrgt(1:ii,1:jj) .ne. atile(1:ii,1:jj) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe1: arrays are identical'
        else
          write(lp,*) 'pe1: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
      if     (mnproc.eq.2) then
        kbad = count( atrgt(1:ii,1:jj) .ne. atile(1:ii,1:jj) )
        if     (kbad.eq.0) then
          write(lp,*) 'pe2: arrays are identical'
        else
          write(lp,*) 'pe2: arrays have ',kbad,' differing elements'
        endif
      endif
      call xcsync(flush_lp)
c
      call xcstop('(normal)')
             stop '(normal)'
      end
