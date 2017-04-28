      program testxc
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      logical, parameter :: lregion = .true.   ! input regional.depth?
c
c     test xcsum.
c
      integer i,j,ksea,nrecl
      integer morig(itdm,jtdm)
      integer mtile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real*4  depth(itdm,jtdm)
      real    aorig(itdm,jtdm)
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    amino,amint,amaxo,amaxt
      real*8  sumo,sumt,sumjo(jtdm),sumjt(jtdm)
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
      amino =  1.e30
      amaxo = -1.e30
      do j= 1,jtdm
        do i= 1,itdm
          if     (depth(i,j).gt.0.0 .and. depth(i,j).lt.2.0**99) then
            morig(i,j) = 1
            if     (mod(j,2).eq.0) then
              aorig(i,j) =  i + (j-1)*100
            else
              aorig(i,j) = -i - (j-1)*100
            endif
            amino = min( amino, aorig(i,j) )
            amaxo = max( amaxo, aorig(i,j) )
          else
            morig(i,j) = 0
            aorig(i,j) = 0.0
          endif
        enddo
      enddo
c
      amint =  1.e30
      amaxt = -1.e30
      do j= 1,jj
        do i= 1,ii
          atile(i,j) = aorig(i+i0,j+j0)
          amint = min( amint, atile(i,j) )
          amaxt = max( amaxt, atile(i,j) )
        enddo
      enddo
c
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
          if     (i+i0.ge.1 .and. i+i0.le.itdm .and.
     &            j+j0.ge.1 .and. j+j0.le.jtdm      ) then
            mtile(i,j) = morig(i+i0,j+j0)
          else
            mtile(i,j) = 0
          endif
        enddo
      enddo
      if     (mnproc.eq.1) then
      write(lp,*) 'itdm,jtdm = ',itdm,jtdm
      write(lp,*) 'ii,  jj   = ',ii,  jj
      ksea = count( aorig(:,:).ne.0.0 )
      write(lp,*) 'sea, land = ',ksea,itdm*jtdm-ksea
      write(lp,*) 'aorig = ',aorig(1,1),aorig(1,2),
     +                       aorig(2,1),aorig(2,2)
      write(lp,*) 'atile = ',atile(1,1),atile(1,2),
     +                       atile(2,1),atile(2,2)
      write(lp,*) 'mtile = ',mtile(1,1),mtile(1,2),
     +                       mtile(2,1),mtile(2,2)
      write(lp,*)
      endif
      call xcsync(flush_lp)
c
c     test xcminr and xcmaxr
c
      call xcminr(amint)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xcminr(amint)'
      if     (abs((amint-amino)/amino).le.1.e-5) then
        write(lp,*) 'amint,amino = ',amint,amino
        write(lp,*) 'mins are identical  (',(amint-amino)/amino,')'
      else
        write(lp,*) 'amint,amino = ',amint,amino
        write(lp,*) 'min error = ',amint-amino,(amint-amino)/amino
      endif
      endif
      call xcsync(flush_lp)
c
      call xcmaxr(amaxt)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xcmaxr(amaxt)'
      if     (abs((amaxt-amaxo)/amaxo).le.1.e-5) then
        write(lp,*) 'amaxt,amaxo = ',amaxt,amaxo
        write(lp,*) 'maxs are identical  (',(amaxt-amaxo)/amaxo,')'
      else
        write(lp,*) 'amaxt,amaxo = ',amaxt,amaxo
        write(lp,*) 'max error = ',amaxt-amaxo,(amaxt-amaxo)/amaxo
      endif
      endif
      call xcsync(flush_lp)
c
c     test xcsum and xcsumj.
c
      call xcsum(sumt, atile,mtile)
      call xxsum(sumo, aorig,morig)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xcsum(sumt, atile,mtile)'
      if     (abs((sumt-sumo)/sumo).le.1.e-5) then
        write(lp,*) 'sumt,sumo = ',sumt,sumo
        write(lp,*) 'sums are identical  (',(sumt-sumo)/sumo,')'
      else
        write(lp,*) 'sumt,sumo = ',sumt,sumo
        write(lp,*) 'sum error = ',sumt-sumo,(sumt-sumo)/sumo
      endif
      endif
      call xcsync(flush_lp)
c
      call xcsumj(sumjt, atile,mtile)
      call xxsumj(sumjo, aorig,morig)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xcsumj(sumtj, atile,mtile)'
      i = 0
      do j= 1,jtdm
        if     (sumjo(j).eq.0.0) then
          if     (sumjt(j).ne.0.0) then
            write(lp,*) 'j,sumt,sumo = ',j,sumjt(j),sumjo(j)
            write(lp,*) 'j,sum error = ',j,sumjt(j)-sumjo(j)
            i = i + 1
          endif
        else
          if     (abs((sumjt(j)-sumjo(j))/sumjo(j)).gt.1.e-5) then
            write(lp,*) 'j,sumt,sumo = ',j,sumjt(j),sumjo(j)
            write(lp,*) 'j,sum error = ',j,sumjt(j)-sumjo(j)
            i = i + 1
          endif
        endif
      enddo
      if     (i.eq.0) then
        write(lp,*) 'j sums are identical'
      else
        write(lp,*) 'j sums differ in ',i,' rows'
      endif
      endif
      call xcsync(flush_lp)
c
c     new mask
c
      do j= 1,jtdm
        do i= 1,itdm
          morig(i,j) = mod(i+j,2)
        enddo
      enddo
c
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
          if     (i+i0.ge.1 .and. i+i0.le.itdm .and.
     &            j+j0.ge.1 .and. j+j0.le.jtdm      ) then
            mtile(i,j) = morig(i+i0,j+j0)
          else
            mtile(i,j) = 0
          endif
        enddo
      enddo
c
c     test.
c
      call xcsum(sumt, atile,mtile)
      call xxsum(sumo, aorig,morig)
      if     (mnproc.eq.1) then
      write(lp,*) 'call xcsum(sumt, atile,mtile)'
      if     (abs((sumt-sumo)/sumo).le.1.e-5) then
        write(lp,*) 'sumt,sumo = ',sumt,sumo
        write(lp,*) 'sums are identical  (',(sumt-sumo)/sumo,')'
      else
        write(lp,*) 'sumt,sumo = ',sumt,sumo
        write(lp,*) 'sum error = ',sumt-sumo,(sumt-sumo)/sumo
      endif
      endif
      call xcsync(flush_lp)
c
      call xcstop('(normal)')
             stop '(normal)'
      end
      subroutine xxsum(sum, a,mask)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      real*8  sum
      real    a(   itdm,jtdm)
      integer mask(itdm,jtdm)
c
c**********
c*
c  1) sum a 2-d array, where mask==1
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    sum             real*8         output    sum of a
c    a               real           input     source array
c    mask            integer        input     mask array
c*
c**********
c
      real*8     zero8
      parameter (zero8=0.0)
c
      real*8  sum8,sum8p,sum8j(jtdm)
      integer i,i1,j
c
c     row sums in 2*nbdy+1 wide strips.
c
!$OMP PARALLEL DO PRIVATE(j,i1,i,sum8,sum8p)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jtdm
        sum8 = zero8
        do i1=1,itdm,2*nbdy+1
          sum8p = zero8
          do i= i1,min(i1+2*nbdy,itdm)
            if     (mask(i,j).eq.1) then
              sum8p = sum8p + a(i,j)
            endif
          enddo
          sum8 = sum8 + sum8p
        enddo
        sum8j(j) = sum8  ! use of sum8 minimizes false sharing of sum8j
      enddo
!$OMP END PARALLEL DO
c
c     serial sum of rwo-sum loop.
c
      sum8 = sum8j(1)
      do j=2,jtdm
        sum8 = sum8 + sum8j(j)
      enddo
      sum = sum8
      return
c     end of xxsum.
      end
      subroutine xxsumj(sumj, a,mask)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      real*8  sumj(     jtdm)
      real    a(   itdm,jtdm)
      integer mask(itdm,jtdm)
c
c**********
c*
c  1) row sum a 2-d array, where mask==1
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    sumj            real*8         output    rwo sums of a
c    a               real           input     source array
c    mask            integer        input     mask array
c*
c**********
c
      real*8     zero8
      parameter (zero8=0.0)
c
      real*8  sum8,sum8p
      integer i,i1,j
c
c     row sums in 2*nbdy+1 wide strips.
c
!$OMP PARALLEL DO PRIVATE(j,i1,i,sum8,sum8p)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jtdm
        sum8 = zero8
        do i1=1,itdm,2*nbdy+1
          sum8p = zero8
          do i= i1,min(i1+2*nbdy,itdm)
            if     (mask(i,j).eq.1) then
              sum8p = sum8p + a(i,j)
            endif
          enddo
          sum8 = sum8 + sum8p
        enddo
        sumj(j) = sum8  ! use of sum8 minimizes false sharing of sumj
      enddo
!$OMP END PARALLEL DO
      return
c     end of xxsumj.
      end
