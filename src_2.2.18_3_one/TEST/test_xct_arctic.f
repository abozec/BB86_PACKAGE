      program testxc
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      logical, parameter :: lregion = .true.   ! input regional.depth?
c
c     test xctilr, for arctic regions only.
c
      integer i,j,k,ksea,nrecl
      real*4  depth(itdm,jtdm)
      real    aorig(itdm,jtdm)
      real    ahalo(1-nbdy:itdm+nbdy,1-nbdy:jtdm+nbdy)
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
c
      integer halo_type,l
      integer halo_t(8)
      data halo_t / 1, 2, 3, 4, 11, 12, 13, 14 /
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
      do l= 1,8
        halo_type = halo_t(l)
c
c     initialize for halo type.
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
      call arctic_fix( aorig,1, halo_type)
      call arctic_halo(aorig,1, halo_type, ahalo)
c
      do k= 1,kdm
        do j= 1,jj
          do i= 1,ii
            if     (aorig(i+i0,j+j0).gt.0.0) then
              atile(i,j,k) = aorig(i+i0,j+j0) + (k-1)*10000
            elseif (aorig(i+i0,j+j0).lt.0.0) then
              atile(i,j,k) = aorig(i+i0,j+j0) - (k-1)*10000
            else
              atile(i,j,k) = 0.0
            endif
          enddo
        enddo
      enddo
      if     (mnproc.eq.ijpr) then
      write(lp,*) 
      write(lp,*) 'halo_type = ',halo_type
      write(lp,*) 'itdm,jtdm = ',itdm,jtdm
      write(lp,*) 'ii,  jj   = ',ii,  jj  
      ksea = count( aorig(:,:).ne.0.0 )
      write(lp,*) 'sea, land = ',ksea,itdm*jtdm-ksea
      write(lp,'(a,4f12.2)') 'aorig1 = ',aorig(1, 1),aorig(ii, 1),
     +                                   aorig(1,jj),aorig(ii,jj)
      write(lp,'(a,4f12.2)') 'atile1 = ',atile(1, 1,1),atile(ii, 1,1),
     +                                   atile(1,jj,1),atile(ii,jj,1)
      write(lp,'(a,4f12.2)') 'aorig9 = ',aorig( 1, 1)+8*10000,
     +                                   aorig(ii, 1)+8*10000,
     +                                   aorig( 1,jj)+8*10000,
     +                                   aorig(ii,jj)+8*10000
      write(lp,'(a,4f12.2)') 'atile9 = ',atile(1, 1,9),atile(ii, 1,9),
     +                                   atile(1,jj,9),atile(ii,jj,9)
      write(lp,*)
      endif
      call xcsync(flush_lp)
c
c     test.
c
      call zztile(atile,1,1)
      call xctilr(atile,1,1, 1,1, halo_type)
      call yytile(atile,1,1, 1,1, halo_type, ahalo)
c
      call zztile(atile,1,1)
      call xctilr(atile,1,1, 0,1, halo_type)
      call yytile(atile,1,1, 0,1, halo_type, ahalo)
c
      call zztile(atile,1,1)
      call xctilr(atile,1,1, 1,0, halo_type)
      call yytile(atile,1,1, 1,0, halo_type, ahalo)
c
      call zztile(atile,1,1)
      call xctilr(atile,1,1, nbdy,nbdy, halo_type)
      call yytile(atile,1,1, nbdy,nbdy, halo_type, ahalo)
c
      call zztile(atile,3,9)
      call xctilr(atile,3,9, nbdy,nbdy, halo_type)
      call yytile(atile,3,9, nbdy,nbdy, halo_type, ahalo)
c
      call zztile(atile,1,kdm)
      call xctilr(atile,1,kdm, nbdy,nbdy, halo_type)
      call yytile(atile,1,kdm, nbdy,nbdy, halo_type, ahalo)
c
      call zztile(atile,1,kdm)
      call xctilr(atile,1,kdm, 2,   nbdy, halo_type)
      call yytile(atile,1,kdm, 2,   nbdy, halo_type, ahalo)
c
      call zztile(atile,1,kdm)
      call xctilr(atile,1,kdm, nbdy,   2, halo_type)
      call yytile(atile,1,kdm, nbdy,   2, halo_type, ahalo)
c
      enddo !l
c
      call xcstop('(normal)')
             stop '(normal)'
      end
      subroutine arctic_fix(aorig,ld, halo_type)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer ld,halo_type
      real    aorig(itdm,jtdm,ld)
c
c     make top edge consistent for arctic patch
c
      real    s
      integer i,io,j,jo,k
c
      if     (halo_type.lt.10) then
        s =  1.0  !scalar
      else
        s = -1.0  !vector
      endif
c
      do k= 1,ld
        j = jtdm
          if     (halo_type.eq.1 .or. halo_type.eq.11) then !p-grid
            jo = jtdm-1
            do i= 1,itdm
              io = itdm-mod(i-1,itdm)
              aorig(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          elseif (halo_type.eq.2 .or. halo_type.eq.12) then !q-grid
            jo = jtdm
            do i= 1,itdm/2
              io = mod(itdm-(i-1),itdm)+1
              aorig(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          elseif (halo_type.eq.3 .or. halo_type.eq.13) then !u-grid
            jo = jtdm-1
            do i= 1,itdm
              io = mod(itdm-(i-1),itdm)+1
              aorig(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          elseif (halo_type.eq.4 .or. halo_type.eq.14) then !v-grid
            jo = jtdm
            do i= 1,itdm/2
              io = itdm-mod(i-1,itdm)
              aorig(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          endif !halo_type
      enddo !k
      return
      end
      subroutine arctic_halo(aorig,ld, halo_type, ahalo) 
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer ld,halo_type
      real    aorig(itdm,jtdm,ld)
      real    ahalo(1-nbdy:itdm+nbdy,1-nbdy:jtdm+nbdy,ld)
c
c     copy aorig into ahalo including the halo.
c
      real    s
      integer i,io,j,jo,k
c
      if     (halo_type.lt.10) then
        s =  1.0  !scalar
      else
        s = -1.0  !vector
      endif
c
      do k= 1,ld
        do j= 1,jtdm
          do i= 1,itdm
            ahalo(i,j,k) = aorig(i,j,k)
          enddo
        enddo
        do j= 1-nbdy,0
          do i= 1,itdm
            ahalo(i,j,k) = 0.0  !southern boundary is closed
          enddo
        enddo
c
        do j= jtdm+1,jtdm+nbdy
          if     (halo_type.eq.1 .or. halo_type.eq.11) then !p-grid
            jo = jtdm-1-(j-jtdm)
            do i= 1,itdm
              io = itdm-mod(i-1,itdm)
              ahalo(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          elseif (halo_type.eq.2 .or. halo_type.eq.12) then !q-grid
            jo = jtdm-(j-jtdm)
            do i= 1,itdm
              io = mod(itdm-(i-1),itdm)+1
              ahalo(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          elseif (halo_type.eq.3 .or. halo_type.eq.13) then !u-grid
            jo = jtdm-1-(j-jtdm)
            do i= 1,itdm
              io = mod(itdm-(i-1),itdm)+1
              ahalo(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          elseif (halo_type.eq.4 .or. halo_type.eq.14) then !v-grid
            jo = jtdm-(j-jtdm)
            do i= 1,itdm
              io = itdm-mod(i-1,itdm)
              ahalo(i,j,k) = s*aorig(io,jo,k)
            enddo !i
          endif !halo_type
        enddo !j
c
        do j= 1-nbdy,jtdm+nbdy
          do i= 1-nbdy,0
            io = itdm+i
            ahalo(i,j,k) = ahalo(io,j,k)  !periodic
          enddo
          do i= itdm+1,itdm+nbdy
            io = i-itdm
            ahalo(i,j,k) = ahalo(io,j,k)  !periodic
          enddo
        enddo
      enddo !k
      return
      end
      subroutine yytile(atile,l1,ld,mh,nh, ht, ahalo)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer l1,ld,mh,nh,ht
      real    atile(1-nbdy:idm +nbdy,1-nbdy:jdm +nbdy,ld)
      real    ahalo(1-nbdy:itdm+nbdy,1-nbdy:jtdm+nbdy)
c
c     check that atile's halo is up to date.
c     based on ahalo being correct in interior and in the halo.
c
      real    o
      integer i,io,j,jo,k,mn
      integer kbad,ksea,ke,kn,ks,kw
c
      if     (mnproc.eq.1) then
        write(lp,'(a,4i5,i3)') 'call xctilr - l1,ld,mh,nh,ht = ',
     &                                        l1,ld,mh,nh,ht
      endif
      call xcsync(flush_lp)
c
      do mn= 1,ijpr
      if     (mn.eq.mnproc) then
c
      kbad = 0
      do k= l1,ld
        o = (k-1)*10000
        do i= 1,ii
          io = i0+i
          do j= 1,jj
            jo = j0+j
            if     (ahalo(io,jo).eq.0.0) then
              if     (atile(i,j,k).ne.0.0) then
                kbad = kbad + 1
              endif
            elseif (ahalo(io,jo).gt.0.0 .and.
     &              atile(i,j,k).ne.ahalo(io,jo)+o) then
              kbad = kbad + 1
            elseif (ahalo(io,jo).lt.0.0 .and.
     &              atile(i,j,k).ne.ahalo(io,jo)-o) then
              kbad = kbad + 1
            endif
          enddo !j
        enddo !i
      enddo !k
      if     (kbad.ne.0) then
        write(6,'(a,2i3,a)') 'mp,np =',mproc,nproc,
     &                       '   tile interior incorrect'
      endif
c
      kbad = 0
      do k= l1,ld
        o  = (k-1)*10000
        ke = 0
        kw = 0
        do i= 1,mh
          io = i0+(1-i)
          do j= 1-nh,jj+nh
            jo = j0+j
            if     (ahalo(io,jo).eq.0.0) then
              if     (atile(1-i,j,k).ne.0.0) then
                kw = kw + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,2f12.2)') 
     &              'at.W = ',1-i,j,k,atile(1-i,j,k),0.0
                  write(lp,'(a,3i5,2f12.2)') 
     &              'ao.W = ',io,jo,k,0.0,atile(1-i,j,k)
                endif
              endif
            elseif (ahalo(io,jo).gt.0.0 .and.
     &              atile(1-i,j,k).ne.ahalo(io,jo)+o) then
              kw = kw + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.W = ',1-i,j,k,atile(1-i,j,k),ahalo(io,jo)+o
                write(lp,'(a,3i5,2f12.2)') 
     &            'ao.W = ',io,jo,k,ahalo(io,jo)+o,atile(1-i,j,k)
              endif
            elseif (ahalo(io,jo).lt.0.0 .and.
     &              atile(1-i,j,k).ne.ahalo(io,jo)-o) then
              kw = kw + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.W = ',1-i,j,k,atile(1-i,j,k),ahalo(io,jo)-o
                write(lp,'(a,3i5,2f12.2)') 
     &            'ao.W = ',io,jo,k,ahalo(io,jo)-o,atile(1-i,j,k)
              endif
            endif
          enddo !j
          io = i0+(ii+i)
          do j= 1-nh,jj+nh
            jo = j0+j
            if     (ahalo(io,jo).eq.0.0) then
              if     (atile(ii+i,j,k).ne.0.0) then
                ke = ke + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,2f12.2)') 
     &              'at.E = ',ii+i,j,k,atile(ii+i,j,k),0.0
                  write(lp,'(a,3i5,2f12.2)') 
     &              'ao.E = ',io,jo,k,0.0,atile(ii+i,j,k)
                endif
              endif
            elseif (ahalo(io,jo).gt.0.0 .and.
     &              atile(ii+i,j,k).ne.ahalo(io,jo)+o) then
              ke = ke + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.E = ',ii+i,j,k,atile(ii+i,j,k),ahalo(io,jo)+o
                write(lp,'(a,3i5,2f12.2)') 
     &            'ao.E = ',io,jo,k,ahalo(io,jo)+o,atile(ii+i,j,k)
              endif
            elseif (ahalo(io,jo).lt.0.0 .and.
     &              atile(ii+i,j,k).ne.ahalo(io,jo)-o) then
              ke = ke + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.E = ',ii+i,j,k,atile(ii+i,j,k),ahalo(io,jo)-o
                write(lp,'(a,3i5,2f12.2)') 
     &            'ao.E = ',io,jo,k,ahalo(io,jo)-o,atile(ii+i,j,k)
              endif
            endif
          enddo !j
        enddo !i
c
        kn = 0
        ks = 0
        do j= 1,nh
          jo = j0+(1-j)
          do i= 1-mh,ii+mh
            io = i0+i
            if     (ahalo(io,jo).eq.0.0) then
              if     (atile(i,1-j,k).ne.0.0) then
                ks = ks + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,2f12.2)') 
     &              'at.S = ',i,1-j,k,atile(i,1-j,k),0.0
                  write(lp,'(a,3i5,2f12.2)') 
     &              'ao.S = ',io,jo,k,0.0,atile(i,1-j,k)
                endif
              endif
            elseif (ahalo(io,jo).gt.0.0 .and.
     &              atile(i,1-j,k).ne.ahalo(io,jo)+o) then
              ks = ks + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.S = ',i,1-j,k,atile(i,1-j,k),ahalo(io,jo)+o
                write(lp,'(a,3i5,2f12.2)') 
     &            'ao.S = ',io,jo,k,ahalo(io,jo)+o,atile(i,1-j,k)
              endif
            elseif (ahalo(io,jo).lt.0.0 .and.
     &              atile(i,1-j,k).ne.ahalo(io,jo)-o) then
              ks = ks + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.S = ',i,1-j,k,atile(i,1-j,k),ahalo(io,jo)-o
                write(lp,'(a,3i5,2f12.2)') 
     &            'ao.S = ',io,jo,k,ahalo(io,jo)-o,atile(i,1-j,k)
              endif
            endif
          enddo
          jo = j0+jj+j
          do i= 1-mh,ii+mh
            io = i0+i
            if     (ahalo(io,jo).eq.0.0) then
              if     (atile(i,jj+j,k).ne.0.0) then
                kn = kn + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,2f12.2)') 
     &              'at.N = ',i,jj+j,k,atile(i,jj+j,k),0.0
                  write(lp,'(a,3i5,2f12.2)')
     &              'ao.N = ',io,jo,k,0.0,atile(i,jj+j,k)
                endif
              endif
            elseif (ahalo(io,jo).gt.0.0 .and.
     &              atile(i,jj+j,k).ne.ahalo(io,jo)+o) then
              kn = kn + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.N = ',i,jj+j,k,atile(i,jj+j,k),ahalo(io,jo)+o
                write(lp,'(a,3i5,2f12.2)')
     &            'ao.N = ',io,jo,k,ahalo(io,jo)+o,atile(i,jj+j,k)
              endif
            elseif (ahalo(io,jo).lt.0.0 .and.
     &              atile(i,jj+j,k).ne.ahalo(io,jo)-o) then
              kn = kn + 1
              if     (k.eq.l1) then
                write(lp,'(a,3i5,2f12.2)') 
     &            'at.N = ',i,jj+j,k,atile(i,jj+j,k),ahalo(io,jo)-o
                write(lp,'(a,3i5,2f12.2)')
     &            'ao.N = ',io,jo,k,ahalo(io,jo)-o,atile(i,jj+j,k)
              endif
            endif
          enddo
        enddo
        if     (kn+ks+ke+kw.ne.0) then
          kbad = kbad + 1
          write(lp,6000) mproc,nproc,k,ks,kn,kw,ke
        endif
      enddo
c
      endif  ! mn.eq.mnproc
      call xcsync(flush_lp)
      enddo  ! mn=1,ijpr
c
      do mn= 1,ijpr
        if     (mn.eq.mnproc) then
          if     (kbad.eq.0) then
            write(lp,6100) mproc,nproc,ld-l1+1
          else
            write(lp,6150) mproc,nproc,kbad,ld-l1+1
          endif
        endif
        call xcsync(flush_lp)
      enddo
      if     (mnproc.eq.1) then
        write(lp,*)
      endif
      call xcsync(flush_lp)
      return
 6000 format('mp,np =',2i3,'   k =',i3,'   ks,kn,kw,ke = ',4i4)
 6100 format('mp,np =',2i3,'   halo correct for all',i3,' levels')
 6150 format('mp,np =',2i3,'   halo incorrect for',i3,' of',
     +   i3,' levels')
      end
      subroutine zztile(atile,l1,ld)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer l1,ld
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
c
c     set the tiles halo to -1.
c
      integer i,j,k
c
      do k= l1,ld
        do j= 1,jj
          do i= 1,nbdy
            atile(ii+i,j,k) = -1.0
            atile( 1-i,j,k) = -1.0
          enddo !i
        enddo !j
        do j= 1,nbdy
          do i= 1-nbdy,ii+nbdy
            atile(i,jj+j,k) = -1.0
            atile(i, 1-j,k) = -1.0
          enddo !i
        enddo !j
      enddo !k
      return
      end
