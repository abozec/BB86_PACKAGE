      program testxc
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      logical, parameter :: lregion = .true.   ! input regional.depth?
c
c     test xctilr, for non-arctic tiles only (itype ignored).
c
      integer i,j,k,ksea,nrecl
      real*4  depth(itdm,jtdm)
      real    aorig(itdm,jtdm,kdm)
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
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
            do k= 1,kdm
              aorig(i,j,k) = i + (j-1)*100 + (k-1)*10000
            enddo
          else
            do k= 1,kdm
              aorig(i,j,k) = 0.0
            enddo
          endif
        enddo
      enddo
c
      do k= 1,kdm
        do j= 1,jj
          do i= 1,ii
            atile(i,j,k) = aorig(i+i0,j+j0,k)
          enddo
        enddo
      enddo
      if     (mnproc.eq.1) then
      write(lp,*) 'itdm,jtdm = ',itdm,jtdm
      write(lp,*) 'ii,  jj   = ',ii,  jj  
      ksea = count( aorig(:,:,1).ne.0.0 )
      write(lp,*) 'sea, land = ',ksea,itdm*jtdm-ksea
      write(lp,*) 'aorig1 = ',aorig(1,1,1),aorig(1,2,1),
     +                        aorig(2,1,1),aorig(2,2,1)
      write(lp,*) 'atile1 = ',atile(1,1,1),atile(1,2,1),
     +                        atile(2,1,1),atile(2,2,1)
      write(lp,*) 'aorig9 = ',aorig(1,1,9),aorig(1,2,9),
     +                        aorig(2,1,9),aorig(2,2,9)
      write(lp,*) 'atile9 = ',atile(1,1,9),atile(1,2,9),
     +                        atile(2,1,9),atile(2,2,9)
      write(lp,*)
      endif
      call xcsync(flush_lp)
c
c     test.
c
      call xctilr(atile,1,1, 1,1, halo_ps)
      call yytile(atile,1,1, 1,1, aorig)
c
      call xctilr(atile,1,1, 0,1, halo_ps)
      call yytile(atile,1,1, 0,1, aorig)
c
      call xctilr(atile,1,1, 1,0, halo_ps)
      call yytile(atile,1,1, 1,0, aorig)
c
      call xctilr(atile,3,9, nbdy,nbdy, halo_ps)
      call yytile(atile,3,9, nbdy,nbdy, aorig)
c
      call xctilr(atile,1,kdm, nbdy,nbdy, halo_ps)
      call yytile(atile,1,kdm, nbdy,nbdy, aorig)
c
      call xctilr(atile,1,kdm, 2,   nbdy, halo_ps)
      call yytile(atile,1,kdm, 2,   nbdy, aorig)
c
      call xctilr(atile,1,kdm, nbdy,   2, halo_ps)
      call yytile(atile,1,kdm, nbdy,   2, aorig)
c
      call xcstop('(normal)')
             stop '(normal)'
      end
      subroutine yytile(atile,l1,ld,mh,nh, aorig)
      use mod_xc    ! HYCOM communication interface
      implicit none
c
      integer l1,ld,mh,nh
      real    atile(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
      real    aorig(itdm,jtdm,ld)
c
c     check that atile's halo is up to date.
c
      integer i,io,j,jo,k,mn
      integer kbad,ksea,ke,kn,ks,kw
c
      if     (mnproc.eq.1) then
        write(lp,'(a,4i5)') 'call xctilr - l1,ld,mh,nh = ',
     &                                     l1,ld,mh,nh
      endif
      call xcsync(flush_lp)
c
      do mn= 1,ijpr
      if     (mn.eq.mnproc) then
c
      kbad = 0
      do k= l1,ld
        ke = 0
        kw = 0
        do i= 1,mh
          if     (nreg.eq.0) then
            io =     i0+1-i
          else
            io = mod(i0+1-i+itdm-1,itdm)+1
          endif
          do j= 1-nh,jj+nh
            jo = j0+j
            if     ((1.le.jo .and. jo.le.jtdm) .and.
     &              (1.le.io .and. io.le.itdm)      ) then
              if     (atile(1-i,j,k).ne.aorig(io,jo,k)) then
                kw = kw + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.W = ',1-i,j,k,atile(1-i,j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'ao.W = ',io,jo,k,aorig(io,jo,k)
                endif
              endif
            else
              if     (atile(1-i,j,k).ne.vland) then
                kw = kw + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.W = ',1-i,j,k,atile(1-i,j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'land = ',io,jo,k,vland
                endif
              endif
            endif
          enddo
          if     (nreg.eq.0) then
            io =     i0+ii+i
          else
            io = mod(i0+ii+i-1,itdm)+1
          endif
          do j= 1-nh,jj+nh
            jo = j0+j
            if     ((1.le.jo .and. jo.le.jtdm) .and.
     &              (1.le.io .and. io.le.itdm)      ) then
              if     (atile(ii+i,j,k).ne.aorig(io,jo,k)) then
                ke = ke + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.E = ',ii+i,j,k,atile(ii+i,j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'ao.E = ',io,jo,k,aorig(io,jo,k)
                endif
              endif
            else
              if     (atile(ii+i,j,k).ne.vland) then
                ke = ke + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.E = ',ii+i,j,k,atile(ii+i,j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'land = ',io,jo,k,vland
                endif
              endif
            endif
          enddo
        enddo
c
        kn = 0
        ks = 0
        do j= 1,nh
          jo = j0+1-j
          do i= 1-mh,ii+mh
            if     (nreg.eq.0) then
              io =     i0+i
            else
              io = mod(i0+i+itdm-1,itdm)+1
            endif
            if     ((1.le.jo .and. jo.le.jtdm) .and.
     &              (1.le.io .and. io.le.itdm)      ) then
              if     (atile(i,1-j,k).ne.aorig(io,jo,k)) then
                ks = ks + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.S = ',i,1-j,k,atile(i,1-j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'ao.S = ',io,jo,k,aorig(io,jo,k)
                endif
              endif
            else
              if     (atile(i,1-j,k).ne.vland) then
                ks = ks + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.S = ',i,1-j,k,atile(i,1-j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'land = ',io,jo,k,vland
                endif
              endif
            endif
          enddo
          jo = j0+jj+j
          do i= 1-mh,ii+mh
            if     (nreg.eq.0) then
              io =     i0+i
            else
              io = mod(i0+i+itdm-1,itdm)+1
            endif
            if     ((1.le.jo .and. jo.le.jtdm) .and.
     &              (1.le.io .and. io.le.itdm)      ) then
              if     (atile(i,jj+j,k).ne.aorig(io,jo,k)) then
                kn = kn + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.N = ',i,jj+j,k,atile(i,jj+j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'ao.N = ',io,jo,k,aorig(io,jo,k)
                endif
              endif
            else
              if     (atile(i,jj+j,k).ne.vland) then
                kn = kn + 1
                if     (k.eq.l1) then
                  write(lp,'(a,3i5,f12.2)') 
     &              'at.N = ',i,jj+j,k,atile(i,jj+j,k)
                  write(lp,'(a,3i5,f12.2)') 
     &              'land = ',io,jo,k,vland
                endif
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
      call xcsync(flush_lp)
      endif
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
