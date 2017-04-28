      subroutine overtn(dtime,dyear)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      real*8  dtime,dyear
c
c --- diagnose meridional heat flux in basin model
c
      real*8     dsmall
      parameter (dsmall=0.001d0)
c
      integer i,j,k,l,noo
c
      logical lfirst
      save    lfirst
      data    lfirst / .true. /
c
      real*8 anum(jtdm),heatf(jtdm),heatfl(jtdm),hfxzon(jtdm),
     &       zonavt(jtdm,kdm),zonavf(jtdm,kdm)
      save   anum,heatfl,hfxzon,zonavt,zonavf
c
c --- integrate meridional heat fluxes vertically and in zonal direction
      do j=1,jtdm
        hfxzon(j)=0.
        heatfl(j)=0.
      enddo
c
      do k=1,kk
        do j=1,jj
          do i=1,ii
            if     (iv(i,j).ne.0) then
              util1(i,j) = (temp(i,j,k,1)+temp(i,j-1,k,1))
              util2(i,j) = vflx(i,j,k)
              util3(i,j) = (temp(i,j,k,1)+temp(i,j-1,k,1))*vflx(i,j,k)
            endif
          enddo
        enddo
        call xcsumj(zonavt(1,k), util1,iv)
        call xcsumj(zonavf(1,k), util2,iv)
        call xcsumj(heatf,       util3,iv)
        if     (lfirst) then
          util4 = 1.0
          call xcsumj(anum, util4,iv)
        endif
        if     (mnproc.eq.1) then
          do j=1,jtdm
            if     (anum(j).ne.0.0) then
              heatfl(j)=heatfl(j)+heatf(j)
              hfxzon(j)=hfxzon(j)+zonavt(j,k)*zonavf(j,k)/anum(j)
            endif
          enddo
        endif
      enddo
      if     (mnproc.eq.1) then
        do j= 1,jtdm
          hfxzon(j)=.5*hfxzon(j)*spcifh/g * 1.e-15
          heatfl(j)=.5*heatfl(j)*spcifh/g * 1.e-15
        enddo
cdiag      print 999, nstep,vflx(31,11,11)
cdiag 999  format(' overtn - nstep,vflx=',i10,d20.12)
c       
c ---   save everything in a special file
        noo=26
c
        if     (lfirst) then
          open (unit=noo,file=flnmovr,status='new',form='formatted')
        endif
c
        write (noo,'(a,f10.2,i6,f7.2)')
     &      ' time,year,day =',dtime,int((dtime+dsmall)/dyear),
     &                                mod(dtime+dsmall, dyear)
        write (noo,'(a/(11f7.3))')
     &     ' northward heat flux (petawatts):', (heatfl(j),j=1,jtdm-1)
        write (noo,'(a/(11f7.3))')
     &     ' meridional overturning component:',(hfxzon(j),j=1,jtdm-1)
        call flush(noo)
c
        write (lp, '(a/(11f7.3))')
     &     ' northward heat flux (petawatts):', (heatfl(j),j=1,jtdm-1)
        call flush(lp)
      endif
      call xcsync(flush_lp)
c
      lfirst = .false.
c
      return
      end subroutine overtn


