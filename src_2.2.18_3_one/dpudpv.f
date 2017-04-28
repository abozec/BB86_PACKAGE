      subroutine dpudpv(dpu,dpv, p,depthu,depthv, margin_dpudpv)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer, intent(in)    :: margin_dpudpv
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &         intent(out)   :: dpu,dpv
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1),
     &         intent(inout) :: p
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: depthu,depthv
c
c --- -----------------------------------------------------------------
c --- define layer depth at  u,v  points with halo out to margin_dpudpv
c --- -----------------------------------------------------------------
c
      interface
          subroutine dpudpvj(dpu,dpv, p,depthu,depthv, margin_dpudpv, j)
          use mod_xc
          integer, intent(in)    :: margin_dpudpv,j
          real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &             intent(out)   :: dpu,dpv
          real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1),
     &             intent(in)    :: p
          real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &             intent(in)    :: depthu,depthv
          end subroutine dpudpvj
      end interface
c
      integer j
c
      if     (margin_dpudpv.lt.0    .or.
     &        margin_dpudpv.ge.nbdy     ) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i2 /)')
     &      'error: dpudpv called with margin_dpudpv = ',margin_dpudpv
        endif
        call xcstop('dpudpv')
               stop 'dpudpv'
      endif
c
c --- p's halo is valid out to margin, is this far enough?
c
      if     (margin.lt.margin_dpudpv+1) then
        call xctilr(p(1-nbdy,1-nbdy,2),1,kk,
     &              margin_dpudpv+1,margin_dpudpv+1, halo_ps)
      endif
c
c --- using single row routine fixes SGI OpenMP bug.
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin_dpudpv,jj+margin_dpudpv
        call dpudpvj(dpu,dpv, p,depthu,depthv, margin_dpudpv, j)
      enddo
      return
      end
      subroutine dpudpvj(dpu,dpv, p,depthu,depthv, margin_dpudpv, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer, intent(in)    :: margin_dpudpv,j
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &         intent(out)   :: dpu,dpv
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1),
     &         intent(in)    :: p
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: depthu,depthv
c
c --- -----------------------------------------------
c --- define layer depth at  u,v  points,  single row
c --- -----------------------------------------------
c
      integer i,k,l
c
      do k=1,kk
c
        do l=1,isu(j)
          do i=max( 1-margin_dpudpv,ifu(j,l)),
     &         min(ii+margin_dpudpv,ilu(j,l))
            dpu(i,j,k)=max(0.,
     &           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     &           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
        enddo
c
        do l=1,isv(j)
          do i=max( 1-margin_dpudpv,ifv(j,l)),
     &         min(ii+margin_dpudpv,ilv(j,l))
            dpv(i,j,k)=max(0.,
     &           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,j-1,k+1)))-
     &           min(depthv(i,j),.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
        enddo
c
      enddo
      return
      end subroutine dpudpvj
