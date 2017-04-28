      subroutine psmooth(a,margin_smooth)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer margin_smooth
      real    a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth data carried at -p- points.
c
c --- see also psmooth_ice and psmooth_dif.
c
c --- util1 used as workspace, so array a can't be util1
c
      integer i,ismth,j,jsmth,msmth
      real    qc,sh
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      msmth = min(margin_smooth,nbdy-1)
c
      if     (margin.lt.msmth+1) then
c ---   update the halo
        call xctilr(a,1,1, msmth+1,msmth+1, halo_ps)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC)
      do j=0-msmth,jj+msmth+1
        do i=0-msmth,ii+msmth+1
          util1(i,j) = a(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,i,sh,jsmth,ismth)
!$OMP&         SCHEDULE(STATIC)
      do j=1-msmth,jj+msmth
        do i=1-msmth,ii+msmth
          if     (ip(i,j).eq.1) then
            sh = 0.0
            do jsmth= -1,1
              do ismth= -1,1
                if     (ip(i+ismth,j+jsmth).eq.1) then
                  sh = sh + c(ismth,jsmth)*util1(i+ismth,j+jsmth)
                else
                  sh = sh + c(ismth,jsmth)*util1(i,      j)
                endif
              enddo
            enddo
            a(i,j) = sh*qc
          endif  !ip.eq.1
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine psmooth

      subroutine psmooth_new(a,b,margin_smooth)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer margin_smooth
      real    a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &        b(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth data carried at -p- points.
c --- input in a, output in b.
c
c --- see also psmooth.
c
c --- a and b must not be the same array.
c
      integer i,ismth,j,jsmth,msmth
      real    qc,sh
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      msmth = min(margin_smooth,nbdy-1)
c
      if     (margin.lt.msmth+1) then
c ---   update the halo
        call xctilr(a,1,1, msmth+1,msmth+1, halo_ps)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i,sh,jsmth,ismth)
!$OMP&         SCHEDULE(STATIC)
      do j=1-msmth,jj+msmth
        do i=1-msmth,ii+msmth
          if     (ip(i,j).eq.1) then
            sh = 0.0
            do jsmth= -1,1
              do ismth= -1,1
                if     (ip(i+ismth,j+jsmth).eq.1) then
                  sh = sh + c(ismth,jsmth)*a(i+ismth,j+jsmth)
                else
                  sh = sh + c(ismth,jsmth)*a(i,      j)
                endif
              enddo
            enddo
            b(i,j) = sh*qc
          endif  !ip.eq.1
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine psmooth_new

      subroutine psmooth_ice(a,margin_smooth)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer margin_smooth
      real    a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth data carried at -p- points.
c --- it also smooths covice=1.0 and covice=0.0 regions separately,
c --- leaving areas with fractional covice untouched.  note that
c --- covice must be valid in the halo out to margin_smooth.
c
c --- see also psmooth and psmooth_dif
c
c --- util1 used as workspace, so array a can't be util1 or covice
c
      integer i,ismth,j,jsmth,msmth
      real    qc,sh,ci
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      msmth = min(margin_smooth,nbdy-1)
c
      if     (margin.lt.msmth+1) then
c ---   update the halo
        call xctilr(a,1,1, msmth+1,msmth+1, halo_ps)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC)
      do j=0-msmth,jj+msmth+1
        do i=0-msmth,ii+msmth+1
          util1(i,j) = a(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,i,sh,ci,jsmth,ismth)
!$OMP&         SCHEDULE(STATIC)
      do j=1-msmth,jj+msmth
        do i=1-msmth,ii+msmth
          if     (ip(i,j).eq.1) then
            ci = covice(i,j)
            if     (ci.eq.0.0 .or.
     &              ci.eq.1.0     ) then  !full sea or full ice
              sh = 0.0
              do jsmth= -1,1
                do ismth= -1,1
                  if     (    ip(i+ismth,j+jsmth).eq.1 .and.
     &                    covice(i+ismth,j+jsmth).eq.ci     ) then
                    sh = sh + c(ismth,jsmth)*util1(i+ismth,j+jsmth)
                  else
                    sh = sh + c(ismth,jsmth)*util1(i,      j)
                  endif
                enddo
              enddo
              a(i,j) = sh*qc
            endif  !full sea or full ice
          endif  !ip.eq.1
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine psmooth_ice

      subroutine psmooth_dif(a,aklist,k,margin_smooth)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer k,margin_smooth
      real    a(     1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &        aklist(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth vcty carried at -p- points.
c --- it return the maximum of the original and smoothed value and
c --- ignores locations where k > aklist(i,j).
c
c --- see also psmooth and psmooth_ice.
c
c --- util1 used as workspace, so array a can't be util1
c --- assumes that aklist's halo is valid out to msmth+1.
c
      integer i,ismth,j,jsmth,msmth
      real    qc,sh
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      msmth = min(margin_smooth,nbdy-1)
c
      if     (margin.lt.msmth+1) then
c ---   update the halo
        call xctilr(a,1,1, msmth+1,msmth+1, halo_ps)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC)
      do j=0-msmth,jj+msmth+1
        do i=0-msmth,ii+msmth+1
          util1(i,j) = a(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,i,sh,jsmth,ismth)
!$OMP&         SCHEDULE(STATIC)
      do j=1-msmth,jj+msmth
        do i=1-msmth,ii+msmth
          if     (    ip(i,j).eq.1 .and.
     &            aklist(i,j).ge.k      ) then
            sh = 0.0
            do jsmth= -1,1
              do ismth= -1,1
                if     (    ip(i+ismth,j+jsmth).eq.1 .and.
     &                  aklist(i+ismth,j+jsmth).ge.k      ) then
                  sh = sh + c(ismth,jsmth)*util1(i+ismth,j+jsmth)
                else
                  sh = sh + c(ismth,jsmth)*util1(i,      j)
                endif
              enddo
            enddo
            a(i,j) = max( a(i,j), sh*qc )
          endif  !ip.eq.1
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine psmooth_dif
