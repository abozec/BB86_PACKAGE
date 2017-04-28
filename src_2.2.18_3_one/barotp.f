      subroutine barotp(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- micom version 2.8
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations.
c ---   on entry: -n- is time t-dt, -m- is time t
c ---   on exit:                    -m- is time t, -n- is time t+dt
c ---   time level 3 is only used internally (n and m are always 1 or 2).
c
c --- LeapFrog version based on:
c ---   Y. Morel, Baraille, R., Pichon A. (2007) "Time splitting and
c ---   linear stability of the slow part of the barotropic component", 
c ---   Ocean Modeling (submitted)
c --- ------------------------------------------------------------------------
c
      logical    lpipe_barotp
      parameter (lpipe_barotp=.false.)
      logical    ldebug_barotp
      parameter (ldebug_barotp=.false.)
c
      real    q,pbudel,pbvdel,utndcy,vtndcy
      real*8  sump
      integer i,j,l,lll,ml,nl,mn,lstep1,mbdy
c
      mbdy = 6
c
      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_uv)
      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_vv)
c
      if     (lpipe .and. lpipe_barotp) then
c ---   compare two model runs.
        call pipe_compare_sym2(utotn, iu,'barotp:utotn',
     &                         vtotn, iv,'barotp:vtotn')
        call pipe_compare_sym1(pvtrop,iq,'barotp:pvtrp')
      endif
c
c --- explicit time integration of barotropic flow (forward-backward scheme)
c --- in order to combine forward-backward scheme with leapfrog treatment of
c --- coriolis term, v-eqn must be solved before u-eqn every other time step
c
      if     (btrlfr .and. delt1.ne.baclin) then  !not on very 1st time step
C ---   start at time level t-dt and go to t+dt.
        lstep1 = lstep + lstep  !more stable, but also more expensive
      else
C ---   start at time level t    and go to t+dt.
        lstep1 = lstep          !original, less stable, method
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
         do j=1,jj
           do i=1,ii
             pbavg(i,j,n) = pbavg(i,j,m)
             ubavg(i,j,n) = ubavg(i,j,m)
             vbavg(i,j,n) = vbavg(i,j,m)
           enddo !i
         enddo !j
      endif !btrlfr
c
      do 840 lll=1,lstep1,2
c
      call xctilr(pbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_ps)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_vv)
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn',
     &    vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm',
     &    vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
c
c --- odd minor time step.
c
      ml=n
      nl=3
c
c --- continuity equation
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,pbudel,pbvdel)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j))
     &               -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1))
     &               -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))
            pbavg(i,j,nl)=
     &        ((1.-wbaro)*pbavg(i,j,ml)+
     &             wbaro *pbavg(i,j,nl) )-
     &         (1.+wbaro)*dlt*(pbudel + pbvdel)*scp2i(i,j)
          enddo
        enddo
      enddo
c
      mn=ml
c
c --- u momentum equation, 1st
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)+
     &       ((vbavg(i  ,j,  mn)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+
     &        (vbavg(i-1,j,  mn)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
c
            ubavg(i,j,nl)=
     &        ((1.-wbaro)*ubavg(i,j,ml)+
     &             wbaro *ubavg(i,j,nl))+
     &         (1.+wbaro)*dlt*(utndcy+utotn(i,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,5f7.3)')
*    &          nstep,i+i0,j+j0,lll,
*    &          'u_old,u_new,p_grad,corio,u_star =',
*    &          ubavg(i,j,ml),ubavg(i,j,nl),
*    &           -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
*    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
*    &          *(pvtrop(i,j)+pvtrop(i,j+1))
*    &          *.125 * dlt,utotn(i,j) * dlt
*           endif
          enddo
        enddo
      enddo
c
      mn = nl
c
c --- v momentum equation, 2nd
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)-
     &       ((ubavg(i,  j  ,mn)*depthu(i,  j  )
     &        +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+
     &        (ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &        +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
c
            vbavg(i,j,nl)=
     &        ((1.-wbaro)*vbavg(i,j,ml)+
     &             wbaro *vbavg(i,j,nl))+
     &         (1.+wbaro)*dlt*(vtndcy+vtotn(i,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,5f7.3)')
*    &          nstep,i+i0,j+j0,lll,
*    &          'v_old,v_new,p_grad,corio,v_star =',
*    &          vbavg(i,j,ml),vbavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
*    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
*    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
*    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
*    &          *(pvtrop(i,j)+pvtrop(i+1,j))
*    &          *.125 * dlt, vtotn(i,j) * dlt
*           endif
          enddo
        enddo
      enddo
c
      if     (ldebug_barotp) then
        call xcsync(flush_lp)
      endif
c
      if     (lbflag.eq.1) then
        call latbdp(nl)
      elseif (lbflag.eq.2) then
        call latbdt(nl,lll)
      elseif (lbflag.eq.3) then
        call latbdf(nl,lll)
      endif
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn',
     &    vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm',
     &    vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
c
c --- even minor time step.
c
      ml=3
      nl=n
c
c --- continuity equation
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,pbudel,pbvdel)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j))
     &               -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1))
     &               -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))
            pbavg(i,j,nl)=
     &        ((1.-wbaro)*pbavg(i,j,ml)+
     &             wbaro *pbavg(i,j,nl) )-
     &         (1.+wbaro)*dlt*(pbudel + pbvdel)*scp2i(i,j)
          enddo
        enddo
      enddo
c
      mn=ml
c
c --- v momentum equation, 1st
c
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)-
     &       ((ubavg(i,  j  ,mn)*depthu(i,  j  )
     &        +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+
     &        (ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &        +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
c
            vbavg(i,j,nl)=
     &        ((1.-wbaro)*vbavg(i,j,ml)+
     &             wbaro *vbavg(i,j,nl))+
     &         (1.+wbaro)*dlt*(vtndcy+vtotn(i,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,5f7.3)')
*    &          nstep,i+i0,j+j0,lll+1,
*    &          'v_old,v_new,p_grad,corio,v_star =',
*    &          vbavg(i,j,ml),vbavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
*    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
*    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
*    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
*    &          *(pvtrop(i,j)+pvtrop(i+1,j))
*    &          *.125 * dlt, vtotn(i,j) * dlt
*           endif
          enddo
        enddo
      enddo
c
      mn=nl
c
c --- u momentum equation, 2nd
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)+
     &       ((vbavg(i  ,j,  mn)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+
     &        (vbavg(i-1,j,  mn)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
c
            ubavg(i,j,nl)=
     &        ((1.-wbaro)*ubavg(i,j,ml)+
     &             wbaro *ubavg(i,j,nl))+
     &         (1.+wbaro)*dlt*(utndcy+utotn(i,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,5f7.3)')
*    &          nstep,i+i0,j+j0,lll+1,
*    &          'u_old,u_new,p_grad,corio,u_star =',
*    &          ubavg(i,j,ml),ubavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
*    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
*    &          *(pvtrop(i,j)+pvtrop(i,j+1))
*    &          *.125 * dlt,utotn(i,j) * dlt
*           endif
          enddo
        enddo
      enddo
c
      if     (ldebug_barotp) then
        call xcsync(flush_lp)
      endif
c
      if     (lbflag.eq.1) then
        call latbdp(nl)
      elseif (lbflag.eq.2) then
        call latbdt(nl,lll+1)
      elseif (lbflag.eq.3) then
        call latbdf(nl,lll+1)
      endif
c
 840  continue  ! lll=1,lstep1,2
c
      if     (lbflag.eq.1) then
c
c ---   correct mean height.
c ---   this should not be required - so there may be a bug in the bc.
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,sump)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j) = pbavg(i,j,nl)*scp2(i,j)
            enddo
          enddo
        enddo
        call xcsum(sump, util1,ip)
        q = sump/area
c
c ---   rhs: pbavg
c ---   lhs: pbavg
c
        margin = 0
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              pbavg(i,j,1) = pbavg(i,j,1) - q
              pbavg(i,j,2) = pbavg(i,j,2) - q
              pbavg(i,j,3) = pbavg(i,j,3) - q
            enddo
          enddo
        enddo
      endif
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,1), ip,'barotp:pbav1')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,2), ip,'barotp:pbav2')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,3), ip,'barotp:pbav3')
      endif
c
      return
      end subroutine barotp
c
c
c> Revision history:
c>
c> Mar. 1995 - changed vertical velocity averaging interval from 10 cm to 1 m
c>             (loops 33,35)
c> Mar. 1995 - changed order of loop nesting in loop 842
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - transferred loops preceding loop 840 to momeq2.f
c> Jan. 2000 - added latbdp for lateral boundary ports
c> Aug. 2001 - two barotropic time steps per loop, for halo efficiency
c> Nov. 2006 - added lbflag==3 (latbdf) and thref_bt (mod_tides)
c> Nov. 2006 - removed thref_bt (and mod_tides)
c> Apr. 2007 - added btrlfr: leapfrog time step; see also momtum
