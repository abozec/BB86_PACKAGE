      subroutine latbdf(n,lll)
      use mod_xc     ! HYCOM communication interface
      use mod_tides  ! HYCOM tides
      implicit none
      include 'common_blocks.h'
c     
      integer n,lll
c     
c     --- apply lateral boundary conditions to   barotropic  flow field
c     
c     --- port flow version:
c     --- NOT similar to the standard 'Browning and Kreiss' MICOM/HYCOM open
c     --- boundary condition. This version uses algorithms based on a 
c     --- 1 invariant Flather boundary condition (setting the gradient
c     --- of the incoming characteristic to zero).
c
c     --- The tangential velocity is not constrained.
c
c     --- see also: latbdp
c     
c     --- the code is as similar as possible to that for the standard case.
c     --- so for example, 'speed' is in fact 1/SQRT(gH) which represents
c     --- c1/g in the notation of (Bleck and Sun, Open boundary conditions
c     --- for MICOM).  The 1/g allows for the use of pressure fields.
c
c     --- Note that East, West, North and South refers to the grid 
c     --- (i.e i,j points) and NOT geographic East, West, North and South
c     
c     --- the first call is made during initialization.
c     
c     --- Iris Lohmann, Carlos Lozano, NCEP, April 2006
c     
      logical, parameter :: ldebug_latbdf=.false.
c
      integer, parameter :: mports=9  !maximum number of ports
c
      integer, parameter :: nchar =120
c     
      logical    lfatal,lfatalp
      integer    i,j,isec,ifrst,ilast,l
      real       aline(nchar),          
     &    dline(itdm+jtdm),xline(itdm+jtdm),
     &    pline(itdm+jtdm),uline(itdm+jtdm)

      real        sum,svspin,fatal
      character*3 char3
c 
      integer nports    
      integer lnport(mports),kdport(mports)
      integer jfport(mports),jlport(mports),
     &        ifport(mports),ilport(mports)
      real    svpnow(mports),svport(mports)

      save lnport
      save svpnow,svport
      save nports,kdport,ifport,ilport,jfport,jlport
c     
      real uportw(jtdm,mports),speedw(jtdm,mports),rspedw(jtdm,mports),
     &     uporte(jtdm,mports),speede(jtdm,mports),rspede(jtdm,mports),
     &     vportn(itdm,mports),speedn(itdm,mports),rspedn(itdm,mports),
     &     vports(itdm,mports),speeds(itdm,mports),rspeds(itdm,mports)
      save uportw,speedw,rspedw,uporte,speede,rspede,
     &     vportn,speedn,rspedn,vports,speeds,rspeds

c     tides stuff

      integer npts_p,kdpt_p(mports),
     &               ifpt_p(mports),ilpt_p(mports),
     &               jfpt_p(mports),jlpt_p(mports),lnpt_p(mports)
      integer npts_v,kdpt_v(mports),
     &               ifpt_v(mports),ilpt_v(mports),
     &               jfpt_v(mports),jlpt_v(mports),lnpt_v(mports)
      integer npts_u,kdpt_u(mports),
     &               ifpt_u(mports),ilpt_u(mports),
     &               jfpt_u(mports),jlpt_u(mports),lnpt_u(mports)

c     the max of itdm and jtdm is ok for any port
c     number of tidal consituents (ncon) from mod_tides
      real z1r_p(mports,max(jtdm,itdm),ncon)
      real z1i_p(mports,max(jtdm,itdm),ncon)
      real z1r_u(mports,max(jtdm,itdm),ncon)
      real z1i_u(mports,max(jtdm,itdm),ncon)
      real z1r_v(mports,max(jtdm,itdm),ncon)
      real z1i_v(mports,max(jtdm,itdm),ncon)
      real tmp1, tmp2

      real tmpr(max(itdm,jtdm),ncon), tmpi(max(itdm,jtdm),ncon)

      real  upred(max(itdm,jtdm)),zpred(max(itdm,jtdm))
      real  udpred(max(itdm,jtdm))
      real  vpred(max(itdm,jtdm))
      real  ulow(max(itdm,jtdm),mports),plow(max(itdm,jtdm),mports)
      real  uu(max(itdm,jtdm)),vv(max(itdm,jtdm))
      integer bnd_init(mports)

      real*8 d_time
      real*8 timermp,frmp
      logical astroflag
      integer jn,in,ic

      save z1r_p, z1i_p
      save z1r_u, z1i_u
      save bnd_init
      save ulow,plow
c     
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
c
c     USER-INPUT: optimization coefficients for the 1 inv algorithm.
      real w_1,w_1c
      save w_1 
      data w_1 / 0.1 / 
c     
      integer lcount
      save    lcount
      data    lcount / 0 /

c     
c     Comment out the next line to run without 
c     boundary tides
      d_time = time_8 + lll*dlt/86400.d0

c     add 15384 for obtaining mjd
      lcount = lcount + 1
c
c     set 1-invariant coefficient
      w_1c=1.0-w_1
c
c     
c     --- the first call just initializes data structures.
c     
      if     (lcount.eq.1) then
c     
        open(unit=uoff+99,file=trim(flnminp)//'ports.input')
c     
c     ---   'nports' = number of boundary port sections.
        call blkini(nports,'nports')
        if     (mnproc.eq.1) then
          write(lp,*)
        endif
        if     (nports.lt.0 .or. nports.gt.mports) then
          if     (mnproc.eq.1) then
            write(lp,*) 
            write(lp,*) 'error in latbdf - illegal nports value'
            if     (nports.gt.mports) then
              write(lp,*) 'increase parameter mports to',nports
            endif
            write(lp,*) 
            call flush(lp)
          endif
          call xcstop('(latbdp)')
          stop '(latbdp)'
        endif
c     
c     ---   read in the ports one at a time
c     
        do l= 1,nports
c     
c     ---     port location is w.r.t. u (EW) or v (NS) grid
c     ---     and identifies the sea at the port
c     ---     the minimum index is 0
c     
c     ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W, 5=S(Fundy))
c     ---     'ifport' = first i-index
c     ---     'ilport' = last  i-index (=ifport for N or S orientation)
c     ---     'jfport' = first j-index
c     ---     'jlport' = last  j-index (=jfport for E or W orientation)
c     ---     'svpnow' = existing port transport in Sv (+ve towards E or S)
c     ---     'svport' = target   port transport in Sv (+ve towards E or S)
c     ---     'lnport' = port length (calculated, not input)
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          call blkinr(svpnow(l),'svpnow','(a6," =",f10.4," Sv")')
          call blkinr(svport(l),'svport','(a6," =",f10.4," Sv")')
          if     (mnproc.eq.1) then
            write(lp,*)
          endif
c     
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
c     
c     ---     sanity check.
c     
          if     (kdport(l).eq.3.or.kdport(l).eq.4) then
            if     (ifport(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'error in latbdp - port direction',
     &               ' and orientation are not consistent'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'error in latbdp - port direction',
     &               ' and orientation are not consistent'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &         jfport(l).gt.jlport(l)     ) then
            if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdp - port',
     &             ' location is not consistent'
              write(lp,*) 
              call flush(lp)
            endif
            call xcstop('(latbdp)')
            stop '(latbdp)'
          endif
        enddo
c     
        close(unit=uoff+99)

c **********************************************

c------ read low frequency input-file

c       initialize the low-frequency u and eta to zero
          
        do i = 1,max(itdm,jtdm)
          do j = 1,mports
            ulow(i,j) =0.0
            plow(i,j) =0.0
          enddo
        enddo

        open(unit=uoff+99,file=trim(flnminp)//'lowfreq.input')

c       the file is read in order west, east, south, north

        do l= 1,nports

          if     (kdport(l).eq.4) then     
c         western port                     
            do j= jfport(l),jlport(l)     
              read(uoff+99,*) ulow(j,l),plow(j,l)
            enddo
          elseif (kdport(l).eq.3) then    
c         eastern port
            do j= jfport(l),jlport(l)
              read(uoff+99,*) ulow(j,l),plow(j,l)
            enddo
          elseif (kdport(l).eq.2) then
c         southern port  
            do i= ifport(l),ilport(l)
              read(uoff+99,*) ulow(i,l),plow(i,l)
            enddo   
          elseif (kdport(l).eq.1) then
c         northern port  
            do i= ifport(l),ilport(l)
              read(uoff+99,*) ulow(i,l),plow(i,l)
            enddo    
          endif  !kdport

        enddo !l=1,nports

        close(uoff+99)

c************************************************
c ****** READ TIDAL CONSTITUENTS ****************

        if (tidflg.ge.1) then


c     initialize the tidal constituents to zero
          
          do i = 1,mports
            do j = 1,max(itdm,jtdm)
              do ic  = 1,ncon
                z1r_p(i,j,ic)  = 0.
                z1i_p(i,j,ic)  = 0.
                z1r_u(i,j,ic)  = 0.
                z1i_u(i,j,ic)  = 0.
                z1r_v(i,j,ic)  = 0.
                z1i_v(i,j,ic)  = 0.
              enddo
            enddo
          enddo


          do j = 1,max(itdm,jtdm)
            uu(j) = 0.
            vv(j) = 0.
            upred(j) = 0.
            vpred(j) = 0.
            zpred(j) = 0.
          enddo

c  --------------------------------------------------------

c     Now the P points

          open(unit=uoff+99,file=trim(flnminp)//'tidalports_p.input')

c     ---   'nports' = number of boundary port sections.
          call blkini(npts_p,'npts_p')
          if     (mnproc.eq.1) then
            write(lp,*)
          endif 

          if     (npts_p.ne.nports) then
            if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error number of ports needs to be same'
              write(lp,*) 
              call flush(lp)
            endif
            call xcstop('(latbdp)')
            stop '(latbdp)'
          endif

   
          
          do l= 1,nports
c     
c     ---     port location is w.r.t. u (EW) or v (NS) grid
c     ---     and identifies the sea at the port
c     ---     the minimum index is 0
c     
c     ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W, 5=Fundy(S))
c     ---     'ifport' = first i-index
c     ---     'ilport' = last  i-index (=ifport for N or S orientation)
c     ---     'jfport' = first j-index
c     ---     'jlport' = last  j-index (=jfport for E or W orientation)
c     ---     'lnport' = port length (calculated, not input)


            call blkini(kdpt_p(l),'kdpt_p')
            
            if     (kdpt_p(l).ne.kdport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the kdport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif
            

            call blkini(ifpt_p(l),'ifpt_p')

            if     (ifpt_p(l).ne.ifport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the ifport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif


            call blkini(ilpt_p(l),'ilpt_p')
            if     (ilpt_p(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the ilport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif


            call blkini(jfpt_p(l),'jfpt_p')

            if     (jfpt_p(l).ne.jfport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the jfport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif



            call blkini(jlpt_p(l),'jlpt_p')
            
            if     (jlpt_p(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the jlport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif
           
            if     (kdport(l).eq.4) then
c     
c     western port
c                 
              do j= jfport(l),jlport(l)        
                do ic = 1, ncon
                  read(uoff+99,*)  tmp1, tmp2
                  z1r_p(l,j,ic) = tmp1
                  z1i_p(l,j,ic) = tmp2
                enddo 
              enddo             
c     
            elseif (kdport(l).eq.3) then
c     
c     eastern port
              
              do j= jfport(l),jlport(l)
                do ic = 1, ncon
                  read(uoff+99,*)  tmp1, tmp2
                  z1r_p(l,j,ic) = tmp1
                  z1i_p(l,j,ic) = tmp2
                enddo 
              enddo
c     
            elseif (kdport(l).eq.1) then
c     
c     northern port
c                  
              do i= ifport(l),ilport(l)
                do ic = 1, ncon
                  read(uoff+99,*)  tmp1, tmp2
                  z1r_p(l,i,ic) = tmp1
                  z1i_p(l,i,ic) = tmp2
                enddo 
              enddo 
c     
            elseif (kdport(l).eq.2) then
c     
c     southern port
c     
              do i= ifport(l),ilport(l)
                do ic = 1, ncon
                   read(uoff+99,*)  tmp1, tmp2
                  z1r_p(l,i,ic) = tmp1
                  z1i_p(l,i,ic) = tmp2
                enddo 
              enddo
              
            endif !kdport=

c     Close the l = 1, nports loop
          enddo !nports

          close(uoff+99)

c -------------------------------------------------------------          
c     Now the normal-velocity points

          
          open(unit=uoff+99,file=trim(flnminp)//'tidalports_v.input')

c     ---   'nports' = number of boundary port sections.
          call blkini(npts_u,'npts_u')
          if     (mnproc.eq.1) then
            write(lp,*)
          endif

          if     (npts_u.ne.nports) then
            if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error number of ports needs to be same'
              write(lp,*) 
              call flush(lp)
            endif
            call xcstop('(latbdp)')
            stop '(latbdp)'
          endif          
          
          do l= 1,nports
c     
c     ---     port location is w.r.t. u (EW) or v (NS) grid
c     ---     and identifies the sea at the port
c     ---     the minimum index is 0
c     
c     ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W, 5=Fundy(S))
c     ---     'ifport' = first i-index
c     ---     'ilport' = last  i-index (=ifport for N or S orientation)
c     ---     'jfport' = first j-index
c     ---     'jlport' = last  j-index (=jfport for E or W orientation)
c     ---     'svpnow' = existing port transport in Sv (+ve towards E or S)
c     ---     'svport' = target   port transport in Sv (+ve towards E or S)
c     ---     'lnport' = port length (calculated, not input)



            call blkini(kdpt_u(l),'kdpt_u')
            
            if     (kdpt_u(l).ne.kdport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the kdport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif
            

            call blkini(ifpt_u(l),'ifpt_u')

            if     (ifpt_u(l).ne.ifport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the ifport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif


            call blkini(ilpt_u(l),'ilpt_u')
            if     (ilpt_u(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the ilport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif


            call blkini(jfpt_u(l),'jfpt_u')

            if     (jfpt_u(l).ne.jfport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the jfport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif


            call blkini(jlpt_u(l),'jlpt_u')
            
            if     (jlpt_u(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
                write(lp,*) 
                write(lp,*) 'Mismatch of the jlport'
                write(lp,*) 
                call flush(lp)
              endif
              call xcstop('(latbdp)')
              stop '(latbdp)'
            endif
            

            if     (kdport(l).eq.4) then
c     
c     western port
c     
              
              do j= jfport(l),jlport(l)  
                do ic = 1, ncon
                  read(uoff+99,*) tmp1, tmp2
                  z1r_u(l,j,ic) = tmp1
                  z1i_u(l,j,ic) = tmp2 
                enddo 
              enddo                   
c     
            elseif (kdport(l).eq.3) then
c     
c     eastern port
              
              do j= jfport(l),jlport(l)
                do ic = 1, ncon
                  read(uoff+99,*) tmp1, tmp2
                  z1r_u(l,j,ic) = tmp1
                  z1i_u(l,j,ic) = tmp2 
                enddo 
              enddo
c     
            elseif (kdport(l).eq.1) then
c     
c     northern port
c     
              do i= ifport(l),ilport(l)
                do ic = 1, ncon
                  read(uoff+99,*) tmp1, tmp2
                  z1r_u(l,i,ic) = tmp1
                  z1i_u(l,i,ic) = tmp2 
                enddo 
              enddo
c     
            elseif (kdport(l).eq.2) then
c     
c     southern port
c                  
              do i= ifport(l),ilport(l)
                do ic = 1, ncon
                  read(uoff+99,*) tmp1, tmp2
                  z1r_u(l,i,ic) = tmp1
                  z1i_u(l,i,ic) = tmp2 
                enddo 
              enddo

            endif !kdport=

c         Close the l = 1, nports loop
          enddo !nports

          close(uoff+99)
          
c*****************END OF READING THE TIDAL CONSTITUENTS
        endif !tidflg.ge.1

c     
c ---   check ports against masks,
c ---   mark the port locations on masks and print them out.
c     
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c     
          if     (kdport(l).eq.4) then
c     
c         western port
c     
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.itdm-2 .or.
     &               j.lt.1 .or. j.gt.jtdm       ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                 j.le.j0 .or. j.gt.j0+jj     ) then
                cycle         ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9 !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0+1,j-j0).ne.1 .or.
     &               iu(i-i0+2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7 !indicate an error
              endif
            enddo
c     
          elseif (kdport(l).eq.3) then
c     
c         eastern port
c     
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.3 .or. i.gt.itdm .or.
     &               j.lt.1 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                 j.le.j0 .or. j.gt.j0+jj     ) then
                cycle         ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9 !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0-1,j-j0).ne.1 .or.
     &               iu(i-i0-2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7 !indicate an error
              endif
            enddo
c     
          elseif (kdport(l).eq.1) then
c     
c         northern port
c     
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm .or.
     &               j.lt.3 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                 j.le.j0 .or. j.gt.j0+jj     ) then
                cycle         ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9 !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0-1).ne.1 .or.
     &               iv(i-i0,j-j0-2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7 !indicate an error
              endif
            enddo
c     
          elseif (kdport(l).eq.2.or.kdport(l).eq.5) then
c     
c         southern port
c     
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm   .or.
     &               j.lt.1 .or. j.gt.jtdm-2     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                 j.le.j0 .or. j.gt.j0+jj     ) then
                cycle         ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9 !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0+1).ne.1 .or.
     &               iv(i-i0,j-j0+2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7 !indicate an error
              endif
            enddo
c     
          endif !kdport=
c     
          if     (lfatalp) then
            write(lp,*) 
            write(lp,*) 'error in latbdp - port ',l,' mislocated',
     &             '  (mnproc = ',mnproc,')'
            write(lp,*) 
            call flush(lp)
          endif
          lfatal = lfatal .or. lfatalp

        enddo !nports
c     
c       local lfatal to global lfatal
c     
        if     (lfatal) then
          fatal = 1.0
        else
          fatal = 0.0
        endif
        call xcmaxr(fatal)
        lfatal = fatal.gt.0.5
c     
c ---   write out  -iu-  and -iv- arrays, if they are not too big
c ---   data are written in strips nchar points wide
        if     (lfatal .or. max(itdm,jtdm).le.2*nchar) then
          util1(1:ii,1:jj) = iu(1:ii,1:jj) ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
              write (lp,'(a,i5,a,i5)')
     &               'iu array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
                write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
            write (lp,*)
          endif
          call xcsync(flush_lp)
c     
          util1(1:ii,1:jj) = iv(1:ii,1:jj) ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
              write (lp,'(a,i5,a,i5)')
     &               'iv array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
                write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
            write (lp,*)
          endif
          call xcsync(flush_lp)
        endif                 ! small region
c     
        if     (lfatal) then
          write(lp,*) 
          write(lp,*) 'error in latbdp - bad port(s)'
          write(lp,*) 
          call flush(lp)
          call xchalt('(latbdp)')
          stop '(latbdp)'
        endif
c     
c ---   restore iu and iv, and zero iuopn and ivopn.
c     
!     $OMP PARALLEL DO PRIVATE(j,i)
!     $OMP&         SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            iu(i,j) = max( iu(i,j), 0 )
            iv(i,j) = max( iv(i,j), 0 )
          enddo
        enddo
!     $OMP PARALLEL DO PRIVATE(j,i)
!     $OMP&         SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            iuopn(i,j) = 0
            ivopn(i,j) = 0
          enddo
        enddo

c     
c  ---  initialize the ports
c     
        do l= 1,nports
          if     (kdport(l).eq.4) then
c     
c         western port
c     
            sum = 0.0
            i = ifport(l)
            j = jfport(l) 
            call xclget(dline(j),lnport(l), depths,i,j,0,1, 0)
            call xclget(xline(j),lnport(l), scuy,  i,  j,0,1, 0)
            do j= jfport(l),jlport(l)
              sum = sum + dline(j)*xline(j)
            enddo
            sum = 1.e6/sum
            do j= jfport(l),jlport(l)
              uportw(j,l) = sum
              speedw(j,l) = sqrt(1.0*thref/(onem*dline(j)))
              rspedw(j,l) = 1.0/speedw(j,l)
              if     (i.ge.i0+ 1-nbdy .and.
     &               i.le.i0+ii+nbdy .and.
     &               j.ge.j0+ 1-nbdy .and.
     &               j.le.j0+jj+nbdy      ) then
                iuopn(i-i0,j-j0) = 1
              endif
            enddo
c     
          elseif (kdport(l).eq.3) then
c     
c         eastern port
c     
            sum = 0.0
            i = ifport(l)-1
            j = jfport(l)
            call xclget(dline(j),lnport(l), depths,i,  j,0,1, 0)
            call xclget(xline(j),lnport(l), scuy,  i+1,j,0,1, 0)
            do j= jfport(l),jlport(l)
              sum = sum + dline(j)*xline(j)
            enddo
            sum = 1.e6/sum
            do j= jfport(l),jlport(l)
              uporte(j,l) = sum
              speede(j,l) = sqrt(1.0*thref/(onem*dline(j)))
              rspede(j,l) = 1.0/speede(j,l)
              if     (i+1.ge.i0+ 1-nbdy .and.
     &               i+1.le.i0+ii+nbdy .and.
     &               j  .ge.j0+ 1-nbdy .and.
     &               j  .le.j0+jj+nbdy      ) then
                iuopn(i-i0+1,j-j0) = 1
              endif
            enddo
c     
          elseif (kdport(l).eq.1) then
c     
c         northern port
c     
            sum = 0.0
            j = jfport(l)-1
            i = ifport(l)
            call xclget(dline(i),lnport(l), depths,i,j,  1,0, 0)
            call xclget(xline(i),lnport(l), scuy,  i,j+1,1,0, 0)
            do i= ifport(l),ilport(l)
              sum = sum + dline(i)*xline(i)
            enddo
            sum = 1.e6/sum
            do i= ifport(l),ilport(l)
              vportn(i,l) = sum
              speedn(i,l) = sqrt(1.0*thref/(onem*dline(i)))
              rspedn(i,l) = 1.0/speedn(i,l)
              if     (i  .ge.i0+ 1-nbdy .and.
     &               i  .le.i0+ii+nbdy .and.
     &               j+1.ge.j0+ 1-nbdy .and.
     &               j+1.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0+1) = 1
              endif
            enddo
c     
          elseif (kdport(l).eq.2.or.kdport(l).eq.5) then
c     
c         southern port
c     
            sum = 0.0
            j = jfport(l)
            i = ifport(l)
            call xclget(dline(i),lnport(l), depths,i,j,1,0, 0)
            call xclget(xline(i),lnport(l), scuy,  i,j,  1,0, 0)
            do i= ifport(l),ilport(l)
              sum = sum + dline(i)*xline(i)
            enddo
            sum = 1.e6/sum
            do i= ifport(l),ilport(l)
              vports(i,l) = sum
              speeds(i,l) = sqrt(1.0*thref/(onem*dline(i)))
              rspeds(i,l) = 1.0/speeds(i,l)
              if     (i.ge.i0+ 1-nbdy .and.
     &               i.le.i0+ii+nbdy .and.
     &               j.ge.j0+ 1-nbdy .and.
     &               j.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0) = 1
              endif
            enddo
c     
          endif
c     
        enddo !nports

        if     (mnproc.eq.1) then
          write(lp,*) 
          call flush(lp)
        endif
c     
c       end of initialization
c     
        call xcsync(flush_lp)
        return
      endif !lcount=1
c     
c --- 'wellposed' treatment of pressure and normal velocity fields
c --- not in fact wellposed with this exterior data
c
c --- set ramping factor 
c     when ramp_time=0:no ramping of tide (=full tide)
c     when ramp_time>0: ramping of tide, if...  
c     ...ramp_orig<=timermp<=(ramp_orig+ramp_time)
      if(ramp_time.gt.0.0 ) then
        timermp=d_time 
        if(timermp.ge.ramp_orig)then
          timermp=(timermp-ramp_orig)/ramp_time
c         frmp=(1-exp(-10*timermp))    
          frmp=(1-exp(-5*timermp))
        else      
          frmp=0.0
        endif   
      else   
        frmp=1.0  
      endif   

c     
      do l=1,nports

        if     (kdport(l).eq.4) then
c     
c       western port
c     
          i = ifport(l)
          j = jfport(l)
          call xclget(dline(j),  lnport(l),
     &           depthu,                 i,j,0,1, 0)            
          call xclget(pline(j),  lnport(l),
     &           pbavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)             
          call xclget(uline(j),lnport(l),
     &           ubavg(1-nbdy,1-nbdy,n), i+1,j,0,1, 0)

            
          if (tidflg.ge.1) then  !tide

            do jn= jfport(l),jlport(l)
              do ic = 1, ncon
                tmpr(jn,ic) = z1r_u(l,jn,ic)
                tmpi(jn,ic) = z1i_u(l,jn,ic)
              enddo
            enddo
                            
            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,uu,j,max(itdm,jtdm),lnport(l)) !normal 
            
c           Note!! uu and vv from tpx are transports
            do jn= jfport(l),jlport(l)
               upred(jn)=uu(jn)*frmp*onem/dline(jn)
            enddo
              
            do jn= jfport(l),jlport(l)
              do ic = 1, ncon
                tmpr(jn,ic) = z1r_p(l,jn,ic)
                tmpi(jn,ic) = z1i_p(l,jn,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,zpred,j,max(itdm,jtdm),lnport(l)) 
                                     
            do j= jfport(l),jlport(l)
              zpred(j)=zpred(j)*onem*frmp
            enddo


            do j= jfport(l),jlport(l) 
c         ----set both u and eta at boundary; 1 invariant weighted:
              uline(j)=ulow(j,l)+upred(j)
     &            +w_1*speedw(j,l)*((plow(j,l)+zpred(j))-pline(j))
              pline(j)=w_1c*(plow(j,l)+zpred(j))+w_1*pline(j)

            enddo
      

          else  !no bnd-tide
 

            do j= jfport(l),jlport(l)
c       ----  set both u and eta; 1 invariant weighted:
              uline(j)=ulow(j,l)+
     &                     w_1*speedw(j,l)*(plow(j,l)-pline(j))
              pline(j)=w_1c*plow(j,l)+w_1*pline(j)


            enddo 

          endif   !tide/no tide            
            
          j = jfport(l)                        
          call xclput(pline(j),  lnport(l),
     &           pbavg(1-nbdy,1-nbdy,n), i,  j,0,1)
          call xclput(uline(j),lnport(l),
     &           ubavg(1-nbdy,1-nbdy,n), i,  j,0,1) 
c     

        elseif (kdport(l).eq.3) then

c         eastern port
c     
          i = ifport(l)-1
          j = jfport(l)
          call xclget(dline(j),  lnport(l),
     &           depthu,                 i+1,j,0,1, 0)
          call xclget(pline(j),  lnport(l), 
     &            pbavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j),lnport(l),
     &           ubavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)


          if (tidflg.ge.1) then !tide

            do jn= jfport(l),jlport(l)
              do ic = 1, ncon
                tmpr(jn,ic) = z1r_u(l,jn,ic)
                tmpi(jn,ic) = z1i_u(l,jn,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,uu,j,max(itdm,jtdm),lnport(l)) !normal 

c           Note!! uu and vv from tpx are transports
            do jn= jfport(l),jlport(l)
              upred(jn)=uu(jn)*frmp*onem/dline(jn)
            enddo

            do jn= jfport(l),jlport(l)
              do ic = 1, ncon
                tmpr(jn,ic) = z1r_p(l,jn,ic)
                tmpi(jn,ic) = z1i_p(l,jn,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,zpred,j,max(itdm,jtdm),lnport(l)) 
            

            do j= jfport(l),jlport(l)               
              zpred(j) = zpred(j)*onem*frmp
            enddo

            do j= jfport(l),jlport(l)
c         ----set u and eta on boundary; 1 invariant weighted:
              uline(j)=ulow(j,l)+upred(j)
     &             +w_1*speede(j,l)*(pline(j)-(plow(j,l)+zpred(j)))
              pline(j)=w_1c*(plow(j,l)+zpred(j))+w_1*pline(j)
            enddo


          else !no bnd-tide


            do j= jfport(l),jlport(l)
c         ----set u and eta on boundary; 1 invariant weighted:
              uline(j)=ulow(j,l)
     &             -w_1*speede(j,l)*(plow(j,l)-pline(j))
              pline(j)=w_1c*plow(j,l)+w_1*pline(j)
            enddo

          endif   !tide/no tide         
                        
          j = jfport(l)
          call xclput(pline(j),  lnport(l), 
     &           pbavg(1-nbdy,1-nbdy,n), i,  j,0,1)
          call xclput(uline(j),lnport(l),
     &           ubavg(1-nbdy,1-nbdy,n), i+1,j,0,1)


        elseif (kdport(l).eq.1) then
c     
c         northern port
c     
          j = jfport(l)-1
          i = ifport(l)
          call xclget(dline(i),  lnport(l),
     &           depthv,                 i,j+1,1,0, 0)
          call xclget(pline(i),  lnport(l), 
     &           pbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i),lnport(l),
     &           vbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)


          if (tidflg.ge.1) then !tide

            do in= ifport(l),ilport(l)
              do ic = 1, ncon
                tmpr(in,ic) = z1r_u(l,in,ic)
                tmpi(in,ic) = z1i_u(l,in,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &           astroflag,uu,i,max(itdm,jtdm),lnport(l)) !normal 

c           Note!! uu and vv from tpx are transports
            do in= ifport(l),ilport(l)
              vpred(in)=uu(in)*frmp*onem/dline(in)
            enddo

            do in= ifport(l),ilport(l)
              do ic = 1, ncon
                tmpr(in,ic) = z1r_p(l,in,ic)
                tmpi(in,ic) = z1i_p(l,in,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,zpred,i,max(itdm,jtdm),lnport(l)) 
           
            do i= ifport(l),ilport(l)               
              zpred(i) =zpred(i)*onem*frmp
            enddo

            do i= ifport(l),ilport(l)   
c         ----set u and eta at boundary; 1 invariant weighted:
              uline(i)=ulow(i,l)+vpred(i)
     &              +w_1*speedn(i,l)*(pline(i)-(plow(i,l)+zpred(i)))
              pline(i)=w_1c*(plow(i,l)+zpred(i))+w_1*pline(i)
            enddo

          else !no bnd-tide

            do i= ifport(l),ilport(l)
c         ----set u and eta at boundary; 1 invariant weighted:
              uline(i)=ulow(i,l)-
     &              w_1*speedn(i,l)*(plow(i,l)-pline(i))
              pline(i)=w_1c*plow(i,l)+w_1*pline(i)
            enddo
          
          endif  !tide/no-tide

          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &           pbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
          call xclput(uline(i),lnport(l),
     &           vbavg(1-nbdy,1-nbdy,n), i,j+1,1,0)
c     

        elseif (kdport(l).eq.2.or.kdport(l).eq.5) then
c     
c         southern port
c     
          j = jfport(l)
          i = ifport(l)
          call xclget(dline(i),  lnport(l),
     &           depthv,                 i,j,  1,0, 0)
          call xclget(pline(i),  lnport(l),
     &           pbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i),lnport(l),
     &           vbavg(1-nbdy,1-nbdy,n), i,j+1,1,0, 0)


          if (tidflg.ge.1) then  !tide

            do in= ifport(l),ilport(l)
              do ic = 1, ncon
                tmpr(in,ic) = z1r_u(l,in,ic)
                tmpi(in,ic) = z1i_u(l,in,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,uu,i,max(itdm,jtdm),lnport(l)) !normal 

c           Note!! uu and vv from tpx are transports
            do in= ifport(l),ilport(l)
              vpred(in)=uu(in)*frmp*onem/dline(in)
            enddo

            do in= ifport(l),ilport(l)
              do ic = 1, ncon
                tmpr(in,ic) = z1r_p(l,in,ic)
                tmpi(in,ic) = z1i_p(l,in,ic)
              enddo
            enddo

            call tides_driver(tmpr,tmpi,d_time,
     &             astroflag,zpred,i,max(itdm,jtdm),lnport(l)) 

            do i= ifport(l),ilport(l)
              zpred(i) =  zpred(i)*onem*frmp 
            enddo

            do i= ifport(l),ilport(l)
c         ----set u and eta at boundary; 1 invariant weighted:
              uline(i)=ulow(i,l)+vpred(i)
     &              +w_1*speeds(i,l)*(pline(i)-(plow(i,l)+zpred(i)))
              pline(i)=w_1c*(plow(i,l)+zpred(i))+w_1*pline(i)
            enddo

          else !no bnd-tide

            do i= ifport(l),ilport(l)
c         ----set u and eta at boundary; 1 invariant weighted:
              uline(i)=ulow(i,l)+
     &              w_1*speeds(i,l)*(plow(i,l)-pline(i))
              pline(i)=w_1c*plow(i,l)+w_1*pline(i)
            enddo

          endif  !tide/no-tide

          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &           pbavg(1-nbdy,1-nbdy,n), i,j,  1,0) 
          call xclput(uline(i),lnport(l),
     &           vbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
c    
        endif !kdport

      enddo !nports
c     
      return
      end
      subroutine latbdp(n)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer n
c
c --- apply lateral boundary conditions to   barotropic  flow field
c
c --- port flow version:
c --- similar to the standard 'Browning and Kreiss' MICOM/HYCOM open
c --- boundary condition, except that the exterior normal velocity
c --- is constant in time and exterior pressure = interior pressure.
c --- tangential velocity is not constrained.
c
c --- see also: latbdp
c
c --- the code is as similar as possible to that for the standard case.
c --- so for example, 'speed' is in fact 1/SQRT(gH) which represents
c --- c1/g in the notation of (Bleck and Sun, Open boundary conditions
c --- for MICOM).  The 1/g allows for the use of pressure fields.
c     
c --- Note that East, West, North and South refers to the grid 
c --- (i.e i,j points) and NOT geographic East, West, North and South
c
c --- the first call is made during initialization.
c
c --- Alan J. Wallcraft,  NRL,  November 1999.
c
      logical, parameter :: ldebug_latbdp=.false.
c
      integer, parameter :: mports=9  !maximum number of ports
c
      integer, parameter :: nchar=120
c
      logical     lfatal,lfatalp
      integer     i,j,isec,ifrst,ilast,l
      real        aline(nchar),
     &            dline(itdm+jtdm),xline(itdm+jtdm),
     &            pline(itdm+jtdm),uline(itdm+jtdm,2)
      real        crs,fin,sum,svspin,uvscl,uvscl2,fatal
      real*8      tstep
      character*3 char3
c
      integer nports,kdport(mports),
     &               ifport(mports),ilport(mports),
     &               jfport(mports),jlport(mports),lnport(mports)
      real    pefold,svpnow(mports),svport(mports)
      real*8  refold
      save    nports,kdport,ifport,ilport,jfport,jlport,lnport
      save    pefold,svpnow,svport,refold
c
      real    uportw(jtdm),speedw(jtdm),rspedw(jtdm),
     &        uporte(jtdm),speede(jtdm),rspede(jtdm),
     &        vportn(itdm),speedn(itdm),rspedn(itdm),
     &        vports(itdm),speeds(itdm),rspeds(itdm)
      save    uportw,speedw,rspedw,uporte,speede,rspede,
     &        vportn,speedn,rspedn,vports,speeds,rspeds
c
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
c
      integer lcount
      save    lcount
      data    lcount / 0 /
c
      lcount = lcount + 1
c
c --- the first call just initializes data structures.
c
      if     (lcount.eq.1) then
c
        open(unit=uoff+99,file=trim(flnminp)//'ports.input')
c
c ---   'nports' = number of boundary port sections.
        call blkini(nports,'nports')
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
        if     (nports.lt.0 .or. nports.gt.mports) then
          if     (mnproc.eq.1) then
          write(lp,*) 
          write(lp,*) 'error in latbdp - illegal nports value'
          if     (nports.gt.mports) then
            write(lp,*) 'increase parameter mports to',nports
          endif
          write(lp,*) 
          call flush(lp)
          endif
          call xcstop('(latbdp)')
                 stop '(latbdp)'
        endif
c
c ---   'pefold' = port transport e-folding time in days
        call blkinr(pefold,'pefold','(a6," =",f10.4," days")')
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
c
c ---   switch units from days to baroclinic time steps
c ---   shift lcount to prevent underflow (lcount*refold.ge.0.001)
c
        tstep  = pefold*(86400.d0/batrop)
        refold = 1.d0/tstep
        lcount = lcount + int(tstep)/1000
c
c ---   read in the ports one at a time
c
        do l= 1,nports
c
c ---     port location is w.r.t. u (EW) or v (NS) grid
c ---     and identifies the sea at the port
c ---     the minimum index is 0
c
c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
c ---     'ifport' = first i-index
c ---     'ilport' = last  i-index (=ifport for N or S orientation)
c ---     'jfport' = first j-index
c ---     'jlport' = last  j-index (=jfport for E or W orientation)
c ---     'svpnow' = existing port transport in Sv (+ve towards E or S)
c ---     'svport' = target   port transport in Sv (+ve towards E or S)
c ---     'lnport' = port length (calculated, not input)
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          call blkinr(svpnow(l),'svpnow','(a6," =",f10.4," Sv")')
          call blkinr(svport(l),'svport','(a6," =",f10.4," Sv")')
          if     (mnproc.eq.1) then
          write(lp,*)
          endif
c
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
c
c ---     sanity check.
c
          if     (kdport(l).gt.2) then
            if     (ifport(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdp - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdp)')
                     stop '(latbdp)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdp - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdp)')
                     stop '(latbdp)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &            jfport(l).gt.jlport(l)     ) then
            if     (mnproc.eq.1) then
            write(lp,*) 
            write(lp,*) 'error in latbdp - port',
     &                   ' location is not consistent'
            write(lp,*) 
            call flush(lp)
            endif
            call xcstop('(latbdp)')
                   stop '(latbdp)'
          endif
        enddo
c
        close(unit=uoff+99)
c
c ---   check ports against masks,
c ---   mark the port locations on masks and print them out.
c
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.itdm-2 .or.
     &                j.lt.1 .or. j.gt.jtdm       ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9  !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0+1,j-j0).ne.1 .or.
     &                iu(i-i0+2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.3 .or. i.gt.itdm .or.
     &                j.lt.1 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9  !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0-1,j-j0).ne.1 .or.
     &                iu(i-i0-2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm .or.
     &                j.lt.3 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9  !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0-1).ne.1 .or.
     &                iv(i-i0,j-j0-2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm   .or.
     &                j.lt.1 .or. j.gt.jtdm-2     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9  !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0+1).ne.1 .or.
     &                iv(i-i0,j-j0+2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          endif
c
          if     (lfatalp) then
            write(lp,*) 
            write(lp,*) 'error in latbdp - port ',l,' mislocated',
     &                  '  (mnproc = ',mnproc,')'
            write(lp,*) 
            call flush(lp)
          endif
          lfatal = lfatal .or. lfatalp
        enddo
c
c       local lfatal to global lfatal
c
        if     (lfatal) then
          fatal = 1.0
        else
          fatal = 0.0
        endif
        call xcmaxr(fatal)
        lfatal = fatal.gt.0.5
c
c ---   write out  -iu-  and -iv- arrays, if they are not too big
c ---   data are written in strips nchar points wide
        if     (lfatal .or. max(itdm,jtdm).le.2*nchar) then
          util1(1:ii,1:jj) = iu(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iu array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
c
          util1(1:ii,1:jj) = iv(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iv array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif  ! small region
c
        if     (lfatal) then
          write(lp,*) 
          write(lp,*) 'error in latbdp - bad port(s)'
          write(lp,*) 
          call flush(lp)
          call xchalt('(latbdp)')
                 stop '(latbdp)'
        endif
c
c ---   restore iu and iv, and zero iuopn and ivopn.
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            iu(i,j) = max( iu(i,j), 0 )
            iv(i,j) = max( iv(i,j), 0 )
          enddo
        enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            iuopn(i,j) = 0
            ivopn(i,j) = 0
          enddo
        enddo
c
c ---   initialize the ports
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            sum = 0.0
            i = ifport(l)
            j = jfport(l)
            call xclget(dline(j),lnport(l), depths,i+1,j,0,1, 0)
            call xclget(xline(j),lnport(l), scuy,  i,  j,0,1, 0)
            do j= jfport(l),jlport(l)
              sum = sum + dline(j)*xline(j)
            enddo
            sum = 1.e6/sum
            do j= jfport(l),jlport(l)
              uportw(j) = sum
              speedw(j) = sqrt(thref/(onem*dline(j)))
              rspedw(j) = 1.0/speedw(j)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1p2e12.5)') 
     &          'w port: ',l,i,j,uportw(j),speedw(j)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                iuopn(i-i0,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            sum = 0.0
            i = ifport(l)-1
            j = jfport(l)
            call xclget(dline(j),lnport(l), depths,i,  j,0,1, 0)
            call xclget(xline(j),lnport(l), scuy,  i+1,j,0,1, 0)
            do j= jfport(l),jlport(l)
              sum = sum + dline(j)*xline(j)
            enddo
            sum = 1.e6/sum
            do j= jfport(l),jlport(l)
              uporte(j) = sum
              speede(j) = sqrt(thref/(onem*dline(j)))
              rspede(j) = 1.0/speede(j)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1p2e12.5)') 
     &          'e port: ',l,i,j,uporte(j),speede(j)
              endif
c
              if     (i+1.ge.i0+ 1-nbdy .and.
     &                i+1.le.i0+ii+nbdy .and.
     &                j  .ge.j0+ 1-nbdy .and.
     &                j  .le.j0+jj+nbdy      ) then
                iuopn(i-i0+1,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            sum = 0.0
            j = jfport(l)-1
            i = ifport(l)
            call xclget(dline(i),lnport(l), depths,i,j,  1,0, 0)
            call xclget(xline(i),lnport(l), scuy,  i,j+1,1,0, 0)
            do i= ifport(l),ilport(l)
              sum = sum + dline(i)*xline(i)
            enddo
            sum = 1.e6/sum
            do i= ifport(l),ilport(l)
              vportn(i) = sum
              speedn(i) = sqrt(thref/(onem*dline(i)))
              rspedn(i) = 1.0/speedn(i)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1p2e12.5)') 
     &          'n port: ',l,i,j,vportn(i),speedn(i)
              endif
c
              if     (i  .ge.i0+ 1-nbdy .and.
     &                i  .le.i0+ii+nbdy .and.
     &                j+1.ge.j0+ 1-nbdy .and.
     &                j+1.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0+1) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            sum = 0.0
            j = jfport(l)
            i = ifport(l)
            call xclget(dline(i),lnport(l), depths,i,j+1,1,0, 0)
            call xclget(xline(i),lnport(l), scuy,  i,j,  1,0, 0)
            do i= ifport(l),ilport(l)
              sum = sum + dline(i)*xline(i)
            enddo
            sum = 1.e6/sum
            do i= ifport(l),ilport(l)
              vports(i) = sum
              speeds(i) = sqrt(thref/(onem*dline(i)))
              rspeds(i) = 1.0/speeds(i)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1p2e12.5)') 
     &          's port: ',l,i,j,vports(i),speeds(i)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0) = 1
              endif
            enddo
c
          endif
c
          if     (mnproc.eq.1) then
          write(lp,*) 'port, now/target velocity = ',
     &                l,svpnow(l)*sum,svport(l)*sum
          call flush(lp)
          endif
        enddo
        if     (mnproc.eq.1) then
        write(lp,*) 
        call flush(lp)
        endif
c
c       end of initialization
c
        call xcsync(flush_lp)
        return
      endif
c
c --- 'wellposed' treatment of pressure and normal velocity fields
c --- not in fact wellposed with this exterior data
c
      tstep  = lcount
      svspin = exp( -tstep*refold )
      do l= 1,nports
        uvscl = svport(l) + svspin*(svpnow(l)-svport(l))
c
        if     (kdport(l).eq.4) then
c
c         western port
c
          i = ifport(l)
          j = jfport(l)
          call xclget(dline(j),  lnport(l),
     &                depthu,                 i+1,j,0,1, 0)
          call xclget(xline(j),  lnport(l),
     &                scuy,                   i,  j,0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+2,j,0,1, 0)
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl*uportw(j)+speedw(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)-speedw(j)*pline(j)
            sum=sum+((crs+fin)-uline(j,1))*dline(j)*xline(j)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl2*uportw(j)+speedw(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)-speedw(j)*pline(j)
            pline(j)  =.5*(crs-fin)*rspedw(j)
            uline(j,1)=(crs+fin)-uline(j,1)
            sum=sum+uline(j,1)*dline(j)*xline(j)
          enddo
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,0,1)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,  j,0,1
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(ub - ',
     &                                      l,lnport(l),i,  j,0,1
                call flush(lp)
              endif
c
        elseif (kdport(l).eq.3) then
c
c         eastern port
c
          i = ifport(l)-1
          j = jfport(l)
          call xclget(dline(j),  lnport(l),
     &                depthu,                 i+1,j,0,1, 0)
          call xclget(xline(j),  lnport(l),
     &                scuy,                   i+1,j,0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i-1,j,0,1, 0)
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl*uporte(j)-speede(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)+speede(j)*pline(j)
            sum=sum+((crs+fin)-uline(j,1))*dline(j)*xline(j)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl2*uporte(j)-speede(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)+speede(j)*pline(j)
            pline(j)  =.5*(fin-crs)*rspede(j)
            uline(j,1)=(fin+crs)-uline(j,1)
            sum=sum+uline(j,1)*dline(j)*xline(j)
*             if     (mnproc.eq.1) then
*             write(lp,'(a,i2,2i5,1p2e12.5)') 
*    &          'e port: ',l,i,j,pline(j),uline(j,1)
*             endif
          enddo
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,0,1)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,  j,0,1
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(ub - ',
     &                                      l,lnport(l),i+1,j,0,1
                call flush(lp)
              endif
c
        elseif (kdport(l).eq.1) then
c
c         northern port
c
          j = jfport(l)-1
          i = ifport(l)
          call xclget(dline(i),  lnport(l),
     &                depthv,                 i,j+1,1,0, 0)
          call xclget(xline(i),  lnport(l),
     &                scux,                   i,j+1,1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j-1,1,0, 0)
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl*vportn(i)-speedn(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)+speedn(i)*pline(i)
            sum=sum+((fin+crs)-uline(i,1))*dline(i)*xline(i)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl2*vportn(i)-speedn(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)+speedn(i)*pline(i)
            pline(i)  =.5*(fin-crs)*rspedn(i)
            uline(i,1)=(fin+crs)-uline(i,1)
            sum=sum+uline(i,1)*dline(i)*xline(i)
          enddo
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j+1,1,0)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,j,  1,0
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(vb - ',
     &                                      l,lnport(l),i,j+1,1,0
                call flush(lp)
              endif
c
        elseif (kdport(l).eq.2) then
c
c         southern port
c
          j = jfport(l)
          i = ifport(l)
          call xclget(dline(i),  lnport(l),
     &                depthv,                 i,j,  1,0, 0)
          call xclget(xline(i),  lnport(l),
     &                scux,                   i,j,  1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j+1,1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j+2,1,0, 0)
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl*vports(i)+speeds(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)-speeds(i)*pline(i)
            sum=sum+((crs+fin)-uline(i,1))*dline(i)*xline(i)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl2*vports(i)+speeds(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)-speeds(i)*pline(i)
            pline(i)  =.5*(crs-fin)*rspeds(i)
            uline(i,1)=(crs+fin)-uline(i,1)
            sum=sum+uline(i,1)*dline(i)*xline(i)
          enddo
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,j,  1,0
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(vb - ',
     &                                      l,lnport(l),i,j,  1,0
                call flush(lp)
              endif
c
        endif
c
*       if     (mod(lcount,512).eq.0) then
*         if     (mnproc.eq.1) then
*         write(lp,*) 'latbdp - l,sv,sum = ',l,uvscl,sum/(onem*1.e6)
*         call flush(lp)
*         endif
*       endif
      enddo
c
      return
      end
      subroutine latbdt(n,lll)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer n,lll
c
c --- apply lateral boundary conditions to   barotropic  flow field
c
c --- nested sub-region version:
c --- Uses the 'Browning and Kreiss' MICOM open boundary condition.
c
c --- Note that 'speed' is in fact thref*SQRT(gH) which represents
c --- c1*thref/g in the notation of (Bleck and Sun, Open boundary conditions
c --- for MICOM).  The thref/g allows for the use of pressure fields.
c     
c --- Note that East, West, North and South refers to the grid 
c --- (i.e i,j points) and NOT geographic East, West, North and South
c
c --- the first call is made during initialization.
c
c --- Alan J. Wallcraft,  NRL,  July, 2001.
c
      logical, parameter :: ldebug_latbdt=.false.
c
      integer, parameter :: mports=9  !maximum number of ports
c
      integer, parameter :: nchar=120
c
      logical     lfatal,lfatalp
      integer     i,j,isec,ifrst,ilast,l,npf,npi,npl
      real        aline(nchar),
     &            pline(itdm+jtdm),uline(itdm+jtdm,2)
      real        crs,fin,fatal
      character*3 char3
c
      integer nports,kdport(mports),
     &               ifport(mports),ilport(mports),
     &               jfport(mports),jlport(mports),lnport(mports)
      save    nports,kdport,ifport,ilport,jfport,jlport,lnport
c
      real    speedw(jtdm),rspedw(jtdm),
     &        speede(jtdm),rspede(jtdm),
     &        speedn(itdm),rspedn(itdm),
     &        speeds(itdm),rspeds(itdm)
      real    plnstw(jtdm),ulnstw(jtdm),vlnstw(jtdm),
     &        plnste(jtdm),ulnste(jtdm),vlnste(jtdm),
     &        plnstn(itdm),ulnstn(itdm),vlnstn(itdm),
     &        plnsts(itdm),ulnsts(itdm),vlnsts(itdm)
      save    speedw,rspedw,plnstw,ulnstw,vlnstw,
     &        speede,rspede,plnste,ulnste,vlnste,
     &        speedn,rspedn,plnstn,ulnstn,vlnstn,
     &        speeds,rspeds,plnsts,ulnsts,vlnsts
c
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
c
      integer lcount
      save    lcount
      data    lcount / 0 /
c
      lcount = lcount + 1
c
c --- the first call just initializes data structures.
c
      if     (lcount.eq.1) then
c
        open(unit=uoff+99,file=trim(flnminp)//'ports.input')
c
c ---   'nports' = number of boundary port sections.
        call blkini(nports,'nports')
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
        if     (nports.lt.0 .or. nports.gt.mports) then
          if     (mnproc.eq.1) then
          write(lp,*) 
          write(lp,*) 'error in latbdt - illegal nports value'
          if     (nports.gt.mports) then
            write(lp,*) 'increase parameter mports to',nports
          endif
          write(lp,*) 
          call flush(lp)
          endif
          call xcstop('(latbdt)')
                 stop '(latbdt)'
        endif
c
c ---   read in the ports one at a time
c
        do l= 1,nports
c
c ---     port location is w.r.t. u (EW) or v (NS) grid
c ---     and identifies the sea at the port
c ---     the minimum index is 0
c
c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
c ---     'ifport' = first i-index
c ---     'ilport' = last  i-index (=ifport for N or S orientation)
c ---     'jfport' = first j-index
c ---     'jlport' = last  j-index (=jfport for E or W orientation)
c ---     'lnport' = port length (calculated, not input)
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          if     (mnproc.eq.1) then
          write(lp,*)
          endif
c
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
c
c ---     sanity check.
c
          if     (kdport(l).gt.2) then
            if     (ifport(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdt - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdt)')
                     stop '(latbdt)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdt - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdt)')
                     stop '(latbdt)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &            jfport(l).gt.jlport(l)     ) then
            if     (mnproc.eq.1) then
            write(lp,*) 
            write(lp,*) 'error in latbdt - port',
     &                   ' location is not consistent'
            write(lp,*) 
            call flush(lp)
            endif
            call xcstop('(latbdt)')
                   stop '(latbdt)'
          endif
        enddo
c
        close(unit=uoff+99)
c
c ---   check ports against masks,
c ---   mark the port locations on masks and print them out.
c
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.itdm-2 .or.
     &                j.lt.1 .or. j.gt.jtdm       ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9  !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0+1,j-j0).ne.1 .or.
     &                iu(i-i0+2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.3 .or. i.gt.itdm .or.
     &                j.lt.1 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9  !indicate an error
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0-1,j-j0).ne.1 .or.
     &                iu(i-i0-2,j-j0).ne.1     ) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm .or.
     &                j.lt.3 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9  !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0-1).ne.1 .or.
     &                iv(i-i0,j-j0-2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm   .or.
     &                j.lt.1 .or. j.gt.jtdm-2     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9  !indicate an error
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0+1).ne.1 .or.
     &                iv(i-i0,j-j0+2).ne.1     ) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  7  !indicate an error
              endif
            enddo
c
          endif
c
          if     (lfatalp) then
            write(lp,*) 
            write(lp,*) 'error in latbdt - port ',l,' mislocated',
     &                  '  (mnproc = ',mnproc,')'
            write(lp,*) 
            call flush(lp)
          endif
          lfatal = lfatal .or. lfatalp
        enddo  !l=1,nports
c
c       local lfatal to global lfatal
c
        if     (lfatal) then
          fatal = 1.0
        else
          fatal = 0.0
        endif
        call xcmaxr(fatal)
        lfatal = fatal.gt.0.5
c
c ---   write out  -iu-  and -iv- arrays, if they are not too big
c ---   data are written in strips nchar points wide
c
        if     (lfatal .or. max(itdm,jtdm).le.2*nchar) then
          util1(1:ii,1:jj) = iu(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iu array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
c
          util1(1:ii,1:jj) = iv(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iv array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif  ! small region
c
        if     (lfatal) then
          write(lp,*) 
          write(lp,*) 'error in latbdt - bad port(s)'
          write(lp,*) 
          call flush(lp)
          call xchalt('(latbdt)')
                 stop '(latbdt)'
        endif
c
c ---   restore iu and iv, and zero iuopn and ivopn.
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            iu(i,j) = max( iu(i,j), 0 )
            iv(i,j) = max( iv(i,j), 0 )
          enddo
        enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            iuopn(i,j) = 0
            ivopn(i,j) = 0
          enddo
        enddo
c
c ---   define the nested boundary input mask.
c
        do j= 1,jj
          do i= 1,ii
            maskbc(i,j) = 0
          enddo
        enddo
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)-1
            do j= jfport(l),jlport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)-1
            do i= ifport(l),ilport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          endif
        enddo  !l=1,nports
c
        if     (ldebug_latbdt) then
          util1(1:ii,1:jj) = maskbc(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'bc array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif !ldebug_latbdt
c
c ---   initialize the ports
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            j = jfport(l)
            call xclget(pline(j),lnport(l), depths,i,j,0,1, 0)
            do j= jfport(l),jlport(l)
              speedw(j) = sqrt(thref/(onem*pline(j)))
              rspedw(j) = 1.0/speedw(j)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1pe12.5)') 
     &          'w port: ',l,i,j,speedw(j)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                iuopn(i-i0,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)-1
            j = jfport(l)
            call xclget(pline(j),lnport(l), depths,i,j,0,1, 0)
            do j= jfport(l),jlport(l)
              speede(j) = sqrt(thref/(onem*pline(j)))
              rspede(j) = 1.0/speede(j)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1pe12.5)') 
     &          'e port: ',l,i,j,speede(j)
              endif
c
              if     (i+1.ge.i0+ 1-nbdy .and.
     &                i+1.le.i0+ii+nbdy .and.
     &                j  .ge.j0+ 1-nbdy .and.
     &                j  .le.j0+jj+nbdy      ) then
                iuopn(i-i0+1,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)-1
            i = ifport(l)
            call xclget(pline(i),lnport(l), depths,i,j,1,0, 0)
            do i= ifport(l),ilport(l)
              speedn(i) = sqrt(thref/(onem*pline(i)))
              rspedn(i) = 1.0/speedn(i)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1pe12.5)') 
     &          'n port: ',l,i,j,speedn(i)
              endif
c
              if     (i  .ge.i0+ 1-nbdy .and.
     &                i  .le.i0+ii+nbdy .and.
     &                j+1.ge.j0+ 1-nbdy .and.
     &                j+1.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0+1) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            i = ifport(l)
            call xclget(pline(i),lnport(l), depths,i,j,1,0, 0)
            do i= ifport(l),ilport(l)
              speeds(i) = sqrt(thref/(onem*pline(i)))
              rspeds(i) = 1.0/speeds(i)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i5,1pe12.5)') 
     &          's port: ',l,i,j,speeds(i)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0) = 1
              endif
            enddo
c
          endif
        enddo  !l=1,nports
        if     (ldebug_latbdt .and. mnproc.eq.1) then
        write(lp,*) 
        call flush(lp)
        endif
c
c       end of initialization
c
        call xcsync(flush_lp)
        return
      endif
c
c --- nested input only required on first barotropic time step.
c
      if     (lll.eq.1) then
        do j= 1,jj
          do i= 1,ii
            if     (maskbc(i,j).eq.1) then
              util1(i,j) = ubnest(i,j,ln0)*wb0+ubnest(i,j,ln1)*wb1
              util2(i,j) = vbnest(i,j,ln0)*wb0+vbnest(i,j,ln1)*wb1
              util3(i,j) = pbnest(i,j,ln0)*wb0+pbnest(i,j,ln1)*wb1
              util4(i,j) = ubpnst(i,j,ln0)*wb0+ubpnst(i,j,ln1)*wb1
              util5(i,j) = vbpnst(i,j,ln0)*wb0+vbpnst(i,j,ln1)*wb1
            endif
          enddo
        enddo
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            j = jfport(l)
            call xclget(plnstw(j),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! pbnest
            call xclget(ulnstw(j),  lnport(l),
     &                  util4(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! ubpnst
            call xclget(vlnstw(j+1),lnport(l)-1,
     &                  util2(1-nbdy,1-nbdy),   i,  j+1,0,1, 0)  ! vbnest
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)-1
            j = jfport(l)
            call xclget(plnste(j),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! pbnest
            call xclget(ulnste(j),  lnport(l),
     &                  util4(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! ubpnst
            call xclget(vlnste(j+1),lnport(l)-1,
     &                  util2(1-nbdy,1-nbdy),   i,  j+1,0,1, 0)  ! vbnest
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)-1
            i = ifport(l)
            call xclget(plnstn(i),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! pbnest
            call xclget(vlnstn(i),  lnport(l),
     &                  util5(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! vbpnst
            call xclget(ulnstn(i+1),lnport(l)-1,
     &                  util1(1-nbdy,1-nbdy),   i+1,j,  1,0, 0)  ! ubnest
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            i = ifport(l)
            call xclget(plnsts(i),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! pbnest
            call xclget(vlnsts(i),  lnport(l),
     &                  util5(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! vbpnst
            call xclget(ulnsts(i+1),lnport(l)-1,
     &                  util1(1-nbdy,1-nbdy),   i+1,j,  1,0, 0)  ! ubnest
          endif
        enddo  !l=1,nports
      endif !lll.eq.1
c
c --- 'wellposed' treatment of pressure and velocity fields.
c --- alternate order of ports in case corners are open.
c
      if     (mod(lll,2).eq.1) then
        npf =  1
        npl =  nports
        npi =  1
      else
        npf =  nports
        npl =  1
        npi = -1
      endif
      do l= npf,npl,npi
c
        if     (kdport(l).eq.4) then
c
c         western port
c
          i = ifport(l)
          j = jfport(l)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+2,j,  0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
          do j= jfport(l),jlport(l)
            crs=                   ulnstw(j) +speedw(j)*plnstw(j)
            fin=1.5*uline(j,1)-0.5*uline(j,2)-speedw(j)*pline( j)
            pline(j)  =0.5*(crs-fin)*rspedw(j)
            uline(j,1)=    (crs+fin)-uline(j,1)
          enddo
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,  0,1)  ! normal
          call xclput(vlnstw(j+1),lnport(l)-1,
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,0,1)  ! tangential
c
        elseif (kdport(l).eq.3) then
c
c         eastern port
c
          i = ifport(l)-1
          j = jfport(l)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i-1,j,  0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            j=jlport(l)
            write(lp,'(a,i4,1p2e12.5)') 'e port, uline:',j,uline(j,1:2)
            write(lp,'(a,i4,1p1e12.5)') 'e port, pline:',j,pline(j)
            write(lp,'(a,i4,1p1e12.5)') 'e port, plnst:',j,plnste(j)
            write(lp,'(a,i4,1p1e12.5)') 'e port, ulnst:',j,ulnste(j)
            endif
          do j= jfport(l),jlport(l)
            crs=                   ulnste(j) -speede(j)*plnste(j)
            fin=1.5*uline(j,1)-0.5*uline(j,2)+speede(j)*pline( j)
            pline(j)  =0.5*(fin-crs)*rspede(j)
            uline(j,1)=    (fin+crs)-uline(j,1)
          enddo
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            j=jlport(l)
            write(lp,'(a,i4,1p2e12.5)') 'e port,   crs:',j,crs,fin
            write(lp,'(a,i4,1p1e12.5)') 'e port, pbavg:',j,pline(j)
            write(lp,'(a,i4,1p1e12.5)') 'e port, ubavg:',j,uline(j,1)
            write(lp,'(a,i4,1p1e12.5)') 'e port, vbavg:',j,vlnste(j)
            write(lp,*)
            call flush(lp)
            endif
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  0,1)  ! normal
          call xclput(vlnste(j+1),lnport(l)-1,
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,0,1)  ! tangential
c
        elseif (kdport(l).eq.1) then
c
c         northern port
c
          j = jfport(l)-1
          i = ifport(l)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j-1,1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            i=ilport(l)
            write(lp,'(a,i4,1p2e12.5)') 'n port, uline:',i,uline(i,1:2)
            write(lp,'(a,i4,1p1e12.5)') 'n port, pline:',i,pline(i)
            write(lp,'(a,i4,1p1e12.5)') 'n port, plnst:',i,plnstn(i)
            write(lp,'(a,i4,1p1e12.5)') 'n port, vlnst:',i,vlnstn(i)
            endif
          do i= ifport(l),ilport(l)
            crs=                   vlnstn(i) -speedn(i)*plnstn(i)
            fin=1.5*uline(i,1)-0.5*uline(i,2)+speedn(i)*pline( i)
            pline(i)  =0.5*(fin-crs)*rspedn(i)
            uline(i,1)=    (fin+crs)-uline(i,1)
          enddo
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            i=ilport(l)
            write(lp,'(a,i4,1p2e12.5)') 'n port,   crs:',i,crs,fin
            write(lp,'(a,i4,1p1e12.5)') 'n port, pbavg:',i,pline(i)
            write(lp,'(a,i4,1p1e12.5)') 'n port, vbavg:',i,uline(i,1)
            write(lp,'(a,i4,1p1e12.5)') 'n port, ubavg:',i,ulnstn(i)
            write(lp,*)
            call flush(lp)
            endif
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,1,0)  ! normal
          call xclput(ulnstn(i+1),lnport(l)-1,
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  1,0)  ! tangential
c
        elseif (kdport(l).eq.2) then
c
c         southern port
c
          j = jfport(l)
          i = ifport(l)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+2,1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
          do i= ifport(l),ilport(l)
            crs=                   vlnsts(i) +speeds(i)*plnsts(i)
            fin=1.5*uline(i,1)-0.5*uline(i,2)-speeds(i)*pline( i)
            pline(i)  =0.5*(crs-fin)*rspeds(i)
            uline(i,1)=    (crs+fin)-uline(i,1)
          enddo
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j,  1,0)  ! normal
          call xclput(ulnsts(i+1),lnport(l)-1,
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  1,0)  ! tangential
c
        endif
c
      enddo  !l=1,nports
c
      return
      end
c
c
c> Revision history:
c>
c> Mar. 2004 -- fixed bug in latbdp's speed calculation 
c> Nov. 2006 -- added latbdf
