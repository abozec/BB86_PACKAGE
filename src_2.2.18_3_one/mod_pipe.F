      module mod_pipe
      use mod_xc  ! HYCOM communication interface
c
c --- HYCOM (named pipe based) debugging interface
c
      logical, save, public  :: lpipe
c
      integer, save, private :: ipunit,lpunit,ishift,jshift,nsym
      logical, save, private :: ldebug,ldebugssh,
     &                          lmaster,lpipeio,lshift,lslave,
     &                          lsym,ltracer
c
      real,    allocatable, dimension(:,:),
     &         save, private :: field1,field2,tmask,vmask,amask

      contains
c
c --- this set of routines facilitates output comparison from two HYCOM
c --- versions running side by side. one model, the 'slave', writes its
c --- output into a named pipe. the other model, the 'master', reads
c --- from the pipe and compares.
c --- differences are recorded in 'PIPE_base.out'.
c
c --- call 'pipe_init' at start of main program.
c
c ---   if the file 'PIPE_MASTER' exists then this is the master,
c ---   if the file 'PIPE_SLAVE'  exists then this is the slave,
c ---   if the file 'PIPE_SYM'    exists then this is master and slave,
c ---   if the file 'PIPE_TRACER' exists then this is master and slave,
c ---   otherwise there is no comparison made.
c
c ---   if the file 'PIPE_SHIFT' exists for the slave, then it
c ---   is a single-line plain text file containing two integers
c ---   specifiying how much to periodically shift the slave arrays
c ---   before sending them to the master.  it is an error for
c ---   'PIPE_SHIFT' to exist (a) on the master and (b) when not
c ---   making a comparison.
c
c ---   if the file 'PIPE_SYM' exists, there is no slave and the
c ---   master compares its own fields for various symmetries.
c ---   it is a single-line plain text file containing an integer
c ---   specifiying what kind of symmetries to test for (0=constant,
c ---   1=transpose, 2=constant-in-j, -2=arctic, 4=4-way, 8=8-way).  
c ---   it is an error for 'PIPE_SYM' to exist when making a 
c ---   master/slave comparison.
c
c ---   if the file 'PIPE_TRACER' exists, there is no slave and the
c ---   master checks that all appropriate tracers are non-negative
c ---   and compares temperature to any "temperature" tracer.
c
c ---   if the file 'PIPE_DEBUG' exists then debugging printout
c ---   is produced for point itest,jtest (>0, see blkdat.f).  
c ---   if itest=-1 the min/max/iospycnal of th3d are printed.
c ---   if jtest=-1 the basin-wide means are printed.
c ---   this works with or without a pipe for comparison.
c
c ---   if the file 'PIPE_DEBUG_SSH' exists then debugging printout
c ---   is produced for SSH at point itest,jtest.
c ---   this works with or without a pipe for comparison.
c
c ---   the 'PIPE_MASTER' and 'PIPE_SLAVE' files contain the location
c ---   of an existing named-pipe.  the 'PIPE_SHIFT' file contains the
c ---   periodic shift to apply on the slave.  the 'PIPE_SYM' file
c ---   contains the kind of symmetries to test for.  the contents of
c ---   the 'PIPE_DEBUG' and 'PIPE_DEBUG_SSH' files are ignored.
c
c --- call 'pipe_compare'      (from master and slave) anywhere in the code
c ---   to check whether data stored in a single array are identical
c
c --- call 'pipe_compare_sym1' anywhere in the code to check whether
c ---   data stored in a single p-grid array are symmetrical.
c ---   note that this can be used in place of 'pipe_compare', since
c ---   it will call the latter in master/slave mode.
c
c --- call 'pipe_compare_sym2' anywhere in the code to check whether
c ---   data stored in vector u and v grid arrays are symmetrical.
c ---   note that this can be used in place of 'pipe_compare', since
c ---   it will call the latter twice (for u and v) in master/slave mode.
c
c --- call 'pipe_comparall'   (from master and slave) after major routines
c ---   to check whether data stored in all major arrays are identical or
c ---   symmetric.
c
      subroutine pipe_init
      implicit none
c
      character*256 cpipe
c
      character*12  cinfo
      integer       irecl
c
      inquire(file='PIPE_MASTER',    exist=lmaster)
      inquire(file='PIPE_SLAVE',     exist=lslave)
      inquire(file='PIPE_SHIFT',     exist=lshift)
      inquire(file='PIPE_SYM',       exist=lsym)
      inquire(file='PIPE_TRACER',    exist=ltracer)
      inquire(file='PIPE_DEBUG',     exist=ldebug)
      inquire(file='PIPE_DEBUG_SSH', exist=ldebugssh)
c
      if     (lmaster .and. lslave) then
        call xchalt('pipe_init: (master/slave ambiguity)')
               stop 'pipe_init: (master/slave ambiguity)'
      endif
      if     (lsym .and. (lmaster .or. lslave)) then
        call xchalt('pipe_init: (sym/master/slave ambiguity)')
               stop 'pipe_init: (sym/master/slave ambiguity)'
      endif
      lpipe   = lmaster .or. lslave .or. lsym
      lpipeio = lmaster .or. lslave
      if     (lshift .and. .not.lslave) then
        call xchalt('pipe_init: (shift ambiguity)')
               stop 'pipe_init: (shift ambiguity)'
      endif
c
      if     (lpipe .or. ltracer) then
c
c ---   allocate arrays for comparison
c
        allocate( field1(itdm,jtdm) )
        allocate( field2(itdm,jtdm) )
        if     (.not.lslave) then
          allocate( tmask( itdm,jtdm) )
          allocate( amask( 1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        endif
        if     (lsym) then
          allocate( vmask( itdm,jtdm) )
        endif
      endif  ! pipe
c
      if     (lpipeio) then
c
c ---   open the pipe and some output files
c
        ipunit=18
        lpunit=19
c
        if     (mnproc.eq.1) then
          if     (lmaster) then
            open (unit=17,file='PIPE_MASTER',status='old',
     &            form='formatted')
            read (     17,'(a)') cpipe
            close(unit=17)
            write(lp,'(a,a)') 'master opening pipe for reading: ',
     &                        cpipe(1:len_trim(cpipe))
            call flush(lp)
#if defined(ALPHA)
c ---       work-around a compiler bug by skipping irecl
            open (unit=ipunit,file=cpipe,status='old',
     &            action='read',
     &            form='unformatted')
#else
            cinfo=' '  !removes spurious compiler warning message
            inquire( iolength=irecl ) cinfo,field1
            open (unit=ipunit,file=cpipe,status='old',
     &            action='read',recl=irecl,
     &            form='unformatted')
#endif
            open (unit=lpunit,file='PIPE_base.out',status='unknown')
          else  ! slave
            open (unit=17,file='PIPE_SLAVE', status='old',
     &            form='formatted')
            read (     17,'(a)') cpipe
            close(unit=17)
            write(lp,'(a,a)') 'slave opening pipe for writing: ',
     &                        cpipe(1:len_trim(cpipe))
            call flush(lp)
#if defined(ALPHA)
c ---       work-around a compiler bug by skipping irecl
            open (unit=ipunit,file=cpipe,status='old',
     &            action='write',
     &            form='unformatted')
#else
            cinfo=' '  !removes spurious compiler warning message
            inquire( iolength=irecl ) cinfo,field1
            open (unit=ipunit,file=cpipe,status='old',
     &            action='write',recl=irecl,
     &            form='unformatted')
#endif
            open (unit=lpunit,file='PIPE_test.out',status='unknown')
c
            if     (lshift) then
              open (unit=17,file='PIPE_SHIFT', status='old',
     &              form='formatted')
              read (     17,*) ishift,jshift
              close(unit=17)
              write(lp,'(a,2i5)') 'slave periodic shift is:',
     &                            ishift,jshift
              call flush(lp)
            endif ! shift
          endif ! master/slave
        endif  !1st tile only.
        call xcsync(flush_lp)
      endif  ! pipeio
c
      if     (lsym) then
        open (unit=17,file='PIPE_SYM', status='old',
     &        form='formatted')
        read (     17,*) nsym
        close(unit=17)
        if     (mnproc.eq.1) then
        lpunit=19
        open (unit=lpunit,file='PIPE_base.out',status='unknown')
        write(lpunit,'(a,i2)') 'symmetry type is:',nsym
        write(lp,    '(a,i2)') 'symmetry type is:',nsym
        call flush(lpunit)
        endif
        if     (nsym.ne. 0 .and.
     &          nsym.ne. 1 .and.
     &          nsym.ne. 2 .and.
     &          nsym.ne.-2      ) then
          if     (mnproc.eq.1) then
          write(lp,'(a)') 'symmetry type is not supported'
          endif
          call xcstop('(pipe_init)')
                 stop '(pipe_init)'
        endif
        call xcsync(flush_lp)
      endif ! sym
c
      return
      end subroutine pipe_init
c
      subroutine pipe_compare(field,mask,what)
      implicit none
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: mask
      character*12,
     &         intent(in) :: what
c
c --- call this routine from anywhere in the code (from both versions, of
c --- course) to check whether data stored in 'field' are identical
c
      integer      i,i1,j,j1
      logical      fail
      character*12 which
c
      if (lpipeio) then
        if (lmaster) then
          do j=1,jj
            do i=1,ii
              amask(i,j) = mask(i,j)
            enddo
          enddo
          call xcaget(tmask,  amask, 1)
        endif !master
        call xcaget(field2, field, 1)
c
        if     (mnproc.eq.1) then
          if (lslave) then
            if (.not.lshift) then
              write (lpunit,'(2a)') 'writing for comparison: ',what
              call flush(lpunit)
              write (ipunit) what, field2
            else  ! shift slave array by ishift,jshift
              do j=1,jtdm
                j1 = mod( j-1+jshift+jtdm, jtdm ) + 1
                do i=1,itdm
                  i1 = mod( i-1+ishift+itdm, itdm ) + 1
                  field1(i1,j1) = field2(i,j)
                enddo
              enddo
              write (lpunit,'(2a)') 'writing for comparison: ',what
              call flush(lpunit)
              write (ipunit) what, field1
            endif
          else  ! master
            read  (ipunit) which,field1
            write (lpunit,'(2a)') 'reading for comparison: ',which
            call flush(lpunit)
            if (what.ne.which) then
              write (lpunit,'(4a)') 'out of sync -- trying to compare ',
     &           what,'  to  ',which
              call xchalt('(pipe_compare)')
                     stop '(pipe_compare)'
            endif
c
            fail=.false.
            do j=1,jtdm
              do i=1,itdm
                if (tmask(i,j).gt.0.0 .and.
     &              field2(i,j).ne.field1(i,j)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  master:',field2(i,j),
     &              '  error:', field2(i,j)-field1(i,j),what
                  fail=.true.
                endif
              enddo
            enddo
            if (fail) then  ! optional
              call xchalt('(pipe_compare)')
                     stop '(pipe_compare)'
            endif
          endif !slave:master
        endif !1st tile
        call xcsync(no_flush) ! wait for 1st tile
      endif !lpipeio
      return
      end subroutine pipe_compare

      subroutine pipe_compare_sym1(field,mask,what)
      implicit none
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: mask
      character*12,
     &         intent(in) :: what
c
c --- call this routine from anywhere in the code 
c --- to check whether data stored in 'field' is symmetric.
c
c --- pass through to pipe_compare when in master/slave mode.
c
      integer      i,io,j,jo
      logical      fail
c
      if     (lpipeio) then
        call pipe_compare(field,mask,what)
      elseif (lsym) then
        do j=1,jj
          do i=1,ii
            amask(i,j) = mask(i,j)
          enddo
        enddo
        call xcaget(tmask,  amask, 1)
        call xcaget(field1, field, 1)
        if     (mnproc.eq.1) then
          write (lpunit,'(2a)') 'comparing: ',what
          call flush(lpunit)
          fail=.false.
          if     (nsym.eq.-2) then  !arctic
            j  = jtdm
            jo = jtdm-1
            do i=1,itdm
              io = itdm-mod(i-1,itdm)
              if (tmask(i,j).gt.0.0 .and.
     &            field1(i,j).ne.field1(io,jo)) then
                write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &            'i,j=',i,j,
     &            '  orig :',field1(i,j),
     &            '  error:',field1(i,j)-field1(io,jo),what
                fail=.true.
              endif
            enddo !i
          else  !standard symteries
            do j=1,jtdm
            do i=1,itdm
              if     (nsym.eq.0) then  ! constant field
                if (field1(i,j).ne.field1(1,1)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j)-field1(1,1),what
                  fail=.true.
                endif
              elseif (nsym.eq.2) then  ! constant field in j direction
                if (field1(i,j).ne.field1(i,1)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j)-field1(i,1),what
                  fail=.true.
                endif
              elseif (nsym.eq.1) then  ! p=p.transpose
                if (tmask(i,j).gt.0.0 .and.
     &              field1(i,j).ne.field1(j,i)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j)-field1(j,i),what
                  fail=.true.
                endif
              endif
            enddo !i
            enddo !j
          endif !nsym==-2:else
          if (fail) then  ! optional
            call flush(lpunit)
            call xchalt('(pipe_compare_sym1)')
                   stop '(pipe_compare_sym1)'
          endif
        endif !1st tile
        call xcsync(no_flush) ! wait for 1st tile
      endif !lpipeio:sym
      return
      end subroutine pipe_compare_sym1

      subroutine pipe_compare_sym2(field_u,mask_u,what_u,
     &                             field_v,mask_v,what_v)
      implicit none
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: field_u,field_v
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: mask_u,mask_v
      character*12,
     &         intent(in) :: what_u,what_v
c
c --- call this routine from anywhere in the code 
c --- to check whether data stored in 'field_[uv]' is symmetric.
c
c --- pass through to pipe_compare when in master/slave mode.
c
      integer      i,io,j,jo
      logical      fail
c
      if     (lpipeio) then
        call pipe_compare(field_u,mask_u,what_u)
        call pipe_compare(field_v,mask_v,what_v)
      elseif (lsym) then
        do j=1,jj
          do i=1,ii
            amask(i,j) = mask_u(i,j)
          enddo
        enddo
        call xcaget(tmask,  amask,   1)
        call xcaget(field1, field_u, 1)
        call xcaget(field2, field_v, 1)
        if     (nsym.eq.-2) then  !arctic
          do j=1,jj
            do i=1,ii
              amask(i,j) = mask_v(i,j)
            enddo
          enddo
          call xcaget(vmask,  amask,   1)
        endif
        if     (mnproc.eq.1) then
          write (lpunit,'(4a)') 'comparing: ',what_u,' and ',what_v
          call flush(lpunit)
          fail=.false.
          if     (nsym.eq.-2) then  !arctic
            j  = jtdm
            jo = jtdm-1
            do i=1,itdm
              io = mod(itdm-(i-1),itdm)+1
              if (tmask(i,j).gt.0.0 .and.
     &            field1(i,j).ne.-field1(io,jo)) then
                write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &            'i,j=',i,j,
     &            '  orig :',field1(i,j),
     &            '  error:',field1(i,j)+field1(io,jo),what_u
                fail=.true.
              endif
            enddo !i
            j  = jtdm
            jo = jtdm
            do i=1,itdm
              io = itdm-mod(i-1,itdm)
              if (vmask(i,j).gt.0.0 .and.
     &            field2(i,j).ne.-field2(io,jo)) then
                write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &            'i,j=',i,j,
     &            '  orig :',field2(i,j),
     &            '  error:',field2(i,j)+field2(io,jo),what_v
                fail=.true.
              endif
            enddo !i
          else  !standard symteries
            do j=1,jtdm
            do i=1,itdm
              if     (nsym.eq.0) then  ! constant field
                if (field1(i,j).ne.field1(1,1)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j)-field1(1,1),what_u
                  fail=.true.
                endif
                if (field2(i,j).ne.field2(1,1)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field2(i,j),
     &              '  error:',field2(i,j)-field2(1,1),what_v
                  fail=.true.
                endif
              elseif (nsym.eq.2) then  ! constant field in j direction
                if (field1(i,j).ne.field1(i,1)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j)-field1(i,1),what_u
                  fail=.true.
                endif
                if (field2(i,j).ne.field2(i,1)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field2(i,j),
     &              '  error:',field2(i,j)-field2(i,1),what_v
                  fail=.true.
                endif
              elseif (nsym.eq.1) then  ! u==v.transpose
                if (tmask(i,j).gt.0.0 .and.
     &              field1(i,j).ne.field2(j,i)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  uvel :',field1(i,j),
     &              '  error:',field1(i,j)-field2(j,i),what_u
                  fail=.true.
                endif
              endif
            enddo !i
            enddo !j
          endif
          if (fail) then  ! optional
            call xchalt('(pipe_compare_sym2)')
                   stop '(pipe_compare_sym2)'
          endif
        endif !1st tile
        call xcsync(no_flush) ! wait for 1st tile
      endif !lpipeio:sym
      return
      end subroutine pipe_compare_sym2

      subroutine pipe_compare_same(fielda,fieldb,mask,what)
      implicit none
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: fielda,fieldb
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: mask
      character*12,
     &         intent(in) :: what
c
c --- call this routine from anywhere in the code 
c --- to check whether data stored in 'fielda' and 'fieldb'
c --- are identical.
c
c --- only active in PIPE_TRACER mode.
c --- typically fielda is temp and fieldb is a "temperature" tracer.
c
      integer      i,j
      logical      fail
c
      if     (ltracer) then
        do j=1,jj
          do i=1,ii
            amask(i,j) = mask(i,j)
          enddo
        enddo
        call xcaget(tmask,  amask,  1)
        call xcaget(field1, fielda, 1)
        call xcaget(field2, fieldb, 1)
        if     (mnproc.eq.1) then
          write (lpunit,'(2a)') 'comparing: ',what
          call flush(lpunit)
          fail=.false.
          do j=1,jtdm
            do i=1,itdm
              if (tmask(i,j).gt.0.0 .and.
     &            field1(i,j).ne.field2(i,j)) then
                write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &            'i,j=',i,j,
     &            '  orig :',field1(i,j),
     &            '  error:',field1(i,j)-field2(i,j),what
                fail=.true.
              endif
            enddo
          enddo
          if (fail) then  ! optional
            call xchalt('(pipe_compare_same)')
                   stop '(pipe_compare_same)'
          endif
        endif !1st tile
        call xcsync(no_flush) ! wait for 1st tile
      endif !ltracer
      return
      end subroutine pipe_compare_same

      subroutine pipe_compare_notneg(field,mask,what)
      implicit none
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: mask
      character*12,
     &         intent(in) :: what
c
c --- call this routine from anywhere in the code 
c --- to check whether data stored in 'field' is non-negative.
c
c --- only active in PIPE_TRACER mode.
c --- typically field is a tracer.
c
      integer      i,j
      logical      fail
c
      if     (ltracer) then
        do j=1,jj
          do i=1,ii
            amask(i,j) = mask(i,j)
          enddo
        enddo
        call xcaget(tmask,  amask,  1)
        call xcaget(field1, field, 1)
        if     (mnproc.eq.1) then
          write (lpunit,'(2a)') 'comparing: ',what
          call flush(lpunit)
          fail=.false.
          do j=1,jtdm
            do i=1,itdm
              if (tmask(i,j).gt.0.0) then
                if (field1(i,j).lt.0.0) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j),what
                  fail=.true.
                elseif (field1(i,j).ne.field1(i,j)) then !NaN
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  orig :',field1(i,j),
     &              '  error:',field1(i,j),what
                  fail=.true.
                endif !errors
              endif !tmask.gt.0
            enddo
          enddo
          if (fail) then  ! optional
            call xchalt('(pipe_compare_notneg)')
                   stop '(pipe_compare_notneg)'
          endif
        endif !1st tile
        call xcsync(no_flush) ! wait for 1st tile
      endif !ltracer
      return
      end subroutine pipe_compare_notneg

      subroutine pipe_comparall(m,n, cinfo)
      implicit none
c
      include 'common_blocks.h'
c
      integer, intent(in) :: m,n
      character*12,
     &         intent(in) :: cinfo
c
c --- write out a standard menu of arrays for testing
c
      logical hycom_isnaninf  !function to detect NaN and Inf
c
      character*99 cformat
      character*12 txt1,txt2
      integer      i,imax,imin,j,jja,jmax,jmin,k,ktr,l,mnp
      real         diso,dmax,dmin,damax,damin
      real*8       tmean,smean,pmean,rmean
      real*8       d1,d2,d3,d4
      real         r1,r2,r3,r4
c
      real*8       tmean0,smean0,rmean0,
     &             tmean1,smean1,rmean1
      save         tmean0,smean0,rmean0,
     &             tmean1,smean1,rmean1
      data         tmean0,smean0,rmean0 / 3*0.0d0 /
c
cdiag if     (mnproc.eq.1) then
cdiag write(lp,'(a,i10)') cinfo,nstep
cdiag call flush(lp)
cdiag endif
c
      if     (ldebugssh) then
c
c ---   printout SSH in cm at point itest,jtest.
c
        if     (min(ittest,jttest).le.0) then
          call xcstop('(pipe_comparall: debug_ssh ambiguity)')
                 stop '(pipe_comparall: debug_ssh ambiguity)'
        endif
        if     (i0.lt.ittest .and. i0+ii.ge.ittest .and.
     &          j0.lt.jttest .and. j0+jj.ge.jttest      ) then
c         ssh,montg,thref*pbavg (cm)
          write (lp,"(i8,i5,i4,1x,a,a,3f15.8)")
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),':',
     &       (100.0/g)*srfhgt(itest,jtest),
     &       (100.0/g)*montg1(itest,jtest),
     &       (100.0/g)*srfhgt(itest,jtest)-
     &       (100.0/g)*montg1(itest,jtest)
        endif  ! ittest,jttest tile
        call xcsync(flush_lp)
      endif !ldebugssh
c
      if     (ldebug .and. ittest.ne.-1 .and. jttest.ne.-1) then
c
c ---   printout at point itest,jtest.
c
        if     (min(ittest,jttest).le.0) then
          call xcstop('(pipe_comparall: debug ambiguity)')
                 stop '(pipe_comparall: debug ambiguity)'
        endif
        if     (i0.lt.ittest .and. i0+ii.ge.ittest .and.
     &          j0.lt.jttest .and. j0+jj.ge.jttest      ) then
        if     (ntracr.eq.0) then
          write(cformat,'(a,a)')
     &      '(i8,i5,i4,1x,a,a/',
     &      '(i8,5x,i4,1x,a,a,2f7.3,2f7.3,f8.4,f9.3,f9.2))'
        else
          write(cformat,'(a,i2,a,a,i2,a)')
     &      '(i8,i5,i4,1x,a,a,',ntracr,'a / ',
     &      '(i8,5x,i4,1x,a,a,2f7.3,2f7.3,f8.4,f9.3,f9.2,',
     &      ntracr,'f8.4))'
        endif
*       write(lp,'(3a)') '"',trim(cformat),'"'
        if     (.not.mxlkrt) then
          write (lp,cformat) 
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),
     &       ':   utot   vtot   temp   saln    dens    thkns     dpth',
     &      ('  tracer',ktr=1,ntracr),
     &      (nstep,k,          cinfo(1:6),':',
     &          u(itest,jtest,k,n)+ubavg(itest,jtest,n),
     &          v(itest,jtest,k,n)+vbavg(itest,jtest,n),
     &       temp(itest,jtest,k,n),
     &       saln(itest,jtest,k,n),
     &       th3d(itest,jtest,k,n)+thbase,
     &         dp(itest,jtest,k,n)/onem,
     &          p(itest,jtest,k+1)/onem,
     &       (tracer(itest,jtest,k,n,ktr),ktr=1,ntracr),
     &       k=1,kk)
        else
c ---     include KT mixed layer values.
          write (lp,cformat) 
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),
     &       ':   utot   vtot   temp   saln    dens    thkns     dpth',
     &      ('  tracer',ktr=1,ntracr),
     &       nstep,0,          cinfo(1:6),':',
     &          0.0,
     &          0.0,
     &       tmix(itest,jtest),     
     &       smix(itest,jtest),     
     &      thmix(itest,jtest)+thbase,
     &     dpmixl(itest,jtest,n)/onem,
     &     dpmixl(itest,jtest,n)/onem,
     &       (0.0,ktr=1,ntracr),
     &      (nstep,k,          cinfo(1:6),':',
     &          u(itest,jtest,k,n)+ubavg(itest,jtest,n),
     &          v(itest,jtest,k,n)+vbavg(itest,jtest,n),
     &       temp(itest,jtest,k,n),
     &       saln(itest,jtest,k,n),
     &       th3d(itest,jtest,k,n)+thbase,
     &         dp(itest,jtest,k,n)/onem,
     &          p(itest,jtest,k+1)/onem,
     &       (tracer(itest,jtest,k,n,ktr),ktr=1,ntracr),
     &       k=1,kk)
        endif
        if     (mxlmy) then
          write(cformat,'(a,a)')
     &      '(i8,i5,i4,1x,a,a/',
     &      '(i8,5x,i4,1x,a,a,g15.5,g15.5,f9.3,f9.2))'
          write (lp,cformat) 
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),
     &       ':             q2            q2l    thkns     dpth',
     &      (nstep,k,          cinfo(1:6),':',
     &         q2(itest,jtest,k,n),
     &        q2l(itest,jtest,k,n),
     &         dp(itest,jtest,k,n)/onem,
     &          p(itest,jtest,k+1)/onem,
     &       k=1,kk)
        endif  !'mxlmy'
        if     (cinfo(1:6).eq.'mxkprf' .and. .not.mxlkrt) then
          write(cformat,'(a,a)')
     &      '(i8,i5,i4,1x,a,a/',
     &      '(i8,5x,i4,1x,a,a,f7.3,f8.2,f7.3,f8.2,f9.3,f9.2))'
          write (lp,cformat) 
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),
     &       ':   temp  t-diff   saln  s-diff    thkns     dpth',
     &      (nstep,k,          cinfo(1:6),':',
     &       temp(itest,jtest,k,n),
     &       dift(itest,jtest,k+1)*1.e4,
     &       saln(itest,jtest,k,n),
     &       difs(itest,jtest,k+1)*1.e4,
     &         dp(itest,jtest,k,n)/onem,
     &          p(itest,jtest,k+1)/onem,
     &       k=1,klist(itest,jtest))
        endif  !'mxkprf'
        endif  ! ittest,jttest tile
        call xcsync(flush_lp)
      endif
c
      if     (ldebug .and. ittest.eq.-1) then
c
c ---   printout min/max/iospycnal th3d
c
 104    format (i8,a3,1x,a,a)
 105    format (i8,i3,1x,a,a,2i5,f9.5,f7.3,f9.5,2i5,i7)
        if     (mnproc.eq.1) then
        write(lp,104)
     &       nstep,'  k',cinfo(1:6),
     &       ': imin jmin  denamin deniso  denamax imax jmax mnproc'
        endif
        call xcsync(flush_lp)
        do k= 1,kk
          diso=sigma(k)-thbase
          dmin= huge
          dmax=-huge
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                if     (th3d(i,j,k,n).lt.dmin) then
                  dmin=th3d(i,j,k,n)
                  imin=i
                  jmin=j
                endif
                if     (th3d(i,j,k,n).gt.dmax) then
                  dmax=th3d(i,j,k,n)
                  imax=i
                  jmax=j
                endif
              endif
            enddo
          enddo
          damin=dmin
          call xcminr(damin)
          damax=dmax
          call xcmaxr(damax)
          do mnp= 1,ijpr
            if     (mnp.eq.mnproc) then
              if     (dmin.eq.damin .or. dmax.eq.damax) then
                write (lp,105) 
     &            nstep,k,cinfo(1:6),
     &            ':',imin,jmin,dmin-diso,
     &                          diso+thbase,
     &                          dmax-diso,imax,jmax,mnproc
              endif
            endif
            call xcsync(flush_lp)
          enddo
        enddo
        call flush(lp)
      endif
c
      if     (ldebug .and. jttest.eq.-1) then
c
c ---   printout basin-wide means.
c
#if defined(ARCTIC)
c ---   Arctic (tripole) domain, top row is replicated (ignore it)
        jja = min( jj, jtdm-1-j0 )
        if     (jja.ne.jj) then
          do i=1,ii
            util5(i,jj)=0.0
            util6(i,jj)=0.0
            util3(i,jj)=0.0
            util4(i,jj)=0.0
          enddo
        endif
#else
        jja = jj
#endif
!$OMP   PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jja
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util5(i,j)=              dp(i,j,1,n)*scp2(i,j)
                util6(i,j)=temp(i,j,1,n)*dp(i,j,1,n)*scp2(i,j)
                util3(i,j)=saln(i,j,1,n)*dp(i,j,1,n)*scp2(i,j)
                util4(i,j)=th3d(i,j,1,n)*dp(i,j,1,n)*scp2(i,j)
              do k=2,kk
                util5(i,j)=util5(i,j)+
     &                                   dp(i,j,k,n)*scp2(i,j)
                util6(i,j)=util6(i,j)+
     &                     temp(i,j,k,n)*dp(i,j,k,n)*scp2(i,j)
                util3(i,j)=util3(i,j)+
     &                     saln(i,j,k,n)*dp(i,j,k,n)*scp2(i,j)
                util4(i,j)=util4(i,j)+
     &                     th3d(i,j,k,n)*dp(i,j,k,n)*scp2(i,j)
              enddo
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcsum(d1, util5,ip)
        call xcsum(d2, util6,ip)
        call xcsum(d3, util3,ip)
        call xcsum(d4, util4,ip)
        pmean=d1
        tmean=d2/pmean
        smean=d3/pmean
        rmean=d4/pmean
c
 106    format (i8,3x,1x,a,a,3f8.4,1p3e10.2)
        if     (mnproc.eq.1) then
        write (lp,106)
     &    nstep,cinfo(1:6),
     &    ': t,s,th',
     &    tmean,smean,rmean+thbase,
     &    tmean-tmean0,smean-smean0,rmean-rmean0
        call flush(lp)
        endif
c
c ---   NaN detection.
        r1 = d1
        r2 = d2
        r3 = d3
        r4 = d4
        if     (hycom_isnaninf(r1) .or.
     &          hycom_isnaninf(r2)  .or.
     &          hycom_isnaninf(r3)  .or.
     &          hycom_isnaninf(r4)      ) then
          if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error - NaN or Inf detected'
          write(lp,*)
          call flush(lp)
          endif !1st tile
        endif !NaN
!$OMP   PARALLEL DO PRIVATE(j,i,k)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if     (iu(i,j).eq.1) then
              util5(i,j)=u(i,j,1,n)
              do k=2,kk
                util5(i,j)=util5(i,j)+u(i,j,k,n)
              enddo
            endif !iu
            if     (iv(i,j).eq.1) then
              util6(i,j)=v(i,j,1,n)
              do k=2,kk
                util6(i,j)=util6(i,j)+v(i,j,k,n)
              enddo
            endif !iv
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcsum(d1, util5,iu)
        call xcsum(d2, util6,iv)
        r1 = d1
        r2 = d2
        if     (hycom_isnaninf(r1) .or.
     &          hycom_isnaninf(r2)      ) then
          if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error - u or v NaN or Inf detected'
          write(lp,*)
          call flush(lp)
          endif !1st tile
        endif !NaN
c
        if     (cinfo(1:6).eq.'ENTER ') then
          tmean1=tmean
          smean1=smean
          rmean1=rmean
        elseif (cinfo(1:6).eq.'tsadvc') then
          if     (mnproc.eq.1) then
          write (lp,106)
     &      nstep,'cn+tsa',
     &      ': t,s,th',
     &      tmean,smean,rmean+thbase,
     &      tmean-tmean1,smean-smean1,rmean-rmean1
          call flush(lp)
          endif
        elseif (cinfo(1:6).eq.'hybgen') then
          if     (mnproc.eq.1) then
          write (lp,106)
     &      nstep,'EXIT  ',
     &      ': t,s,th',
     &      tmean,smean,rmean+thbase,
     &      tmean-tmean1,smean-smean1,rmean-rmean1
          call flush(lp)
          endif
        endif
c
        tmean0=tmean
        smean0=smean
        rmean0=rmean
      endif
c
      if     (ltracer .and. cinfo(1:1).ne.'i') then
        do ktr= 1,ntracr
          if     (mnproc.eq.1) then
          write (lpunit,'(a,i10)') cinfo,nstep
          endif
          call xcsync(flush_lp)
          if     (trcflg(ktr).eq.2) then
c
c ---       compare temp and this temperature tracer.
c
            do k=1,kk
              write (txt1,'(a9,i3)') 'temp(kn) ',k
              call pipe_compare_same(  temp(1-nbdy,1-nbdy,k,n),
     &                               tracer(1-nbdy,1-nbdy,k,n,ktr),
     $                               ip,txt1)
            enddo
          else
c
c ---       check that tracer is non-negative.
c
            do k=1,kk
              write (txt1,'(a6,i3,i3)') 'tracer',ktr,k
              call pipe_compare_notneg(tracer(1-nbdy,1-nbdy,k,n,ktr),
     $                                 ip,txt1)
            enddo
          endif
        enddo !ktr
        if     (mnproc.eq.1) then
        write (lpunit,'(a,i10,a)') cinfo,nstep,' -- OK'
        endif
        call xcsync(flush_lp)
      endif !ltracer
c
      if     (lpipe) then
c
c ---   pipe_compare_sym[12] works for both lsym and lpipeio.
c
        if     (mnproc.eq.1) then
        write (lpunit,'(a,i10)') cinfo,nstep
        endif
        call xcsync(flush_lp)
        txt1='ubavg(n)    '
        txt2='vbavg(n)    '
        call pipe_compare_sym2(ubavg(1-nbdy,1-nbdy,n),iu,txt1,
     &                         vbavg(1-nbdy,1-nbdy,n),iv,txt2)
        txt1='pbavg(n)    '
        call pipe_compare_sym1(pbavg(1-nbdy,1-nbdy,n),ip,txt1)
        txt1='montg(1)    '
        call pipe_compare_sym1(montg(1-nbdy,1-nbdy,1,1),ip,txt1)
        if     (cinfo(1:6).eq.'icloan' .or.
     &          cinfo(1:6).eq.'icecpl'     ) then  !ice fields
          txt1='covice      '
          call pipe_compare_sym1(covice,ip,txt1)
          txt1='flxice      '
          call pipe_compare_sym1(flxice,ip,txt1)
          txt1='fswice      '
          call pipe_compare_sym1(fswice,ip,txt1)
          txt1='sflice      '
          call pipe_compare_sym1(sflice,ip,txt1)
        endif
        if     (cinfo(1:6).eq.'icloan' .or.
     &          cinfo(1:6).eq.'icecpl' .or.
     &          cinfo(1:6).eq.'thermf'     ) then  !surface fields
          txt1='surflx      '
          call pipe_compare_sym1(surflx,ip,txt1)
          txt1='sswflx      '
          call pipe_compare_sym1(sswflx,ip,txt1)
          txt1='salflx      '
          call pipe_compare_sym1(salflx,ip,txt1)
        endif
        do k=1,kk
          write(txt1(10:12),'(i3)') k
          write(txt2(10:12),'(i3)') k
c
          txt1(1:9) = 'u(kn)  k='
          txt2(1:9) = 'v(kn)  k='
          call pipe_compare_sym2(   u(1-nbdy,1-nbdy,k,n),iu,txt1,
     &                              v(1-nbdy,1-nbdy,k,n),iv,txt2)
          txt1(1:9) = 'dp(kn) k='
          call pipe_compare_sym1(  dp(1-nbdy,1-nbdy,k,n),ip,txt1)
          txt1(1:9) = 'temp(kn) '
          call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,txt1)
          txt1(1:9) = 'saln(kn) '
          call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,txt1)
          txt1(1:9) = 'th3d(kn) '
          call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,txt1)
*         write(lpunit,'(a,a12,a,i3)') '***',txt1,'*** k=',k
        enddo
        if     (mnproc.eq.1) then
        write (lpunit,'(a,i10,a)') cinfo,nstep,' -- OK'
        endif
        call xcsync(flush_lp)
      endif !lpipe
c
      return
      end subroutine pipe_comparall
c
      end module mod_pipe
c
c
c> Revision history:
c>
c> Oct  2000 - added PIPE_DEBUG     for debugging printout
c> Aug  2001 - added PIPE_SHIFT     for shifted   comparision
c> Aug  2001 - added PIPE_SYM       for symmetric comparision
c> Feb  2004 - added PIPE_SYM       arctic option for arctic bipolar patch
c> Jun  2005 - added PIPE_DEBUG_SSH for debugging SSH printout
