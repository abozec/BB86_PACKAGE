c
c-----------------------------------------------------------------------
c
c     machine dependent I/O routines.
c     message passing version, with all I/O from first processor.
c     contained in module mod_za.
c
c     author:  Alan J. Wallcraft,  NRL.
c
c-----------------------------------------------------------------------
c
      subroutine zagetc(cline,ios, iunit)
      implicit none
c
      character*80, intent(out)   :: cline
      integer,      intent(out)   :: ios
      integer,      intent(in)    :: iunit
c
c**********
c*
c  1) machine specific routine for reading one text line from a file.
c
c  2) The read is performed on the first processor only.
c*
c**********
c
      integer        iline,ibuf
      common/czgetc/ iline(81,0:1),ibuf
      save  /czgetc/
c
      integer i
c
c --- I/O from first processor only
c
      ibuf = mod(ibuf+1,2)
c
      if     (mnproc.eq.1) then
        read(iunit,'(a)',iostat=ios) cline
        do i= 1,80
          iline(i,ibuf) = ichar(cline(i:i))
        enddo
        iline(81,ibuf) = ios
      endif
c
c     broadcast to all other processors
c
      call xcgetc(iline(:,ibuf))
      do i= 1,80
        cline(i:i) = char(iline(i,ibuf))
      enddo
      ios = iline(81,ibuf)  ! iostat value
      return
      end subroutine zagetc

      subroutine zaiopn(cstat, iaunit)
      implicit none
c
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cstat
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zaiost before first call to zaiopn.
c     see also 'zaiope' and 'zaiopf'.
c
c  2) the filename is taken from the environment variable FORxxxA,
c       where xxx = iaunit, with default fort.xxxa.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zaiord and zaiowr.
c     the file should be closed using zaiocl.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w,wminy,wmaxy
      common/czioxw/ w(n2drec),wminy(jtdm),wmaxy(jtdm)
      save  /czioxw/
c
      integer   ios,nrecl
      character cfile*256,cenv*7
      character cact*9
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiopn)')
               stop '(zaiopn)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiopn - iaunit = ',iaunit
*     call flush(lp)
c
c     get filename.
c
      write(cenv,"('FOR',i3.3,'A')") iaunit
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(cfile,"('fort.',i3.3,'a')") iaunit
      endif
c
c     dummy I/O?
c
      if     (cfile.eq.'/dev/null' .or.
     &        cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
        call xcsync(no_flush)
#if defined(TIMER)
        call xctmr1(16)
#endif
        return
      endif
c
c     open file.
c
      inquire(iolength=nrecl) w
c
      if     (cstat.eq.'OLD' .or.
     &        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     &        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
#if defined(X1)
      call asnunit(iaunit+uaoff,'-F event,cachea:4096:4:2 -B on',ios)
      if     (ios.ne.0) then
        write(lp,9050) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif !ios
#endif
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+uaoff,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+uaoff,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif !ios
#endif
      if     (cstat.eq.'scratch' .or.
     &        cstat.eq.'SCRATCH'     ) then
        open(unit=iaunit+uaoff,
     &       form='unformatted', status='scratch',
     &       access='direct', recl=nrecl, action=cact, iostat=ios)
      else
        open(unit=iaunit+uaoff, file=cfile,
     &       form='unformatted', status=cstat,
     &       access='direct', recl=nrecl, action=cact, iostat=ios)
      endif
      if     (ios.ne.0) then
        write(lp,9100) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        if     (cstat.eq.'scratch' .or.
     &          cstat.eq.'SCRATCH'     ) then
          write(6,'(a)')  "status='SCRATCH'"
        else
          write(6,'(3a)') 'FILENAME="',trim(cfile),'"'
        endif
        call flush(lp)
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif !ios
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiopn -  array I/O unit ',
     &   i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in zaiopn -  can''t asnunit ',i3,
     &   ', for array I/O.'/ /)
#endif
 9100 format(/ /10x,'error in zaiopn -  can''t open unit ',i3,
     &   ', for array I/O.'/ /)
      end subroutine zaiopn

      subroutine zaiope(cenv,cstat, iaunit)
      implicit none
c
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cenv,cstat
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zaiost before first call to zaiope.
c     see also 'zaiopn' and 'zaiopf'.
c
c  2) the filename is taken from environment variable 'cenv'.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zaiord and zaiowr.
c      arrays passed to these routines must conform to 'h'.
c     the file should be closed using zaiocl.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w,wminy,wmaxy
      common/czioxw/ w(n2drec),wminy(jtdm),wmaxy(jtdm)
      save  /czioxw/
c
      integer   ios,nrecl
      character cfile*256
      character cact*9
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiope)')
               stop '(zaiope)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiope - iaunit = ',iaunit
*     call flush(lp)
c
c     get filename.
c
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(lp,9300) trim(cenv)
        write(lp,*) 'iaunit = ',iaunit
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif
c
c     dummy I/O?
c
      if     (cfile.eq.'/dev/null' .or.
     &        cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
        call xcsync(no_flush)
#if defined(TIMER)
        call xctmr1(16)
#endif
        return
      endif
c
c     open file.
c
      inquire(iolength=nrecl) w
c
      if     (cstat.eq.'OLD' .or.
     &        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     &        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
#if defined(X1)
      call asnunit(iaunit+uaoff,'-F event,cachea:4096:4:2 -B on',ios)
      if     (ios.ne.0) then
        write(lp,9050) iaunit,trim(cfile)
        write(lp,*) 'ios = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif !ios
#endif
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+uaoff,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+uaoff,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit,trim(cfile)
        write(lp,*) 'ios = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile,
     &     form='unformatted', status=cstat,
     &     access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif !ios
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiope -  array I/O unit ',
     &   i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in zaiope -  can''t asnunit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in zaiope -  can''t open unit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
 9300 format(/ /10x,'error in zaiope -  environment variable ',a,
     &   ' not defined'/ /)
      end subroutine zaiope

      subroutine zaiopf(cfile,cstat, iaunit)
      implicit none
c
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cfile,cstat
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zaiost before first call to zaiopf.
c     see also 'zaiopn' and 'zaiope'.
c
c  2) the filename is taken from 'cfile'.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zaiord and zaiowr.
c      arrays passed to these routines must conform to 'h'.
c     the file should be closed using zaiocl.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w,wminy,wmaxy
      common/czioxw/ w(n2drec),wminy(jtdm),wmaxy(jtdm)
      save  /czioxw/
c
      integer   ios,nrecl
      character cact*9
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiopf)')
               stop '(zaiopf)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiopf - iaunit = ',iaunit
*     call flush(lp)
c
c     dummy I/O?
c
      if     (cfile.eq.'/dev/null' .or.
     &        cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
        call xcsync(no_flush)
#if defined(TIMER)
        call xctmr1(16)
#endif
        return
      endif
c
c     open file.
c
      inquire(iolength=nrecl) w
c
      if     (cstat.eq.'OLD' .or.
     &        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     &        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
#if defined(X1)
      call asnunit(iaunit+uaoff,'-F event,cachea:4096:4:2 -B on',ios)
      if     (ios.ne.0) then
        write(lp,9050) iaunit,trim(cfile)
        write(lp,*) 'ios   = ',ios
        call flush(lp)
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif !ios
#endif
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+uaoff,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+uaoff,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit,trim(cfile)
        write(lp,*) 'ios   = ',ios
        call flush(lp)
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile,
     &     form='unformatted', status=cstat,
     &     access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        call flush(lp)
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif !ios
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiopf -  array I/O unit ',
     &   i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in zaiopf -  can''t asnunit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in zaiopf -  can''t open unit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
      end subroutine zaiopf

      subroutine zaiopi(lopen, iaunit)
      implicit none
c
      logical, intent(out)   :: lopen
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) is an array i/o unit open?
c
c  2) must call zaiost before first call to zaiopi.
c*
c**********
c
      lopen = iarec(iaunit).ne.-1
      return
      end subroutine zaiopi

      subroutine zaiost
      implicit none
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for initializing array i/o.
c
c  2) see also zaiopn, zaiord, zaiowr, and zaiocl.
c*
c**********
c
      integer i
c
      if     (mnproc.eq.1) then
        write(lp,'(/a/)')
     &    'zaiost - Array I/O is Fortran DA I/O from the 1st task'
      endif
c
      do i= 1,999
        iarec(i) = -1
      enddo
#if defined(TIMER)
c
c     initialize timers.
c
      call xctmrn(16,'zaio**')
      call xctmrn(17,'zaiord')
      call xctmrn(18,'zaiowr')
      call xctmrn(19,'zaioIO')
#endif
      return
      end subroutine zaiost

      subroutine zaiocl(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o file closing.
c
c     must call zaiopn for this array unit before calling zaiocl.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c*
c**********
c
      integer ios
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiocl)')
               stop '(zaiocl)'
      endif
c
      iarec(iaunit) = -1
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiocl - iaunit = ',iaunit
*     call flush(lp)
c
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        close(unit=iaunit+uaoff, status='keep')
#if defined(T3E) || defined(YMP) || defined(X1)
        call asnunit(iaunit+uaoff,'-R',ios)
#endif
      endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiocl -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiocl

      subroutine zaiofl(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o buffer flushing.
c
c     must call zaiopn for this array unit before calling zaiocl.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c*
c**********
c
      integer   irlen
      character cfile*256
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiofl)')
               stop '(zaiofl)'
      endif
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
c
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        inquire(unit=iaunit+uaoff, name=cfile, recl=irlen)
        close(  unit=iaunit+uaoff, status='keep')
        open(   unit=iaunit+uaoff, file=cfile, form='unformatted',
     &          access='direct', recl=irlen)
      endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiofl -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiofl

      subroutine zaioiq(iaunit, irec)
      implicit none
c
      integer, intent(in)    :: iaunit
      integer, intent(out)   :: irec
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o inquiry.
c
c  2) returns the number of records processed, or -1 for a closed file.
c*
c**********
c
      irec = iarec(iaunit)
      return
      end subroutine zaioiq

      subroutine zaiorw(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o file rewinding.
c
c     must call zaiopn for this array unit before calling zaiocl.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c*
c**********
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiorw)')
               stop '(zaiorw)'
      endif
c
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = 0
      endif
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiorw - iaunit,rec = ',iaunit,iarec(iaunit)
*     call flush(lp)
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiorw -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiorw

      subroutine zaiord3(h, l, mask,lmask, hmin,hmax,  iaunit)
      implicit none
c
      logical, intent(in)    :: lmask
      integer, intent(in)    :: l,iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin(l),hmax(l)
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(out)   :: h
#else
      real,    intent(out)   :: hmin(l),hmax(l)
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(out)   :: h
#endif
c
c**********
c*
c  1) machine specific routine for 3-d array reading.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the
c     array, ignoring array elements set to 2.0**100.
c     if lmask==.true. the range is calculated only where mask.ne.0,
c     with all other values unchanged in h on exit.  It is then an
c     error if mask.ne.0 anywhere the input is 2.0**100.
c*
c**********
c
c     this version just calls zaiord l times.
c
      integer k
c
      do k= 1,l
        call zaiord(h(1-nbdy,1-nbdy,k), mask,lmask,
     &              hmin(k),hmax(k), iaunit)
      enddo
c
      return
      end subroutine zaiord3

      subroutine zaiord(h, mask,lmask, hmin,hmax,  iaunit)
      implicit none
c
      logical, intent(in)    :: lmask
      integer, intent(in)    :: iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin,hmax
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(out)   :: h
#else
      real,    intent(out)   :: hmin,hmax
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(out)   :: h
#endif
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array reading.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the
c     array, ignoring array elements set to 2.0**100.
c     if lmask==.true. the range is calculated only where mask.ne.0,
c     with all other values unchanged in h on exit.  It is then an
c     error if mask.ne.0 anywhere the input is 2.0**100.
c
c  5) Optimized by Dan Moore, Planning Systems Inc., August 2005.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w,wminy,wmaxy
      common/czioxw/ w(n2drec),wminy(jtdm),wmaxy(jtdm)
      save  /czioxw/
c
      real*4         awtmp(itdm,jtdm)
      equivalence   (w(1),awtmp(1,1))  !saves a data copy
c
      real*4         htmp
      common/czioxr/ htmp(idm*jdm)
      save  /czioxr/
c
      character cfile*256
      integer   ios, i,j
      real*4    wmin,wmax
      real      rmin(1),rmax(1)
#if defined(TIMER)
c
      call xctmr0(17)
#endif
c
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiord)')
               stop '(zaiord)'
      endif
c
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
c
      wmin =  spval
      wmax = -spval
c
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*)'In initial proc 1 section'
*     write(lp,*) 'zaiord - iaunit,rec = ',iaunit,iarec(iaunit)
*     write(lp,*) 'zaiord - mask.1,1    = ',amsk(1,1)
*     write(lp,*) 'zaiord - h.1,1       = ',atmp(1,1)
*     call flush(lp)
c
      if     (iarec(iaunit).eq.-99) then  !dummy I/O
        w(1:n2drec) = 0.0
      else !standard I/O
        call zaiordd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(lp,9100) iarec(iaunit),iaunit
          write(lp,*) 'ios = ',ios
          call flush(lp)
          cfile = ' '
          inquire(unit=iaunit+uaoff,name=cfile)
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
          call flush(lp)
          call xchalt('(zaiord)')
                 stop '(zaiord)'
        endif !ios
       endif ! dummy I/O
*      write(lp,*)'I/O on proc 1 finished?'
*      write(lp,*)'lmask = ',lmask
c
       if     (.not.lmask)then
c        Get global min and max on first processor.
c        must be done here because tiles need not cover the full domain.
!$OMP    PARALLEL DO PRIVATE(j,i)
!$OMP&            SCHEDULE(STATIC,jblk)
         do j= 1,jtdm
#if defined(ENDIAN_IO)
           call zaio_endian(w(1+(j-1)*itdm),itdm)  !swap to big-endian
#endif
           wminy(j) =  spval
           wmaxy(j) = -spval
           do i= 1,itdm
             if     (w(i+(j-1)*itdm).ne.spval) then
               wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
               wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
             endif
           enddo !i
         enddo !j
         wmin = minval(wminy(1:jtdm))
         wmax = maxval(wmaxy(1:jtdm))
       endif !Not Lmask (global min,max on processor #1)
c
       endif !end I/O from first processor only
c
c --- put field from 1st processor to all tiles
      call xcaput4(awtmp,htmp, 1)
c
c --- Each processor loads h from htmp (where mask = 1)
c --- Each processor does local min max if lmask is true
c
      if     (lmask) then
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j= 1,jj
#if defined(ENDIAN_IO)
          call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          wminy(j) =  spval
          wmaxy(j) = -spval
          do i= 1,ii
            if     (mask(i,j).ne.0) then
              h(i,j) =                  htmp(i+(j-1)*ii)
              wminy(j) = min( wminy(j), htmp(i+(j-1)*ii) )
              wmaxy(j) = max( wmaxy(j), htmp(i+(j-1)*ii) )
            endif
          enddo !i
        enddo !j
        wmin = minval(wminy(1:jj))
        wmax = maxval(wmaxy(1:jj))
      else
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            h(i,j) = htmp(i+(j-1)*ii)
          enddo !i
        enddo !j
      endif !lmask:else
c
c --- Min/Max broadcast/gather
c
      rmin(1) = wmin
      rmax(1) = wmax
      call xcminr(rmin)
      call xcmaxr(rmax)
      hmin = rmin(1)
      hmax = rmax(1)
c
      if     (lmask .and. hmax.eq.spval) then
        if     (mnproc.eq.1) then
        write(lp,9200) iarec(iaunit),iaunit
        call flush(lp)
        cfile = ' '
        inquire(unit=iaunit+uaoff,name=cfile)
        write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
        call flush(lp)
        endif !master
        call xcstop('(zaiord)')
               stop '(zaiord)'
      endif
#if defined(TIMER)
c
      call xctmr1(17)
#endif
      return
c
 9000 format(/ /10x,'error in zaiord -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiord -  can''t read record',
     &   i4,' on array I/O unit ',i3,'.'/ /)
 9200 format(/ /10x,'error in zaiord -  record',
     &   i4,' on array I/O unit ',i3,
     &   ' has 2.0**100 outside masked region.'/ /)
      end subroutine zaiord

      subroutine zaiordd(a,n, iunit,irec,ios)
      implicit none
c
      integer, intent(in)    :: n,iunit,irec
      integer, intent(out)   :: ios
      real*4,  intent(out)   :: a(n)
c
c**********
c*
c 1)  direct access read a single record.
#if defined(ENDIAN_IO)
c     on exit, a is still little-endian
#endif
c
c 2)  expressed as a subroutine because i/o with
c     implied do loops can be slow on some machines.
c*
c**********
c
#if defined(TIMER)
c
      call xctmr0(19)
#endif
      read(unit=iunit, rec=irec, iostat=ios) a
#if defined(TIMER)
c
      call xctmr1(19)
#endif
      return
      end subroutine zaiordd

      subroutine zaiosk(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for skipping an array read.
c
c     must call zaiopn for this array unit before calling zaiosk.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c*
c**********
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiosk)')
               stop '(zaiosk)'
      endif
c
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiosk - iaunit,rec = ',iaunit,iarec(iaunit)
*     call flush(lp)
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiosk -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiosk

      subroutine zaiowr3(h, l, mask,lmask, hmin,hmax, iaunit, lreal4)
      implicit none
c
      logical, intent(in)    :: lmask,lreal4
      integer, intent(in)    :: l,iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin(l),hmax(l)
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(inout) :: h
#else
      real,    intent(out)   :: hmin(l),hmax(l)
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(inout) :: h
#endif
c
c**********
c*
c  1) machine specific routine for 3-d array writing.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the array.
c     if lmask==.true. the range is only where mask.ne.0, with all other
c     values output as 2.0**100.
c
c  5) If lreal4==.true. then h is overwritten on exit with real*4 version
c     of the same array.  This is typically used for reproducability on
c     restart.
c*
c**********
c
c     this version just calls zaiowr l times.
c
      integer k
c
      do k= 1,l
        call zaiowr(h(1-nbdy,1-nbdy,k), mask,lmask,
     &              hmin(k),hmax(k), iaunit, lreal4)
      enddo
      return
      end subroutine zaiowr3

      subroutine zaiowr(h, mask,lmask, hmin,hmax,  iaunit, lreal4)
      implicit none
c
      logical, intent(in)    :: lmask,lreal4
      integer, intent(in)    :: iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin,hmax
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(inout) :: h
#else
      real,    intent(out)   :: hmin,hmax
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(inout) :: h
#endif
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array writing.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
c
c  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+uaoff refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the array.
c     if lmask==.true. the range is only where mask.ne.0, with all other
c     values output as 2.0**100.
c
c  5) If lreal4==.true. then h is overwritten on exit with real*4 version
c     of the same array.  This is typically used for reproducability on
c     restart.
c
c  6) Optimized by Dan Moore, Planning Systems Inc., August 2005.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w,wminy,wmaxy
      common/czioxw/ w(n2drec),wminy(jtdm),wmaxy(jtdm)
      save  /czioxw/
c
      real*4         awtmp(itdm,jtdm)
      equivalence   (w(1),awtmp(1,1))  !saves a data copy
c
      real*4         htmp
      common/czioxr/ htmp(idm*jdm)
      save  /czioxr/
c
      character cfile*256
      integer   ios, i,j
      real      rmin(1),rmax(1)
      real*4    data_void(1),vsave4
#if defined(TIMER)
c
      call xctmr0(18)
#endif
c
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiowr)')
               stop '(zaiowr)'
      endif
c
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
c
      data_void(1) = spval
#if defined(ENDIAN_IO)
      call zaio_endian(data_void,1) !swap to big-endian
#endif
c
c --- Copy into real*4 buffer, and find Min,Max
c
      if     (lreal4) then
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              if     (mask(i,j).ne.0) then
                htmp(i+(j-1)*ii) = h(i,j)
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              else
                htmp(i+(j-1)*ii) = spval
              endif
#if defined(REAL4)
! ---         h(i,j) = htmp(i+(j-1)*ii)  ! h is already real*4
#else
              h(i,j) = htmp(i+(j-1)*ii)  ! h is not real*4, so update it
#endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        else !.not.lmask
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              htmp(i+(j-1)*ii) = h(i,j)
              if     (htmp(i+(j-1)*ii).ne.spval) then
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              endif
#if defined(REAL4)
! ---         h(i,j) = htmp(i+(j-1)*ii)  ! h is already real*4
#else
              h(i,j) = htmp(i+(j-1)*ii)  ! h is not real*4, so update it
#endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        endif !lmask:else
      else !.not.lreal4
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              if     (mask(i,j).ne.0) then
                htmp(i+(j-1)*ii) = h(i,j)
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              else
                htmp(i+(j-1)*ii) = spval
              endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        else !.not.lmask
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              htmp(i+(j-1)*ii) = h(i,j)
              if     (htmp(i+(j-1)*ii).ne.spval) then
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        endif !lmask:else
      endif !lreal4:else
c
      rmin(1) = minval(wminy(1:jj))
      rmax(1) = maxval(wmaxy(1:jj))
      call xcminr(rmin)
      call xcmaxr(rmax)
      hmin = rmin(1)
      hmax = rmax(1)
c
c --- I/O from first processor only
#if defined(ENDIAN_IO)
c --- htmp and data_void are already big-endian
#endif
c
      vsave4 = vland4
      vland4 = data_void(1)
      call xcaget4(awtmp,htmp, 1)  !htmp to w (awtmp) on 1st processor
      vland4 = vsave4
c
      if     (mnproc.eq.1) then
        do i= itdm*jtdm+1,n2drec
          w(i) = data_void(1)
        enddo
        if     (iarec(iaunit).ne.-99) then  !standard I/O
          call zaiowrd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
          if     (ios.ne.0) then
            write(lp,9100) iarec(iaunit),iaunit
            call flush(lp)
            cfile = ' '
            inquire(unit=iaunit+uaoff,name=cfile)
            write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
            call flush(lp)
            call xchalt('(zaiowr)')
                   stop '(zaiowr)'
          endif !ios
        endif !standard I/O
      endif !I/O from first processor only
c
#if defined(TIMER)
c
      call xctmr1(18)
#endif
      return
c
 9000 format(/ /10x,'error in zaiowr -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiowr -  can''t write record',
     &   i4,' on array I/O unit ',i3,'.'/ /)
      end subroutine zaiowr

      subroutine zaiowrd(a,n, iunit,irec,ios)
      implicit none
c
      integer, intent(in)    :: n,iunit,irec
      integer, intent(out)   :: ios
      real*4,  intent(in)    :: a(n)
c
c**********
c*
c 1)  direct access write a single record.
#if defined(ENDIAN_IO)
c     a is already big-endian.
#endif
c
c 2)  expressed as a subroutine because i/o with
c     implied do loops can be slow on some machines.
c*
c**********
c
#if defined(TIMER)
c
      call xctmr0(19)
#endif
      write(unit=iunit, rec=irec, iostat=ios) a
#if defined(TIMER)
c
      call xctmr1(19)
#endif
      return
      end subroutine zaiowrd
