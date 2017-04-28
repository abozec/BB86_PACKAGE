      logical function hycom_isnaninf(a)
      implicit none
c
      real a
c
c**********
c*
c 1)  return .true. if a is NaN or +Inf or -Inf.
c*
c**********
c
      hycom_isnaninf = .not. (a.ge.-huge(a) .and. a.le.huge(a))
      end
