subroutine solution(Ap, Ai, Ax, b, x)

!program main
use comvar
use mUMFPACK
implicit none

! real components
integer :: Ap(0:Nz)
integer :: Ai(SM )
real(8) :: Ax(SM )
real(8) :: b(Nz),x(Nz)

! UMFPACK constant
integer :: sys=UMFPACK_A
! ---------- SELECT ONE CHOICE: ----------
! C pointers
  type(c_ptr) :: symbolic,numeric
! integer pointers
! integer(c_intptr_t) :: symbolic,numeric
! ----------------------------------------
! zero-based arrays
real(8) :: control(0:UMFPACK_CONTROL-1),info(0:UMFPACK_INFO-1)
integer :: i

print '(a)',"Basic Fortran interface (umf4* calls)"

! double-int (di)
call umf4def (control)
call umf4sym (Nz, Nz, Ap, Ai, Ax, symbolic, control, info)
print '(a,f8.3,a,i0)',"symbolic phase: time",Info(UMFPACK_SYMBOLIC_TIME)
call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info)
print '(a,f8.3,a,i0)',"numeric phase:  time",Info(UMFPACK_NUMERIC_TIME),", flop ",int(Info(UMFPACK_FLOPS),8)
call umf4fsym (symbolic)
call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info)
print '(a,f8.3,a,i0)',"solve phase:    time",Info(UMFPACK_SOLVE_TIME),", flop ",int(Info(UMFPACK_SOLVE_FLOPS),8)
call umf4fnum (numeric)
do i=lbound(x,1),ubound(x,1)
  print '(a,i0,a,f0.10)',"x (",i,") = ",x(i)
enddo



end subroutine
