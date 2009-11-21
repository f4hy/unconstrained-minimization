  ! Methods for the final project in math 693A
  !
  ! This final contains UMICK CHOlSOLVE MACHINEPS UMSTOP UMSTOP0 and
  ! MODELHESS.
  !
  ! Writen by Brendan Fahy. Adpated from psuedocode in Appendix A "A
  ! Modula system of Algorithems for Unconstrained Minimization and
  ! Nonlinear Equations"
  !
  ! I am unsure of the copyright and lisence status of the code from
  ! the Appendix.


module optimization
  integer n
  real :: macheps

  
end module optimization


! subroutine UMINICK(x0)
!   use optimization
!   real, intent(in) :: x0(n)

! end subroutine UMINICK


real function macheps()
  macheps = 1.0e0
  do while(1.0e0+macheps .NE. 1.0e0)
     macheps = macheps / 2.0e0
  end do
  macheps = macheps * 2.0e0
  print *, macheps
  return 
end function macheps
