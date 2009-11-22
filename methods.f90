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
  implicit none
  integer n
  real :: typf
  real :: gradtol
  integer ::termcode
  integer :: retcode, itncount, itnlimit
  real ::  steptol
  

  contains
    real function norm(v) 
      real :: v(:)
      norm = sqrt(sum(v*v))
    end function norm
end module optimization


! subroutine UMINICK(x0)
!   use optimization
!   real, intent(in) :: x0(n)

! end subroutine UMINICK


real function macheps()
  ! Determine the smallest possible real epsilon that makes an
  ! additive change for the given precision
  implicit none
  macheps = 1.0e0
  do while(1.0e0+macheps .NE. 1.0e0)
     macheps = macheps / 2.0e0
  end do
  macheps = macheps * 2.0e0
  print *, macheps
  return 
end function macheps

subroutine UMSTOP0(x0,func,grad,Sx,consecmax)
  ! Decide weather to terminate after iteraction zero because x0 is
  ! too close to a critical point
  
  ! Uses termcode, typf, gradtol and n from the module

  use optimization
  implicit none
  real, intent(in) :: x0(n)
  real, intent(in) :: func
  real, intent(in) :: grad(n)
  real, intent(in) :: Sx(n)
  integer, intent(out) :: consecmax
  real :: temp(n)               !Temp storage while looking for max
                                !(The compiler should optimize this
                                !away if it is smart)
  integer :: i

  consecmax = 0
  
  do i=0,n
     temp(i) = abs(grad(i)) * (max(abs(x0(i)),1/Sx(i)) / max(abs(func),typf))
  end do
  if (maxval(temp) .le. gradtol * 1.0e-3 ) then
     termcode = 1 
     print *, "UMSTOP0 failed"
  end if
  termcode = 0
end subroutine UMSTOP0

subroutine UMSTOP(xplus,xc,func,grad,Sx,consecmax,maxtaken)
  ! Decide weather to terminate after iteraction zero because x0 is
  ! too close to a critical point

  ! Uses n, typf, retcode, gradtol, steptol, itncount, itnlimit from module
  use optimization
  implicit none
  real, intent(in) :: xplus(n)
  real, intent(in) :: xc(n)
  real, intent(in) :: func
  real, intent(in) :: grad(n)
  real, intent(in) :: Sx(n)
  integer, intent(inout) :: consecmax
  logical, intent(out) :: maxtaken
  real :: temp(n) 
  integer :: i

  termcode = 0

  if (retcode .EQ. 1) then 
     termcode = 3
     return
  else
     
     do i=0,n
        temp(i) = abs(grad(i)) * (max(abs(xplus(i)),1/Sx(i)) / max(abs(func),typf))
     end do
     if (maxval(temp) .le. gradtol) then
        termcode = 1 

        return
     end if
     
     do i=0,n
        temp(i) = abs(xplus(i) - xc(i)) / max(abs(xplus(i)),1/Sx(i))
     end do
     
     if(maxval(temp) .le. steptol) then
        termcode = 2
        return
     else if (itncount .ge. itnlimit) then
        termcode = 4
        return

     else if(maxtaken) then
        consecmax = consecmax +1
        if (consecmax .eq. 5) then 
           termcode = 5
           return
        end if
        
     else
        consecmax = 0
        return
     end if
  end if
  
end subroutine UMSTOP

