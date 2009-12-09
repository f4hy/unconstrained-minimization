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
  use size
  implicit none
  ! integer :: n
  real :: typf = 1
  integer ::termcode
  integer :: retcode, iterations

  ! real :: maxstep = 1.0e8

  real :: gradtol = 1.0e-12
  real ::  steptol = 1.0e-12
  real :: maxoff
  real :: macheps
  logical :: maxtaken
  ! logical :: newttaken
  logical :: hook = .FALSE.
  integer :: method

  real :: typef = 0.0

  logical :: nohessian = .FALSE.
  logical :: linesearch = .FALSE.
  logical :: dogleg = .FALSE.
  logical :: analyticgrad = .TRUE.
  logical :: analytichessian = .TRUE.

  integer, parameter :: linsearch = 1
  integer, parameter :: maxiterations = 10000
  

  contains
    real function norm(v) 
      real :: v(:)
      norm = sqrt(dot_product(v,v))
      ! norm = sqrt(sum(v*v))
    end function norm
    subroutine computemacheps()
      ! Determine the smallest possible real epsilon that makes an
      ! additive change for the given precision
      implicit none
      macheps = 1.0e0
      do while(1.0e0+macheps .NE. 1.0e0)
         macheps = macheps / 2.0e0
      end do
      macheps = macheps * 2.0e0
      return 
    end subroutine computemacheps

end module optimization

subroutine initalize()
  use optimization

  integer :: input =0
  
  open(unit=7,file="testfile")


  alpha = 1.0
  beta = 1.0

  if (n.lt. 1) then
     termcode = -1
     return
  end if

  call computemacheps()
  


  print *, "Machine precision set down to", macheps

  do while(input .lt. 1 .or.  input .gt. 6)
     print *, "To what accuracy do you want to minimize the problem?"
     print *, "1. 10^-4"
     print *, "2. 10^-8"
     print *, "3. 10^-10"
     print *, "4. 10^-12"
     print *, "5. Macheps*10",macheps*10
     print *, "6. Macheps",macheps
     read(5,*) input
     print *, "input",input
  end do
  if(input .eq. 1) then
     steptol = 1.0e-4
  else if(input .eq. 2) then
     steptol = 1.0e-8
  else if(input .eq. 3) then
     steptol = 1.0e-10
  else if(input .eq. 4) then
     steptol = 1.0e-12
  else if(input .eq. 5) then
     steptol = macheps*10
  else if(input .eq. 6) then
     steptol = macheps
  end if
  gradtol = steptol


  input = 0

  do while(input .lt. 1 .or.  input .gt. 5)
     print *, "Which analytics would you like to use?"
     print *, "1. All analytic"
     print *, "2. Only analytic gradient"
     print *, "3. No analytic hessian or gradient"
     print *, "4. Analytic gradient NO HESSIAN AT ALL"
     print *, "5. no analytic gradient and NO HESSIAN AT ALL!"
     read(5,*) input
     print *, "input",input
  end do

  if(input .eq. 1) then
     analyticgrad = .TRUE.
     analytichessian = .TRUE.
  else if(input .eq. 2) then
     analyticgrad = .TRUE.
     analytichessian = .FALSE.
  else if(input .eq. 3) then
     analyticgrad = .FALSE.
     analytichessian = .FALSE.
  else if(input .eq. 4) then
     analyticgrad = .TRUE.
     analytichessian = .FALSE.
     nohessian = .TRUE.
  else if(input .eq. 5) then
     analyticgrad = .FALSE.
     analytichessian = .FALSE.
     nohessian = .TRUE.
  end if
  
  input = 0
  if(nohessian) input = 2


  do while(input .ne. 1 .AND.  input.ne. 2)
     print *, "Select a method"
     print *, "1. Dogleg"
     print *, "2. linesearch"
     read(5,*) input
     print *, "input",input
  end do

  if(input .eq. 1) then
     dogleg = .TRUE.
  else
     linesearch = .TRUE.
  end if

  input = 0

end subroutine INITALIZE





subroutine UMSTOP0(x0,func,grad,Sx)
  ! Decide weather to terminate after iteraction zero because x0 is
  ! too close to a critical point
  
  ! Uses termcode, typf, gradtol and n from the module

  use optimization
  implicit none
  real, intent(in) :: x0(n)
  real, intent(in) :: func
  real, intent(in) :: grad(n)
  real, intent(in) :: Sx(n)
  real :: temp(n)               !Temp storage while looking for max
                                !(The compiler should optimize this
                                !away if it is smart)
  integer :: i

  typef = func
  
  do i=1,n
     temp(i) = abs(grad(i)) * (max(abs(x0(i)),1/Sx(i)) / max(abs(func),typf))
  end do
  if (maxval(temp) .le. gradtol * 1.0e-3 ) then
     termcode = 1 
     print *, "UMSTOP0 failed"
  end if
  termcode = 0
end subroutine UMSTOP0

subroutine UMSTOP(xplus,xc,func,grad,Sx)
  ! Decide weather to terminate after iteraction zero because x0 is
  ! too close to a critical point

  ! Uses n, typf, retcode, gradtol, steptol, iterations, maxiterations from module
  use optimization
  implicit none
  real, intent(in) :: xplus(n)
  real, intent(in) :: xc(n)
  real, intent(in) :: grad(n)
  real, intent(in) :: Sx(n)
  real, intent(in) :: func
  real :: temp(n) 
  integer :: i
  integer ::consecmax

  termcode = 0
  
  if (retcode .EQ. 1) then 
     termcode = 3
     print *, "Terminating because method failed"
     return
  else
     
     if(abs(typef - func) .gt. 10*typf) then
        typef = func
     end if
     do i=1,n
        temp(i) = max(abs(grad(i)),  max(abs(xplus(i)),1.0) / max(abs(func),typf))
     end do
     print *, "temp",temp
     ! temp = abs(grad)
     if (maxval(temp) .le. gradtol) then
        print *, "gradient small, local minimum found"
        print *, temp
      termcode = 1 
        return
     end if
     
     ! do i=1,n
     !    temp(i) = abs(xplus(i) - xc(i)) / max(abs(xplus(i)),1/Sx(i))
     ! end do
     do i=1,n
        temp(i) = abs(xplus(i) - xc(i)) / max(abs(xplus(i)),1.0)
     end do
     if(maxval(temp) .le. steptol) then
        print *, "Converged within tolerance after",iterations, "iterations"
        termcode = 2
        return
     else if (iterations .ge. maxiterations) then
        print *, "Did not converge after",maxiterations, "steps"
        termcode = 4
        return

     else if(maxtaken) then
        consecmax = consecmax +1
        if (consecmax .eq. 5) then 
           print *, "Five steps of max lenght were taken"
           termcode = 5
           return
        end if
        
     else
        consecmax = 0
        return
     end if
  end if
  
end subroutine UMSTOP

subroutine CHOLSOLVE(grad,L,s)
  ! Solves (LL^T)s = -g for s
  use optimization
  real, intent(in) :: grad(n)
  real, intent(in) :: L(n,n)    !should be lower triangular
  real, intent(out) :: s(n)

  CALL LSOLVE(grad,L,s)
  CALL LTSOLVE(s,L,s)


  s = -s
end subroutine CHOLSOLVE

subroutine LSOLVE(b,L,y)
  ! Solve Ly = b for y
  use optimization
  implicit none
  real, intent(in) :: L(n,n)
  real, intent(in) :: b(n)
  real, intent(out) :: y(n)
  
  integer :: i,j

  y(1) = b(1)/L(1,1)
  do i = 2, n
     y(i) = b(i)
     do j=1,i-1
        y(i) = y(i) - (L(i,j) * y(j))
     end do
     y(i) = y(i) / L(i,i)
  end do

end subroutine LSOLVE

subroutine LTSOLVE(y,L,x)
  ! Solve L^Tx = y for x
  use optimization
  implicit none
  real, intent(in) :: L(n,n)
  real, intent(in) :: y(n)
  real, intent(out) :: x(n)
  
  integer :: i,j

  x(n) = y(n) / L(n,n)
  
  do i = n-1, 1,-1
     x(i) = y(i)
     do j=i+1,n
        x(i) = x(i) - (L(j,i) * x(j))
     end do
     x(i) = x(i) / L(i,i)
  end do

end subroutine LTSOLVE


subroutine CHOLDECOMP(H,maxoffl,L,maxadd)
  ! Find the LL^T decomposition of H+D. D is a diangonal to force
  ! positive definateness.

  ! Modified decomposition from Gill Murra and Wright
  !
  ! H is the input matrix (only upper tri used) (in) 
  ! L is the LL^T decomposition (only lower tri used) (out)
  ! maxadd is the largest component of D (out)

  use optimization
  implicit none
  real, intent(in) :: H(n,n)
  real, intent(out) :: L(n,n)
  real, intent(out) :: maxadd
  real :: maxoffl


  real :: minl,minl2,minljj
  integer :: i,j,k

  minl=0.0
  minl2=0.0
  minljj=0.0
  L = 0

  minl = sqrt(sqrt(macheps)) * maxoffl
  
  if ( maxoffl .eq. 0) then
     maxoffl = sqrt(maxval(diag(H)))
  end if

  minl2 = sqrt(macheps) * maxoffl

  maxadd = 0

  do j=1,n
     L(j,j) = H(j,j)
     do i=1,j-1
        L(j,j) = L(j,j) - L(j,i)**2
     end do
     minljj = 0
     do i=j+1,n
        L(i,j) = H(j,i)
        do k=1,j-1
           L(i,j) = L(i,j) - (L(i,k) * L(j,k))
        end do
        minljj = max(abs(L(i,j)),minljj)
     end do
     minljj = max(minljj/maxoffl,minl)
     if (L(j,j) .gt. minljj**2) then !Normal Cholesky
        L(j,j) = sqrt(L(j,j))
     else
        if (minljj .gt. minl2) then
           minljj = minl2
           maxadd = max(maxadd,minljj**2-L(j,j))
           L(j,j) = minljj
        end if
     end if
     do i = j+1,n
        L(i,j) = L(i,j)/L(j,j)
     end do
  end do

  contains
    function diag(A) result(d)
      
      real :: A(n,n)
      real :: d(n)
      integer :: i
      do i=1,n
         d(i) = A(i,i)
      end do
    end function diag
end subroutine CHOLDECOMP

subroutine modelhess(Sx,H,L)
  ! Force H to be positive definite by adding to the diagonal.
  ! Then calculate the decomposition of H into LL^T
  !
  ! Needs macheps from the module
  ! Calls CHOLDECOMP
  !
  ! Sx is the scaling factors as input
  ! H is input as the hessian and output to be positive definiate
  ! L is the lower triangular decompposition of H (after modified)
  ! 
  use optimization
  
  implicit none
  real, intent(in) :: Sx(n)
  real, intent(inout) :: H(n,n)
  real, intent(out) :: L(n,n)

  integer :: i,j

  real:: sqrteps
  real :: mu
  real :: maxdiag, mindiag, maxposdiag
  real :: maxadd
  real :: maxev,minev
  real :: offrow
  real :: sdd
  real :: maxoffl

  do i=1,n
     do j=i,n
        H(i,j) = H(i,j)/(Sx(i) * Sx(j))
     end do
  end do


  sqrteps = sqrt(macheps)
  

  maxdiag = maxval(diag(H))
  
  mindiag = minval(diag(H))

  maxposdiag = max(0.0,maxdiag)
  
  if (mindiag .le. sqrteps * maxposdiag) then
     mu = 2 * (maxposdiag - mindiag) * sqrteps - mindiag
     maxdiag = maxdiag + mu
  else
     mu =0
  end if
  maxoff = maxval(H-(maxdiag*ident())) !Clever trick, but probably
                                     !inefficent. There should be a
                                     !better way of doing this using
                                     !where statments
  if (maxoff *(1+2*sqrteps) .gt. maxdiag) then
     mu = mu+(maxoff-maxdiag) + 2*sqrteps *maxoff
     maxdiag = maxoff *(1+2*sqrteps)
  end if

  if (maxdiag .eq.0) then
     mu = 1
     maxdiag = 1
  end if

  if (mu .gt. 0) then
     H = H + (mu*ident())
  end if

  maxoffl = sqrt(max(maxdiag,(maxoff/n)))


  call choldecomp(H,maxoffl,L,maxadd)


  if (maxadd .gt. 0) then       !Not possitive deff
     maxev = H(1,1)
     minev = H(1,1)
     ! REMOVED A SUM (YAY VECTOR OPS)
     offrow = sum(H) - sum(diag(H))
     maxev = maxval(diag(H))+offrow
     minev = minval(diag(H))-offrow
     !
     sdd = max((maxev-minev)*sqrteps-minev,0.0)
     mu = min(maxadd,sdd)
     
     H = H+(mu*ident())
  end if

  maxoffl = 0.0
  call choldecomp(H,maxoffl,L,maxadd)

  
  ! Undo scaling
  do i=1,n
     do j=i,n
        H(i,j) = H(i,j)*Sx(i)*Sx(j)
     end do
     do j=1,i
        L(i,j) = L(i,j) * Sx(i)
     end do
  end do
  
contains
  function diag(A) result(d)
    
    real :: A(n,n)
    real :: d(n)
    integer :: i
    do i=1,n
       d(i) = A(i,i)
    end do
  end function diag
  
  function ident() result(I)
    real :: I(n,n)
    integer :: j
    I =0
    do j=1,n
       I(j,j) = 1.0
    end do
  end function ident

  function lowermask() result(L)
    real :: L(n,n)
    integer :: i,j
    L = 0
    do i=1,n
       do j=1,i
          L(i,j) = 1
       end do
    end do
  end function lowermask
  

end subroutine modelhess

