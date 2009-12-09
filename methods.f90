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
  real :: gradtol = 0
  integer ::termcode
  integer :: retcode, iterations

  ! real :: maxstep = 1.0e8

  real,parameter ::  steptol = 1.0e-8
  real :: maxoff
  real :: macheps
  logical :: maxtaken
  ! logical :: newttaken
  logical :: hook = .FALSE.
  integer :: method

  real :: globtol = 1.0e-8

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
      norm = sqrt(sum(v*v))
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
     return
  else
     
     do i=1,n
        temp(i) = max(abs(grad(i)),  abs(xplus(i)) / max(abs(func),1.0))
     end do
     print *, "temp",temp
     ! temp = abs(grad)
     if (maxval(temp) .le. gradtol) then
        print *, "gradient small, local minimum found"
        print *, temp
        termcode = 1 
        return
     end if
     
     do i=1,n
        temp(i) = abs(xplus(i) - xc(i)) / max(abs(xplus(i)),1/Sx(i))
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



subroutine FDJAC(x0,f0,grad,Jacob)
  use optimization
  implicit none
  
  real, intent(in) :: x0(n)
  real, intent(in) :: f0(n)
  interface
     function grad(p) Result(del)
       use size
       real :: p(n)
       real :: del(n)
     end function grad
  end interface
  real, intent(out) :: Jacob(n,n)

  real :: sqrteta
  real :: tempj
  real :: stepsizej
  real :: xtemp(n)
  real :: tempgrad(n)
  integer :: i,j

  xtemp = x0

  sqrteta = sqrt(macheps)
  
  do j =1, n
     stepsizej = sqrteta * xtemp(j)
     tempj = xtemp(j)
     xtemp(j) = xtemp(j) + stepsizej
     stepsizej = xtemp(j) - tempj
     tempgrad = grad(xtemp)
     do i=1,n
        Jacob(i,j) = (tempgrad(i) - f0(i)) /stepsizej
     end do
     xtemp(j) = tempj
  end do

end subroutine FDJAC


  
function ffdhessf(x0) result(H)
  use optimization
  use analytics
  implicit none
  
  real, intent(in) :: x0(n)
  ! real, intent(in) :: f0
  ! interface
  !    real function fn(p)
  !      use size
  !      real :: p(n)
  !    end function fn
  ! end interface
  
  real :: H(n,n)
  
  real :: cuberteta
  integer :: i,j
  real :: stepsize(n)
  real :: tempi,tempj,tempx(n)
  real :: fneighbor(n),fii,fij
  real :: f0

  f0 = fn(x0)

  H = 0

  
  tempx = x0
  cuberteta = macheps**(1.0/3.0)
  do i =1,n
     stepsize(i) = cuberteta*tempx(i)
     tempi = tempx(i)
     tempx(i) = tempx(i) + stepsize(i)
     stepsize(i) = tempx(i) - tempi
     fneighbor(i) = fn(tempx)
     tempx(i) = tempi
  end do

  do i =1,n
     tempi = tempx(i)
     tempx(i) = tempx(i) + 2*stepsize(i)
     fii = fn(tempx)
     H(i,i) = (f0-fneighbor(i) + (fii-fneighbor(i)))/(stepsize(i)*stepsize(i))
     tempx(i) = tempi + stepsize(i)
     do j =i+1,n
        tempj = tempx(j)
        tempx(j) = tempx(j) + stepsize(j)
        fij = fn(tempx)
        H(i,j) = ((f0-fneighbor(i)) + (fij-fneighbor(j)))/(stepsize(i)*stepsize(j))
        tempx(j) = tempj
     end do
     tempx(i) = tempi
  end do
end function ffdhessf

function FFDGRAD(x0) result(g)
  use optimization
  use analytics
  implicit none
  
  real, intent(in) :: x0(n)
  real :: g(n)
  
  real :: sqrteta
  real :: stepsizej
  real :: tempj
  real :: tempx(n)
  integer :: i,j
  real :: f0

  f0 = fn(x0)

  tempx = x0
  sqrteta = sqrt(macheps)


  do j=1,n
     stepsizej = sqrteta * tempx(j)
     tempj = tempx(j)
     tempx(j) = tempx(j) + stepsizej
     stepsizej = tempx(j) - tempj
     g(j) = ( fn(tempx)-f0)/stepsizej
     tempx(j) = tempj
  end do
end function FFDGRAD

  
function FFDHESSG(x0) result(H)
  use optimization
  use analytics
  implicit none
  
  real, intent(in) :: x0(n)
  
  real :: H(n,n)
  real :: g0(n)
  integer :: i,j
  g0 = grad(x0)

  call FDJAC(x0,g0,grad,H)
  do i =1,n-1
     do j =i+1,n
        H(i,j) = (H(i,j)+H(j,i))/2
     end do
  end do
end function FFDHESSG








! subroutine trustregionupdate(x,fc,grad,funct,L,step,delta,xprev,fprev,nextx,nextf,hess)
!   use optimization
!   implicit none
!   real, intent(in) :: x(n)
!   real, intent(in) :: fc     !f(x)
!   real, intent(in) :: grad(n)  
!   real, intent(in) :: L(n,n)  
!   ! optional :: hess
!   real, intent(in), optional :: hess(n,n)    !Model hessian
!   real, intent(in) :: step(n)  
!   real, intent(inout) :: delta
!   real, intent(inout) :: xprev(n)
!   real, intent(inout) :: fprev

!   interface
!      real function funct(p)
!        use size
!        real :: p(n)
!      end function funct
!   end interface
  
!   real, intent(out) :: nextx(n)
!   real, intent(out) :: nextf
  
!   integer :: i,j
!   real :: dF,dfpred
!   real :: deltatemp
!   ! real :: alpha = 1.0e-4
!   real :: initslope
!   real :: rellength
!   real :: temp

!   maxtaken = .FALSE.

  
!   nextx = x+step
!   nextf = funct(nextx)
!   dF = nextf - fc

!   initslope = dot_product(grad,step)

!   if (retcode .ne. 3) fprev = 0
!   if (retcode .eq. 3 .and. nextf .ge. fprev .or. dF .gt. 1.0e-4 * initslope ) then
!      nextx = xprev
!      nextf = fprev
!      delta = delta/2.0
!      retcode =0
!      return
!   end if

  

!   if(df .ge. alpha * initslope) then
!      ! f(nextx) too large

!      rellength = maxval(abs(step)/abs(nextx))
!      if (rellength .lt. steptol) then
!         !  nextx-x is too small
!         retcode = 1
!         print *, "RET 1!"

!         nextx = x
!      else
!         retcode = 2
!         deltatemp = -initslope * norm(step) / (2 * (df - initslope) )
!         if (deltatemp .lt. 0.1*delta) then
!            delta = 0.1 * delta
!         else if (deltatemp .gt. 0.5*delta ) then
!            delta = 0.5 * delta
!         else
!            delta = deltatemp
!         end if
!      end if
!   else
!      dfpred = initslope
!      if (present(hess) .AND. HOOK) then  !hook step
!         do i = 1,n
!            temp = 0.5 * hess(i,i) * step(i)**2
!            do j = i+1,n
!               temp = temp + hess(i,j) * step(i) * step(j)
!            end do
!            dfpred = dfpred + temp
!         end do
!      else                       !Dog leg
!         do i=1,n
!            do j=i,n
!               temp = L(j,i) * step(j)
!            end do
!            dfpred = dfpred + (temp * temp/2)
!         end do
!      end if
!      if (retcode .ne. 2 .and. abs(dfpred - df) .le. 0.1 * abs(df) .or. df .le. initslope &
!           & .and. newttaken .eqv. .FALSE. .and. delta .le. .99 *maxstep) then
!         retcode = 3
!         print *, "RET 3!"
!         xprev = nextx
!         fprev = nextf
!         delta = min(2 * delta, maxstep)
!      else
!         retcode = 0
!         print *, "RET 0!"
!         if (norm(step) .ge. 0.99 * maxstep) then
!            maxtaken = .TRUE.
!         end if
!         if (df .ge. 0.1 * dfpred ) then
!            print *, "oops"
!            delta = delta/2
           
!         else if(df .le. 0.75 * dfpred) then
!            delta = min(2*delta,maxstep)
!         end if
        
!      end if
!   end if


! end subroutine trustregionupdate

! subroutine dogstep(grad,L,Sn,newtlen,delta,firstdog,cauchylen,eta,shat,vhat,step)
!   use optimization
!   implicit none
  
!   real, intent(in) :: grad(n)
!   real, intent(in) :: L(n,n)
!   real, intent(in) :: Sn(n)
!   real, intent(in) :: newtlen
!   real, intent(inout) :: delta
!   logical, intent(inout) :: firstdog
!   real, intent(inout) :: cauchylen
!   real, intent(inout) :: eta
!   real, intent(inout) :: shat(n)
!   real, intent(inout) :: vhat(n)
!   real, intent(out) :: step(n)

!   integer :: i,j
!   real :: lambda
!   real :: temp
!   real :: tempv

!   if (newtlen .lt. delta) then
!      newttaken = .TRUE.
!      step = Sn
!      delta = newtlen
!      print *, "dogstep3",step
!      return
!   else
!      newttaken = .FALSE.
!      if (firstdog) then
!         print *, "FIRST DOGGED!"
!         firstdog = .FALSE.
!         alpha = norm(grad)**2
!         beta = 0
!         do i = 1,n
!            temp = 0
!            do j=1,n
!               temp = temp +L(j,i) *grad(j)
!            end do
!            beta = beta + temp*temp
!         end do
!         shat = (alpha/beta) * grad
!         cauchylen = alpha * sqrt(alpha) / beta
!         eta = 0.2+(0.8*alpha**2/(beta * abs(dot_product(grad,Sn))))
!         vhat = eta*Sn-shat
!         if (delta .eq. -1) then
!            delta = min(cauchylen,maxstep)
!         end if
!      end if
!      if(eta * newtlen .lt. delta) then
!         step = (delta/newtlen) * Sn
!         print *, "dogstep1",step
  
!         return
!      else if (cauchylen .ge. delta) then
!         step = (delta/ cauchylen) * shat
!         print *, "dogstep2",step
!         return
!      else
!         temp = dot_product(vhat,shat)
!         tempv = dot_product(vhat,vhat)
!         lambda = (-temp + sqrt( temp**2 - tempv * (cauchylen**2-delta**2))) * tempv
!         step = shat + lambda*vhat
        
!      end if
!   end if


! end subroutine dogstep


! subroutine dogdriver(x0,f0,grad,funct,L,Sn,delta,nextx,nextf)
!   use optimization
!   implicit none
  
!   real, intent(in) :: x0(n)
!   real, intent(in) :: f0
!   real, intent(in) :: grad(n)
  
!   interface
!      real function funct(p)
!        use size
!        real :: p(n)
!      end function funct
!   end interface

!   real, intent(in) :: L(n,n)
!   real, intent(in) :: Sn(n)
  
!   real, intent(inout) :: delta
  
!   real, intent(out) :: nextx(n)
!   real, intent(out) :: nextf


!   logical :: firstdog = .TRUE.
  
!   real :: newtlen

!   real :: cauchylen
!   real :: eta
!   real :: shat(n)
!   real :: vhat(n)
!   real :: xprev(n)
!   real :: fprev
!   real :: step(n)
!   real :: H(n,n)

!   integer :: count = 0

!   alpha = -1.0
!   beta = -1.0
!   delta = -1.0
!   cauchylen = -1.0

!   firstdog = .TRUE.

!   eta = macheps

!   H = 1

!   newttaken = .FALSE.

!   retcode = 4

!   newtlen = norm(Sn)
  
!   count = 0

!   do while(retcode .GT. 1)
!      call dogstep(grad,L,Sn,newtlen,delta,firstdog,cauchylen,eta, &
!           & shat,vhat,step)
!      call trustregionupdate(x0,f0,grad,funct,L,step,delta,xprev,fprev,nextx,nextf)
!   end do
  
! end subroutine dogdriver


