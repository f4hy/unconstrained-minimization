

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
     stepsizej = sqrteta * sign(max(abs(xtemp(j)),1.0),xtemp(j))
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
     stepsize(i) = cuberteta*sign(max(abs(tempx(i)),1.0),tempx(i))
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
     stepsizej = sqrteta *  sign(max(abs(tempx(j)),1.0),tempx(j))
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



