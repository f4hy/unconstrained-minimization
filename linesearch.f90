module linesearch
  ! integer, parameter :: Max_iterations =10000
  ! real, parameter :: tol = 1.0e-12

  real :: a = 1.0e0
  real :: contraction = 0.5e0
  real :: c = 0.5e0
end module linesearch

  
subroutine backtrackinglinesearch(x0,p,x,f,grad,hessian)
  use linesearch
  implicit none
  real, intent(in) :: x0(:),p(:)
  real, intent(out) :: x(:)
  optional :: hessian
  interface
     real function f(p)
       use size
       real:: p(n)
     end function f
     function grad(p) Result(del)
       use size
       real :: p(n)
       real :: del(n)
     end function grad
     function hessian(p) result(hess)
       use size
       real :: hess(n,n)
       real  :: p(n)
     end function hessian
     subroutine steepestdecent(grad,point,direction)
       use size
       real :: point(:)
       real :: direction(:)
       interface
          function grad(p) Result(del)
            use size
            real :: p(n)
            real :: del(n)
          end function grad
       end interface
     end subroutine steepestdecent
     subroutine newtondirection(grad,hessian,point,direction)
       real :: point(:)
       real :: direction(:)
       interface
          function grad(p) Result(del)
            use size
            real :: p(n)
            real :: del(n)
          end function grad
          function hessian(p) result(hess)
            use size
            real :: hess(n,n)
            real  :: p(n)
          end function hessian
       end interface
     end subroutine newtondirection
  end interface

  real :: previous
  ! integer :: iterations = 0


  ! iterations = 0
  x = x0
  ! print 30, iterations, x , a, f(x)
  
  
  previous = f(x)
  a = 1.0e0
  
  if(present(hessian))  then
     print *, "newtoning"
     call newtondirection(grad,hessian,x,p)
  else
     call steepestdecent(grad,x,p)
  end if

  ! Could be more sophisticated here on how to find step size.
  do while (f(x+a*p) .GT. f(x) + c*a*sum(p * grad(x)))
     a = contraction*a
  end do
  x = x+a*p
  ! print 30, iterations, x , a, f(x), (abs(f(x) - previous))
  ! if(abs(f(x) - previous) .LT. tol) then
  !    exit
  ! end if
  ! iterations = iterations+1


! 30 format(I6,2X,4(F15.12,2X),G20.10)
end subroutine backtrackinglinesearch


subroutine steepestdecent(grad,point,direction)
  real, intent(in) :: point(:)
  real, intent(out) :: direction(:)
  real, dimension(size(point)) :: g
  interface
     function grad(p) Result(del)
       use size
       real :: p(n)
       real :: del(n)
     end function grad
  end interface
  g = grad(point)
  direction = - g / norm(g)
  contains
    real function norm(v) 
      real :: v(:)
      norm = sqrt(sum(v*v))
    end function norm
end subroutine steepestdecent

subroutine newtondirection(grad,hessian,point,direction)
  use size
  real, intent(in) :: point(:)
  real, intent(out) :: direction(:)
  real :: D
  interface
     function grad(p) Result(del)
       use size
       real :: p(n)
       real :: del(n)
     end function grad
     function hessian(p) result(hess)
       use size
       real :: hess(n,n)
       real :: p(n)
     end function hessian
  end interface
  real :: g(n)
  real :: ipv(n)
  real :: info
  real :: H(n,n)

  g= -grad(point)
  H = hessian(point)

  if(H(2,1) .eq. 0.0) then
     do i=1,n
        do j=i+1,n
           H(j,i) = H(i,j)
        end do
     end do
  end if
!   print *, H
! 10 format(4(g20.8))
!   call exit(1)

  call DGESV( n , 1, hessian(point), n, ipv, g, n, info )
  direction = g
  
end subroutine newtondirection





