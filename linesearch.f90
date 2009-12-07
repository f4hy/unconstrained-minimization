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
       real, dimension(size(x)) :: p
     end function f
     function grad(p) Result(del)
       real, dimension(size(x)) :: p
       real, dimension(size(x)) :: del
     end function grad
     function hessian(p) result(hess)
       real, dimension(size(x),size(x)) :: hess
       real, dimension(size(x))  :: p
     end function hessian
     subroutine steepestdecent(grad,point,direction)
       real :: point(:)
       real :: direction(:)
       interface
          function grad(p) Result(del)
            real :: p(size(x))
            real, dimension(size(x)) :: del
          end function grad
       end interface
     end subroutine steepestdecent
     subroutine newtondirection(grad,hessian,point,direction)
       real :: point(:)
       real :: direction(:)
       interface
          function grad(p) Result(del)
            real, dimension(1:size(x)) :: p
            real, dimension(1:size(x)) :: del
          end function grad
          function hessian(p) result(hess)
            real, dimension(size(x),size(x)) :: hess
            real, dimension(size(x))  :: p
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
       real, dimension(1:size(point)) :: p
       real, dimension(1:size(point)) :: del
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
  real, intent(in) :: point(:)
  real, intent(out) :: direction(:)
  real :: D
  interface
     function grad(p) Result(del)
       real, dimension(size(point)) :: p
       real, dimension(size(point)) :: del
     end function grad
     function hessian(p) result(hess)
       real, dimension(size(point),size(point)) :: hess
       real, dimension(size(point))  :: p
     end function hessian
  end interface
  real, dimension(size(point)) :: g
  real :: ipv(size(point))
  real :: info
  D = size(point)
  g= -grad(point)
  call SGESV( size(point), 1, hessian(point), size(point), ipv, g, size(point), info )
  direction = g
  
end subroutine newtondirection





