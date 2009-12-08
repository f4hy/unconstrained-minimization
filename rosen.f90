program rosen
  real :: x0(2)
  interface
     real function rosenbock(p)
       real :: p(2)
     end function rosenbock

     function grad(p) Result(del)
       real :: p(1:2)
       real :: del(1:2)
     end function grad
     function hessian(p) result(hess)
       real :: hess(1:2,1:2)
       real :: p(2)
     end function hessian

     subroutine backtrackinglineserach(x0,p,x,f,grad,hessian)
       real :: p(:),x0(:)
       real :: x(:)
       optional :: hessian
       interface
          real function f(p)
            real, dimension(2) :: p
          end function f
          function grad(p) Result(del)
            real, dimension(2) :: p
            real, dimension(2) :: del
          end function grad
          function hessian(p) result(hess)
            real, dimension(2,2) :: hess
            real, dimension(2)  :: p
          end function hessian
       end interface
     end subroutine backtrackinglineserach

  end interface

  x0 = [-1.2,1.0]

  print *, rosenbock(x0)
  print *, grad(x0)
  
  call initalize(2)

  call minimize(x0,rosenbock,grad,hessian)

  !  print *, "Newton"
  ! call backtrackinglinesearch(x0,p,x,f,grad,hessian)
    
  ! ! call backtrackinglineserach(x0,p,x,rosenbock,grad,hessian)
  ! print *, "steepest"
  ! ! call backtrackinglineserach(x0,p,x,rosenbock,grad)
  ! call test()
  ! 
end program rosen

real function rosenbock(p)
  real :: p(2)
  real :: x1,x2
  x1 = p(1)
  x2 = p(2)
  rosenbock = 100*(x2-x1**2)**2+(1-x1)**2
end function rosenbock


function grad(p) Result(del)
  real :: p(1:2)
  real :: del(1:2)
  x1 = p(1)
  x2 = p(2)
  del(1) = -400*x1*((x2-x1**2))+2*x1-2.0
  del(2) = 200*(x2-x1**2)
end function grad

function hessian(p) result(hess)
  real :: hess(1:2,1:2)
  real :: p(2)
  real :: x1,x2
  x1 = p(1)
  x2 = p(2)
  hess(1,1) = 2 + 800*x1**2 - 400*(-x1**2 + x2)
  hess(1,2) = -400*x1
  hess(2,1) = -400*x1
  hess(2,2) = 200
end function hessian

