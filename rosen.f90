module size
  integer, parameter :: n = 2
end module size

program rosen
  integer, parameter :: n = 2
  
  real :: x0(n)
  interface
     real function rosenbock(p)
       use size
       real :: p(n)
     end function rosenbock

     function grad(p) Result(del)
       use size
       real :: p(n)
       real :: del(n)
     end function grad
     function hessian(p) result(hess)
       use size
       real :: p(n)
       real :: hess(n,n)
     end function hessian


  end interface

  print *, "start?"
  x0 = [-1.2,1.0]
  
  print *, rosenbock(x0)
  print *, grad(x0)

  call initalize(n)

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
  use size
  real :: p(n)

  real :: f(n)
  integer :: i

  if(mod(n,2) .ne. 0) then
     print *, "N must be multipule of 2"
     call exit(n)
  end if

  do i = 1,n/2
     f(2*i-1) = 10*(p(2*i) - p(2*i-1)**2)
     f(2*i) = 1-p(2*i-1)
  end do

  f = f**2
  rosenbock = sum(f)
  
  ! real :: x1,xn
  ! x1 = p(1)
  ! x2 = p(2)
  ! rosenbock = 10*(x2-x1**2)**2+(1-x1)**2
end function rosenbock


function grad(p) Result(del)
  use size
  real :: p(n)
  real :: del(1:2)
  real :: x1,x2
  integer :: i

  do i = 1,n/2
     x1 = p(2*i-1)
     x2 = p(2*i)
     del(2*i-1) =  -400*x1*((x2-x1**2))+2*x1-2.0
     del(2*i) = 200*(x2-x1**2)
  end do
  
  ! del = del
     
end function grad

function hessian(p) result(hess)
  use size
  real :: hess(1:2,1:2)
  real :: p(2)
  real :: x1,x2
  integer :: one,two

  hess = 0

  do i = 1,n/2
     one = 2*i-1
     two = 2*i
     x1 = p(one)
     x2 = p(two)

     hess(one,one) = 2 + 800*x1**2 - 400*(-x1**2 + x2)
     hess(one,two) = -400*x1
     hess(two,one) = -400*x1
     hess(two,two) = 200

  end do

end function hessian

