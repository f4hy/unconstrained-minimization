module size
  integer, parameter :: n = 4
end module size

module analytics
contains 
  real function fn(p)
    use size
    real :: p(n)
    real :: x1,x2,x3,x4

    x1 = p(1)
    x2 = p(2)
    x3 = p(3)
    x4 = p(4)

    fn = 100*(x1**2-x2)**2+(1-x1)**2 + 90*(x3**2-x4)**2
    fn = fn + (1-x3)**2 + 10.1*((1-x2)**2+(1-x4)**2) + 19.8*(1-x2)*(1-x4)
  end function fn


  function grad(p) Result(del)
    use size
    real :: p(n)
    real :: del(n)
    real :: x1,x2,x3,x4

    x1 = p(1)
    x2 = p(2)
    x3 = p(3)
    x4 = p(4)

    del(1) = -2*(1-x1) + 400*x1*(x1**2-x2)
    del(2) = -20.2*(1-x2)-200*(x1**2-x2)-19.8*(1-x4)
    del(3) = -2*(1-x3)+360*x3*(x3**2-x4)
    del(4) = -19.8*(1-x2)-20.2*(1-x4)-180*(x3**2-x4)

  end function grad

  function hessian(p) result(hess)
    use size
    real :: hess(n,n)
    real :: p(n)
    real :: x1,x2,x3,x4
    integer :: i

    hess = 0

    do i = 1,n/4
       x1 = p(1)
       x2 = p(2)
       x3 = p(3)
       x4 = p(4)

       hess(1,1) = 2+800*x1**2+400*(x1**2-x2)
       hess(1,2) = -400*x1
       hess(1,3) = 0
       hess(1,4) = 0

       hess(2,1) = -400*x1
       hess(2,2) = 220.2
       hess(2,3) = 0
       hess(2,4) = 19.8

       hess(3,1) = 0
       hess(3,2) = 0
       hess(3,3) = 2+720*x3**2+360*(x3**2-x4)
       hess(3,4) = -360*x3

       hess(4,1) = 0
       hess(4,2) = 19.8
       hess(4,3) = -360*x3
       hess(4,4) = 200.2



    end do

  end function hessian

end module analytics


program wood
  use size
  use analytics
  real :: x0(n)

 
  x0(1) = -3
  x0(2) = -1
  x0(3) = -3
  x0(4) = -1
  
  print *, x0
  print *, fn(x0)
  print *, grad(x0)

  call minimize(x0,fn,grad,hessian)

  !  print *, "Newton"
  ! call backtrackinglinesearch(x0,p,x,f,grad,hessian)
    
  ! ! call backtrackinglineserach(x0,p,x,powell,grad,hessian)
  ! print *, "steepest"
  ! ! call backtrackinglineserach(x0,p,x,powell,grad)
  ! call test()
  ! 
end program wood

