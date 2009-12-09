module size
  integer, parameter :: n = 2
end module size

module analytics
  
contains
  real function fn(p)
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
    fn = sum(f)
  end function fn


  function grad(p) Result(del)
    use size
    real :: p(n)
    real :: del(n)
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
    real :: hess(n,n)
    real :: p(n)
    real :: x1,x2
    integer :: one,two
    integer :: i

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
end module analytics


program rosen
  use size
  use analytics
  real :: x0(n)
  integer :: i

 
  do i = 1,n/2
     x0(2*i-1) = -1.2
     x0(2*i) = 1.0
  end do
  
  print *, fn(x0)
  print *, grad(x0)

  call minimize(x0,fn,grad,hessian)

end program rosen


