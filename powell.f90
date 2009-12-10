module size
  integer, parameter :: n = 4
end module size

module analytics
contains 
  real function fn(p)
    use size
    real :: p(n)

    real :: f(n)
    integer :: i

    if(mod(n,4) .ne. 0) then
       print *, "N must be multipule of 2"
       call exit(n)
    end if

    do i = 1,n/4
       f(4*i-3) = p(4*1-3) + 10* p(4*i-2)
       f(4*i-2) = sqrt(5.0) * (p(4*i-1 ) - p(4*i))
       f(4*i-1) = (p(4*i-2) - 2*p(4*i-1))**2
       f(4*i) = sqrt(10.0)*(p(4*i-3) - p(4*i))**2
    end do

    f = f**2
    fn = sum(f)

  end function fn


  function grad(p) Result(del)
    use size
    real :: p(n)
    real :: del(n)
    real :: x1,x2,x3,x4
    integer :: i

    do i = 1,n/4
       x1 = p(4*i-3)
       x2 = p(4*i-2)
       x3 = p(4*i-1)
       x4 = p(4*i)  

       del(4*i-3) = 2*(x1+10*x2) + 40*(x1-x4)**3
       del(4*i-2) = 10*(x1 + 10*x2) + 4*(x2-2*x3)**3
       del(4*i-1) = -8*(x2-2*x3)**3+10*(x3-x4)
       del(4*i)   = -40*(x1-x4)**3-10*(x3-x4)
    end do

    ! del = del

  end function grad

  function hessian(p) result(hess)
    use size
    real :: hess(n,n)
    real :: p(n)
    real :: x1,x2,x3,x4
    integer :: one,two,thr,for
    integer :: i

    hess = 0

    do i = 1,n/4
       one = (4*i-3)
       two = (4*i-2)
       thr = (4*i-1)
       for = (4*i)  
       x1 = p(one)
       x2 = p(two)
       x3 = p(thr)
       x4 = p(for)

       hess(one,one) = 2+120*(x1-x4)**2
       hess(one,two) = 20.0
       hess(one,thr) = 0.0
       hess(one,for) = -120*(x1-x4)**2

       hess(two,one) = 20.0
       hess(two,two) = 200 - 12*(x2-2*x3)**2
       hess(two,thr) = -24*(x2-2*x3)**2
       hess(two,for) = 0.0

       hess(thr,one) = 0
       hess(thr,two) = -24*(x2-2*x3)**2
       hess(thr,thr) = 10+48*(x2-2*x3)**2
       hess(thr,for) = -10

       hess(for,one) = -120*(x1-x4)**2
       hess(for,two) = 0
       hess(for,thr) = -10
       hess(for,for) = 10+120*(x1-x4)**2



    end do

  end function hessian

end module analytics


program powell
  use size
  use analytics
  real :: x0(n)
  integer :: i

  print *, "start?"
 
  do i = 1,n/4
     x0(4*i-3) = 3
     x0(4*i-2) = -1
     x0(4*i-1) = 0
     x0(4*i)   = 1
  end do
  
  print *, x0
  print *, fn(x0)
  print *, grad(x0)

  call minimize(x0,fn,grad,hessian)
end program powell

