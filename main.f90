subroutine minimize(x0,fn,grad,hessian)
  use optimization
  real, intent(in) :: x0(n)

  optional :: hessian

  interface
     real function fn(p)
       real, dimension(size(x0)) :: p
     end function fn
     function grad(p) Result(del)
       real, dimension(size(x0)) :: p
       real, dimension(size(x0)) :: del
     end function grad
     function hessian(p) result(hess)
       real, dimension(size(x0),size(x0)) :: hess
       real, dimension(size(x0))  :: p
     end function hessian


     subroutine backtrackinglinesearch(x0,p,x,f,grad,hessian)
       real :: p(:),x0(:)
       real :: x(:)
       optional :: hessian
       interface
          real function f(p)
            real, dimension(size(x0)) :: p
          end function f
          function grad(p) Result(del)
            real, dimension(size(x0)) :: p
            real, dimension(size(x0)) :: del
          end function grad
          function hessian(p) result(hess)
            real, dimension(size(x0),size(x0)) :: hess
            real, dimension(size(x0))  :: p
          end function hessian
       end interface
     end subroutine backtrackinglinesearch
  end interface

  real :: Scaling(n)
  real :: x(n)
  real :: p(n)
  real :: xstep(n)
  

  ! logical :: DONE = .FALSE.


  Scaling =1 

  call initalize(n)

  call UMSTOP0(x0,fn(x0),grad(x0),Scaling)

  x = x0


    
  ! call UMSTOP0(x0,fn(x0),grad(x0),Sx,consecmax)

  do iteration=0,maxiterations
     call takestep()

     call UMSTOP(xstep,x,fn(xstep),grad(xstep),Scaling)
     x = xstep
     if(termcode .gt. 0) then
        print *, "terminating with code",termcode
        exit
     end if
  end do

  print *, "reached max iterations"
  print *, iteration

  contains
    subroutine takestep()

      call backtrackinglinesearch(x,p,xstep,fn,grad,hessian)
      ! call backtrackinglinesearch(x,p,xstep,fn,grad)
    end subroutine takestep

  
end subroutine minimize

