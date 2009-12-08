subroutine minimize(x0,fn,grad,hessian)
  use optimization
  implicit none
  real, intent(in) :: x0(n)

  optional :: hessian

  interface
     real function fn(p)
       use size
       real :: p(n)
     end function fn
     function grad(p) Result(del)
       use size
       real :: p(n)
       real :: del(n)
     end function grad
     function hessian(p) result(hess)
       use size
       real:: hess(n,n)
       real :: p(n)
     end function hessian


     subroutine backtrackinglinesearch(x0,p,x,fn,grad,hessian)
       real :: p(:),x0(:)
       real :: x(:)
       optional :: hessian
       interface
          real function fn(p)
            use size
            real :: p(n)
          end function fn
          function grad(p) Result(del)
            use size
            real :: p(n)
            real :: del(n)
          end function grad
          function hessian(p) result(hess)
            use size
            real:: hess(n,n)
            real :: p(n)
          end function hessian
       end interface
     end subroutine backtrackinglinesearch
     function ffdhessf(x0) result(H)
       use optimization
       use analytics
       implicit none
       
       real, intent(in) :: x0(n)
       
       real :: H(n,n)
     end function ffdhessf
     function FFDGRAD(x0) result(g)
       use optimization
       use analytics
       implicit none
       
       real, intent(in) :: x0(n)
       real :: g(n)
     end function FFDGRAD
     function FFDHESSG(x0) result(H)
       use optimization
       use analytics
       implicit none
       
       real, intent(in) :: x0(n)
       
       real :: H(n,n)
     end function FFDHESSG
  end interface

  real :: Scaling(n)
  real :: x(n)
  real :: p(n)
  real :: nextx(n)
  real :: nextf
  real :: Sn(n)
  real :: L(n,n)
  real :: delta=-1
  real :: temphess(n,n)

  logical :: dogleg = .TRUE.

  ! logical :: DONE = .FALSE.


  Scaling =1 

  call initalize(n)

  print *, n
  print *, x0
  print 10, fn(x0)
  print *, "grad"
  print 10, grad(x0)
  print *, "hess"
  print 10, hessian(x0)

  print *, "ffdgrad"
  print 10, FFDGRAD(x0)
  print *, "fdhessg"
  print 10, FFDHESSG(x0)
  print *, "FFDHESSF"
  print 10, ffdhessf(x0)
  
  call exit(1)
  

  call UMSTOP0(x0,fn(x0),grad(x0),Scaling)

  x = x0


    
  ! call UMSTOP0(x0,fn(x0),grad(x0),Sx,consecmax)

  do iterations=0,maxiterations
     
     call takestep()

     call UMSTOP(nextx,x,grad(nextx),Scaling)
     x = nextx
     if(termcode .gt. 0) then
        print *, "terminating with code",termcode
        exit
     end if
  end do

  print *, "reached max iterations",iterations
  print *, nextx,x

10  format(4(g20.8))
  contains
    subroutine takestep()

      if(dogleg) then
         temphess = hessian(x)
         call modelhess(Scaling,temphess,L)
         call cholsolve(grad(x),L,Sn)
         call dogdriver(x,fn(x),grad(x),fn,L,Sn,delta,nextx,nextf)
      else
         call backtrackinglinesearch(x,p,nextx,fn,grad,hessian)
         ! call backtrackinglinesearch(x,p,nextx,fn,grad)
      end if
      ! call backtrackinglinesearch(x,p,xstep,fn,grad)
    end subroutine takestep


end subroutine minimize

