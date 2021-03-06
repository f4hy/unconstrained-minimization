subroutine minimize(x0,fn,grad,hessian)
  use optimization
  use size
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
            real :: p(n)
            real:: hess(n,n)
          end function hessian
       end interface
     end subroutine backtrackinglinesearch
     function ffdhessf(x0) result(H)
       use optimization
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
  real :: fc
  real :: gc(n)
  real :: maxstep = 1.0

  Scaling =1 

  call initalize()

  ! print *, n
  ! print *, x0
  ! print 10, fn(x0)
  ! print *, "grad"
  ! print 10, grad(x0)
  ! print *, "hess"
  ! print 10, hessian(x0)

  ! print *, "ffdgrad"
  ! print 10, FFDGRAD(x0)
  ! print *, "fdhessg"
  ! print 10, FFDHESSG(x0)
  ! print *, "FFDHESSF"
  ! print 10, ffdhessf(x0)

  ! call exit(1)

  
  call UMSTOP0(x0,fn(x0),grad(x0))

  if(termcode .gt. 0) then
     print *, "terminating before starting with code",termcode
     call exit(1)
  end if

  x = x0

  write(10,*) x(1),",",x(2),",", fn(x)
  do iterations=0,maxiterations
 
     call takestep()

     call UMSTOP(nextx,x,fn(x),grad(nextx))

     x = nextx

     write(7,*) fn(x)
     write(10,*) x(1),",",x(2),",", fn(x)
     if(termcode .gt. 0) then
        ! print *, "terminating with code",termcode
        exit
     end if
  end do

  print *, "stopped after iterations",iterations
  print *, "final point"
  print 10, x
  print *, "Function value at final point"
  print 10, fn(x)

open(unit=11,file="answer.out")
  write(11,*) "terminated with code",termcode
  write(11,*) "stopped after iterations",iterations
  write(11,*) "final point"
  write(11,10)  x
  write(11,*) "Function value at final point"
  write(11,10) fn(x)
  

  close(7)
  close(8)
  close(9)
  close(10)
  close(11)
  
10 format(4(g15.8))
contains
  subroutine takestep()
10  format(4(g20.8))

    integer :: i,j

    if(nohessian) then
       if(analyticgrad) then
          ! print *, "lineseach stepest decent"
          call backtrackinglinesearch(x,p,nextx,fn,grad)
          return
       else
          ! print *, "lineseach stepest decent fake grad"
          call backtrackinglinesearch(x,p,nextx,fn,ffdgrad)
       end if
    else if(linesearch) then
       if(analytichessian) then
          ! print *, "lineseach newton"
          call backtrackinglinesearch(x,p,nextx,fn,grad,hessian)
       else if (analyticgrad) then
          ! print *, "lineseach newton fake hess"
          call backtrackinglinesearch(x,p,nextx,fn,grad,ffdhessg)
       else
          ! print *, "lineseach newton fake grad+hess"
          call backtrackinglinesearch(x,p,nextx,fn,ffdgrad,ffdhessf)
       end if

    else if (dogleg) then
       fc = fn(x)
       if(analytichessian) then
          ! print *, "dog leg"
          temphess = hessian(x)
          gc = grad(x)
       else if(analyticgrad) then
          ! print *, "dog leg fake hess"
          temphess = ffdhessg(x)
          gc = grad(x)
       else 
          ! print *, "dog leg fake grad+hess"
          temphess = ffdhessf(x)
          gc = ffdgrad(x)
       end if
          call modelhess(Scaling,temphess,L)
          call cholsolve(-gc,L,Sn)
          do i=1,n
             do j=i+1,n
                L(i,j) = L(j,i)
             end do
          end do
          ! print *, "grad"
          ! print 10, gc
          ! print *, "fullhess"
          ! print 10, hessian(x)
          ! print *, "hess"
          ! print 10, temphess
          ! print *, "L"
          ! print 10, L
          ! print *, "Sn"
          ! print 10, Sn
          ! call exit(1)
          call dogdriver(x ,fc,fn,gc,L,-Sn,maxstep,delta,nextx,nextf)
    else
       print *, "ERROR: NO METHOD SELECTED!"
       call exit(1)
    end if
    ! Newton
  end subroutine takestep
end subroutine minimize



