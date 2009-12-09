subroutine dogdriver(xc,fc,fn,grad,L,Sn,maxstep,delta,nextx,nextf)
  use optimization
  implicit none
  
  real,intent(in) :: xc(n)
  real,intent(in) :: fc
  interface
     real function fn(p)
       use size
       real :: p(n)
     end function fn
  end interface
  
  real,intent(in) :: grad(n)
  real,intent(in) :: L(n,n)
  real,intent(in) :: sn(n)
  real,intent(in) ::maxstep
  real,intent(inout) :: delta
  
  real,intent(out) :: nextx(n)
  real,intent(out) :: nextf

  logical :: firstdog

  real :: newtlen

  real :: cauchylen = -1

  real :: eta
  real :: ssd(n)
  real :: v(n)
  real :: step(n)

  real :: xprev
  real :: fprev

  logical :: newttaken
  logical :: done

  maxtaken = .FALSE.
  retcode = 0

  firstdog = .TRUE.
  
  newtlen = norm(Sn)
  
  do while( .not. done)
     call dogstep(grad,L,Sn,newtlen,maxstep,delta,firstdog,cauchylen,eta,ssd,v,step,newttaken)
     print *, "finished step", step
     call trustregup(xc,fc,fn,grad,L,step,newttaken,maxstep,delta,xprev,fprev,nextx,nextf)
     if(retcode .lt. 2) then
        done = .TRUE.
     end if

  end do

  print *, "RETCODE AFTER DOG",retcode

end subroutine dogdriver

subroutine dogstep(grad,L,Sn,newtlen,maxstep,delta,firstdog,cauchylen,eta,ssd,v,step,newttaken)
  use optimization
  implicit none
  
  real,intent(in) :: grad(n)
  real,intent(in) :: L(n,n)
  real,intent(in) :: Sn(n)
  real,intent(in) :: newtlen
  real,intent(in) :: maxstep
  
  real,intent(inout) :: delta
  logical,intent(inout) :: firstdog
  real, intent(inout) :: cauchylen
  real, intent(inout) :: eta
  real, intent(inout) :: ssd(n)
  real, intent(inout) :: v(n)
  real,intent(out) :: step(n)
  logical,intent(out) :: newttaken

  real :: alpha,beta
  real :: temp,tempv
  real :: lambda
  integer :: i,j



  if(newtlen .le. delta) then
     ! NEWTONSTEP
     print *, "NEWTONSTEP"
     newttaken = .TRUE.
     step = Sn
     delta = newtlen
     return
  end if

  ! Not newton step
  newttaken = .FALSE.
  if(firstdog) then
     ! First dog
     ! COMPUTE DOG LEG
     print *, "computing dog leg curve"
10   format(4(g20.8))

     firstdog = .FALSE.
     alpha = norm(grad)
     print *, "alpha",alpha
     print 10, "grad",grad
     print *, "L"
     print 10, L
     beta = 0.0
     do i =1,n
        temp = 0.0
        do j=i,n
           temp = temp + L(j,i)*grad(j)
        end do
        beta = beta + temp*temp
        print *, "beta..",beta,temp
     end do
     ssd = -(alpha/beta)*grad
     cauchylen = alpha*sqrt(alpha) / beta
     eta = 0.2+(0.8* alpha**2 / (beta * abs(dot_product(grad,Sn)))) 
     print *, "eta should be <1",eta
     print *, "beta",beta
     print *, "grad",grad
     print *, "alpha",alpha
     print *, "Sn",Sn
     print *, "grad*Sn",abs(dot_product(grad,Sn))*beta
     print *, 0.2+(0.8* alpha**2 / (beta * abs(dot_product(grad,Sn)))) 
     if (eta .gt. 1) call exit(1)
     v = eta*Sn-ssd
     if(delta .eq. -1.0) then
        delta = min(cauchylen,maxstep)
     end if
  end if
  
  if (eta * newtlen .le. delta) then
     ! Partial newton
     print *, "partial newton"
     step = (delta/newtlen) * Sn
     return
  else if (cauchylen .ge. delta) then
     ! Steepest descent
     print *, "stepest decent"
     print *, delta
     print *, cauchylen
     print *, ssd
     step = (delta/cauchylen) * ssd
     return
  else
     ! Combo
     print *, "combostep"
     temp = dot_product(v,ssd)
     tempv = dot_product(v,v)
     lambda = (-temp + sqrt(temp**2-tempv*(cauchylen**2-delta**2)) )* tempv
     step = ssd + lambda*v
     return
  end if
  print *, "shouldnt happen"

end subroutine dogstep

subroutine trustregup(xc,fc,fn,grad,L,step,newttaken,maxstep,delta,xprev,fprev,nextx,nextf)
  use optimization
  implicit none
  
  real,intent(in) :: xc(n)
  real,intent(in) :: fc
  interface
     real function fn(p)
       use size
       real :: p(n)
     end function fn
  end interface
  real,intent(in)  :: grad(n)
  real,intent(in)  :: L(n,n)
  real,intent(in)  :: step(n)
  logical,intent(in) :: newttaken
  real,intent(in)  :: maxstep
  ! real,intent(in) :: steptol
  
  real,intent(inout) :: delta
  real,intent(inout) :: xprev(n)
  real,intent(inout) :: fprev
  
  real,intent(out) :: nextx(n)
  real,intent(out) :: nextf
  ! logical,intent(out) :: maxtaken

  real :: alpha
  real :: initslope
  real :: temp
  real :: dfpred
  real :: df
  real :: steplen,rellength
  real :: deltatemp
  integer :: i,j


  maxtaken = .FALSE.
  alpha = 1.0e-4

  steplen = norm(step)

  nextx = xc+step
  print *, "stuff"
  print *, fc
  print *, fn(xc)
  print *, "xc",xc
  print *, step
  print *, xc+step
  print *, "nextx",nextx

  print *, fn(nextx)
  print *, "funccal?", fc, fn(xc-step)
  nextf = fn(nextx)
  print *, "funcall!",nextf
  df = nextf-fc
  print *, "df should be negatve",df

  initslope = dot_product(grad,step)
  
  if(retcode .ne. 3) then
     fprev = 0
  end if

  if(retcode .eq.3 .and. (nextf .ge. fprev .OR. df .gt. alpha * initslope)) then
     ! RESET and terminate step
     print *, "reset and terminate step"
     retcode = 0
     nextx = xprev
     nextf = fprev
     delta = delta/2
     return
  
  else if (df .gt. alpha*initslope) then
     rellength = maxval( (abs(step)) / abs(nextx) )
     
     if (rellength .lt. steptol) then
        ! nextx-xc too small, terminate
        print *, " nextx-xc too small, terminate"
        retcode = 1
        nextx = xc
        return
     else
        ! reduce detla, continue
        print *, "reduce detla, continue"
        retcode = 2
        deltatemp = (-initslope * steplen) / (2.0*(df-initslope))
        
        if (deltatemp .lt. 0.1 * delta) then
           delta = 0.1 * delta
        else if (deltatemp .gt. 0.5 * delta) then
           delta = 0.5 *delta
        else
           delta = deltatemp
        end if
           
     end if
  else
     ! f(nextx) small enough
     print *, " f(nextx) small enough"

     dfpred = initslope
     
     temp = 0.0
     do i=1,n
        do j=1,n
           temp = temp + L(j,i)*step(j)
        end do
        dfpred = dfpred + (temp * temp/2.0)
     end do

     
     if (retcode .ne. 2 .and. ( abs(dfpred-df) .ge. 0.1 * abs(df) .or. df .le. initslope) .and. newttaken .eqv. .FALSE. .and. delta .le. 0.99 * maxstep) then
        ! Double delta and continue
        print *, " Double delta and continue"
        
        retcode = 3
        xprev = nextx
        fprev = nextf
        delta = (min(2*delta,maxstep))
        return
     else
        ! accept new iterate choose new delta
        print *, " accept new iterate choose new delta"
        
        retcode = 0
        if (steplen .gt. 0.99 * maxstep) then
           maxtaken = .true.
        end if
        if(df .ge. 0.1 * dfpred) then
           ! decrease delta
           print *, "decrease delta"
           delta = delta /2
        else if (df .le. 0.75 * dfpred) then
           ! increase delta
           print *, " increase delta"
           delta = min(2*delta,maxstep)
        else
           print *, "delta unchanged"
        end if
     end if
  end if

end subroutine trustregup
