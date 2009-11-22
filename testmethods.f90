  ! Unit tests for methods.f90
  
program unittests

  call testmacheps()

  call testcholsolve()

  ! call testlsolve()

  ! call testltsolve()
  ! 

  call testcholdecomp()

  print *, "tests ran succuessfully"
end program unittests

subroutine testmacheps()
  use optimization

  call UMINICK(1)

  if (1 + macheps/2 .NE. 1) then
     print *,  "macheps failed", macheps
     call exit(1)
  end if
  print *, "macheps passed", macheps
end subroutine testmacheps



subroutine testcholsolve()
  real :: A(2,2)
  real :: b(2)
  real :: y(2)

  real :: C(3,3)
  real :: d(3)
  real :: w(3)

  A = reshape( (/ 1,0,3,4 /), (/2 , 2/) )
  b = (/4,4/)

  call UMINICK(2)

  call CHOLSOLVE(b,A,y)


  if(y(1)  .ne. -4.0 .or. y(2) .ne. -1.0/4.0) then
     print *, "cholsolve failed"
     call exit(1)
  end if

  C = reshape( (/ 1,0,0,2,3,0,4,5,6 /), (/3,3/) )
  d = (/18,18,18/)

  call UMINICK(3)

  call CHOLSOLVE(d,C,w)


  if(w(1) .ne. -18 .or. w(2) .ne. -2 .or. w(3) .ne. -1.0/2.0) then
     print *, "cholsolve failed"
     call exit(1)
  end if

  
  print *, "cholsolve passed"

end subroutine testcholsolve

subroutine testcholdecomp()
  real :: A(2,2)
  real :: L(2,2)
  real :: maxadd
  A = reshape( (/ 1,0,0,1 /), (/2,2/) )
  
  call UMINICK(2)

  call choldecomp(A,L,maxadd)
  
  if(L(1,1) .ne. 1 .or. L(1,2).ne. 0 .or. L(2,1) .ne. 0 .or. L(2,2) .ne. 1) then
     print *, "choldecomp failed"
     call exit(1)
  end if


  print *, "choldecomp passed"
end subroutine testcholdecomp












! subroutine testlsolve()
!   real :: A(2,2)
!   real :: b(2)
!   real :: y(2)

!   real :: C(3,3)
!   real :: d(3)
!   real :: w(3)


!   A = reshape( (/ 1,0,3,4 /), (/2 , 2/) )
!   b = (/2,2/)



!   call UMINICK(2)

!   call Lsolve(b,A,y)

!   if(y(1) .ne. 2 .or. y(2) .ne. -1) then
!      print *, "lsolve failed"
!      call exit(1)
!   end if

!   C = reshape( (/ 1,0,0,2,3,0,4,5,6 /), (/3,3/) )
!   d = (/9,9,9/)

!   call UMINICK(3)

!   call Lsolve(d,C,w)

!   if(w(1) .ne. 9 .or. w(2) .ne. -3 .or. w(3) .ne. -2) then
!      print *, "ltsolve failed"
!      call exit(1)
!   end if

  
!   print *, "lsolve passed"

! end subroutine testlsolve

! subroutine testltsolve()
!   real :: A(2,2)
!   real :: b(2)
!   real :: y(2)

!   real :: C(3,3)
!   real :: d(3)
!   real :: w(3)

!   A = reshape( (/ 1,0,3,4 /), (/2 , 2/) )
!   b = (/4,4/)

!   call UMINICK(2)

!   call Ltsolve(b,A,y)


!   if(y(1) .ne. 1 .or. y(2) .ne. 1) then
!      print *, "ltsolve failed"
!      call exit(1)
!   end if
  

!   C = reshape( (/ 1,0,0,2,3,0,4,5,6 /), (/3,3/) )
!   d = 18

!   call UMINICK(3)

!   call Ltsolve(d,C,w)

!   if(w(1) .ne. 4 .or. w(2) .ne. 1 .or. w(3) .ne. 3) then
!      print *, "ltsolve failed"
!      call exit(1)
!   end if

  
!   print *, "ltsolve passed"

! end subroutine testltsolve

