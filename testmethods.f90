  ! Unit tests for methods.f90
  
program unittests

  call testmacheps()

  print *, "tests ran succuessfully"
end program unittests

subroutine testmacheps()
  real, external :: macheps
  if (1 + macheps()/2 .NE. 1) then
     print *,  "macheps failed", macheps()
     call exit(1)
  end if
  print *, "macheps passed"
end subroutine testmacheps

