program rand_test
  use,intrinsic :: ISO_Fortran_env
  real(kind=4)  :: r(50)
  integer       :: i(50), n, j(1)

  call random_list()
  call random_list()
  call random_list()


contains
subroutine random_list()

  call random_number(r)

  ! Uniform distribution requires floor: Thanks to @francescalus 
  i = floor( r*300.0 )

  print *, i
  print *, sum(i)/float(50)
end subroutine
end program
