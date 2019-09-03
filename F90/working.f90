program sim
  implicit none
  integer :: j
  real :: i
  ! real, dimension(6) :: y
  real, dimension(2) :: y, f
  real, parameter :: PI = 3.1415927
  ! y = (/0.0,PI/2,0.0,9.0,0.0,0.0/)
  y = (/PI, PI/)
  ! do j = 1, 6
  !    print *, y(j)
  ! end do
  f = statespace(y)
  print *, f

contains
! Writing the statespace function
function statespace(y)
  implicit none
  real, dimension(2), intent(in) :: y
  real, dimension(2) :: statespace
  ! real, dimension(2), intent(out) :: statespace
  print *, y
  statespace = (/cos(y(1)), -cos(y(1))/)
  ! statespace = 1.00
  ! return 1
  ! statespace = (/cos(y(1)), cos(y(1)), cos(y(1)), cos(y(1)), cos(y(1)), cos(y(1))/)
end function statespace
end program sim
