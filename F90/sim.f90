program sim
  ! use csv_file
  implicit none
  real, dimension(6) :: y, f, k1, k2, k3, k4, yt, yout
  real, parameter :: PI = 3.1415927
  real :: h, hh, h6
  integer :: n, j
  h = 1.0/100.0
  hh = h*0.5
  h6 = h/6.0
  y = (/0.0,PI/2,0.0,9.0,0.0,0.0/)
  n = 6
  j = 1000
  ! Opening a csv file
  ! open(unit=1,file='data.csv', status='unknown')
  ! open (unit = 2, file = "data.txt")
  ! Beginning the for loop
  do j = 1, 1000
    ! print *, y
    ! call csv_write(1, y, .true.)
    ! write(2, *) y(6)
    k1 = statespace(y)
    yt = y+hh*k1
    k2 = statespace(yt)
    yt = y+hh*k2
    k3 = statespace(yt)
    yt = y+h*k3
    k3 = k3+k2
    k4 = statespace(yt)
    yout = y+h6*(k1+2*k3+k4)
    y = yout
    ! print *, yout
  end do

contains
! Writing the statespace function
function statespace(y)
  implicit none
  real, dimension(6), intent(in) :: y
  real, dimension(6) :: statespace
  real :: cp, sp, cs, ss, dt, dp, ds
  cp = cos(y(2))
  sp = sin(y(2))
  cs = cos(y(3))
  ss = sin(y(3))
  dt = y(4)
  dp = y(5)
  ds = y(6)

  statespace(1) = dt
  statespace(2) = dp
  statespace(3) = ds
  statespace(4) = (-5*cs*dt*(ds + cp*dt)*ss)/3. - (dp*(1/sp)*(ds*(1 + 5*cs**2 - 5*ss**2) + cp*dt*(13 + 5*cs**2 - 5*ss**2)))/6.
  statespace(5) = (ds*(10*cs*dp*ss + dt*sp*(1 - 5*cs**2 + 5*ss**2)) + cp*dt*(10*cs*dp*ss + dt*sp*(7 - 5*cs**2 + 5*ss**2)))/6.
  statespace(6) = ((1/sp)*(-60*cs*sp*ss*dp**2 + 10*cs*sp*ss*dt**2*(5 + sp**2) + 2*cp*ds*(dp + &
  10*cs*dt*sp*ss + 5*dp*cs**2 - 5*dp*ss**2) + dp*dt*(19 + 5*cs**2*(7 + 5*sp**2) - 35*ss**2 -&
  sp**2*(7 + 25*ss**2)) + dt*cp**2*(-30*cs*dt*sp*ss + dp*(7 - 25*cs**2 + 25*ss**2))))/12.

end function statespace

end program sim
