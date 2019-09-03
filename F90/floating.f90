program test
REAL A
DATA A /.01/
print *, A

DATA X/77777/, Y/7/
Y1 = 1 / Y
Z = X / Y
Z1 = X * Y1
print *, Z
print *, Z1

! REAL B, C
! DATA C /1000.2/
! B = C - 1000.0
! PRINT *, B
end program test
