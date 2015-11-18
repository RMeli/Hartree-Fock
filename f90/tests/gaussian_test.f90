PROGRAM gaussian_test

    USE GAUSSIAN

    IMPLICIT NONE

    REAL*8 :: aa = 1.0D0, bb = 2.0D0 ! Gaussian exponential coefficients
    REAL*8, dimension(3) :: Ra = (/0.0D0,0.0D0,0.0D0/), Rb = (/1.0D0,1.0D0,1.0D0/)

    REAL :: Na, Nb

    Na = norm(0,0,0,aa)
    Nb = norm(0,0,0,bb)

    WRITE(*,*) "Na = ", Na
    WRITE(*,*) "Nb = ", Nb

END PROGRAM gaussian_test
