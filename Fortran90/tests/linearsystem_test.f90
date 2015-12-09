PROGRAM ls
    ! --------------------------------------
    ! Solve linear system of equations AX=B.
    ! --------------------------------------
    !
    ! Results compared with MATLAB(R)
    !
    ! --------------------------------------

    USE LA
    USE OUTPUT

    IMPLICIT NONE

    INTEGER, PARAMETER :: d = 3

    REAL*8, dimension(d,d) :: A
    REAL*8, dimension(d) :: b
    REAL*8, dimension(d) :: x

    A(1,:) = (/1.0D0, 0.0D0, 0.0D0/)
    A(2,:) = (/0.0D0, 2.0D0, 0.0D0/)
    A(3,:) = (/0.0D0, 0.0D0, 3.0D0/)

    b = (/1.0D0,2.0D0,3.0D0/)

    CALL LINEAR_SYSTEM(d,A,b,x)

    CALL print_real_matrix(1,d,x)

    A(1,:) = (/1.0D0, 2.0D0, 3.0D0/)
    A(2,:) = (/-2.0D0, 3.0D0, -5.0D0/)
    A(3,:) = (/0.0D0, 7.0D0, -4.0D0/)

    b = (/-2.0D0,5.0D0,6.0D0/)

    CALL LINEAR_SYSTEM(d,A,b,x)

    CALL print_real_matrix(1,d,x)

    A(1,:) = (/-9.0D0, 2.0D0, 3.0D0/)
    A(2,:) = (/-2.0D0, -3.0D0, -5.0D0/)
    A(3,:) = (/0.0D0, 7.0D0, -4.0D0/)

    b = (/-2.0D0,-4.0D0,6.0D0/)

    CALL LINEAR_SYSTEM(d,A,b,x)

    CALL print_real_matrix(1,d,x)

END PROGRAM ls
