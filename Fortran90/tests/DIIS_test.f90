PROGRAM DIIS_test

    USE DIIS
    USE OUTPUT

    REAL*8, dimension(3,3) :: A
    REAL*8, dimension(3*3) :: B

    A(1,:) = (/1.0D0, 2.0D0, 3.0D0/)
    A(2,:) = (/4.0D0, 5.0D0, 6.0D0/)
    A(3,:) = (/7.0D0, 8.0D0, 9.0D0/)

    CALL ravel(3,A,B)

    CALL print_real_matrix(1,3*3,B)

END PROGRAM DIIS_test
