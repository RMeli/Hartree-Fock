PROGRAM eigs_test

    USE UTILS, only: eigs, print_real_matrix

    INTEGER, PARAMETER :: d = 3

    REAL*8, dimension(d,d) :: A
    REAL*8, dimension(d,d) :: V
    REAL*8, dimension(d) :: l

    A(1,:) = (/3.0D0, 2.0D0, 1.5D0/)
    A(2,:) = (/2.0D0, -6.0D0, 3.0D0/)
    A(3,:) = (/1.5D0, 5.0D0, 5.0D0/)

    CALL EIGS(d,A,V,l)

    WRITE(*,*) "Matrix A"
    CALL print_real_matrix(d,d,A)

    WRITE(*,*) "Matrix V"
    CALL print_real_matrix(d,d,V)

    WRITE(*,*) "Vector l"
    CALL print_real_matrix(d,1,l)

END PROGRAM eigs_test
