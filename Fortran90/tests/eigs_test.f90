PROGRAM eigs_test

    USE UTILS, only: eigs, print_real_matrix

    INTEGER, PARAMETER :: d = 2

    REAL*8, dimension(d,d) :: A
    REAL*8, dimension(d,d) :: V
    REAL*8, dimension(d) :: l

    A(1,:) = (/-1.53866, -0.515838/)
    A(2,:) = (/-0.515838, -2.43973/)

    CALL EIGS(d,A,V,l)

    WRITE(*,*) "Matrix A"
    CALL print_real_matrix(d,d,A)

    WRITE(*,*) "Matrix V"
    CALL print_real_matrix(d,d,V)

    WRITE(*,*) "Vector l"
    CALL print_real_matrix(d,1,l)

END PROGRAM eigs_test
