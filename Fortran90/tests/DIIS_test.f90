PROGRAM DIIS_test

    USE DIIS
    USE OUTPUT

    IMPLICIT NONE

    REAL*8, dimension(3,3) :: A
    REAL*8, dimension(3,3) :: D
    REAL*8, dimension(3*3) :: B
    REAL*8, dimension(3*3) :: E
    REAL*8, dimension(3) :: w
    REAL*8 :: c
    REAL*8, allocatable, dimension(:,:,:) :: Mlist
    REAL*8, allocatable, dimension(:,:) :: Elist

    A(1,:) = (/1.0D0, 2.0D0, 3.0D0/)
    A(2,:) = (/4.0D0, 5.0D0, 6.0D0/)
    A(3,:) = (/7.0D0, 8.0D0, 9.0D0/)

    CALL ravel(3,A,B)

    CALL print_real_matrix(1,3*3,B)

    D(1,:) = (/3.0D0, 2.0D0, 7.0D0/)
    D(2,:) = (/4.0D0, 4.0D0, 6.0D0/)
    D(3,:) = (/6.0D0, 8.0D0, 9.0D0/)

    CALL DIIS_error(3,D,A,A,A,B,c)

    CALL print_real_matrix(1,3*3,B)

    WRITE(*,*) c

    CALL addFock(3,1,Mlist,A)
    CALL addFock(3,2,Mlist,D)

    CALL print_real_matrix(3,3,Mlist(1,:,:))
    CALL print_real_matrix(3,3,Mlist(2,:,:))

    CALL ravel(3,D,E)

    CALL addError(3,1,Elist,B/1000)
    CALL addError(3,2,Elist,E)

    CALL print_real_matrix(1,3*3,Elist(1,:))
    CALL print_real_matrix(1,3*3,Elist(2,:))

    CALL DIIS_B(Elist,2,B)

    CALL print_real_matrix(3,3,B)

    CALL DIIS_weigts(3,B,w)

    CALL print_real_matrix(1,3,w)



END PROGRAM DIIS_test
