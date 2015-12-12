PROGRAM DIIS_test

    USE DIIS
    USE OUTPUT

    IMPLICIT NONE

    REAL*8, dimension(3,3) :: A
    REAL*8, dimension(3*3) :: B
    REAL*8, dimension(3,3) :: C
    REAL*8, dimension(3,3) :: D
    REAL*8, dimension(3*3) :: E
    REAL*8, allocatable, dimension(:,:) :: F
    REAL*8 :: maxerr
    REAL*8, allocatable, dimension(:) :: w
    REAL*8, allocatable, dimension(:,:) :: Elist
    REAL*8, allocatable, dimension(:,:,:) :: Mlist

    INTEGER :: info, dim, i, j

    A(1,:) = (/1.0D0, 2.0D0, 3.0D0/)
    A(2,:) = (/4.0D0, 5.0D0, 6.0D0/)
    A(3,:) = (/7.0D0, 8.0D0, 9.0D0/)

    WRITE(*,*) "Matrix A:"
    CALL print_real_matrix(3,3,A)

    CALL ravel(3,A,B)

    WRITE(*,*)
    WRITE(*,*) "Vector B (ravel(A)):"

    CALL print_real_vector(3*3,B)

    A(1,:) = (/0.0D0, 1.0D0, 0.0D0/)
    A(2,:) = (/-1.0D0, 0.0D0, 0.0D0/)
    A(3,:) = (/0.0D0, 0.0D0, 1.0D0/)

    C(1,:) = (/1.0D0, 0.0D0, 0.0D0/)
    C(2,:) = (/0.0D0, 0.0D0, 1.0D0/)
    C(3,:) = (/0.0D0, -1.0D0, 0.0D0/)

    D(1,:) = (/1.0D0, 0.0D0, 0.0D0/)
    D(2,:) = (/0.0D0, 1.0D0, 0.0D0/)
    D(3,:) = (/0.0D0, 0.0D0, 1.0D0/)

    CALL DIIS_error(3,A,C,D,D,E,maxerr)

    WRITE(*,*)
    WRITE(*,*) "Error matrix:"

    CALL print_real_matrix(3,3,MATMUL(A,MATMUL(C,D)) - MATMUL(D,MATMUL(C,A)))

    WRITE(*,*)
    WRITE(*,*) "Error E:"

    CALL print_real_vector(3*3,E)

    WRITE(*,*)
    WRITE(*,*) "Maxerr: ", maxerr

    CALL addError(3,1,Elist,E)
    CALL addError(3,2,Elist,E)
    CALL addError(3,3,Elist,E)

    WRITE(*,*)
    WRITE(*,*) "Elist(1):"

    CALL print_real_vector(3*3,Elist(1,:))

    WRITE(*,*)
    WRITE(*,*) "Elist(2):"

    CALL print_real_vector(3*3,Elist(2,:))

    WRITE(*,*)
    WRITE(*,*) "Elist(3):"

    CALL print_real_vector(3*3,Elist(3,:))

    WRITE(*,*)
    WRITE(*,*) "Elist:"

    CALL print_real_matrix(3,3*3,Elist)

    ALLOCATE(F(4,4),w(3))

    CALL DIIS_B(3,elist,3,F)

    WRITE(*,*)
    WRITE(*,*) "Matrix F:"

    CALL print_real_matrix(4,4,F)

    dim = 4
    info = -1

    DO WHILE(info .NE. 0)

        CALL DIIS_weigts(dim,F,w,info)

        IF (info .NE. 0) THEN
            WRITE(*,*)
            WRITE(*,*) "IMPOSSIBLE TO SOLVE THE LINEAR SYSTEM: REDUCING MATRIX F"

            CALL DIIS_reduce_B(dim,F)

            dim = dim - 1

            WRITE(*,*)
            WRITE(*,*) "Reduced matrix F:"

            CALL print_real_matrix(dim,dim,F)

            DEALLOCATE(w)
            ALLOCATE(w(dim-1))
        END IF

    END DO

    WRITE(*,*)
    WRITE(*,*) "Weights:"

    CALL print_real_vector(dim-1,w)

    CALL addFock(3,1,Mlist,A)
    CALL addFock(3,2,Mlist,2*A)
    CALL addFock(3,3,Mlist,3*A)

    j = 3

    A(:,:) = 0.0D0

    DO i = 1,dim-1

        A = A + w(i) * Mlist(j,:,:)

        j = j - 1

    END DO

    WRITE(*,*)
    WRITE(*,*) "New operator A:"

    CALL print_real_matrix(3,3,A)

END PROGRAM DIIS_test
