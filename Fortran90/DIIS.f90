! --------------------------------------------------------------------
!
! Copyright (C) 2015 Rocco Meli
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! ---------------------------------------------------------------------

MODULE DIIS
    ! ------------------------------------------------------
    ! DIIS algorithm.
    ! ------------------------------------------------------
    !
    ! Implement DIIS algorithm utilities in order to speedup
    ! SCF convergence
    !
    ! ------------------------------------------------------

    USE LA, only: LINEAR_SYSTEM
    USE OUTPUT

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE ravel(n,M,Mr)
        ! -------------------------------------
        ! Reshape n x n matrix M to n**2 array.
        ! -------------------------------------

        ! INPUT
        INTEGER, intent(in) :: n
        REAL*8, dimension(n,n) :: M

        ! INTERMEDIATE VARIABLES
        INTEGER :: i, j ! loop indices

        ! OUTPUT
        REAL*8, dimension(n*n) :: Mr

        DO i = 1, n
            DO j = 1, n
                Mr(j+n*(i-1)) = M(i,j)
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE DIIS_error(Kf,F,P,S,X,error,maxerror)
        ! ---------------------------------------
        ! Compute error vector
        ! ---------------------------------------
        !
        ! Source:
        !   P. Pulay
        !   Improved SCF Convergence Acceleration
        !   Journal of Computatuonal Chemistry
        !   1982
        !
        ! ---------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                       ! Basis set size
        REAL*8, dimension(Kf,Kf), intent(in) :: F       ! Fock matrix
        REAL*8, dimension(Kf,Kf), intent(in) :: P       ! Density matrix
        REAL*8, dimension(Kf,Kf), intent(in) :: S       ! Overlap matrix
        REAL*8, dimension(Kf,Kf), intent(in) :: X       ! Transformation matrix

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(Kf,Kf):: err                  ! Error matrix

        ! OUTPUT
        REAL*8, dimension(Kf*Kf), intent(out) :: error  ! Error vector
        REAL*8, intent(out) :: maxerror                 ! Maximal error

        err = MATMUL(F,MATMUL(P,S)) - MATMUL(S,MATMUL(P,F)) ! Compute error matrix

        err = MATMUL(TRANSPOSE(X),MATMUL(err,X)) ! Balanced error matrix

        CALL ravel(Kf,err,error) ! Create error vector

        maxerror = MAXVAL(ABS(error)) ! Maximal absolute error

    END SUBROUTINE DIIS_error



    SUBROUTINE DIIS_B(Kf,elist,steps,B)
        ! ----------------------------------------
        ! Compute matrix B of DIIS algorithm.
        ! ----------------------------------------
        !
        ! Source:
        !   T. Helgaker, P. Jørgensen and J. Olsen
        !   Molecular Electronic-Structure Theory
        !   Wiley
        !   2000
        !
        ! ----------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                               ! Basis set size
        INTEGER, intent(in) :: steps                            ! Current SCF iteration
        REAL*8, dimension(steps,Kf*Kf), intent(in) :: elist     ! List of errors

        ! INTERMEDIATE VARIABLES
        INTEGER :: i, j

        ! OUTPUT
        REAL*8, dimension(steps+1,steps+1), intent(out) :: B     ! Matrix B

        !CALL print_real_matrix(steps,Kf*Kf,elist)

        DO i = 1, steps ! Rows of B
            DO j = 1, steps ! Cols of B

                B(i,j) = DOT_PRODUCT(elist(i,:),elist(j,:))

            END DO
        END DO

        B(:,steps+1) = -1.0D0

        B(steps+1,:) = -1.0D0

        B(steps+1,steps+1) = 0.0D0

    END SUBROUTINE DIIS_B

    SUBROUTINE DIIS_weigts(dim,B,w,info)
        ! ---------------------------------------
        ! Compute DIIS weights.
        ! ---------------------------------------
        !
        ! Source:
        !   T. Helgaker, P. Jørgensen and J. Olsen
        !   Molecular Electronic-Structure Theory
        !   Wiley
        !   2000
        !
        ! ----------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: dim                          ! Dimension of b
        REAL*8, dimension(dim,dim), intent(in) :: B         ! Matrix B

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(dim) :: sol                       ! Solution of B*SOL=RHS (weight and Lagrande multiplier)
        REAL*8, dimension(dim) :: rhs                       ! Right hand side of the linear system
        REAL*8, dimension(dim,dim) :: BB        ! Copy of matrix B (side effects) (TODO: avoid copy)

        ! OUTPUT
        REAL*8, dimension(dim-1), intent(out) :: w          ! Weights (SOL without the Lagrange multiplier)
        INTEGER, intent(out) :: info                        ! Information about the solution of the linear system

        rhs(:) = 0.0D0

        rhs(dim) = -1.0D0

        !CALL print_real_matrix(step+1,step+1,B)

        BB = B

        CALL LINEAR_SYSTEM(dim,BB,rhs,sol,info) ! Solve linear system B*SOL=RHS

        w = sol(1:dim-1) ! Discard Lagrange multiplier

    END SUBROUTINE DIIS_weigts

    SUBROUTINE DIIS_reduce_B(dim,B)
        ! ---------------------------------------------------------------
        ! Reduce the matrix B by one in order to eliminate ill behaviour.
        ! ---------------------------------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: dim
        REAL*8, allocatable, dimension(:,:), intent(inout) :: B

        ! INTERMEDIATE VARIABLES
        REAL*8, allocatable, dimension(:,:) :: BB
        INTEGER :: i, j

        ALLOCATE(BB(dim-1,dim-1))

        DO i = 2, dim
            DO j = 2, dim
                BB(i-1,j-1) = B(i,j)
            END DO
        END DO

        DEALLOCATE(B)

        CALL MOVE_ALLOC(BB,B)

    END SUBROUTINE DIIS_reduce_B

    SUBROUTINE DIIS_Fock(Kf,step,F,P,S,X,Flist,elist)
        ! ---------------------------------------
        ! Compute the DIIS Fock matrix.
        ! ---------------------------------------
        !
        ! Source:
        !   P. Pulay
        !   Improved SCF Convergence Acceleration
        !   Journal of Computatuonal Chemistry
        !   1982
        !
        ! ---------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                           ! Basis set size
        INTEGER, intent(in) :: step                         ! Current SCF step
        REAL*8, dimension(Kf,Kf), intent(in) :: P           ! Current density matrix
        REAL*8, dimension(Kf,Kf), intent(in) :: S           ! Overlap matrix
        REAL*8, dimension(Kf,Kf), intent(in) :: X           ! Transformation matrix

        ! INPUT/OUTPUT
        REAL*8, allocatable, dimension(:,:,:), intent(inout) :: Flist   ! List of Fock matrices
        REAL*8, allocatable, dimension(:,:), intent(inout) :: elist     ! List of error vectors
        REAL*8, dimension(Kf,Kf), intent(inout) :: F                    ! Fock matrix

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(Kf*Kf) :: error                               ! Error vector at current SCF step
        REAL*8 :: maxerror                                              ! Maximal error at current iteration
        REAL*8, allocatable, dimension(:,:) :: B                        ! Matrix B
        REAL*8, allocatable, dimension(:) :: w                          ! Weights
        INTEGER :: i, j                                                 ! Loop index
        INTEGER :: info                                                 ! Information about the solution of the linear system
        INTEGER :: dim                                                  ! Actual dimension of the linear system

        ! Compute error at current iteration
        CALL DIIS_error(Kf,F,P,S,X,error,maxerror)

        WRITE(*,*) "#########"
        WRITE(*,*) "MAXERROR:", maxerror
        WRITE(*,*) "#########"

        ! Store Fock matrix and error vector
        CALL addFock(Kf,step,Flist,F)
        CALL addError(Kf,step,elist,error)

        info = -1
        dim = step+1

        ALLOCATE(B(dim,dim))
        ALLOCATE(w(dim-1))

        CALL DIIS_B(Kf,elist,step,B)

        DO WHILE (info .NE. 0)

            CALL DIIS_weigts(dim,B,w,info)

            IF (info .NE. 0) THEN ! Check the solution of the system
                WRITE(*,*)
                WRITE(*,*) "IMPOSSIBLE TO SOLVE THE LINEAR SYSTEM: REDUCING MATRIX B"

                CALL DIIS_reduce_B(dim,B) ! Reduce B in order to eliminate ill behaviour

                dim = dim - 1 ! Reduce dimensionality of the system

                ! Deallocate/allocate weight vector (TODO: useless/performance)
                DEALLOCATE(w)
                ALLOCATE(w(dim-1))
            END IF

        END DO

        WRITE(*,*)
        WRITE(*,*) "SOLVED DIIS SYSTEM."

        w = w(dim-1:1:-1)

        F(:,:) = 0.0D0 ! Erase current Fock matrix F (already stored in FLIST)

        ! Create a new Fock matrix according to the DIIS algorithm
        j = step

        DO i = 1,dim-1

            F = F + w(i) * Flist(j,:,:)

            j = j - 1

        END DO

        DEALLOCATE(B,w)

    END SUBROUTINE

    SUBROUTINE addFock(Kf,step,Flist,F)
        ! -------------------------------
        ! Add new Fock matrix F to FLIST.
        ! -------------------------------

        ! INPUT
        INTEGER, intent(in) :: Kf                                       ! Basis set size
        INTEGER, intent(in) :: step                                     ! Current SCF step
        REAL*8, dimension(Kf,Kf), intent(in) :: F                       ! Current Fock matrix

        ! INTERMEDIATE VARIABLES
        INTEGER :: i                                                    ! Loop index
        REAL*8, allocatable, dimension(:,:,:) :: cFlist                 ! Temporary copy of Fock list

        ! INPUT/OUTPUT
        REAL*8, allocatable, dimension(:,:,:), intent(inout) :: Flist   ! Fock matrix list

        IF (ALLOCATED(Flist)) THEN
            ALLOCATE(cFlist(step,Kf,Kf)) ! Allocate temporary list

            DO i = 1, step - 1
                cFlist(i,:,:) = Flist(i,:,:) ! Compy old Fock operators into new list (TODO: better solution?)
            END DO

            cFlist(step,:,:) = F ! Add the new element

            DEALLOCATE(Flist)

            CALL MOVE_ALLOC(cFlist,Flist)

        ELSE ! Flist not allocated

            ALLOCATE(Flist(1,Kf,Kf))

            Flist(1,:,:) = F

        END IF

    END SUBROUTINE

    SUBROUTINE addError(Kf,step,elist,e)
        ! -----------------------------------
        ! Append new error vector E to ELIST.
        ! -----------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                   ! Basis set size
        INTEGER, intent(in) :: step                 ! Current SCF step
        REAL*8, dimension(Kf*Kf), intent(in) :: e   ! Current error

        ! INTERMEDIATE VARIABLES
        INTEGER :: i                                                    ! Loop index
        REAL*8, allocatable, dimension(:,:) :: celist                 ! Temporary copy of error list

        ! INPUT/OUTPUT
        REAL*8, allocatable, dimension(:,:), intent(inout) :: elist   ! Error list

        IF (ALLOCATED(elist)) THEN
            ALLOCATE(celist(step,Kf*Kf)) ! Allocate temporary list

            DO i = 1, step - 1
                celist(i,:) = elist(i,:) ! Compy old Fock operators into new list (TODO: better solution?)
            END DO

            celist(step,:) = e ! Add the new element

            DEALLOCATE(elist)

            CALL MOVE_ALLOC(celist,elist)

        ELSE ! Flist not allocated

            ALLOCATE(elist(1,Kf*Kf))

            elist(1,:) = e

        END IF

    END SUBROUTINE

END MODULE DIIS
