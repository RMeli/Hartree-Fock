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

MODULE LA
    ! --------------
    ! LINEAR ALGEBRA
    ! --------------

    IMPLICIT NONE

    CONTAINS

    ! ------------------
    ! EIGENVALUE PROBLEM
    ! ------------------
    SUBROUTINE EIGS(d,M,V,lambda)
        ! -------------------------------------
        ! Solve eigenvalue problem using LAPACK
        ! -------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: d                        ! Dimension of M (d x d)
        REAL*8, dimension(d,d),intent(in) :: M          ! Matrix M

        ! INTERMEDIATE VARIABLES
        INTEGER, PARAMETER :: LWMAX = 10000             ! Maximal workspace size
        INTEGER :: LWORK                                ! Query optimal workspace size ()
        INTEGER :: INFO                                 ! Information flag for DSYEV
        REAL*8, dimension(LWMAX) :: WORK

        ! OUTPUT
        REAL*8, dimension(d,d), intent(out) :: V        ! Eigenvectors
        REAL*8, dimension(d), intent(out) :: lambda     ! Eigenvalues

        V = M

        LWORK = -1
        CALL  DSYEV('V','U', d, V, d, lambda, WORK, LWORK, INFO )       ! Query optimal workspace size
        LWORK = MIN(LWMAX, INT(WORK(1))) * 2 * d                        ! Optimal workspace size
                                                                        ! See (http://www.netlib.org/lapack/lug/node120.html)

        CALL  DSYEV('V','U', d, V, d, lambda, WORK, LWORK, INFO )       ! Solve the eigenvalue problem

        IF (INFO .NE. 0) THEN
            WRITE(*,*) "ERROR: IMPOSSIBLE TO SOLVE THE EIGENVALUE PROBLEM!"
        END IF

    END SUBROUTINE EIGS


END MODULE LA
