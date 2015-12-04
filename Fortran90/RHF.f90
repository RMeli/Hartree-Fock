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

MODULE RHF
    ! ----------------------------------------------
    ! RESTRICTED HARTREE-FOCK
    ! ----------------------------------------------
    !
    ! Self-consistent solution of Roothaan equations
    ! (Restricted Hartree-Fock)
    !
    ! ----------------------------------------------

    USE UTILS, only: EIGS, print_real_matrix, print_ee_list
    USE DENSITY, only: P_density, delta_P
    USE ELECTRONIC, only: G_ee, EE_list
    USE ENERGY, only: E_tot
    USE CONSTANTS
    USE OVERLAP, only: S_overlap, X_transform
    USE CORE, only: H_core

    IMPLICIT NONE

    CONTAINS

        ! --------
        ! SCF STEP
        ! --------
        SUBROUTINE RHF_step(Kf,Ne,H,X,ee,Pold,Pnew,F,orbitalE,verbose)
            ! ----------------------------------------------------------------------
            ! Perform a single step of the SCF procedure to solve Roothan equations.
            ! ----------------------------------------------------------------------
            !
            ! Source:

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                           ! Number of basis functions
            !INTEGER, intent(in) :: Nn                           ! Number of nuclei
            INTEGER, intent(in) :: Ne                           ! Number of electrons
            !REAL*8, dimension(Kf,3), intent(in) :: basis_D      ! Contraction coefficients
            !REAL*8, dimension(Kf,3), intent(in) :: basis_A      ! Contraction exponential coefficients
            !INTEGER, dimension(Kf,3), intent(in) :: basis_L     ! Basis function angular momenta
            !REAL*8, dimension(Kf,3), intent(in) :: basis_R      ! Basis set centers
            !REAL*8, dimension(Nn,3), intent(in) :: Rn           ! Nuclear coordinates
            !INTEGER, dimension(Nn), intent(in) :: Zn            ! Nuclear charges
            REAL*8, dimension(Kf,Kf), intent(in) :: H           ! Core Hamiltonian
            REAL*8, dimension(Kf,Kf), intent(in) :: X           ! Transformation maxrix
            REAL*8, dimension(Kf,Kf), intent(in) :: Pold        ! Old density matrix
            REAL*8, dimension(Kf,Kf,Kf,Kf), intent(in) :: ee    ! Electron-electron list
            LOGICAL, intent(in) :: verbose                      ! Flag to print matrices

            ! INTERMEDIATE VARIABLES
            REAL*8, dimension(Kf,Kf) :: G                       ! Electron-electron repulsion matrix
            REAL*8, dimension(Kf,Kf) :: Fx                      ! Fock matrix in the orthogonal basis set
            REAL*8, dimension(Kf,Kf) :: Cx                      ! Coefficient matrix in the orthogonal basis set
            REAL*8, dimension(Kf,Kf) :: C                       ! Coefficient matrix in the original basis set

            ! INPUT / OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: Pnew       ! New density matrix
            REAL*8, dimension(Kf,Kf), intent(out) :: F          ! Fock operator
            REAL*8, dimension(Kf), intent(out) :: orbitalE      ! Orbital energies

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Density matrix P:"
                CALL print_real_matrix(Kf,Kf,Pold)
            END IF

            CALL G_ee(Kf,ee,Pold,G) ! Compute new electron-electron repulsion matrix G

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Electron-electron repulsion matrix G:"
                CALL print_real_matrix(Kf,Kf,G)
            END IF

            F = H + G ! Compute new Fock operator

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix F:"
                CALL print_real_matrix(Kf,Kf,F)
            END IF

            Fx = MATMUL(TRANSPOSE(X),MATMUL(F,X))

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix in orthogonal orbital basis Fx:"
                CALL print_real_matrix(Kf,Kf,Fx)
            END IF

            ! ------------------------------------
            ! Solve orthogonalized Roothan eqution
            ! ------------------------------------
            CALL EIGS(Kf,Fx,Cx,orbitalE) ! Compute coefficients (in orthonormal basis) and orbital energies

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients in orthogonal orbital basis Cx:"
                CALL print_real_matrix(Kf,Kf,Cx)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Orbital energies:"
                CALL print_real_matrix(Kf,1,orbitalE)
            END IF

            C = MATMUL(X,Cx) ! Compute coefficients in the original basis

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients:"
                CALL print_real_matrix(Kf,Kf,C)
            END IF

            CALL P_density(Kf,Ne,C,Pnew) ! Compute new density

        END SUBROUTINE RHF_step


        ! ---------
        ! SCF Cycle
        ! ---------
        SUBROUTINE RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,verbose)
            ! --------------------------------
            ! Compute total energy (SCF cycle)
            ! --------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                                   ! Basis set size
            INTEGER, intent(in) :: c                                    ! Number of contractions
            INTEGER, intent(in) :: Ne                                   ! Number of electrons
            INTEGER, intent(in) :: Nn                                   ! Number of nuclei
            INTEGER, dimension(Kf,3), intent(in) :: basis_L             ! Angular momenta of basis set Gaussians
            REAL*8, dimension(Kf,3), intent(in) :: basis_R              ! Centers of basis set Gaussians
            REAL*8, dimension(Kf,c), intent(in) :: basis_D, basis_A     ! Basis set coefficients
            REAL*8, dimension(Nn,3), intent(in) :: Rn                   ! Nuclear positions
            INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
            LOGICAL, intent(in) :: verbose                              ! Verbose flag

            ! OUTPUT
            REAL*8,intent(out) :: final_E ! Converged total energy

            ! --------
            ! Matrices
            ! --------

            REAL*8, dimension(Kf,Kf) :: S     ! Overlap matrix
            REAL*8, dimension(Kf,Kf) :: X     ! Transformation matrix

            REAL*8, dimension(Kf,Kf) :: Hc    ! Core Hamiltonian

            REAL*8, dimension(Kf,Kf) :: Pold  ! Old density matrix
            REAL*8, dimension(Kf,Kf) :: Pnew  ! New density matrix

            REAl*8, dimension(Kf,Kf,Kf,Kf) :: ee ! List of electron-electron integrals

            REAL*8, dimension(Kf,Kf) :: F     ! Fock matrix
            REAL*8, dimension(Kf) :: E           ! Orbital energies

            ! --------------
            ! SCF parameters
            ! --------------

            LOGICAL :: converged                    ! Convergence parameter
            INTEGER, PARAMETER :: maxiter = 250     ! Maximal number of iterations (TODO: user defined maxiter)
            INTEGER :: step                         ! SCF steps counter

            !!!
            !!! NEVER INITIALIZE ON DECLARATION
            !!! "A local variable that is initialized when declared has an implicit save attribute.
            !!! The variable is initialized only the first time the unction is called.
            !!! On subsequent calls the old value is retained."
            !!!

            converged = .FALSE.
            step = 0

            ! -----------------
            ! HF initialization
            ! -----------------

            CALL S_overlap(Kf,c,basis_D,basis_A,basis_L,basis_R,S) ! Compute overlap matrix

            IF (verbose) THEN
                WRITE(*,*) "Overlap matrix S:"
                CALL print_real_matrix(Kf,Kf,S)
            END IF

            CALL X_transform(Kf,S,X) ! Compute transformation matrix

            IF (verbose) THEN
                WRITE(*,*) "Transformation matrix X:"
                CALL print_real_matrix(Kf,Kf,X)
            END IF

            CALL H_core(Kf,c,Nn,basis_D,basis_A,basis_L,basis_R,Rn,Zn,Hc,verbose)

            IF (verbose) THEN
                WRITE(*,*) "Core Hamiltonian Hc:"
                CALL print_real_matrix(Kf,Kf,Hc)
            END IF

            CALL EE_list(Kf,c,basis_D,basis_A,basis_L,basis_R,ee)

            IF (verbose) THEN
                WRITE(*,*) "Two-electron integrals:"
                CALL print_ee_list(Kf,ee)
            END IF

            Pold(:,:) = 0.0D0
            Pnew(:,:) = 0.0D0

            ! ---------
            ! SCF cycle
            ! ---------

            DO WHILE ((converged .EQV. .FALSE.) .AND. step .LT. maxiter)
                step = step + 1

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*)
                    WRITE(*,*)
                    WRITE(*,*) "--------"
                    WRITE(*,*) "SCF step #", step
                    WRITE(*,*) "--------"
                END IF

                CALL RHF_step(Kf,Ne,Hc,X,ee,Pold,Pnew,F,E,verbose)

                IF (verbose) THEN
                    WRITE(*,*   )
                    WRITE(*,*) "Total energy:", E_tot(Kf,Nn,Rn,Zn,Pold,F,Hc)
                END IF

                IF ( delta_P(Kf,Pold,Pnew) < 1.0e-12) THEN
                    converged = .TRUE.

                    final_E = E_tot(Kf,Nn,Rn,Zn,Pold,F,Hc)

                    IF (verbose) THEN
                        WRITE(*,*)
                        WRITE(*,*) "SCF cycle converged!"
                        WRITE(*,*)
                        WRITE(*,*)
                        WRITE(*,*)
                        WRITE(*,*) "TOTAL ENERGY:", final_E
                    END IF
                END IF

                Pold = Pnew

            END DO ! SCF

            IF (converged .EQV. .FALSE.) THEN
                WRITE(*,*)
                WRITE(*,*) "SCF NOT CONVERGED!"
                CALL EXIT(-1)
            END IF

        END SUBROUTINE RHF_SCF


END MODULE RHF
