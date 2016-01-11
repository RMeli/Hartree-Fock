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

MODULE UHF
    ! --------------------------------------------------
    ! UNRESTRICTED HARTREE-FOCK
    ! --------------------------------------------------
    !
    ! Self-consistent solution of Pople-Nesbet equations
    ! (Unrestricted Hartree-Fock)
    !
    ! --------------------------------------------------

    USE LA, only: EIGS
    USE OUTPUT, only: print_real_matrix, print_ee_list
    USE DENSITY, only: P_density_spin, delta_P
    USE ELECTRONIC, only: G_ee_spin, EE_list
    USE ENERGY, only: E_tot_spin
    USE OVERLAP, only: S_overlap, X_transform
    USE CORE, only: H_core
    USE CONSTANTS
    USE INIT
    USE DIIS

    IMPLICIT NONE

    CONTAINS

        ! --------
        ! SCF STEP
        ! --------
        SUBROUTINE UHF_step(Kf,Ne,H,X,ee,Paold,Panew,Pbold,Pbnew,Fa,Fb,orbitalEa,orbitalEb,step,verbose)
            ! ---------------------------------------------------------------------------
            ! Perform a single step of the SCF procedure to solve Pople-Nesbet equations.
            ! ---------------------------------------------------------------------------
            !
            ! Source:
            !   A. Szabo and N. S. Ostlund
            !   Modern Quantum Chemistry
            !   Dover
            !   1996
            !
            ! ---------------------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                           ! Number of basis functions
            INTEGER, intent(in) :: Ne                           ! Number of electrons
            REAL*8, dimension(Kf,Kf), intent(in) :: H           ! Core Hamiltonian
            REAL*8, dimension(Kf,Kf), intent(in) :: X           ! Transformation maxrix
            REAL*8, dimension(Kf,Kf), intent(in) :: Paold       ! Old density matrix for alpha spin
            REAL*8, dimension(Kf,Kf), intent(in) :: Pbold       ! Old density matrix for beta spin
            REAL*8, dimension(Kf,Kf,Kf,Kf), intent(in) :: ee    ! Electron-electron list
            INTEGER, intent(in) :: step                         ! Current step
            LOGICAL, intent(in) :: verbose                      ! Flag to print matrices

            ! INTERMEDIATE VARIABLES
            REAL*8, dimension(Kf,Kf) :: Ga                       ! Electron-electron repulsion matrix for alpha spin
            REAL*8, dimension(Kf,Kf) :: Fax                      ! Fock matrix in the orthogonal basis set for alpha spin
            REAL*8, dimension(Kf,Kf) :: Cax                      ! Coefficient matrix in the orthogonal basis set for alpha spin
            REAL*8, dimension(Kf,Kf) :: Ca                       ! Coefficient matrix in the original basis set for alpha spin
            REAL*8, dimension(Kf,Kf) :: Gb                       ! Electron-electron repulsion matrix for beta spin
            REAL*8, dimension(Kf,Kf) :: Fbx                      ! Fock matrix in the orthogonal basis set for beta spin
            REAL*8, dimension(Kf,Kf) :: Cbx                      ! Coefficient matrix in the orthogonal basis set for beta spin
            REAL*8, dimension(Kf,Kf) :: Cb                       ! Coefficient matrix in the original basis set for beta spin
            INTEGER :: Nea                                       ! Alpha spin electrons
            INTEGER :: Neb                                       ! Beta spin electrons
            REAL*8 :: mix                                        ! Mixing coefficient for alpha HOMO/LUMO
            REAL*8, dimension(Kf) :: mix_HOMO                    ! Temporary storage of old HOMO coefficients

            ! INPUT / OUTPUT
            REAL*8, dimension(Kf,Kf), intent(inout) :: Fa        ! Fock operator for alpha spin
            REAL*8, dimension(Kf,Kf), intent(inout) :: Fb        ! Fock operator for beta spin

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: Panew       ! New density matrix for alpha spin
            REAL*8, dimension(Kf,Kf), intent(out) :: Pbnew       ! New density matrix for beta spin
            REAL*8, dimension(Kf), intent(out) :: orbitalEa      ! Orbital energies for alpha spin
            REAL*8, dimension(Kf), intent(out) :: orbitalEb      ! Orbital energies for beta spin

            ! COMPUTE SPIN POPULATION

            IF (MODULO(Ne,2) .EQ. 0) THEN ! Even number of electrons
                Nea = Ne / 2
                Neb = Ne / 2
            ELSE ! Odd number of electrons
                Nea = (Ne + 1) / 2
                Neb = Ne - Nea
            END IF

            IF (step .NE. 1) THEN

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Density matrix P (alpha spin):"
                    CALL print_real_matrix(Kf,Kf,Paold)
                END IF

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Density matrix P (beta spin):"
                    CALL print_real_matrix(Kf,Kf,Pbold)
                END IF

                CALL G_ee_spin(Kf,ee,Paold+Pbold,Paold,Ga) ! Compute new electron-electron repulsion matrix Ga for alpha spin
                CALL G_ee_spin(Kf,ee,Paold+Pbold,Pbold,Gb) ! Compute new electron-electron repulsion matrix Gb for beta spin

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Electron-electron repulsion matrix G (alpha spin):"
                    CALL print_real_matrix(Kf,Kf,Ga)
                END IF

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Electron-electron repulsion matrix G (beta spin):"
                    CALL print_real_matrix(Kf,Kf,Gb)
                END IF

                Fa = H + Ga ! Compute new Fock operator (alpha spin)
                Fb = H + Gb ! Compute new Fock operator (alpha spin)

            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix F (alpha spin):"
                CALL print_real_matrix(Kf,Kf,Fa)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix F (beta spin):"
                CALL print_real_matrix(Kf,Kf,Fb)
            END IF

            Fax = MATMUL(TRANSPOSE(X),MATMUL(Fa,X))
            Fbx = MATMUL(TRANSPOSE(X),MATMUL(Fb,X))

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix in orthogonal orbital basis Fx (alpha spin)"
                CALL print_real_matrix(Kf,Kf,Fax)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix in orthogonal orbital basis Fx (beta spin):"
                CALL print_real_matrix(Kf,Kf,Fbx)
            END IF

            ! ------------------------------------------
            ! Solve orthogonalized Pople-Nesbet equtions
            ! ------------------------------------------

            CALL EIGS(Kf,Fax,Cax,orbitalEa) ! Compute coefficients (in orthonormal basis) and orbital energies (alpha spin)
            CALL EIGS(Kf,Fbx,Cbx,orbitalEb) ! Compute coefficients (in orthonormal basis) and orbital energies (beta spin)

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients in orthogonal orbital basis Cx (alpha spin):"
                CALL print_real_matrix(Kf,Kf,Cax)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients in orthogonal orbital basis Cx (beta spin):"
                CALL print_real_matrix(Kf,Kf,Cbx)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Orbital energies (alpha spin):"
                CALL print_real_matrix(Kf,1,orbitalEa)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Orbital energies (beta spin):"
                CALL print_real_matrix(Kf,1,orbitalEb)
            END IF

            Ca = MATMUL(X,Cax) ! Compute coefficients in the original basis
            Cb = MATMUL(X,Cbx) ! Compute coefficients in the original basis

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients (alpha spin):"
                CALL print_real_matrix(Kf,Kf,Ca)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients (beta spin):"
                CALL print_real_matrix(Kf,Kf,Cb)
            END IF

            ! Mix alpha LUMO/HOMO orbitals in order to destroy spatial simmetry

            mix = 1.0D0

            IF ((step .EQ. 1) .AND. (Nea .EQ. Neb)) THEN
                mix_HOMO = Ca(Nea,:)

                Ca(Nea,:) = 1.0D0 / DSQRT(1.0D0 + mix**2) * (Ca(Nea,:) + mix * Ca(Nea+1,:))
                Ca(Nea+1,:) = 1.0D0 / DSQRT(1.0D0 + mix**2) * (mix_HOMO - mix * Ca(Nea+1,:))
            END IF

            ! Update density

            CALL P_density_spin(Kf,Nea,Ca,Panew) ! Compute new density matrix for alpha spin
            CALL P_density_spin(Kf,Neb,Cb,Pbnew) ! Compute new density matrix for beta spin

        END SUBROUTINE UHF_step


        ! ---------
        ! SCF Cycle
        ! ---------
        SUBROUTINE UHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,Ptot,verbose)
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
            REAL*8,intent(out) :: final_E                               ! Converged total energy
            REAl*8, dimension(Kf,Kf), intent(out) :: Ptot               ! Total density matrix

            ! --------
            ! Matrices
            ! --------

            REAL*8, dimension(Kf,Kf) :: S     ! Overlap matrix
            REAL*8, dimension(Kf,Kf) :: X     ! Transformation matrix

            REAL*8, dimension(Kf,Kf) :: Hc    ! Core Hamiltonian

            REAL*8, dimension(Kf,Kf) :: Paold  ! Old density matrix (alpha spin)
            REAL*8, dimension(Kf,Kf) :: Panew  ! New density matrix (alpha spin)
            REAL*8, dimension(Kf,Kf) :: Pbold  ! Old density matrix (beta spin)
            REAL*8, dimension(Kf,Kf) :: Pbnew  ! New density matrix (beta spin)

            REAl*8, dimension(Kf,Kf,Kf,Kf) :: ee ! List of electron-electron integrals

            REAL*8, dimension(Kf,Kf) :: Fa     ! Fock matrix (alpha spin)
            REAL*8, dimension(Kf,Kf) :: Fb     ! Fock matrix (beta spin)
            REAL*8, dimension(Kf) :: Ea        ! Orbital energies (alpha spin)
            REAL*8, dimension(Kf) :: Eb        ! Orbital energies (beta spin)

            ! --------------
            ! SCF parameters
            ! --------------

            LOGICAL :: converged                    ! Convergence parameter
            INTEGER, PARAMETER :: maxiter = 1000    ! Maximal number of iterations (TODO: user defined maxiter)
            INTEGER :: step                         ! SCF steps counter

            !!!
            !!! NEVER INITIALIZE ON DECLARATION VARIABLES OTHER THAN PARAMETERS
            !!! "A local variable that is initialized when declared has an implicit save attribute.
            !!! The variable is initialized only the first time the unction is called.
            !!! On subsequent calls the old value is retained."
            !!!

            converged = .FALSE.
            step = 0

            ! ------------------
            ! RHF initialization
            ! ------------------

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

            Paold(:,:) = 0.0D0
            Panew(:,:) = 0.0D0
            Pbold(:,:) = 0.0D0
            Pbnew(:,:) = 0.0D0

            ! -------------
            ! Initial guess
            ! -------------

            CALL huckel_guess(Kf,Hc,S,Fa,1.45D0)
            CALL huckel_guess(Kf,Hc,S,Fb,1.45D0)

            ! ---------
            ! SCF cycle
            ! ---------

            DO WHILE ((converged .EQV. .FALSE.) .AND. (step .LT. maxiter))
                step = step + 1

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*)
                    WRITE(*,*)
                    WRITE(*,*) "--------"
                    WRITE(*,*) "SCF step #", step
                    WRITE(*,*) "--------"
                END IF

                CALL UHF_step(Kf,Ne,Hc,X,ee,Paold,Panew,Pbold,Pbnew,Fa,Fb,Ea,Eb,step,verbose)

                IF (verbose) THEN
                    WRITE(*,*   )
                    WRITE(*,*) "Total energy:", E_tot_spin(Kf,Nn,Rn,Zn,Paold,Pbold,Fa,Fb,Hc)
                END IF

                IF ( (delta_P(Kf,Paold,Panew) .LT. 1.0e-6) .AND. (delta_P(Kf,Pbold,Pbnew) .LT. 1.0e-6)  ) THEN
                    converged = .TRUE.

                    final_E = E_tot_spin(Kf,Nn,Rn,Zn,Paold,Pbold,Fa,Fb,Hc)

                    IF (verbose) THEN
                        WRITE(*,*)
                        WRITE(*,*) "SCF cycle converged!"
                        WRITE(*,*)
                        WRITE(*,*)
                        WRITE(*,*)
                        WRITE(*,*) "TOTAL ENERGY:", final_E
                    END IF
                END IF

                Paold = Panew
                Pbold = Pbnew

            END DO ! SCF

            IF (converged .EQV. .FALSE.) THEN
                WRITE(*,*)
                WRITE(*,*) "SCF NOT CONVERGED!"
                CALL EXIT(-1)
            END IF

            Ptot = Paold + Pbold

        END SUBROUTINE UHF_SCF


        ! ------------------------
        ! SCF USING DIIS ALGORITHM
        ! ------------------------
        SUBROUTINE UHF_step_DIIS(Kf,Ne,H,S,X,ee,Paold,Panew,Pbold,Pbnew,Fa,Fb,orbitalEa,orbitalEb,step,verbose,&
                                Falist,Ealist,Fblist,Eblist,DIIS_step,DIIS_flag)
            ! ---------------------------------------------------------------------------
            ! Perform a single step of the SCF procedure to solve Pople-Nesbet equations.
            ! ---------------------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                           ! Number of basis functions
            INTEGER, intent(in) :: Ne                           ! Number of electrons
            REAL*8, dimension(Kf,Kf), intent(in) :: H           ! Core Hamiltonian
            REAL*8, dimension(Kf,Kf), intent(in) :: S           ! Overlap matrix
            REAL*8, dimension(Kf,Kf), intent(in) :: X           ! Transformation maxrix
            REAL*8, dimension(Kf,Kf), intent(in) :: Paold       ! Old density matrix for alpha spin
            REAL*8, dimension(Kf,Kf), intent(in) :: Pbold       ! Old density matrix for beta spin
            REAL*8, dimension(Kf,Kf,Kf,Kf), intent(in) :: ee    ! Electron-electron list
            INTEGER, intent(in) :: step                         ! Current step
            LOGICAL, intent(in) :: verbose                      ! Flag to print matrices

            ! INTERMEDIATE VARIABLES
            REAL*8, dimension(Kf,Kf) :: Ga                       ! Electron-electron repulsion matrix for alpha spin
            REAL*8, dimension(Kf,Kf) :: Fax                      ! Fock matrix in the orthogonal basis set for alpha spin
            REAL*8, dimension(Kf,Kf) :: Cax                      ! Coefficient matrix in the orthogonal basis set for alpha spin
            REAL*8, dimension(Kf,Kf) :: Ca                       ! Coefficient matrix in the original basis set for alpha spin
            REAL*8, dimension(Kf,Kf) :: Gb                       ! Electron-electron repulsion matrix for beta spin
            REAL*8, dimension(Kf,Kf) :: Fbx                      ! Fock matrix in the orthogonal basis set for beta spin
            REAL*8, dimension(Kf,Kf) :: Cbx                      ! Coefficient matrix in the orthogonal basis set for beta spin
            REAL*8, dimension(Kf,Kf) :: Cb                       ! Coefficient matrix in the original basis set for beta spin
            INTEGER :: Nea                                       ! Alpha spin electrons
            INTEGER :: Neb                                       ! Beta spin electrons
            REAL*8 :: mix                                        ! Mixing coefficient for alpha HOMO/LUMO
            REAL*8, dimension(Kf) :: mix_HOMO                    ! Temporary storage of old HOMO coefficients
            REAL*8, dimension(Kf*Kf) :: errora
            REAL*8 :: maxerrora
            REAL*8, dimension(Kf*Kf) :: errorb
            REAL*8 :: maxerrorb

            ! INPUT / OUTPUT
            REAL*8, dimension(Kf,Kf), intent(inout) :: Fa                   ! Fock operator for alpha spin
            REAL*8, dimension(Kf,Kf), intent(inout) :: Fb                   ! Fock operator for beta spin
            REAL*8, allocatable, dimension(:,:,:), intent(inout) :: Falist  ! List of Fock operators (alpha spin)
            REAL*8, allocatable, dimension(:,:), intent(inout) :: Ealist    ! Error list (alpha spin)
            REAL*8, allocatable, dimension(:,:,:), intent(inout) :: Fblist  ! List of Fock operators (beta spin)
            REAL*8, allocatable, dimension(:,:), intent(inout) :: Eblist    ! Error list (beta spin)
            INTEGER, intent(inout) :: DIIS_step                             ! DIIS step counter
            LOGICAL, intent(inout) :: DIIS_flag                             ! DIIS algorithm flag

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: Panew       ! New density matrix for alpha spin
            REAL*8, dimension(Kf,Kf), intent(out) :: Pbnew       ! New density matrix for beta spin
            REAL*8, dimension(Kf), intent(out) :: orbitalEa      ! Orbital energies for alpha spin
            REAL*8, dimension(Kf), intent(out) :: orbitalEb      ! Orbital energies for beta spin

            ! COMPUTE SPIN POPULATION

            IF (MODULO(Ne,2) .EQ. 0) THEN ! Even number of electrons
                Nea = Ne / 2
                Neb = Ne / 2
            ELSE ! Odd number of electrons
                Nea = (Ne + 1) / 2
                Neb = Ne - Nea
            END IF

            IF (step .NE. 1) THEN

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Density matrix P (alpha spin):"
                    CALL print_real_matrix(Kf,Kf,Paold)
                END IF

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Density matrix P (beta spin):"
                    CALL print_real_matrix(Kf,Kf,Pbold)
                END IF

                CALL G_ee_spin(Kf,ee,Paold+Pbold,Paold,Ga) ! Compute new electron-electron repulsion matrix Ga for alpha spin
                CALL G_ee_spin(Kf,ee,Paold+Pbold,Pbold,Gb) ! Compute new electron-electron repulsion matrix Gb for beta spin

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Electron-electron repulsion matrix G (alpha spin):"
                    CALL print_real_matrix(Kf,Kf,Ga)
                END IF

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Electron-electron repulsion matrix G (beta spin):"
                    CALL print_real_matrix(Kf,Kf,Gb)
                END IF

                Fa = H + Ga ! Compute new Fock operator (alpha spin)
                Fb = H + Gb ! Compute new Fock operator (alpha spin)

                IF (DIIS_flag .EQV. .FALSE.) THEN
                    CALL DIIS_error(Kf,Fa,Paold,S,X,errora,maxerrora)
                    CALL DIIS_error(Kf,Fb,Pbold,S,X,errorb,maxerrorb)

                    IF ((maxerrora .LT. 1.0D-1) .AND. (maxerrorb .LT. 1.0D-1)) THEN ! Initiate DIIS procedure

                        DIIS_flag = .TRUE.

                        IF (verbose .EQV. .TRUE.) THEN
                            WRITE(*,*)
                            WRITE(*,*) "Using DIIS accelerated SCF algorithm."
                        END IF

                    END IF
                END IF

            END IF

            ! --------------
            ! DIIS ALGORITHM
            ! --------------

            IF (DIIS_flag .EQV. .TRUE.) THEN
                DIIS_step = DIIS_step + 1

                ! Update Fock matrix following DIIS algorithm
                CALL DIIS_Fock(Kf,DIIS_step,Fa,Paold,S,X,Falist,Ealist,verbose)
                CALL DIIS_Fock(Kf,DIIS_step,Fb,Pbold,S,X,Fblist,Eblist,verbose)
            END IF


            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix F (alpha spin):"
                CALL print_real_matrix(Kf,Kf,Fa)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix F (beta spin):"
                CALL print_real_matrix(Kf,Kf,Fb)
            END IF

            Fax = MATMUL(TRANSPOSE(X),MATMUL(Fa,X))
            Fbx = MATMUL(TRANSPOSE(X),MATMUL(Fb,X))

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix in orthogonal orbital basis Fx (alpha spin)"
                CALL print_real_matrix(Kf,Kf,Fax)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix in orthogonal orbital basis Fx (beta spin):"
                CALL print_real_matrix(Kf,Kf,Fbx)
            END IF

            ! ------------------------------------------
            ! Solve orthogonalized Pople-Nesbet equtions
            ! ------------------------------------------

            CALL EIGS(Kf,Fax,Cax,orbitalEa) ! Compute coefficients (in orthonormal basis) and orbital energies (alpha spin)
            CALL EIGS(Kf,Fbx,Cbx,orbitalEb) ! Compute coefficients (in orthonormal basis) and orbital energies (beta spin)

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients in orthogonal orbital basis Cx (alpha spin):"
                CALL print_real_matrix(Kf,Kf,Cax)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients in orthogonal orbital basis Cx (beta spin):"
                CALL print_real_matrix(Kf,Kf,Cbx)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Orbital energies (alpha spin):"
                CALL print_real_matrix(Kf,1,orbitalEa)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Orbital energies (beta spin):"
                CALL print_real_matrix(Kf,1,orbitalEb)
            END IF

            Ca = MATMUL(X,Cax) ! Compute coefficients in the original basis
            Cb = MATMUL(X,Cbx) ! Compute coefficients in the original basis

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients (alpha spin):"
                CALL print_real_matrix(Kf,Kf,Ca)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients (beta spin):"
                CALL print_real_matrix(Kf,Kf,Cb)
            END IF

            ! Mix alpha LUMO/HOMO orbitals in order to destroy spatial simmetry

            !mix = 1.0D0

            !IF ((step .EQ. 1) .AND. (Nea .EQ. Neb)) THEN
            !    mix_HOMO = Ca(Nea,:)

            !    Ca(Nea,:) = 1.0D0 / DSQRT(1.0D0 + mix**2) * (Ca(Nea,:) + mix * Ca(Nea+1,:))
            !    Ca(Nea+1,:) = 1.0D0 / DSQRT(1.0D0 + mix**2) * (mix_HOMO - mix * Ca(Nea+1,:))
            !END IF

            ! Update density

            CALL P_density_spin(Kf,Nea,Ca,Panew) ! Compute new density matrix for alpha spin
            CALL P_density_spin(Kf,Neb,Cb,Pbnew) ! Compute new density matrix for beta spin


        END SUBROUTINE UHF_step_DIIS

        ! ------------------------------
        ! SCF CYCLE USING DIIS ALGORITHM
        ! ------------------------------
        SUBROUTINE UHF_DIIS(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,Ptot,verbose)
            ! --------------------------------------------------------------------------------
            ! Compute total energy (SCF cycle) using DIIS algorithm to accelerate convergence.
            ! --------------------------------------------------------------------------------

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

            ! INTERMEDIATE VARIABLES
            REAL*8, allocatable, dimension(:,:,:) :: Falist             ! List of Fock matrices (alpha spin)
            REAL*8, allocatable, dimension(:,:) :: Ealist               ! List of errors (alpha spin)
            REAL*8, allocatable, dimension(:,:,:) :: Fblist             ! List of Fock matrices (beta spin)
            REAL*8, allocatable, dimension(:,:) :: Eblist               ! List of errors (beta spin)

            ! OUTPUT
            REAL*8,intent(out) :: final_E                               ! Converged total energy
            REAl*8, dimension(Kf,Kf), intent(out) :: Ptot               ! Total density matrix

            ! --------
            ! Matrices
            ! --------

            REAL*8, dimension(Kf,Kf) :: S     ! Overlap matrix
            REAL*8, dimension(Kf,Kf) :: X     ! Transformation matrix

            REAL*8, dimension(Kf,Kf) :: Hc    ! Core Hamiltonian

            REAL*8, dimension(Kf,Kf) :: Paold  ! Old density matrix (alpha spin)
            REAL*8, dimension(Kf,Kf) :: Panew  ! New density matrix (alpha spin)
            REAL*8, dimension(Kf,Kf) :: Pbold  ! Old density matrix (beta spin)
            REAL*8, dimension(Kf,Kf) :: Pbnew  ! New density matrix (beta spin)

            REAl*8, dimension(Kf,Kf,Kf,Kf) :: ee ! List of electron-electron integrals

            REAL*8, dimension(Kf,Kf) :: Fa     ! Fock matrix (alpha spin)
            REAL*8, dimension(Kf,Kf) :: Fb     ! Fock matrix (beta spin)
            REAL*8, dimension(Kf) :: Ea        ! Orbital energies (alpha spin)
            REAL*8, dimension(Kf) :: Eb        ! Orbital energies (beta spin)

            ! --------------
            ! SCF parameters
            ! --------------

            LOGICAL :: converged                    ! Convergence parameter
            INTEGER, PARAMETER :: maxiter = 1000    ! Maximal number of iterations (TODO: user defined maxiter)
            INTEGER :: step                         ! SCF steps counter
            INTEGER :: DIIS_step                    ! DIIS step counter
            LOGICAL :: DIIS_flag                    ! DIIS algorithm flag

            !!!
            !!! NEVER INITIALIZE ON DECLARATION VARIABLES OTHER THAN PARAMETERS
            !!! "A local variable that is initialized when declared has an implicit save attribute.
            !!! The variable is initialized only the first time the unction is called.
            !!! On subsequent calls the old value is retained."
            !!!

            converged = .FALSE.
            step = 0
            DIIS_step = 0
            DIIS_flag = .FALSE.

            ! ------------------
            ! UHF initialization
            ! ------------------

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

            Paold(:,:) = 0.0D0
            Panew(:,:) = 0.0D0
            Pbold(:,:) = 0.0D0
            Pbnew(:,:) = 0.0D0

            ! -------------
            ! Initial guess
            ! -------------

            CALL core_guess(Kf,Hc,Fa)
            CALL core_guess(Kf,Hc,Fb)
            !CALL huckel_guess(Kf,Hc,S,Fa,1.75D0)
            !CALL huckel_guess(Kf,Hc,S,Fb,1.75D0)

            ! ---------
            ! SCF cycle
            ! ---------

            DO WHILE ((converged .EQV. .FALSE.) .AND. (step .LT. maxiter))
                step = step + 1

                IF (verbose) THEN
                    WRITE(*,*)
                    WRITE(*,*)
                    WRITE(*,*)
                    WRITE(*,*) "--------"
                    WRITE(*,*) "SCF step #", step
                    WRITE(*,*) "--------"
                END IF

                !
                ! ------------------------------------------------------------------------
                !
                CALL UHF_step_DIIS(Kf,Ne,Hc,S,X,ee,Paold,Panew,Pbold,Pbnew,Fa,Fb,Ea,Eb,step,verbose,&
                                    Falist,Ealist,Fblist,Eblist,DIIS_step,DIIS_flag)
                !
                ! ------------------------------------------------------------------------
                !

                IF (verbose) THEN
                    WRITE(*,*   )
                    WRITE(*,*) "Total energy:", E_tot_spin(Kf,Nn,Rn,Zn,Paold,Pbold,Fa,Fb,Hc)
                END IF

                IF ( (delta_P(Kf,Paold,Panew) .LT. 1.0e-6) .AND. (delta_P(Kf,Pbold,Pbnew) .LT. 1.0e-6)  ) THEN
                    converged = .TRUE.

                    final_E = E_tot_spin(Kf,Nn,Rn,Zn,Paold,Pbold,Fa,Fb,Hc)

                    IF (verbose) THEN
                        WRITE(*,*)
                        WRITE(*,*) "SCF cycle converged!"
                        WRITE(*,*)
                        WRITE(*,*)
                        WRITE(*,*)
                        WRITE(*,*) "TOTAL ENERGY:", final_E
                    END IF
                END IF

                Paold = Panew
                Pbold = Pbnew

            END DO ! SCF

            IF (converged .EQV. .FALSE.) THEN
                WRITE(*,*)
                WRITE(*,*) "SCF NOT CONVERGED!"
                CALL EXIT(-1)
            END IF

            Ptot = Paold + Pbold

        END SUBROUTINE UHF_DIIS


END MODULE UHF
