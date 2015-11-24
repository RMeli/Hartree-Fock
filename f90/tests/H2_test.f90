PROGRAM HF

    USE RHF
    USE OVERLAP, only: S_overlap, X_transform
    USE CORE, only: H_core
    USE ELECTRONIC, only: EE_list
    USE DENSITY, only: delta_P
    USE ENERGY, only: E_tot
    USE UTILS

    IMPLICIT NONE

    INTEGER, PARAMETER :: c = 3 ! Number of contractions (STO-3G)

    REAL*8, PARAMETER :: zeta_H = 1.24D0    ! STO coefficient correction for H

    ! ------------------
    ! BASIS SET FOR HeH+
    ! ------------------
    INTEGER, PARAMETER :: K = 2     ! Number of basis functions
    INTEGER, PARAMETER :: Nn = 2    ! Number of nuclei

    INTEGER, dimension(K,3) :: basis_L              ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K,3) :: basis_R               ! Centers of basis set Gaussians
    REAL*8, dimension(K,c) :: basis_D, basis_A  ! Basis set coefficients

    REAL*8, dimension(Nn,3) :: Rn   ! Nuclear positions
    INTEGER, dimension(Nn) :: Zn    ! Nuclear charges

    ! -----------------
    ! HF DECLARATION
    ! -----------------

    REAL*8, dimension(K,K) :: S     ! Overlap matrix
    REAL*8, dimension(K,K) :: X     ! Transformation matrix

    REAL*8, dimension(K,K) :: Hc    ! Core Hamiltonian

    REAL*8, dimension(K,K) :: Pold  ! Old density matrix
    REAL*8, dimension(K,K) :: Pnew  ! New density matrix

    REAl*8, dimension(K,K,K,K) :: ee ! List of electron-electron integrals

    REAL*8, dimension(K,K) :: F     ! Fock matrix
    REAL*8, dimension(K) :: E           ! Orbital energies

    LOGICAL :: converged = .FALSE.      ! Convergence parameter
    INTEGER, PARAMETER :: maxiter = 100 ! Maximal number of iterations
    INTEGER :: step = 0                 ! SCF steps counter


    ! -------------
    ! MOLECULE HeH+
    ! -------------

    INTEGER, PARAMETER :: Ne = 2    ! Total number of electrons



    Rn(1,1:3) = (/0.0D0, 0.0D0, 0.0D0/) ! Position of first H atom
    Rn(2,1:3) = (/1.4D0, 0.0D0, 0.0D0/) ! Position of the second H atom

    Zn = (/1, 1/)

    basis_R(1,1:3) = Rn(1,1:3) ! Position of first H atom
    basis_R(2,1:3) = Rn(2,1:3) ! Position of the second H atom

    basis_L(1,1:3) = (/0, 0, 0/) ! Angular momenta for the first basis function
    basis_L(2,1:3) = (/0, 0, 0/) ! Angular momenta for the second basis function

    basis_D(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)

    basis_A(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A(2,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)

    ! -----------------
    ! HF INITIALIZATION
    ! -----------------

    CALL S_overlap(K,basis_D,basis_A,basis_L,basis_R,S) ! Compute overlap matrix

    WRITE(*,*) "Overlap matrix S:"
    CALL print_real_matrix(K,K,S)

    CALL X_transform(K,S,X) ! Compute transformation matrix

    WRITE(*,*) "Transformation matrix X:"
    CALL print_real_matrix(K,K,X)

    CALL H_core(K,Nn,basis_D,basis_A,basis_L,basis_R,Rn,Zn,Hc)
    CALL EE_list(K,basis_D,basis_A,basis_L,basis_R,ee)

    CALL print_ee_list(K,ee)

    Pold(:,:) = 0.0D0
    Pnew(:,:) = 0.0D0

    ! ----------
    ! SCF CYCLES
    ! ----------

    CALL print_real_matrix(K,K,X)

    WRITE(*,*) "SCF step #", step

    DO WHILE ((converged .EQV. .FALSE.) .AND. step .LT. maxiter)
        step = step + 1

        CALL RHF_step(K,Ne,Hc,X,ee,Pold,Pnew,F,E,.TRUE.)

        WRITE(*,*   )
        WRITE(*,*) "Total energy:", E_tot(K,Nn,Rn,Zn,Pold,F,Hc)

        IF ( delta_P(K,Pold,Pnew) < 1.0e-12) THEN
            converged = .TRUE.

            WRITE(*,*)
            WRITE(*,*) "SCF cycle converged!"
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*) "TOTAL ENERGY:", E_tot(K,Nn,Rn,Zn,Pold,F,Hc)
        END IF

        Pold = Pnew

    END DO

END PROGRAM HF
