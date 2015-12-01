PROGRAM HF_H2

    USE RHF
    USE OVERLAP, only: S_overlap, X_transform
    USE CORE, only: H_core
    USE ELECTRONIC, only: EE_list
    USE DENSITY, only: delta_P
    USE ENERGY, only: E_tot
    USE UTILS
    USE CONSTANTS

    IMPLICIT NONE

    REAL*8, PARAMETER :: zeta_H = 1.24D0    ! STO coefficient correction for H
    REAL*8, PARAMETER :: zeta_O_1 = 7.66D0 ! STO coefficient correction for O
    REAL*8, PARAMETER :: zeta_O_2 = 2.25D0 ! STO coefficient correction for O

    ! -----------------
    ! BASIS SET FOR H2O
    ! -----------------
    INTEGER, PARAMETER :: K = 7     ! Number of basis functions
    INTEGER, PARAMETER :: Nn = 3    ! Number of nuclei

    INTEGER, dimension(K,3) :: basis_L              ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K,3) :: basis_R               ! Centers of basis set Gaussians
    REAL*8, dimension(K,c) :: basis_D, basis_A  ! Basis set coefficients

    ! ------------
    ! MOLECULE H2O
    ! ------------

    REAL*8, dimension(Nn,3) :: Rn    ! Nuclear positions
    INTEGER, dimension(Nn) :: Zn    ! Nuclear charges

    INTEGER, PARAMETER :: Ne = 10    ! Total number of electrons

    REAL*8 :: final_E               ! Total converged energy

    Rn(1,1:3) = (/1.809*SIN(104.52/180*PI/2.0D0), 0.0D0, 0.0D0/)        ! Position of first H atom
    Rn(2,1:3) = (/-1.809*SIN(104.52/180*PI/2.0D0), 0.0D0, 0.0D0/)       ! Position of the second H atom
    Rn(3,1:3) = (/0.0D0, 1.809*COS(104.52/180*PI/2.0D0), 0.0D0/)        ! Position of the O atom

    Zn = (/1, 1, 8/)

    basis_R(1,1:3) = Rn(1,1:3)   ! Center of H1 1s
    basis_R(2,1:3) = Rn(2,1:3)   ! Center of H2 1s
    basis_R(3,1:3) = Rn(3,1:3)   ! Center of O 1s
    basis_R(4,1:3) = Rn(3,1:3)   ! Center of O 2s
    basis_R(5,1:3) = Rn(3,1:3)   ! Center of O 2px
    basis_R(6,1:3) = Rn(3,1:3)   ! Center of O 2py
    basis_R(7,1:3) = Rn(3,1:3)   ! Center of O 2pz

    basis_L(1,1:3) = (/0, 0, 0/) ! Angular momenta for H 1s
    basis_L(2,1:3) = (/0, 0, 0/) ! Angular momenta for H 1s
    basis_L(3,1:3) = (/0, 0, 0/) ! Angular momenta for O 1s
    basis_L(4,1:3) = (/0, 0, 0/) ! Angular momenta for O 2s
    basis_L(5,1:3) = (/1, 0, 0/) ! Angular momenta for O 2px
    basis_L(6,1:3) = (/0, 1, 0/) ! Angular momenta for O 2py
    basis_L(7,1:3) = (/0, 0, 1/) ! Angular momenta for O 2pz

    basis_D(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D(3,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D(4,1:c) = (/0.700115D0, 0.399513D0, -0.0999672D0/)
    basis_D(5,1:c) = (/0.391957D0, 0.607684D0, 0.1559163D0/)
    basis_D(6,1:c) = (/0.391957D0, 0.607684D0, 0.1559163D0/)
    basis_D(7,1:c) = (/0.391957D0, 0.607684D0, 0.1559163D0/)

    basis_A(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A(2,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A(3,1:c) = (/0.109818D0 * zeta_O_1**2, 0.405771 * zeta_O_1**2, 2.22766 * zeta_O_1**2/)
    basis_A(4,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)
    basis_A(5,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)
    basis_A(6,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)
    basis_A(7,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)

    ! ------------------------
    ! TOTAL ENERGY CALCULATION
    ! ------------------------

    CALL SCF(K,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,.TRUE.)

END PROGRAM HF_H2
