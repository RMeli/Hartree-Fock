PROGRAM overlap_test

    USE OVERLAP
    USE UTILS, only: print_real_matrix
    USE CONSTANTS

    IMPLICIT NONE

    REAL*8, PARAMETER :: zeta_H = 1.24D0    ! STO coefficient correction for H
    REAL*8, PARAMETER :: zeta_He = 2.0925D0 ! STO coefficient correction for He
    REAL*8, PARAMETER :: zeta_O_1 = 7.66D0 ! STO coefficient correction for He
    REAL*8, PARAMETER :: zeta_O_2 = 2.25D0 ! STO coefficient correction for He

    ! ----------------
    ! BASIS SET FOR H2
    ! ----------------
    INTEGER, PARAMETER :: K_H2 = 2 ! Number of basis functions

    INTEGER, dimension(K_H2,3) :: basis_L_H2            ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K_H2,3) :: basis_R_H2             ! Centers of basis set Gaussians
    REAL*8, dimension(K_H2,c) :: basis_D_H2, basis_A_H2 ! Basis set coefficients

    REAL*8, dimension(K_H2,K_H2) :: S_H2 ! Overlap matrix

    ! ------------------
    ! BASIS SET FOR HeH+
    ! ------------------
    INTEGER, PARAMETER :: K_HeH = 2 ! Number of basis functions

    INTEGER, dimension(K_HeH,3) :: basis_L_HeH              ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K_HeH,3) :: basis_R_HeH               ! Centers of basis set Gaussians
    REAL*8, dimension(K_HeH,c) :: basis_D_HeH, basis_A_HeH  ! Basis set coefficients

    REAL*8, dimension(K_HeH,K_HeH) :: S_HeH ! Overlap matrix

    ! -----------------
    ! BASIS SET FOR H2O
    ! -----------------
    INTEGER, PARAMETER :: K_H2O = 7 ! Number of basis functions

    INTEGER, dimension(K_H2O,3) :: basis_L_H2O              ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K_H2O,3) :: basis_R_H2O               ! Centers of basis set Gaussians
    REAL*8, dimension(K_H2O,c) :: basis_D_H2O, basis_A_H2O  ! Basis set coefficients

    REAL*8, dimension(K_H2O,K_H2O) :: S_H2O ! Overlap matrix

    ! -----------
    ! MOLECULE H2
    ! -----------

    basis_R_H2(1,1:3) = (/0.0D0, 0.0D0, 0.0D0/) ! Position of first H atom
    basis_R_H2(2,1:3) = (/1.4D0, 0.0D0, 0.0D0/) ! Position of the second H atom

    basis_L_H2(1,1:3) = (/0, 0, 0/) ! Angular momenta for the first basis function
    basis_L_H2(2,1:3) = (/0, 0, 0/) ! Angular momenta for the second basis function

    basis_D_H2(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D_H2(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)

    basis_A_H2(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A_H2(2,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)

    CALL S_overlap(K_H2,basis_D_H2,basis_A_H2,basis_L_H2,basis_R_H2,S_H2)

    WRITE(*,*) "Overlap matrix S:"
    CALL print_real_matrix(K_H2,K_H2,S_H2)

    ! ------------
    ! MOLECULE He+
    ! ------------

    basis_R_HeH(1,1:3) = (/0.0D0, 0.0D0, 0.0D0/)    ! Position of the H atom
    basis_R_HeH(2,1:3) = (/1.4632D0, 0.0D0, 0.0D0/) ! Position of the He atom

    basis_L_HeH(1,1:3) = (/0, 0, 0/) ! Angular momenta for the first basis function
    basis_L_HeH(2,1:3) = (/0, 0, 0/) ! Angular momenta for the second basis function

    basis_D_HeH(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D_HeH(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)

    basis_A_HeH(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A_HeH(2,1:c) = (/0.109818D0 * zeta_He**2, 0.405771 * zeta_He**2, 2.22766 * zeta_He**2/)

    CALL S_overlap(K_HeH,basis_D_HeH,basis_A_HeH,basis_L_HeH,basis_R_HeH,S_HeH)

    WRITE(*,*) "Overlap matrix S:"
    CALL print_real_matrix(K_HeH,K_HeH,S_HeH)

    ! ------------
    ! MOLECULE H2O
    ! ------------

    basis_R_H2O(1,1:3) = (/0.0D0, 1.43233673D0, -0.96104039D0/)     ! Center of H1 1s
    basis_R_H2O(2,1:3) = (/0.0D0, -1.43233673D0, -0.96104039D0/)    ! Center of H2 1s
    basis_R_H2O(3,1:3) = (/0.0D0, 0.0D0, 0.24026010D0/)             ! Center of O 1s
    basis_R_H2O(4,1:3) = (/0.0D0, 0.0D0, 0.24026010D0/)             ! Center of O 2s
    basis_R_H2O(5,1:3) = (/0.0D0, 0.0D0, 0.24026010D0/)             ! Center of O 2px
    basis_R_H2O(6,1:3) = (/0.0D0, 0.0D0, 0.24026010D0/)             ! Center of O 2py
    basis_R_H2O(7,1:3) = (/0.0D0, 0.0D0, 0.24026010D0/)             ! Center of O 2pz

    basis_L_H2O(1,1:3) = (/0, 0, 0/) ! Angular momenta for H 1s
    basis_L_H2O(2,1:3) = (/0, 0, 0/) ! Angular momenta for H 1s
    basis_L_H2O(3,1:3) = (/0, 0, 0/) ! Angular momenta for O 1s
    basis_L_H2O(4,1:3) = (/0, 0, 0/) ! Angular momenta for O 2s
    basis_L_H2O(5,1:3) = (/1, 0, 0/) ! Angular momenta for O 2px
    basis_L_H2O(6,1:3) = (/0, 1, 0/) ! Angular momenta for O 2py
    basis_L_H2O(7,1:3) = (/0, 0, 1/) ! Angular momenta for O 2pz

    basis_D_H2O(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D_H2O(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D_H2O(3,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D_H2O(4,1:c) = (/0.700115D0, 0.399513D0, -0.0999672D0/)
    basis_D_H2O(5,1:c) = (/0.391957D0, 0.607684D0, 0.1559163D0/)
    basis_D_H2O(6,1:c) = (/0.391957D0, 0.607684D0, 0.1559163D0/)
    basis_D_H2O(7,1:c) = (/0.391957D0, 0.607684D0, 0.1559163D0/)

    basis_A_H2O(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A_H2O(2,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A_H2O(3,1:c) = (/0.109818D0 * zeta_O_1**2, 0.405771 * zeta_O_1**2, 2.22766 * zeta_O_1**2/)
    basis_A_H2O(4,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)
    basis_A_H2O(5,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)
    basis_A_H2O(6,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)
    basis_A_H2O(7,1:c) = (/0.0751386D0 * zeta_O_2**2, 0.231031 * zeta_O_2**2, 0.994203 * zeta_O_2**2/)

    CALL S_overlap(K_H2O,basis_D_H2O,basis_A_H2O,basis_L_H2O,basis_R_H2O,S_H2O)

    WRITE(*,*) "Overlap matrix S:"
    CALL print_real_matrix(K_H2O,K_H2O,S_H2O)



END PROGRAM overlap_test
