PROGRAM overlap_test

    USE OVERLAP
    USE UTILS, only: print_real_matrix

    IMPLICIT NONE

    INTEGER, PARAMETER :: K = 2 ! Number of basis functions
    INTEGER, PARAMETER :: c = 3 ! Number of contractions

    INTEGER :: i, j

    INTEGER, dimension(K,3) :: basis_L          ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K,3) :: basis_R           ! Centers of basis set Gaussians
    REAL*8, dimension(K,c) :: basis_D, basis_A  ! Basis set coefficients

    REAL*8, dimension(K,K) :: S ! Overlap matrix

    REAL*8, PARAMETER :: zeta_H = 1.24D0    ! STO coefficient correction for H
    REAL*8, PARAMETER :: zeta_He = 2.0925D0 ! STO coefficient correction for He

    ! -----------
    ! MOLECULE H2
    ! -----------

    basis_R(1,1:3) = (/0.0D0, 0.0D0, 0.0D0/) ! Position of first H atom
    basis_R(2,1:3) = (/1.4D0, 0.0D0, 0.0D0/) ! Position of the second H atom

    basis_L(1,1:3) = (/0, 0, 0/) ! Angular momenta for the first basis function
    basis_L(2,1:3) = (/0, 0, 0/) ! Angular momenta for the second basis function

    basis_D(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)

    basis_A(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A(2,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)

    CALL S_overlap(K,basis_D,basis_A,basis_L,basis_R,S)

    WRITE(*,*) "Overla pamtrix S:"
    CALL print_real_matrix(K,K,S)

    ! ------------
    ! MOLECULE He+
    ! ------------

    basis_R(1,1:3) = (/0.0D0, 0.0D0, 0.0D0/)    ! Position of the H atom
    basis_R(2,1:3) = (/1.4632D0, 0.0D0, 0.0D0/) ! Position of the He atom

    basis_L(1,1:3) = (/0, 0, 0/) ! Angular momenta for the first basis function
    basis_L(2,1:3) = (/0, 0, 0/) ! Angular momenta for the second basis function

    basis_D(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)

    basis_A(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A(2,1:c) = (/0.109818D0 * zeta_He**2, 0.405771 * zeta_He**2, 2.22766 * zeta_He**2/)

    CALL S_overlap(K,basis_D,basis_A,basis_L,basis_R,S)

    WRITE(*,*) "Overla pamtrix S:"
    CALL print_real_matrix(K,K,S)

    ! ------------
    ! MOLECULE H2O
    ! ------------

END PROGRAM overlap_test
