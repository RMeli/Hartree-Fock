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

PROGRAM electronic_test

    USE ELECTRONIC
    USE OUTPUT, only: print_ee_list

    IMPLICIT NONE

    ! TODO allows flexibility
    INTEGER, PARAMETER :: c = 3 ! Number of contractions (STO-3G)

    REAL*8, PARAMETER :: zeta_H = 1.24D0    ! STO coefficient correction for H
    REAL*8, PARAMETER :: zeta_He = 2.0925D0 ! STO coefficient correction for He
    REAL*8, PARAMETER :: zeta_O_1 = 7.66D0 ! STO coefficient correction for He
    REAL*8, PARAMETER :: zeta_O_2 = 2.25D0 ! STO coefficient correction for He

    ! ------------------
    ! BASIS SET FOR HeH+
    ! ------------------
    INTEGER, PARAMETER :: K_HeH = 2     ! Number of basis functions
    INTEGER, PARAMETER :: Nn_HeH = 2    ! Number of nuclei

    INTEGER, dimension(K_HeH,3) :: basis_L_HeH              ! Angular momenta of basis set Gaussians
    REAL*8, dimension(K_HeH,3) :: basis_R_HeH               ! Centers of basis set Gaussians
    REAL*8, dimension(K_HeH,c) :: basis_D_HeH, basis_A_HeH  ! Basis set coefficients

    REAL*8, dimension(K_HeH,K_HeH) :: Vn1_HeH   ! Nucleus-electron interaction matrix
    REAL*8, dimension(K_HeH,K_HeH) :: Vn2_HeH   ! Nucleus-electron interaction matrix

    REAL*8, dimension(Nn_HeH,3) :: R_HeH    ! Nuclear positions
    INTEGER, dimension(Nn_HeH) :: Zn_HeH    ! Nuclear charges

    ! -------------
    ! MOLECULE HeH+
    ! -------------

    REAL*8, dimension(K_HeH,K_HeH,K_HeH,K_HeH) :: ee_HeH

    R_HeH(1,1:3) = (/0.0D0, 0.0D0, 0.0D0/)       ! Position of first H atom
    R_HeH(2,1:3) = (/1.4632D0, 0.0D0, 0.0D0/)    ! Position of the second H atom

    Zn_HeH = (/1, 2/)

    basis_R_HeH(1,1:3) = R_HeH(1,1:3)    ! Position of the H atom
    basis_R_HeH(2,1:3) = R_HeH(2,1:3)    ! Position of the He atom

    basis_L_HeH(1,1:3) = (/0, 0, 0/) ! Angular momenta for the first basis function
    basis_L_HeH(2,1:3) = (/0, 0, 0/) ! Angular momenta for the second basis function

    basis_D_HeH(1,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)
    basis_D_HeH(2,1:c) = (/0.444635D0, 0.535328D0, 0.154329D0 /)

    basis_A_HeH(1,1:c) = (/0.109818D0 * zeta_H**2, 0.405771 * zeta_H**2, 2.22766 * zeta_H**2/)
    basis_A_HeH(2,1:c) = (/0.109818D0 * zeta_He**2, 0.405771 * zeta_He**2, 2.22766 * zeta_He**2/)

    WRITE(*,*) "#############"
    WRITE(*,*) "Molecule HeH+"
    WRITE(*,*) "#############"

    CALL ee_list(K_HeH,basis_D_HeH,basis_A_HeH,basis_L_HeH,basis_R_HeH,ee_HeH)

    WRITE(*,*) "HeH+ electron-electron integrals:"
    CALL print_ee_list(K_HeH,ee_HeH)

END PROGRAM electronic_test
