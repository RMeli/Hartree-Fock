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

PROGRAM HF_N2

    USE INPUT
    USE RHF
    USE FORCES

    IMPLICIT NONE

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO)

    ! System and basis set informations
    REAL*8, allocatable, dimension(:,:) :: Rn           ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Zn            ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
    INTEGER, allocatable, dimension(:) :: basis_idx     ! Basis set atomic index
    INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D      ! Conttaction linear coefficients

    REAL*8, allocatable, dimension(:,:) :: F            ! Forces

    CHARACTER (len=4) :: calculation                    ! Type of calculation

    ! ------------------------------------------------
    ! LOAD SYSTEM AND BASIS SET INFORMATIONS FROM FILE
    ! ------------------------------------------------

    CALL load("tests/H2O_f.in",calculation,Ne,Nn,K,c,Rn,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx)

    ! -----
    ! FORCE
    ! -----

    ALLOCATE(F(Nn,3))

    ! Force on atom 1
    CALL force_fd(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,Rn,F,1D-4)

    WRITE(*,*) "######"
    WRITE(*,*) "FORCES"
    WRITE(*,*) "######"
    CALL print_real_matrix(Nn,3,F)

    ! ---
    ! END
    ! ---

    DEALLOCATE(Rn,Zn,basis_R,basis_L,basis_A,basis_D,F) ! Deallocate allocated memory

END PROGRAM HF_N2
