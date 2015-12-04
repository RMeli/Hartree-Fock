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

PROGRAM input_test

    USE INPUT
    USE UTILS

    IMPLICIT NONE

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO: split valence basis set)

    ! System and basis set informations
    REAL*8, allocatable, dimension(:,:) :: R            ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Z             ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
    INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D      ! Conttaction linear coefficients

    CALL load("tests/N2_f.in",Ne,Nn,K,c,R,Z,basis_R,basis_L,basis_A,basis_D)

    CALL print_real_matrix(Nn,3,R)
    CALL print_integer_matrix(Nn,1,Z)
    CALL print_real_matrix(K,3,basis_R)
    CALL print_integer_matrix(K,3,basis_L)
    CALL print_real_matrix(K,c,basis_A)
    CALL print_real_matrix(K,c,basis_D)

    DEALLOCATE(R,Z,basis_R,basis_L,basis_A,basis_D)


END PROGRAM input_test
