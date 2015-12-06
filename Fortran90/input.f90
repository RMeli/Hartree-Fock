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

MODULE INPUT

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE load(fname,calculation,Ne,Nn,K,c,Rn,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx)

            IMPLICIT NONE

            ! OUTPUT
            CHARACTER (len=4), intent(out):: calculation    ! Calculation name
            INTEGER, intent(out) :: Ne                      ! Number of electrons
            INTEGER, intent(out) :: Nn                      ! Number of nuclei
            INTEGER, intent(out) :: K                       ! Basis set size
            INTEGER, intent(out) :: c                       ! Contractions (TODO: split valence basis set)
            CHARACTER(len=*), intent(in) :: fname           ! Input file name

            ! INTERMEDIATE VARIABLES
            INTEGER :: i ! Loop index

            ! ALLOCATABLE
            REAL*8, allocatable, dimension(:,:) :: Rn           ! Atomic potisions
            INTEGER, allocatable, dimension(:) :: Zn            ! Atomic chatges
            REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
            INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
            REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
            REAL*8, allocatable, dimension(:,:) :: basis_D      ! Conttaction linear coefficients
            INTEGER, allocatable, dimension(:) :: basis_idx     ! Basis set atom index

            ! Open file containing system and basis set informations
            OPEN(unit=100,file=fname,form="formatted",status="old",action="read")

            READ(100,*) calculation
            READ(100,*) Ne
            READ(100,*) Nn
            READ(100,*) K
            READ(100,*) C

            ! Allocate arrays once we know system's and basis set's specifications (Ne,Nn,K,c)
            ALLOCATE(Rn(Nn,3))
            ALLOCATE(Zn(Nn))
            ALLOCATE(basis_R(K,3))
            ALLOCATE(basis_L(K,3))
            ALLOCATE(basis_A(K,c))
            ALLOCATE(basis_D(K,c))
            ALLOCATE(basis_idx(K))

            ! Read atomic positions
            DO i = 1, Nn
                READ(100,*) Rn(i,:), Zn(i)
            END DO

            ! Read basis set informations
            DO i = 1, K
                READ(100,*) basis_R(i,:), basis_L(i,:), basis_A(i,:), basis_D(i,:), basis_idx(i)
            END DO

            CLOSE(unit=100) ! Close the file

        END SUBROUTINE load

END MODULE INPUT
