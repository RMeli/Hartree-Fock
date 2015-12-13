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

PROGRAM BOMD_test

    USE INPUT
    USE DYNAMICS
    USE OUTPUT

    IMPLICIT NONE

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO)

    ! System and basis set informations
    REAL*8, allocatable, dimension(:,:) :: pos              ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Zn                ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R          ! Basis functions' centers
    INTEGER, allocatable, dimension(:) :: basis_idx         ! Basis set atomic index
    INTEGER, allocatable, dimension(:,:) :: basis_L         ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A          ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D          ! Conttaction linear coefficients

    REAL*8, allocatable, dimension(:,:) :: F                ! Forces

    CHARACTER (len=5) :: calculation                        ! Type of calculation
    CHARACTER (len=3) :: method                             ! Type of calculation

    REAL*8, allocatable, dimension(:,:) :: vel              ! Nuclear velocities
    REAL*8, allocatable, dimension(:) :: mass             ! Nuclear mass
    CHARACTER(len=2), allocatable, dimension(:) :: atoms    ! Atom name

    ! BOMD parameters

    INTEGER :: i
    INTEGER, PARAMETER :: steps = 1000
    REAL*8, PARAMETER:: dt = 1e-3

    ! ------------------------------------------------
    ! LOAD SYSTEM AND BASIS SET INFORMATIONS FROM FILE
    ! ------------------------------------------------

    CALL load("tests/H2_f_stretched.in",calculation,method,Ne,Nn,K,c,pos,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx)

    ALLOCATE(F(Nn,3),mass(Nn),vel(Nn,3),atoms(Nn))

    vel(1,:) = (/0.0D0, 0.0D0, 0.0D0/)
    vel(2,:) = (/0.0D0, 0.0D0, 0.0D0/)

    mass = (/1.008D0, 1.0008D0/)

    atoms = (/"H","H"/)

    ! ----
    ! BOMD
    ! ----

    ! Initial force
    CALL force_fd(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,pos,F,1D-4,"UHF")

    OPEN(unit=123,file="BOMD.xyz",form="formatted",status="new",action="write") ! Open file 123

    DO i = 1, steps
        WRITE(*,*) "Step", i
        WRITE(*,*) "Distance: ", DSQRT(DOT_PRODUCT(pos(1,:)-pos(2,:),pos(1,:)-pos(2,:)))
        CALL append_xyz(Nn,atoms,pos,123)
        CALL BO_step(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,mass,pos,vel,F,dt,"UHF")
    END DO

    CLOSE(unit=123) ! Close file 123

    ! ---
    ! END
    ! ---

    DEALLOCATE(pos,Zn,basis_R,basis_L,basis_A,basis_D,F,vel,mass) ! Deallocate allocated memory

END PROGRAM BOMD_test
