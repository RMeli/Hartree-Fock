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

PROGRAM HF

    USE INPUT
    USE RHF
    USE UHF

    IMPLICIT NONE

    CHARACTER (len=20) :: fname ! Input filename
    INTEGER :: nargs            ! Number of input arguments

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO)

    ! System and basis set informations
    CHARACTER (len=4) :: calculation                   ! Calculation
    REAL*8, allocatable, dimension(:,:) :: Rn           ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Zn            ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
    INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D      ! Contaction linear coefficients
    INTEGER, allocatable, dimension(:) :: basis_idx     ! Basis set atomic index

    REAL*8 :: final_E ! Total converged energy

    ! --------------------------
    ! READ COMMAND LINE ARGUMENT
    ! --------------------------

    nargs = IARGC()

    IF (nargs .LT. 1) THEN
        WRITE(*,*) "NO INPUT FILE!"
        CALL EXIT(-1)
    END IF

    CALL GETARG(1, fname)

    WRITE(*,*) fname

    ! ------------------------------------------------
    ! LOAD SYSTEM AND BASIS SET INFORMATIONS FROM FILE
    ! ------------------------------------------------

    CALL load(fname,calculation,Ne,Nn,K,c,Rn,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx)

    ! ------------------------
    ! TOTAL ENERGY CALCULATION
    ! ------------------------

    IF (calculation .EQ. "RHF ") THEN

        WRITE(*,*) "----------------------------"
        WRITE(*,*) "RHF TOTAL ENERGY CALCULATION"
        WRITE(*,*) "----------------------------"

        CALL RHF_SCF(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,.TRUE.)

    ELSE IF (calculation .EQ. "UHF ") THEN

        WRITE(*,*) "----------------------------"
        WRITE(*,*) "UHF TOTAL ENERGY CALCULATION"
        WRITE(*,*) "----------------------------"

        CALL UHF_SCF(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,.TRUE.)

    END IF

    DEALLOCATE(Rn,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx) ! Deallocate allocated memory

END PROGRAM HF
