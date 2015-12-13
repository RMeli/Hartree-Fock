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
    USE FORCES
    USE MULTIPOLE

    IMPLICIT NONE

    CHARACTER (len=50) :: fname ! Input filename
    INTEGER :: nargs            ! Number of input arguments

    INTEGER :: i ! Loop index

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO)

    ! System and basis set informations
    CHARACTER (len=5) :: calculation                    ! Calculation
    CHARACTER (len=3) :: method                         ! Method
    REAL*8, allocatable, dimension(:,:) :: Rn           ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Zn            ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
    INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D      ! Contaction linear coefficients
    INTEGER, allocatable, dimension(:) :: basis_idx     ! Basis set atomic index

    REAL*8, allocatable, dimension(:,:) :: P            ! Final density matrix

    REAL*8, dimension(3) :: mu                          ! Dipole
    REAL*8, allocatable, dimension(:,:) :: force        ! Forces acting on nuclei

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

    CALL load(fname,calculation,method,Ne,Nn,K,c,Rn,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx)

    ! ------------------------
    ! TOTAL ENERGY CALCULATION
    ! ------------------------

    IF (calculation .EQ. "SP") THEN

        ALLOCATE(P(K,K))

        ! -----------------------
        ! RESTRICTED HARTREE-FOCK
        ! -----------------------
        IF (method .EQ. "RHF ") THEN

            IF (MODULO(Ne,2) .NE. 0) THEN
                WRITE(*,*) "ERROR: Restricted calculation with an odd number of electrons."
                CALL EXIT(-1)
            END IF

            WRITE(*,*) "----------------------------"
            WRITE(*,*) "RHF TOTAL ENERGY CALCULATION"
            WRITE(*,*) "----------------------------"
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)



            CALL RHF_DIIS(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,P,.TRUE.)

            CALL dipole(K,c,Nn,basis_D,basis_A,basis_L,basis_R,P,Rn,Zn,mu)

            WRITE(*,*)
            WRITE(*,*) "Dipole:", mu

        ! -------------------------
        ! UNRESTRICTED HARTREE-FOCK
        ! -------------------------
        ELSE IF (method .EQ. "UHF ") THEN

            WRITE(*,*) "----------------------------"
            WRITE(*,*) "UHF TOTAL ENERGY CALCULATION"
            WRITE(*,*) "----------------------------"
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)

            CALL UHF_DIIS(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,P,.TRUE.)

            CALL dipole(K,c,Nn,basis_D,basis_A,basis_L,basis_R,P,Rn,Zn,mu)

            WRITE(*,*)
            WRITE(*,*) "Dipole:", mu

        ELSE
            WRITE(*,*) "ERROR: Invalid method for single point (SP) calculation."
            CALL EXIT(-1)
        END IF

        DEALLOCATE(P)

    ! -----------------
    ! FORCE CALCULATION
    ! -----------------

    ELSE IF (calculation .EQ. "FORCE") THEN

        ALLOCATE(force(Nn,3))

        ! -----------------------
        ! RESTRICTED HARTREE-FOCK
        ! -----------------------
        IF (method .EQ. "RHF ") THEN

            IF (MODULO(Ne,2) .NE. 0) THEN
                WRITE(*,*) "ERROR: Restricted calculation with an odd number of electrons."
                CALL EXIT(-1)
            END IF

            WRITE(*,*) "-----------------------"
            WRITE(*,*) "FORCE CALCULATION (RHF)"
            WRITE(*,*) "-----------------------"
            WRITE(*,*)

            CALL force_fd(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,Rn,force,1D-4,method)

            WRITE(*,*) "Forces:"

            DO i = 1, Nn
                WRITE(*,*) "Atom", i, force(i,:)
            END DO


        ! -------------------------
        ! UNRESTRICTED HARTREE-FOCK
        ! -------------------------
        ELSE IF (method .EQ. "UHF") THEN

            WRITE(*,*) "-----------------------"
            WRITE(*,*) "FORCE CALCULATION (UHF)"
            WRITE(*,*) "-----------------------"
            WRITE(*,*)

            CALL force_fd(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,Rn,force,1D-4,method)

            WRITE(*,*) "Forces:"

            DO i = 1, Nn
                WRITE(*,*) "Atom", i, force(i,:)
            END DO

        ELSE
            WRITE(*,*) "ERROR: Invalid method for force (FORCE) calculation."
            CALL EXIT(-1)
        END IF

        DEALLOCATE(force)

    ELSE
        WRITE(*,*) "ERROR: Invalid calculation."
    END IF

    DEALLOCATE(Rn,Zn,basis_R,basis_L,basis_A,basis_D,basis_idx) ! Deallocate allocated memory

END PROGRAM HF
