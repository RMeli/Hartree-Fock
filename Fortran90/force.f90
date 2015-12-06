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

MODULE FORCES
    ! -------------------------------------------
    ! FORCES
    ! -------------------------------------------
    !
    ! Compute electronic forces acting on nuclei.
    !
    !--------------------------------------------

    USE RHF

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE displace_basis(Kf,basis_R,basis_idx,idx,dir,delta,basis_RR)
            ! -----------------------------------------------------
            ! Displace basis functions on atom IDX in direction DIR
            ! -----------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                                   ! Basis set size
            REAL*8, dimension(Kf,3), intent(in) :: basis_R              ! Centers of basis set Gaussians
            INTEGER, dimension(Kf), intent(in):: basis_idx              ! Basis set atom index
            INTEGER, intent(in) :: idx                                  ! Atom index (displace basis functions at atom idx)
            INTEGER, intent(in) :: dir                                  ! Direction of displacement (1=x, 2=y, 3=z)
            REAL*8, intent(in) :: delta                                 ! Displacement

            ! INTERMEDIATE VARIABLES
            INTEGER:: i                                                 ! Loop index

            ! OUTPUT
            REAL*8, dimension(Kf,3) :: basis_RR                         ! Displaced basis set

            basis_RR = basis_R ! Initialize displaced basis centers at original ones

            DO i = 1, Kf
                IF (basis_idx(i) .EQ. idx) THEN
                    ! Displace basis set center along DIR
                    basis_RR(i,dir) = basis_R(i,dir) + delta
                END IF
            END DO

        END SUBROUTINE

        SUBROUTINE force_idx_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,Rn,idx,force,delta)
            ! ------------------------------------------------------------------------
            ! Compute forces on atom IDX by finite difference (computationally costly)
            ! ------------------------------------------------------------------------

            ! TODO: Subroutine that computes force along given direction (x, y or z)

            ! INPUT
            INTEGER, intent(in) :: Kf                                   ! Basis set size
            INTEGER, intent(in) :: Ne                                   ! Number of electrons
            INTEGER, intent(in) :: Nn                                   ! Number of nuclei
            INTEGER, intent(in) :: c                                    ! Number of contractions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L             ! Angular momenta of basis set Gaussians
            REAL*8, dimension(Kf,3), intent(in) :: basis_R              ! Centers of basis set Gaussians
            REAL*8, dimension(Kf,c), intent(in) :: basis_D, basis_A     ! Basis set coefficients
            INTEGER, dimension(Kf), intent(in):: basis_idx            ! Basis set atom index
            REAL*8, dimension(Nn,3), intent(in) :: Rn                   ! Nuclear positions
            INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
            INTEGER, intent(in) :: idx                                  ! Atom index (compute forces for atom idx)
            REAL*8, intent(in) :: delta

            ! OUTPUT
            REAL*8, dimension(Nn,3), intent(out) :: force               ! Vector of forces (on every atom)

            ! INTERMEDIATE VARIABLES
            REAL*8, dimension(3) :: dx
            REAL*8, dimension(3) :: dy
            REAL*8, dimension(3) :: dz
            REAL*8 :: dEp, dEm
            REAL*8, dimension(Nn,3) :: Rnp, Rnm
            REAL*8, dimension(Kf,3) :: basis_RR                         ! Displaced basis set

            ! PARAMETERS
            LOGICAL, PARAMETER :: scf_verbose = .FALSE.

            ! --------------
            ! INITIALIZATION
            ! --------------

            dx = (/delta,0.0D0,0.0D0/)
            dy = (/0.0D0,delta,0.0D0/)
            dz = (/0.0D0,0.0D0,delta/)

            ! -------------
            ! FORCE ALONG X
            ! -------------

            dEp = 0.0D0
            dEm = 0.0D0

            ! Displace atom and basis set at new atomic position

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dx

            CALL displace_basis(Kf,basis_R,basis_idx,idx,1,delta,basis_RR)

            ! Total energy calculation at new atomic position

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_RR,Zn,Rnp,dEp,scf_verbose)

            ! Displace atom and basis set at new atomic position

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dx

            CALL displace_basis(Kf,basis_R,basis_idx,idx,1,-delta,basis_RR)

            ! Total energy calculation at new atomic position

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_RR,Zn,Rnm,dEm,scf_verbose)

            ! Compute force along x

            force(idx,1) = - (dEp - dEm) / (2.0D0 * delta) !finite_difference(dEp,dEm,Rnp(idx,1),Rnm(idx,1))

            ! -------------
            ! FORCE ALONG Y
            ! -------------

            dEp = 0.0D0
            dEm = 0.0D0

            ! Displace atom and basis set at new atomic position

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dy

            CALL displace_basis(Kf,basis_R,basis_idx,idx,2,delta,basis_RR)

            ! Total energy calculation at new atomic position

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_RR,Zn,Rnp,dEp,scf_verbose)

            ! Displace atom and basis set at new atomic position

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dy

            CALL displace_basis(Kf,basis_R,basis_idx,idx,2,-delta,basis_RR)

            ! Total energy calculation at new atomic position

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_RR,Zn,Rnm,dEm,scf_verbose)

            ! Compute force along y

            force(idx,2) = - (dEp - dEm) / (2.0D0 * delta) !finite_difference(dEp,dEm,Rnp(idx,2),Rnm(idx,2))

            ! -------------
            ! FORCE ALONG Z
            ! -------------

            dEp = 0.0D0
            dEm = 0.0D0

            ! Displace atom and basis set at new atomic position

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dz

            CALL displace_basis(Kf,basis_R,basis_idx,idx,3,delta,basis_RR)

            ! Total energy calculation at new atomic position

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_RR,Zn,Rnp,dEp,scf_verbose)

            ! Displace atom and basis set at new atomic position

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dz

            CALL displace_basis(Kf,basis_R,basis_idx,idx,3,-delta,basis_RR)

            ! Total energy calculation at new atomic position

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_RR,Zn,Rnm,dEm,scf_verbose)

            ! Compute force along z

            force(idx,3) = - (dEp - dEm) / (2.0D0 * delta) !finite_difference(dEp,dEm,Rnp(idx,3),Rnm(idx,3))

        END SUBROUTINE

        SUBROUTINE force_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,Rn,force,delta)
            ! -------------------------------------------------------------------------
            ! Compute forces on all atoms by finite difference (computationally costly)
            ! -------------------------------------------------------------------------

            ! INPUT
            INTEGER, intent(in) :: Kf                                   ! Basis set size
            INTEGER, intent(in) :: Ne                                   ! Number of electrons
            INTEGER, intent(in) :: Nn                                   ! Number of nuclei
            INTEGER, intent(in) :: c                                    ! Number of contractions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L             ! Angular momenta of basis set Gaussians
            REAL*8, dimension(Kf,3), intent(in) :: basis_R              ! Centers of basis set Gaussians
            REAL*8, dimension(Kf,c), intent(in) :: basis_D, basis_A     ! Basis set coefficients
            INTEGER, dimension(Kf), intent(in):: basis_idx            ! Basis set atom index
            REAL*8, dimension(Nn,3), intent(in) :: Rn                   ! Nuclear positions
            INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
            REAL*8, intent(in) :: delta

            ! OUTPUT
            REAL*8, dimension(Nn,3), intent(out) :: force               ! Vector of forces (on every atom)

            ! INTERMEDIATE VARIABLES
            INTEGER :: i

            ! Loop on every atom of the system
            DO i = 1, Nn
                CALL force_idx_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,Rn,i,force,delta)
            END DO

        END SUBROUTINE force_fd

END MODULE FORCES
