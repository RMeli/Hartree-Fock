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

        FUNCTION finite_difference(f1,f2,x1,x2) result(df)

            IMPLICIT NONE

            ! INPUT
            REAL*8, intent(in) :: f1, f2, x1, x2

            ! OUTPUT
            REAL*8 :: df

            df = (f2 - f1) / (x2 - x1)

        END FUNCTION finite_difference

        SUBROUTINE force_idx_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,idx,force,delta)
            ! ------------------------------------------------------------------------
            ! Compute forces on atom IDX by finite difference (computationally costly)
            ! ------------------------------------------------------------------------

            ! INPUT
            INTEGER, intent(in) :: Kf                                   ! Basis set size
            INTEGER, intent(in) :: Ne                                   ! Number of electrons
            INTEGER, intent(in) :: Nn                                   ! Number of nuclei
            INTEGER, intent(in) :: c                                    ! Number of contractions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L             ! Angular momenta of basis set Gaussians
            REAL*8, dimension(Kf,3), intent(in) :: basis_R              ! Centers of basis set Gaussians
            REAL*8, dimension(Kf,c), intent(in) :: basis_D, basis_A     ! Basis set coefficients
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

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dx

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnp,dEp,scf_verbose)

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dx

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnm,dEm,scf_verbose)

            force(idx,1) = - (dEp - dEm) / (2.0D0 * delta) !finite_difference(dEp,dEm,Rnp(idx,1),Rnm(idx,1))

            ! -------------
            ! FORCE ALONG Y
            ! -------------

            dEp = 0.0D0
            dEm = 0.0D0

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dy

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnp,dEp,scf_verbose)

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dy

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnm,dEm,scf_verbose)

            force(idx,2) = - (dEp - dEm) / (2.0D0 * delta) !finite_difference(dEp,dEm,Rnp(idx,2),Rnm(idx,2))

            ! -------------
            ! FORCE ALONG Z
            ! -------------

            dEp = 0.0D0
            dEm = 0.0D0

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dz

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnp,dEp,scf_verbose)

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dz

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnm,dEm,scf_verbose)

            force(idx,3) = - (dEp - dEm) / (2.0D0 * delta) !finite_difference(dEp,dEm,Rnp(idx,3),Rnm(idx,3))

        END SUBROUTINE

        SUBROUTINE force_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,force,delta)
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
            REAL*8, dimension(Nn,3), intent(in) :: Rn                   ! Nuclear positions
            INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
            REAL*8, intent(in) :: delta

            ! OUTPUT
            REAL*8, dimension(Nn,3), intent(out) :: force               ! Vector of forces (on every atom)

            ! INTERMEDIATE VARIABLES
            INTEGER :: i

            DO i = 1, Nn
                CALL force_idx_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,i,force,delta)
            END DO

        END SUBROUTINE force_fd

END MODULE FORCES
