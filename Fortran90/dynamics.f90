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

MODULE DYNAMICS
    ! ----------------------------------------------------------
    ! DYNAMICS
    ! ----------------------------------------------------------
    !
    ! Performs classical and Born-oppenheimer molecular dynamics
    ! using the velocity Verlet integration algorithm.
    !
    !-----------------------------------------------------------

    USE FORCES, only: force_fd
    USE RHF, only: RHF_SCF

    IMPLICIT NONE

    CONTAINS

    ! --------------------
    ! POSITION PROPAGATION
    ! --------------------
    SUBROUTINE position_prop(Nn,mass,pos,vel,force_t,dt)
        ! ----------------------------------------------------------------------
        ! Update the position of every particle according to the force at time t
        ! ----------------------------------------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Nn                           ! Number of nuclei (number of moving particles)
        REAL*8, dimension(Nn) :: mass                       ! Nuclear masses
        REAL*8, intent(in) :: dt                            ! Integration time step
        REAL*8, dimension(Nn,3), intent(in) :: force_t      ! Force acting on nuclei (PES gradient) at time t
        REAL*8, dimension(Nn,3), intent(in) :: vel       ! Nuclear velocities at time t

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: pos       ! Nuclear positions

        ! INTERMEDIATE VARIABLES
        INTEGER :: i                                    ! Loop index

        ! Update positions of all atoms
        DO i = 1, Nn
            pos(i,:) = pos(i,:) + vel(i,:) * dt + 0.5D0 * force_t(i,:) / mass(i) * dt**2.0D0
        END DO

    END SUBROUTINE

    ! --------------------
    ! VELOCITY PROPAGATION
    ! --------------------
    SUBROUTINE velocity_prop(Nn,mass,vel,force_t,force_tdt,dt)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Nn                           ! Number of nuclei (number of moving particles)
        REAL*8, dimension(Nn) :: mass                       ! Nuclear masses
        REAL*8, intent(in) :: dt                            ! Integration time step
        REAL*8, dimension(Nn,3), intent(in) :: force_t      ! Force acting on nuclei (PES gradient) at time t
        REAL*8, dimension(Nn,3), intent(in) :: force_tdt    ! Force acting on nuclei (PES gradient) at time t+dt

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: vel       ! Nuclear velocities

        ! INTERMEDIATE VARIABLES
        INTEGER :: i                                        ! Loop index

        ! Update velocities of all atoms
        DO i = 1, Nn
            vel(i,:) = vel(i,:) + 0.5D0 * (force_t(i,:) + force_tdt(i,:)) / mass(i) * dt
        END DO

    END SUBROUTINE


    ! ----------------------------
    ! CLASSICAL MOLECULAR DYNAMICS
    ! ----------------------------
    SUBROUTINE classical_step(Nn,mass,pos,vel,force_t,dU,dt)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Nn                           ! Number of nuclei (number of moving particles)
        REAL*8, dimension(Nn) :: mass                       ! Nuclear masses
        REAL*8, intent(in) :: dt                            ! Integration time step

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: pos       ! Nuclear positions
        REAL*8, dimension(Nn,3), intent(inout) :: vel       ! Nuclear velocities
        REAL*8, dimension(Nn,3), intent(inout) :: force_t   ! Force acting on nuclei (PES gradient) at time t

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(Nn,3) :: force_tdt                ! Force acting on nuclei (PES gradient) at time t+dt

        ! For classical MD we need an analytical description of the PES
        INTERFACE
            FUNCTION dU(Nn,pos)
                INTEGER, intent(in) :: Nn
                REAL*8, dimension(Nn,3), intent(in) :: pos
                REAL*8, dimension(Nn,3) :: dU
            END FUNCTION
        END INTERFACE

        ! Update positions
        CALL position_prop(Nn,mass,pos,vel,force_t,dt) ! Use forces at time t

        ! Compute forces at time t+dt
        force_tdt = (-1) * dU(Nn,pos) ! POS has been updated by POSITION_PROP

        ! Update velocities
        CALL velocity_prop(Nn,mass,vel,force_t,force_tdt,dt)

        force_t = force_tdt ! Store the force at t+dt for the next step

    END SUBROUTINE classical_step

    ! -----------------------------------
    ! BORN-OPPENHEIMER MOLECULAR DYNAMICS
    ! -----------------------------------
    SUBROUTINE BO_step(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,mass,pos,vel,force_t,dt)
        ! ----------------------------------------------------------------------------------------------------------
        ! Perform a step of Born-Oppenheimer Molecular Dynamics using Hartree-Fock solution of Schrödinger equation.
        ! ----------------------------------------------------------------------------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                                   ! Basis set size
        INTEGER, intent(in) :: c                                    ! Number of contractions
        INTEGER, intent(in) :: Ne                                   ! Number of electrons
        INTEGER, intent(in) :: Nn                                   ! Number of nuclei
        INTEGER, dimension(Kf,3), intent(in) :: basis_L             ! Angular momenta of basis set Gaussians
        REAL*8, dimension(Kf,c), intent(in) :: basis_D, basis_A     ! Basis set coefficients
        INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
        INTEGER, dimension(Kf), intent(in) :: basis_idx             ! Basis set atom index (link between basis set function and atoms)
        REAL*8, dimension(Nn) :: mass                               ! Nuclear masses
        REAL*8, intent(in) :: dt                                    ! Integration time step

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: pos               ! Nuclear positions
        REAL*8, dimension(Nn,3), intent(inout) :: vel               ! Nuclear velocities
        REAL*8, dimension(Nn,3), intent(inout) :: force_t           ! Force acting on nuclei (PES gradient) at time t
        REAL*8, dimension(Kf,3), intent(inout) :: basis_R           ! Basis set centers

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(Nn,3) :: force_tdt                        ! Force acting on nuclei (PES gradient) at time t+dt
        INTEGER i, idx

        ! PARAMETERS
        REAL*8, PARAMETER :: delta = 1e-6                           ! Finite difference step for forces

        ! ------------------------
        ! Update nuclear position
        ! ------------------------

        CALL position_prop(Nn,mass,pos,vel,force_t,dt) ! Use forces at time t

        ! -----------------------------------------------------------
        ! Update basis set centers according to new nuclear positions
        ! -----------------------------------------------------------

        DO i = 1, Kf
            idx = basis_idx(i)

            basis_R(i,:) = pos(idx,:) ! Displace basis set center at new atomic position
        END DO

        ! ---------------------------------------
        ! Compute forces at new nuclear positions
        ! ---------------------------------------

        CALL force_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,basis_idx,Zn,pos,force_tdt,delta)

        ! -------------------------
        ! Update nuclear velocities
        ! -------------------------

        CALL velocity_prop(Nn,mass,vel,force_t,force_tdt,dt)

        ! ----------------------------------
        ! Store new forces for the next step
        ! ----------------------------------

        force_t = force_tdt ! Store the force at t+dt for the next step

    END SUBROUTINE BO_step

END MODULE DYNAMICS
