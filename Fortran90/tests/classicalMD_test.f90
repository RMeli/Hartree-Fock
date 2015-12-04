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

PROGRAM classicalMD

    USE DYNAMICS, only: classical_step
    USE OUTPUT
    USE UTILS

    IMPLICIT NONE

    INTEGER, PARAMETER :: steps = 1000
    REAL*8, PARAMETER:: dt = 1e-3

    INTEGER :: i ! Loop index

    INTEGER, PARAMETER :: Nn = 3
    CHARACTER(len=2), dimension(Nn) :: atoms
    REAL*8, dimension(Nn,3) :: pos
    REAL*8, dimension(Nn,3) :: vel
    REAL*8, dimension(Nn,3) :: force
    REAL*8, dimension(Nn) :: mass

    atoms = (/"H","H","H"/)

    pos(1,:) = (/0.0D0, 0.0D0, 0.0D0/)
    pos(2,:) = (/1.2*1.0D0, 0.0D0, 0.0D0/)
    pos(3,:) = (/0.5*1.0D0, 1.3*1.0D0, 0.0D0/)

    vel(1,:) = (/0.0D0, 0.0D0, 0.0D0/)
    vel(2,:) = (/0.0D0, 0.0D0, 0.0D0/)
    vel(3,:) = (/0.0D0, 0.0D0, 0.0D0/)

    force = (-1) * dU(Nn,pos) ! Compute the initial force

    mass = (/1.0D0, 1.0D0, 10.0D0/)

    OPEN(unit=123,file="classicalMD.xyz",form="formatted",status="new",action="write") ! Open file 123

    DO i = 1, steps
        CALL append_xyz(Nn,atoms,pos,123)
        CALL classical_step(Nn,mass,pos,vel,force,dU,dt)
    END DO

    CLOSE(unit=123) ! Close file 123

    CONTAINS
        FUNCTION dU(Nn,pos)

            ! INPUT
            INTEGER, intent(in) :: Nn                   ! Number of particles
            REAL*8, dimension(Nn,3), intent(in) :: pos  ! Particle's position

            ! OUTPUT
            REAL*8, dimension(Nn,3) :: dU               ! Gradient of the PES

            ! PARAMETERS
            REAL*8, PARAMETER :: k = 500 ! Spring constant (TODO: different force constant for different interactions)
            REAL*8, PARAMETER :: d = 1 ! Equilibrium distance (TODO: different equilibrium distance for different interactions)

            ! INTERMEDIATE VARIABLES
            REAL*8 :: r, dx, dy, dz
            INTEGER :: i, j ! Loop indices

            dU(:,:) = 0.0D0

            DO i = 1, Nn-1
                DO j = i+1, Nn
                    dx = pos(i,1) - pos(j,1)
                    dy = pos(i,2) - pos(j,2)
                    dz = pos(i,3) - pos(j,3)

                    r = DSQRT(dx**2 + dy**2 + dz**2)

                    dU(i,1) = dU(i,1) - k * (d - r) / r * dx
                    dU(i,2) = dU(i,2) - k * (d - r) / r * dy
                    dU(i,3) = dU(i,3) - k * (d - r) / r * dz

                    dU(j,1) = dU(j,1) + k * (d - r) / r * dx
                    dU(j,2) = dU(j,2) + k * (d - r) / r * dy
                    dU(j,3) = dU(j,3) + k * (d - r) / r * dz
                END DO
            END DO


        END FUNCTION

END PROGRAM classicalMD
