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

MODULE OUTPUT

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE append_xyz(Nn,atoms,pos,unit)

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Nn                               ! Number of particles
            CHARACTER(len=2), dimension(Nn), intent(in) :: atoms    ! Atom type
            REAL*8, dimension(Nn,3), intent(in) :: pos              ! Particles' position
            INTEGER, intent(in) :: unit                             ! Output unit

            ! INTERMEDIATE VARIABLES
            INTEGER :: i                                    ! Loop index

            WRITE(unit,*) Nn
            WRITE(unit,*) ! TODO: add step number
            DO i = 1, Nn
                WRITE(unit,*) atoms(i), pos(i,:)
            END DO

        END SUBROUTINE

END MODULE OUTPUT
