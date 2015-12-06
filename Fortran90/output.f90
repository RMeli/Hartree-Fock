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

        ! --------------
        ! PRINT MATRICES
        ! --------------

        SUBROUTINE print_real_matrix(r,c,M)
            ! ----------------------------------------------
            ! Print RxC matrix M containing REAL*8 elements.
            ! ----------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: r, c
            REAL*8, dimension(r,c), intent(in) :: M

            ! VARIABLES
            INTEGER :: i

            DO i = 1, r
                WRITE(*,'(20G16.6)') M(i,1:c)
            END DO

        END SUBROUTINE print_real_matrix



        SUBROUTINE print_integer_matrix(r,c,M)
            ! -----------------------------------------------
            ! Print RxC matrix M containing INTEGER elements.
            ! -----------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: r, c
            INTEGER, dimension(r,c), intent(in) :: M

            ! VARIABLES
            INTEGER :: i

            DO i = 1, r
                WRITE(*,'(20I16)') M(i,1:c)
            END DO

        END SUBROUTINE print_integer_matrix




        ! ---------------------------
        ! ELECTRON-ELECTRON INTEGRALS
        ! ---------------------------
        SUBROUTINE print_ee_list(Kf,ee)
            ! --------------------------------------
            ! Print electron-electron integral list.
            ! --------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf               ! Basis set size
            REAL*8, dimension(Kf,Kf,Kf,Kf) :: ee    ! List of electron-electron integrals

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j, k, l   ! Loop indices

            DO i = 1,Kf
                DO j = 1,Kf
                    DO k = 1,Kf
                        DO l = 1,Kf

                            WRITE(*,*) "(", i, j, k, l, ") ", ee(i,j,k,l)

                        END DO ! l
                    END DO ! k
                END DO ! j
            END DO ! i

        END SUBROUTINE print_ee_list



        ! -------------
        ! APPEND TO XYZ
        ! -------------
        SUBROUTINE append_xyz(Nn,atoms,pos,unit)
            ! -----------------------------------------------
            ! Append nuclear positions to XYZ trajectory file
            ! -----------------------------------------------

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
