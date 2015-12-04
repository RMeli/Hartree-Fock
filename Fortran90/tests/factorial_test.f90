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

PROGRAM factorial_test

    USE FACT

    IMPLICIT NONE

    INTEGER :: i
    INTEGER, PARAMETER :: N = 10

    DO i = 0,N
        WRITE(*,*) i, "! =", factorial(i)
    END DO

    DO i = -1,N
        WRITE(*,*) i, "! =", factorial2(i)
    END DO

    DO i = 0,N
        WRITE(*,*) "binom(10,", i, ") =", binom(10,i)
    END DO
END PROGRAM factorial_test
