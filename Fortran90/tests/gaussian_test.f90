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

PROGRAM gaussian_test

    USE GAUSSIAN

    IMPLICIT NONE

    REAL*8 :: aa = 1.0D0, bb = 2.0D0 ! Gaussian exponential coefficients
    REAL*8, dimension(3) :: Ra = (/0.0D0,0.0D0,0.0D0/), Rb = (/1.0D0,1.0D0,1.0D0/)

    REAL :: Na, Nb

    Na = norm(0,0,0,aa)
    Nb = norm(0,0,0,bb)

    WRITE(*,*) "Na = ", Na
    WRITE(*,*) "Nb = ", Nb

END PROGRAM gaussian_test
