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

PROGRAM boys_test

    USE NUCLEAR, only: boys

    ! ------------------------------------------------------------------------
    !
    ! Source:
    !   On the evaluation of Boys functions using downward  recursion relation
    !   B. A. Mamedov
    !   Journal of Mathematical Chemistry
    !   2004
    !
    ! ------------------------------------------------------------------------

    WRITE(*,*) 8, 16, boys(8,16.0D0)
    WRITE(*,*) 15, 27, boys(15,27.0D0)
    WRITE(*,*) 20, 30, boys(2,30.0D0)
    WRITE(*,*) 25, 13, boys(25,13.0D0)
    WRITE(*,*) 75, 30, boys(75,30.0D0)

    WRITE(*,*)

    WRITE(*,*) 8, 42, boys(8,42.0D0)
    WRITE(*,*) 16, 50, boys(16,50.0D0)
    WRITE(*,*) 21, 56, boys(21,56.0D0)
    WRITE(*,*) 12, 60, boys(12,60.0D0)
    WRITE(*,*) 18, 58, boys(18,58.0D0)

    WRITE(*,*)
    WRITE(*,*) 8, 63, boys(8,63.0D0)
    WRITE(*,*) 33, 85, boys(33,85.0D0)
    WRITE(*,*) 100, 120, boys(100,120.0D0)


END PROGRAM boys_test
