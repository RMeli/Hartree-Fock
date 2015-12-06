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

PROGRAM eigs_test

    USE LA, only: EIGS
    USE OUTPUT, only: print_real_matrix

    INTEGER, PARAMETER :: d = 2

    REAL*8, dimension(d,d) :: A
    REAL*8, dimension(d,d) :: V
    REAL*8, dimension(d) :: l

    A(1,:) = (/-1.53866, -0.515838/)
    A(2,:) = (/-0.515838, -2.43973/)

    CALL EIGS(d,A,V,l)

    WRITE(*,*) "Matrix A"
    CALL print_real_matrix(d,d,A)

    WRITE(*,*) "Matrix V"
    CALL print_real_matrix(d,d,V)

    WRITE(*,*) "Vector l"
    CALL print_real_matrix(d,1,l)

END PROGRAM eigs_test
