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

MODULE FACT

    IMPLICIT NONE

    CONTAINS
        ! ---------
        ! FACTORIAL
        ! ---------
        FUNCTION factorial(n) result(fact)
            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: n

            ! INTERMEDIATE VARIABLE
            INTEGER :: i

            ! OUTPUT
            INTEGER :: fact

            fact = 1

            DO i = 2, n
                fact = fact * i
            END DO

        END FUNCTION factorial



        ! ----------
        ! FACTORIAL2
        ! ----------
        FUNCTION factorial2(n) result(fact2)
            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: n

            ! INTERMEDIATE VARIABLE
            INTEGER :: i

            ! OUTPUT
            INTEGER :: fact2

            IF (n == -1) THEN
                fact2 = 1
            ELSE IF (MOD(n,2) == 0) THEN ! EVEN NUMBER
                fact2 = 1

                DO i = 2, n, 2
                    fact2 = fact2 * i
                END DO
            ELSE ! ODD NUMBER
                fact2 = 1

                DO i = 1, n, 2
                    fact2 = fact2 * i
                END DO
            END IF

        END FUNCTION factorial2



        ! --------------------
        ! BINOMIAL COEFFICIENT
        ! --------------------
        FUNCTION binom(n,k)
            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: n, k

            ! OUTPUT
            INTEGER :: binom

            binom = factorial(n)
            binom = binom / factorial(k)
            binom = binom / factorial(n-k)

        END FUNCTION binom

END MODULE FACT
