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
                fact2 = -1
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

END MODULE FACT
