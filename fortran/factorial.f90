MODULE FACT

    IMPLICIT NONE

CONTAINS

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





    FUNCTION factorial2(n) result(fact2)
        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: n

        ! INTERMEDIATE VARIABLE
        INTEGER :: i

        ! OUTPUT
        INTEGER :: fact2

        IF (n == -1) THEN
            WRITE(*,*) n, "is -1"
            fact2 = 1
        ELSE IF (MOD(n,2) == 0) THEN ! EVEN NUMBER
            WRITE(*,*) n, "is even"
            DO i = 2, n, 2
                fact2 = fact2 * i
            END DO
        ELSE ! ODD NUMBER
            WRITE(*,*) n, "is odd"
            DO i = 1, n, 2
                fact2 = fact2 * i
            END DO
        END IF

    END FUNCTION factorial2

END MODULE FACT




PROGRAM test

    USE FACT, only : factorial, factorial2

    IMPLICIT NONE

    WRITE(*,*) "0!", factorial(0)
    WRITE(*,*) "1!", factorial(1)
    WRITE(*,*) "2!", factorial(2)
    !WRITE(*,*) "3!", factorial(3)
    !WRITE(*,*) "3!", factorial(4)
    !WRITE(*,*) "5!", factorial(5)

    WRITE(*,*) "-1!!", factorial2(-1)
    WRITE(*,*) "0!!", factorial2(0)
    WRITE(*,*) "1!!", factorial2(1)
    WRITE(*,*) "2!!", factorial2(2)
    !WRITE(*,*) "3!!", factorial2(3)
    !WRITE(*,*) "4!!", factorial2(4)
    !WRITE(*,*) "5!!", factorial2(5)

END PROGRAM test
