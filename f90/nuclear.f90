MODULE NUCLEAR

    USE FACT
    USE CONSTANTS
    USE GAUSSIAN

    IMPLICIT NONE

    REAL*8 :: xx
    INTEGER :: nunu

    CONTAINS
        FUNCTION f(j,l,m,a,b)
            !----------------------------------------
            ! Expansion coefficient.
            !----------------------------------------
            !
            ! Source:
            !     Handbook of Computational Chemistry
            !     David Cook
            !     Oxford University Press
            !     1998
            !
            !----------------------------------------

            ! INPUT
            INTEGER, intent(in) :: j, l, m
            REAL*8, intent(in) :: a, b

            ! INTERMEDIATE VARIABLES
            INTEGER :: k    ! Loop index
            REAl*8 :: tmp

            ! OUTPUT
            REAL*8 :: f

            f = 0.0D0

            DO k = MAX(0,j-m), MIN(j,l)
                tmp = 1.0D0
                tmp = tmp * binom(l,k) * binom(m,j-k)
                tmp = tmp * a**(l-k) * b**(m+k-j)

                f = f + tmp
            END DO

        END FUNCTION f


        FUNCTION boys0(x)
            ! ----------------------------
            ! Boys function for s orbitals
            ! ----------------------------
            !
            ! Source:
            !   Szabo and Ostlund
            !   Modern Quantum Chemistry
            !   Doever
            !   1989
            !
            ! ----------------------------

            IMPLICIT NONE

            !INPUT
            REAL*8,intent(in) :: x

            ! OUTPUT
            REAL*8 :: boys0

            boys0 = 0.0D0

            IF (x .LT. 1.0D-12) THEN ! Asymptotic value
                boys0 = 1.0D0 - x / 3.0D0
            ELSE
                boys0 = DSQRT(PI / x) * DERF(DSQRT(x)) / 2.0D0
            END IF

        END FUNCTION boys0

        FUNCTION INT(tt) result(res)

            IMPLICIT NONE

            ! INPUT
            REAL, intent(in) :: tt

            ! OUTPUT
            REAL :: res

            res = tt**(2*nunu) * EXP(-xx*tt**2)

        END FUNCTION INT

        ! -------------
        ! BOYS FUNCTION
        ! -------------
        FUNCTION boys(nu,x)
            ! -------------------------------------------------------------------------
            ! Boys function
            ! -------------------------------------------------------------------------
            !
            ! Sources:
            !
            !   Weiss and Ochsenfeld
            !   A Rigorous and Optimized Strategy for the Evaluation
            !       of the Boys Function Kernel in Molecular Electronic StructureTheory
            !   Journal of Computational Chemistry
            !   2015
            !
            !   Handbook of Computational Chemistry
            !   David Cook
            !   Oxford University Press
            !   1998
            !
            ! -------------------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: nu
            REAL*8, intent(in) :: x

            ! INTERMEDIATE VARIABLES
            !INTEGER :: i
            !REAL*8 :: exp, fact
            !REAL*8 :: b, bnew, sum,  resid
            REAL :: result, abserr
            INTEGER :: neval, ier

            ! PARAMETERS
            !INTEGER, PARAMETER :: trans = 10        ! Transition between small X and large X
            !REAL*8, PARAMETER :: tol = 1.0D-12      ! Boys series convergence

            ! OUTPUT
            REAL*8 :: boys

            nunu = nu   ! Set global MODULE variable to input variable (side effect: change in INT)
            xx = x      ! Set global MODULE variable to input variable (side effect: change in INT)

            ! TODO DOUBLE precision routine
            CALL qags(INT,0.0,1.0,1.0E-5,1.0E-5, result, abserr, neval, ier) ! Call QUADPACK routine

            IF (ier .NE. 0) THEN
                WRITE(*,*) ier
            END IF

            boys = result

            !boys = 0.0D0

            !IF (nu == 0) THEN
            !    boys = boys0(x) ! Call Boys function from Szabo and Ostlund
            !IF (nu .LT. trans) THEN
            !    resid = 1 ! Residual

            !    exp = 0.5D0 * DEXP(-x) ! Multiplicative exponential factor

            !    sum = DGAMMA(nu + 0.5D0) / DGAMMA(nu + 1.5D0) ! First term of the sum

            !    b = exp * sum

            !    i = 1 ! Start from second step
            !    DO WHILE (resid .GT. tol) ! TODO Add maximal numbe of iterations
            !        sum = sum + DGAMMA(nu + 0.5D0) * x**i / DGAMMA(nu + i + 1.5D0)

            !        bnew = exp * sum
            !        resid = ABS(bnew-b)

            !        b = bnew

            !        i = i + 1
            !    END DO

            !    boys = b

            !ELSE
            !    fact = DGAMMA(nu + 0.5D0) / (2.0D0 * x**(nu + 0.5D0))

            !    exp = 0.5D0 * DEXP(-x) / x ! Multiplicative exponential factor

            !    sum = DGAMMA(nu + 0.5D0) / DGAMMA(nu + 1.5D0) ! First term of the sum

            !    b = fact - exp * sum

            !    i = 1 ! Start from second step
            !    DO WHILE (resid .GT. tol) ! TODO Add maximal numbe of iterations
            !        sum = sum + DGAMMA(nu + 0.5D0) * x**(-i) / DGAMMA(nu - i + 1.5D0)

            !        bnew = fact - exp * sum

            !        resid = ABS(bnew-b)

            !        b = bnew

            !        i = i + 1
            !    END DO

            !    boys = b

            !END IF

        END FUNCTION boys

        FUNCTION A(l,r,i,l1,l2,Ra,Rb,Rc,Rp,eps)
            ! -------------------------------------
            ! Factor A
            ! -------------------------------------
            !
            ! Source:
            !   Handbook of Computational Chemistry
            !   David Cook
            !   Oxford University Press
            !   1998
            !
            !--------------------------------------

            ! INPUT
            INTEGER, intent(in) :: l, r, i , l1, l2
            REAL*8, intent(in) :: Ra, Rb, Rc, Rp, eps

            ! OUTPUT
            REAL*8 :: A

            A = 1.0D0
            A = A * (-1)**(l)
            A = A * f(l,l1,l2,Rp-Ra,Rp-Rb)
            A = A * (-1)**i
            A = A * factorial(l)
            A = A * (Rp-Rc)**(l - 2*r - 2*i)
            A = A * eps**(r+i)
            A = A / (factorial(r)*factorial(i)*factorial(l - 2*r - 2*i))

        END FUNCTION A


        FUNCTION nuclear_coeff(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb,Rn,Zn) result(Vn)
            ! -------------------------------------
            ! Compute
            ! -------------------------------------
            !
            ! Source:
            !   Handbook of Computational Chemistry
            !   David Cook
            !   Oxford University Press
            !   1998
            !
            !--------------------------------------

            ! INPUT
            INTEGER, intent(in) :: ax, ay, az, bx, by, bz   ! Angular momentum coefficients
            REAL*8, intent(in) :: aa, bb                    ! Exponential Gaussian coefficients
            REAL*8, dimension(3), intent(in) :: Ra, Rb      ! Gaussian centers
            REAL*8, dimension(3), intent(in) :: Rn          ! Nuclear position
            INTEGER, intent(in) :: Zn                       ! Nuclear charge

            ! INTERMEDIATE VARIABLES
            REAL*8 :: eps
            REAL*8 :: g                                     ! Gaussian produc exponential coefficient
            REAL*8, dimension(3) :: Rp                      ! Gaussian produc center
            REAL*8 :: cp                                    ! Gaussian product multiplicative constant
            INTEGER :: l, r, i, m, s, j, n, t, k            ! Loop indices
            REAL*8 :: AAx, AAy, AAz                            ! Temporary calls to FUNCTION A
            INTEGER :: nu                                   ! Boys function index

            ! OUTPUT
            REAL*8 :: Vn                                    ! Nuclear matrix element


            CALL gaussian_product(aa,bb,Ra,Rb,g,Rp,cp)

            eps = 1.0D0 / (4.0D0 * g)

            Vn = 0.0D0

            DO l = 0, ax+bx
                DO r = 0, FLOOR(l / 2.0)
                    DO i = 0, FLOOR((l - 2*r)/2.0)
                        AAx = A(l,r,i,ax,bx,Ra(1),Rb(1),Rn(1),Rp(1),eps)

                        DO m = 0, ay+by
                            DO s = 0, FLOOR(m / 2.0)
                                DO j = 0, FLOOR((m - 2*s)/2.)
                                    AAy = A(m,s,j,ay,by,Ra(2),Rb(2),Rn(2),Rp(2),eps)

                                    DO n = 0, az+bz
                                        DO t = 0, FLOOR(n / 2.0)
                                            DO k = 0, FLOOR((n - 2*t)/2.0)
                                                AAz =  A(n,t,k,az,bz,Ra(3),Rb(3),Rn(3),Rp(3),eps)

                                                nu = l + m + n - 2 * (r + s + t) - (i + j + k)

                                                Vn = Vn + AAx * AAy * AAz * boys(nu,g*DOT_PRODUCT(Rp-Rn,Rp-Rn))

                                            END DO ! k
                                        END DO ! t
                                    END DO ! n
                                END DO ! j
                            END DO ! s
                        END DO ! m
                    END DO ! i
                END DO ! r
            END DO ! l

            Vn = Vn * (-Zn) * norm(ax,ay,az,aa) * norm(bx,by,bz,bb) * cp * 2.0D0 * PI / g

        END FUNCTION nuclear_coeff


        ! -------------
        ! NUCLEUS-ELECRON POTENTIAL MATRIX
        ! -------------
        SUBROUTINE V_nuclear(Kf,basis_D,basis_A,basis_L,basis_R,Vn,Rn,Zn)
            ! ------------------------------------------
            ! Compute nucleus-electrona potential matrix
            ! ------------------------------------------

            IMPLICIT NONE

            ! TODO Allow flexibility for basis sets other than STO-3G
            ! HARD CODED
            INTEGER, PARAMETER :: c = 3 ! Number of contractions per basis function

            ! INPUT
            INTEGER, intent(in) :: Kf                       ! Number of basis functions
            REAL*8, dimension(Kf,3), intent(in) :: basis_R  ! Basis set niclear positions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L ! Basis set angular momenta
            REAL*8, dimension(Kf,3), intent(in) :: basis_D  ! Basis set contraction coefficients
            REAL*8, dimension(Kf,3), intent(in) :: basis_A  ! Basis set exponential contraction coefficients
            REAL*8, dimension(3), intent(in) :: Rn          ! Nuclear position
            INTEGER, intent(in) :: Zn                       ! Nuclear charge (> 0)

            ! INTERMEDIATE VARIABLES
            INTEGER :: i,j,k,l
            REAL*8 :: tmp

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: Vn     ! Nucleus-electron potential matrix

            Vn(:,:) = 0.0D0

            DO i = 1,Kf
                DO j = 1,Kf
                    DO k = 1,c
                        DO l = 1,c
                            tmp = basis_D(i,k) * basis_D(j,l)
                            tmp = tmp * nuclear_coeff(  basis_L(i,1),&  ! lx for basis function i
                                                        basis_L(i,2),&  ! ly for basis function i
                                                        basis_L(i,3),&  ! lz for basis function i
                                                        basis_L(j,1),&  ! lx for basis function j
                                                        basis_L(j,2),&  ! ly for basis function j
                                                        basis_L(j,3),&  ! lz for basis function j
                                                        basis_A(i,k),&  ! Exponential coefficient for basis function i, contraction k
                                                        basis_A(j,l),&  ! Exponential coefficient for basis function j, contraction l
                                                        basis_R(i,:),&  ! Center of basis function i
                                                        basis_R(j,:),&  ! Center of basis function j
                                                        Rn,Zn)          ! Nuclear position and nuclear charge

                            Vn(i,j) = Vn(i,j) + tmp
                        END DO ! l
                    END DO ! k
                END DO ! j
            END DO ! i

        END SUBROUTINE V_nuclear


END MODULE NUCLEAR
