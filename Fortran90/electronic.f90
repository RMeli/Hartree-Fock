MODULE ELECTRONIC

    USE FACT
    USE GAUSSIAN
    USE CONSTANTS
    USE NUCLEAR, only: f, boys

    IMPLICIT NONE

    CONTAINS

        FUNCTION theta(l,l1,l2,a,b,r,g) result(t)
            !----------------------------------------
            ! Coefficient theta
            !----------------------------------------
            !
            ! Source:
            !     Handbook of Computational Chemistry
            !     David Cook
            !     Oxford University Press
            !     1998
            !
            !----------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: l, l1, l2, r
            REAL*8, intent(in) :: a, b, g

            ! OUTPUT
            REAL*8 :: t

            t = 1.0D0
            t = t * f(l,l1,l2,a,b) * factorial(l) * g**(r-l)
            t = t / ( factorial(r) * factorial(l - 2*r) )

        END FUNCTION theta


        FUNCTION B(l,ll,r,rr,i,l1,l2,Ra,Rb,Rp,g1,l3,l4,Rc,Rd,Rq,g2,delta)
            !----------------------------------------
            ! Coefficient B
            !----------------------------------------
            !
            ! Source:
            !     Handbook of Computational Chemistry
            !     David Cook
            !     Oxford University Press
            !     1998
            !
            !----------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: l, ll, r, rr, i, l1, l2, l3, l4
            REAL*8, intent(in) :: Ra, Rb, Rp, g1, Rc, Rd, Rq, g2, delta

            ! OUTPUT
            REAL*8 :: B

            B = 1.0D0
            B = B * (-1)**(l) * theta(l,l1,l2,Rp-Ra,Rp-Rb,r,g1)
            B = B * theta(ll,l3,l4,Rq-Rc,Rq-Rd,rr,g2) * (-1)**i
            B = B * (2.0D0*delta)**(2*(r+rr)) * factorial(l + ll - 2 * r - 2 * rr)
            B = B * delta**i * (Rp-Rq)**(l+ll - 2*(r+rr+i))
            B = B / ( (4.0D0*delta)**(l+ll) * factorial(i) * factorial(l+ll - 2*(r+rr+i)) )

        END FUNCTION B


        FUNCTION electronic_coeff(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,aa,bb,cc,dd,Ra,Rb,Rc,Rd) result(G)
            !--------------------------------------------------------------------------
            ! Compute electron-electron integral between two cartesia Gaussian function
            !--------------------------------------------------------------------------
            !
            ! Source:
            !     Handbook of Computational Chemistry
            !     David Cook
            !     Oxford University Press
            !     1998
            !
            !--------------------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz   ! Angular momentum coefficients
            REAL*8, intent(in) :: aa, bb, cc, dd                                    ! Exponential Gaussian coefficients
            REAL*8, dimension(3), intent(in) :: Ra, Rb, Rc, Rd                      ! Gaussian centers

            ! INTERMEDIATE VARIABLES
            REAL*8 :: delta
            REAL*8 :: g1, g2                        ! Gaussian produc exponential coefficient
            REAL*8, dimension(3) :: Rp, Rq          ! Gaussian produc center
            REAL*8 :: c1, c2                        ! Gaussian product multiplicative constant
            REAL*8 :: BBx, BBy, BBz                 ! Temporary calls to FUNCTION B
            INTEGER :: nu                           ! Boys function index
            INTEGER :: l, r, i, ll, rr              ! Indices loop 1
            INTEGER :: m, s, j, mm, ss              ! Indices loop 2
            INTEGER :: n, t, k, nn, tt              ! Indices loop 3

            ! OUTPUT
            REAL*8 :: G ! Electron-electron repulsion coefficient

            G = 0.0D0

            CALL gaussian_product(aa,bb,Ra,Rb,g1,Rp,c1) ! Gaussian product of left gaussians
            CALL gaussian_product(cc,dd,Rc,Rd,g2,Rq,c2) ! Gaussian product of right gaussians

            delta = 1.0D0 / (4.0D0 * g1) + 1.0D0 / (4.0D0 * g2)

            DO l = 0, ax + bx
                DO r = 0, FLOOR(l / 2.0)
                    DO ll = 0, cx + dx
                        DO rr = 0, FLOOR(ll / 2.0)
                            DO i = 0, FLOOR((l+ll-2*r-2*rr) / 2.0)

                                 BBx = B(l,ll,r,rr,i,ax,bx,Ra(1),Rb(1),Rp(1),g1,&
                                         cx,dx,Rc(1),Rd(1),Rq(1),g2,delta)

                                 DO m = 0, ay + by
                                     DO s = 0, FLOOR(m / 2.0)
                                         DO mm = 0, cy + dy
                                             DO ss = 0, FLOOR(mm / 2.0)
                                                 DO j = 0, FLOOR((m+mm-2*s-2*ss) / 2.0)
                                                    BBy = B(m,mm,s,ss,j,ay,by,Ra(2),Rb(2),Rp(2),&
                                                            g1,cy,dy,Rc(2),Rd(2),Rq(2),g2,delta)

                                                    DO n = 0, az + bz
                                                        DO t = 0, FLOOR(n / 2.0)
                                                            DO nn = 0, cz + dz
                                                                DO tt = 0, FLOOR(nn / 2.0)
                                                                    DO k = 0, FLOOR((n+nn-2*t-2*tt) / 2.0)

                                                                        BBz = B(n,nn,t,tt,k,az,bz,Ra(3),Rb(3),Rp(3),g1,&
                                                                                cz,dz,Rc(3),Rd(3),Rq(3),g2,delta)

                                                                        nu = l+ll+m+mm+n+nn - 2*(r+rr+s+ss+t+tt) - (i + j + k)

                                                                        G = G + BBx * BBy * BBz *&
                                                                         boys(nu,DOT_PRODUCT(Rp-Rq,Rp-Rq)*g1*g2/(g1+g2))

                                                                    END DO ! tt
                                                                END DO ! nn
                                                            END DO ! k
                                                        END DO ! t
                                                    END DO ! n

                                                END DO ! ss
                                            END DO ! mm
                                        END DO ! j
                                    END DO ! s
                                END DO ! m

                            END DO ! rr
                        END DO ! ll
                    END DO ! i
                END DO ! r
            END DO ! k

            G = G * norm(ax,ay,az,aa) * norm(bx,by,bz,bb) * norm(cx,cy,cz,cc) * norm(dx,dy,dz,dd)
            G = G * c1 * c2 * 2.0D0 * PI**2 / (g1 * g2) * SQRT(PI / (g1 + g2))

        END FUNCTION electronic_coeff


        ! ---------------------------------------
        ! GENERATE LIST OF TWO-ELECTRON INTEGRALS
        ! ---------------------------------------
        SUBROUTINE ee_list(Kf,basis_D,basis_A,basis_L,basis_R,ee)

            IMPLICIT NONE

            ! TODO Allow flexibility for basis sets other than STO-3G
            ! HARD CODED
            INTEGER, PARAMETER :: c = 3 ! Number of contractions per basis function

            ! INPUT
            INTEGER, intent(in) :: Kf                       ! Number of basis functions
            REAL*8, dimension(Kf,3), intent(in) :: basis_R  ! Basis set niclear positions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L ! Basis set angular momenta
            REAL*8, dimension(Kf,c), intent(in) :: basis_D  ! Basis set contraction coefficients
            REAL*8, dimension(Kf,c), intent(in) :: basis_A  ! Basis set exponential contraction coefficients

            ! INTERMEDIATE VARIABLES
            INTEGER :: i,j,k,l  ! Indices running over basis functions
            INTEGER :: m,n,o,p  ! Indices running over contractions
            REAL*8 :: tmp

            ! OUTPUT
            REAL*8, dimension(Kf,Kf,Kf,Kf), intent(out) :: ee   ! List of electron-electron integrals

            ee(:,:,:,:) = 0.0D0

            DO i = 1,Kf
                DO j = 1,Kf
                    DO k = 1,Kf
                        DO l = 1,Kf

                            DO m = 1,c
                                DO n = 1,c
                                    DO o = 1,c
                                        DO p = 1,c

                                            ! Contraction coefficients
                                            tmp = basis_D(i,m) * basis_D(j,n) * basis_D(k,o) * basis_D(l,p)

                                            tmp = tmp * electronic_coeff(   basis_L(i,1),&  ! lx for basis function i
                                                                            basis_L(i,2),&  ! ly for basis function i
                                                                            basis_L(i,3),&  ! lz for basis function i
                                                                            basis_L(j,1),&  ! lx for basis function j
                                                                            basis_L(j,2),&  ! ly for basis function j
                                                                            basis_L(j,3),&  ! lz for basis function j
                                                                            basis_L(k,1),&  ! lx for basis function k
                                                                            basis_L(k,2),&  ! ly for basis function k
                                                                            basis_L(k,3),&  ! lz for basis function k
                                                                            basis_L(l,1),&  ! lx for basis function l
                                                                            basis_L(l,2),&  ! ly for basis function l
                                                                            basis_L(l,3),&  ! lz for basis function l
                                                                            basis_A(i,m),&  ! Exponential coefficient for basis function i, contraction m
                                                                            basis_A(j,n),&  ! Exponential coefficient for basis function j, contraction n
                                                                            basis_A(k,o),&  ! Exponential coefficient for basis function k, contraction o
                                                                            basis_A(l,p),&  ! Exponential coefficient for basis function l, contraction p
                                                                            basis_R(i,:),&  ! Center of basis function i
                                                                            basis_R(j,:),&  ! Center of basis function j
                                                                            basis_R(k,:),&  ! Center of basis function i
                                                                            basis_R(l,:))   ! Center of basis function j


                                            ee(i,j,k,l) = ee(i,j,k,l) + tmp

                                        END DO ! p
                                    END DO ! o
                                END DO ! m
                            END DO ! n

                        END DO ! l
                    END DO ! k
                END DO ! j
            END DO ! i

        END SUBROUTINE ee_list


        ! ----------------------------------
        ! ELECTRON-ELECTRON REPULSION MATRIX
        ! ----------------------------------
        SUBROUTINE G_ee(Kf,ee,P,G)
            ! -----------------------------------------------------------------------
            ! Compute the electron-electron repulsion matrix using the density matrix
            ! -----------------------------------------------------------------------
            !
            ! Source:
            !   Szabo and Ostlund
            !   Modern Quantum Chemistry
            !   Doever
            !   1989
            !
            ! ------------------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                           ! Number of basis functions
            REAL*8, dimension(Kf,Kf,Kf,Kf), intent(in) :: ee    ! Electron-electron list
            REAL*8, dimension(Kf,Kf) :: P                       ! Density matrix

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j, k, l

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: G          ! Electron-electron repulsion matrix

            G(:,:) = 0.0D0

            DO i = 1, Kf
                DO j = 1, Kf
                    DO k = 1, Kf
                        DO l = 1, Kf
                            G(i,j) = G(i,j) + P(k,l) * (ee(i,j,k,l) - 0.5D0 * ee(i,l,k,j))
                        END DO ! l
                    END DO ! k
                END DO ! j
            END DO ! i

        END SUBROUTINE G_ee

END MODULE ELECTRONIC
