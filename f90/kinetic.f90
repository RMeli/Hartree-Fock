MODULE KINETIC

    USE CONSTANTS
    USE GAUSSIAN
    USE OVERLAP, only: Si

    IMPLICIT NONE

    CONTAINS

        FUNCTION Ki(ac,a1,a2,bc,b1,b2,aa,bb,Ra,Rb,Ra1,Rb1,Ra2,Rb2,Rc,R1,R2,cp)
            ! ----------------------------------------------------------------------------------------------
            ! Compute kinetic energy between two unnormalized Cartesian Gaussian functions along direction i
            ! ----------------------------------------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: ac, a1, a2, bc, b1, b2       ! Cartesia Gaussian angular momenta
            REAL*8, intent(in) :: aa, bb                        ! Gaussian exponential coefficients
            REAL*8, intent(in) :: Ra, Rb, Ra1, Rb1, Ra2, Rb2    ! Gaussian centers' coordinates
            REAL*8, intent(in) :: Rc, R1, R2                    ! Gaussian product center coordinates
            REAL*8, intent(in) :: cp                            ! Gaussian product coefficient

            ! INTERMEDIATE VARIABLE
            REAL*8 :: k

            ! OUTPUT
            REAL*8 :: Ki

            k = 0.0D0
            k = k + ac * bc * Si(ac - 1,bc - 1,aa,bb,Ra,Rb,Rc)
            k = k - 2 * aa * bc * Si(ac + 1,bc - 1,aa,bb,Ra,Rb,Rc)
            k = k - 2 * ac * bb * Si(ac - 1,bc + 1,aa,bb,Ra,Rb,Rc)
            k = k + 4 * aa * bb * Si(ac + 1,bc + 1,aa,bb,Ra,Rb,Rc)
            k = k * 0.5D0

            Ki = 1
            Ki = Ki * cp * (PI / (aa + bb))**(3.0D0/2.0D0) * k
            Ki = Ki * Si(a1,b1,aa,bb,Ra1,Rb1,R1)
            Ki = Ki * Si(a2,b2,aa,bb,Ra2,Rb2,R2)

        END FUNCTION Ki



        FUNCTION kinetic_coeff(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb) result(K)
            ! -----------------------------------------------------------------
            ! Compute kinetic integral between two Cartesian Gaussian functions
            ! -----------------------------------------------------------------
            !
            ! Source:
            !   The Mathematica Journal
            !   Evaluation of Gaussian Molecular Integrals
            !   II. Kinetic-Energy Integrals
            !   Minhhuy Hô and Julio Manuel Hernández-Pérez
            !
            ! -----------------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: ax, ay, az, bx, by, bz   ! Angular momentum coefficients
            REAL*8, intent(in) :: aa, bb                    ! Exponential Gaussian coefficients
            REAL*8, dimension(3), intent(in) :: Ra, Rb      ! Gaussian centers

            ! INTERMEDIATE VARIABLE
            REAL*8 :: pp                ! Gaussian produc exponential coefficient
            REAL*8, dimension(3) :: Rp  ! Gaussian produc center
            REAL*8 :: cp                ! Gaussian product multiplicative constant

            ! OUTPUT
            REAL*8 :: K

            CALL gaussian_product(aa,bb,Ra,Rb,pp,Rp,cp) ! Compute PP, RP and CP

            K = 0
            K = K + Ki(ax,ay,az,bx,by,bz,aa,bb,Ra(1),Rb(1),Ra(2),Rb(2),Ra(3),Rb(3),Rp(1),Rp(2),Rp(3),cp)    ! Kx
            K = K + Ki(ay,az,ax,by,bz,bx,aa,bb,Ra(2),Rb(2),Ra(3),Rb(3),Ra(1),Rb(1),Rp(2),Rp(3),Rp(1),cp)    ! Ky
            K = K + Ki(az,ax,ay,bz,bx,by,aa,bb,Ra(3),Rb(3),Ra(1),Rb(1),Ra(2),Rb(2),Rp(3),Rp(1),Rp(2),cp)    ! Kz
            K = K * norm(ax,ay,az,aa) * norm(bx,by,bz,bb)

        END FUNCTION kinetic_coeff

        ! --------------
        ! KINETIC MATRIX
        ! --------------
        SUBROUTINE T_kinetic(Kf,basis_D,basis_A,basis_L,basis_R,T)
            ! ----------------------------------------------
            ! Compute kinetic matrix between basis function.
            ! ----------------------------------------------

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

            ! INTERMEDIATE VARIABLES
            INTEGER :: i,j,k,l
            REAL*8 :: tmp

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: T

            T(:,:) = 0.0D0

            DO i = 1,Kf
                DO j = 1,Kf
                    DO k = 1,c
                        DO l = 1,c
                            tmp = basis_D(i,k) * basis_D(j,l)
                            tmp = tmp * kinetic_coeff(  basis_L(i,1),&  ! lx for basis function i
                                                        basis_L(i,2),&  ! ly for basis function i
                                                        basis_L(i,3),&  ! lz for basis function i
                                                        basis_L(j,1),&  ! lx for basis function j
                                                        basis_L(j,2),&  ! ly for basis function j
                                                        basis_L(j,3),&  ! lz for basis function j
                                                        basis_A(i,k),&  ! Exponential coefficient for basis function i, contraction k
                                                        basis_A(j,l),&  ! Exponential coefficient for basis function j, contraction l
                                                        basis_R(i,:),&  ! Center of basis function i
                                                        basis_R(j,:))   ! Center of basis function j

                            T(i,j) = T(i,j) + tmp
                        END DO ! l
                    END DO ! k
                END DO ! j
            END DO ! i

        END SUBROUTINE T_kinetic


END MODULE KINETIC
