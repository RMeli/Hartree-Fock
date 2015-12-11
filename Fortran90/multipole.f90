MODULE MULTIPOLE

    USE CONSTANTS
    USE FACT
    USE GAUSSIAN
    USE OVERLAP

    IMPLICIT NONE

    CONTAINS


    FUNCTION Sije(e,a,b,aa,bb,Rai,Rbi,Rpi,Rci)
        ! --------------------------------------------------
        ! Compute dipole moment of order E along direction I
        ! --------------------------------------------------
        !
        ! Source:
        !   T. Helgaker, P. Jørgensen and J. Olsen
        !   Molecular Electronic-Structure Theory
        !   Wiley
        !   2000
        !
        ! ---------------------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: e                        ! Dipole order
        INTEGER, intent(in) :: a,b                      ! Angular momentum coefficients of the Gaussians along direction i
        REAL*8, intent(in) :: aa, bb                    ! Exponential coefficients of the Gaussians
        REAL*8, intent(in) :: Rai, Rbi                  ! Centers of the Gaussians
        REAL*8, intent(in) :: Rpi                       ! Center of Gaussian product
        REAL*8, intent(in) :: Rci                       ! Center of the multipole

        ! INTERMEDIATE VARIABLE
        REAL*8, dimension(-1:a+1,-1:b+1,-1:e+1) :: S    ! Vector of overlap coefficients
        INTEGER :: i, j, k
        REAL*8 :: factor

        ! OUTPUT
        REAL*8 :: Sije                                  ! Dipole integral

        S(-1,:,:) = 0.0D0
        S(:,-1,:) = 0.0D0
        S(:,:,-1) = 0.0D0

        S(0,0,0) = S00(aa,bb,Rai,Rbi)

        DO i = 0, a
            DO j = 0, b
                DO k = 0, e
                    factor = 1.0D0 / (2 * (aa+bb)) * (i * S(i-1,j,k) + j * S(i,j-1,k) + k * S(i,j,k-1))

                    S(i+1,j,k) = (Rpi - Rai) * S(i,j,k) + factor
                    S(i,j+1,k) = (Rpi - Rbi) * S(i,j,k) + factor
                    S(i,j,k+1) = (Rpi - Rci) * S(i,j,k) + factor
                END DO
            END DO
        END DO

        Sije = S(a,b,e)

    END FUNCTION Sije

    FUNCTION multipole_coeff_OS(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb,e,f,g,Rc) result(S)
        ! -------------------------------------------------------------------
        ! Compute multipole integral between two Cartesian Gaussian functions
        ! -------------------------------------------------------------------
        !
        ! Source:
        !   T. Helgaker, P. Jørgensen and J. Olsen
        !   Molecular Electronic-Structure Theory
        !   Wiley
        !   2000
        !
        ! -------------------------------------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: ax, ay, az, bx, by, bz   ! Angular momentum coefficients
        REAL*8, intent(in) :: aa, bb                    ! Exponential Gaussian coefficients
        REAL*8, dimension(3), intent(in) :: Ra, Rb      ! Gaussian centers
        INTEGER , intent(in) :: e, f, g                 ! Multipole order along x, y and z
        REAL*8, dimension(3), intent(in) :: Rc          ! Multipole center

        ! INTERMEDIATE VARIABLES
        REAL*8 :: pp ! Gaussian produc exponential coefficient
        REAL*8, dimension(3) :: Rp ! Gaussian produc center
        REAL*8 :: cp ! Gaussian product multiplicative constant

        ! OUTPUT
        REAL*8 :: S

        CALL gaussian_product(aa,bb,Ra,Rb,pp,Rp,cp) ! Compute PP, RP and CP

        S = 1
        S = S * Sije(e,ax,bx,aa,bb,Ra(1),Rb(1),Rp(1),Rc(1))       ! Overlap along x
        S = S * Sije(f,ay,by,aa,bb,Ra(2),Rb(2),Rp(2),Rc(2))       ! Overlap along y
        S = S * Sije(g,az,bz,aa,bb,Ra(3),Rb(3),Rp(3),Rc(3))       ! Overlap along z
        S = S * norm(ax,ay,az,aa) * norm(bx,by,bz,bb)             ! Normalization of Gaussian functions

    END FUNCTION multipole_coeff_OS

    SUBROUTINE S_multipole_i(Kf,c,basis_D,basis_A,basis_L,basis_R,dir,e,R,S)
        ! ----------------------------------------------
        ! Compute overlap matrix between basis function.
        ! ----------------------------------------------

        IMPLICIT NONE

        ! TODO Allow flexibility for basis sets other than STO-3G
        INTEGER, intent(in) :: c                        ! Number of contractions per basis function

        ! INPUT
        INTEGER, intent(in) :: Kf                       ! Number of basis functions
        REAL*8, dimension(Kf,3), intent(in) :: basis_R  ! Basis set niclear positions
        INTEGER, dimension(Kf,3), intent(in) :: basis_L ! Basis set angular momenta
        REAL*8, dimension(Kf,3), intent(in) :: basis_D  ! Basis set contraction coefficients
        REAL*8, dimension(Kf,3), intent(in) :: basis_A  ! Basis set exponential contraction coefficients
        REAL*8, dimension(3), intent(in) :: R           ! Multipole center
        INTEGER,  intent(in) :: dir                     ! Multipole direction
        INTEGER,  intent(in) :: e                       ! Multipole exponent along direction DIR

        ! INTERMEDIATE VARIABLES
        INTEGER :: i,j,k,l
        REAL*8 :: tmp
        INTEGER, dimension(3) :: multi_dir              ! Mulrtipole exponent along direction DIR

        ! OUTPUT
        REAL*8, dimension(Kf,Kf), intent(out) :: S

        IF (dir .EQ. 1) THEN

            multi_dir = (/e, 0, 0/)

        ELSE IF (dir .EQ. 2) THEN

            multi_dir = (/0, e, 0/)

        ELSE IF (dir .EQ. 3) THEN

            multi_dir = (/0, 0, e/)

        END IF

        S(:,:) = 0.0D0

        DO i = 1,Kf
            DO j = 1,Kf
                DO k = 1,c
                    DO l = 1,c
                        tmp = basis_D(i,k) * basis_D(j,l)
                        tmp = tmp * multipole_coeff_OS( basis_L(i,1),&  ! lx for basis function i
                                                        basis_L(i,2),&  ! ly for basis function i
                                                        basis_L(i,3),&  ! lz for basis function i
                                                        basis_L(j,1),&  ! lx for basis function j
                                                        basis_L(j,2),&  ! ly for basis function j
                                                        basis_L(j,3),&  ! lz for basis function j
                                                        basis_A(i,k),&  ! Exponential coefficient for basis function i, contraction k
                                                        basis_A(j,l),&  ! Exponential coefficient for basis function j, contraction l
                                                        basis_R(i,:),&  ! Center of basis function i
                                                        basis_R(j,:),&  ! Center of basis function j
                                                        multi_dir(1),&  ! Multipole exponent along x
                                                        multi_dir(2),&  ! Multipole exponent along y
                                                        multi_dir(3),&  ! Multipole exponent along z
                                                        R)              ! Multipole center

                        S(i,j) = S(i,j) + tmp
                    END DO ! l
                END DO ! k
            END DO ! j
        END DO ! i

    END SUBROUTINE S_multipole_i

    SUBROUTINE dipole(Kf,c,Nn,basis_D,basis_A,basis_L,basis_R,P,Rn,Zn,mu)
        ! ----------------------------
        ! Compute multipole vector
        ! ----------------------------
        !
        ! Source:
        !   A. Szabo and N. S. Ostlund
        !   Modern Quantum Chemistry
        !   Dover
        !   1996
        !
        ! ----------------------------

        ! TODO Allow flexibility for basis sets other than STO-3G
        INTEGER, intent(in) :: c                        ! Number of contractions per basis function

        ! INPUT
        INTEGER, intent(in) :: Kf                       ! Number of basis functions
        INTEGER, intent(in) :: Nn                       ! Number of nuclei
        REAL*8, dimension(Kf,3), intent(in) :: basis_R  ! Basis set niclear positions
        INTEGER, dimension(Kf,3), intent(in) :: basis_L ! Basis set angular momenta
        REAL*8, dimension(Kf,3), intent(in) :: basis_D  ! Basis set contraction coefficients
        REAL*8, dimension(Kf,3), intent(in) :: basis_A  ! Basis set exponential contraction coefficients
        REAL*8, dimension(Kf,Kf), intent(in) :: P       ! Density matrix
        REAL*8, dimension(Nn,3), intent(in) :: Rn       ! Nuclear positions
        INTEGER, dimension(Nn), intent(in) :: Zn

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(Kf,Kf) :: S
        REAL*8, dimension(3) :: R
        INTEGER :: i, j, k

        ! OUTPUT
        REAL*8, dimension(3), intent(out) :: mu

        mu(:) = 0.0D0

        R(:) = 0.0D0

        ! -----------
        ! X component
        ! -----------
        CALL S_multipole_i(Kf,c,basis_D,basis_A,basis_L,basis_R,1,1,R,S)

        DO i = 1, Kf
            DO j = 1, Kf
                mu(1) = mu(1) - P(i,j) * S(j,i)
            END DO
        END DO

        ! -----------
        ! Y component
        ! -----------
        CALL S_multipole_i(Kf,c,basis_D,basis_A,basis_L,basis_R,2,1,R,S)

        DO i = 1, Kf
            DO j = 1, Kf
                mu(2) = mu(2) - P(i,j) * S(j,i)
            END DO
        END DO

        ! -----------
        ! Z component
        ! -----------
        CALL S_multipole_i(Kf,c,basis_D,basis_A,basis_L,basis_R,3,1,R,S)

        DO i = 1, Kf
            DO j = 1, Kf
                mu(3) = mu(3) - P(i,j) * S(j,i)
            END DO
        END DO

        ! Take into account nuclear charges
        DO k = 1, Nn
            mu = mu + Zn(k) * Rn(k,:)
        END DO

    END SUBROUTINE


END MODULE MULTIPOLE
