MODULE FORCE

    IMPLICIT NONE

    ! PARAMETERS
    INTEGER, PARAMETER :: c = 3 ! Number of contractions (STO-3G)

    FUNCTION finite_difference(f1,f2,x1,x2) result(df)

        IMPLICIT NONE

        ! INPUT
        REAL*8, intent(in) :: f1, f2, x1, x2

        ! OUTPUT
        REAL*8 :: df

        df = (f2 - f1) / (x2 - x1)

    END FUNCTION finite_difference

    SUBROUTINE force_fd(Kf,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,Etot,idx,force)
        ! -------------------------------------------------------------------------
        ! Compute forces on atom IDX by finite difference (computationally costrly)
        ! -------------------------------------------------------------------------

        ! INPUT
        INTEGER, intent(in) :: Kf                                   ! Basis set size
        INTEGER, intent(in) :: Ne                                   ! Number of electrons
        INTEGER, intent(in) :: Nn                                   ! Number of nuclei
        INTEGER, dimension(K,3), intent(in) :: basis_L              ! Angular momenta of basis set Gaussians
        REAL*8, dimension(K,3), intent(in) :: basis_R               ! Centers of basis set Gaussians
        REAL*8, dimension(K,c), intent(in) :: basis_D, basis_A      ! Basis set coefficients
        REAL*8, dimension(Nn,3), intent(in) :: Rn                   ! Nuclear positions
        INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
        REAL*8, intent(in) :: Etot                                  ! PES value of current geometry
        REAL*8, intent(in) :: idx                                   ! Atom index (compute forces for atom idx)

        ! OUTPUT
        REAL*8, dimension(3), intent(out) :: force                  ! Force acting on atom IDX

        ! INTERMEDIATE VARIABLES
        REAL*8 :: dx, dy, dz

        ! --------------
        ! SCF PARAMETERS
        ! --------------

        REAL*8, dimension(Kf,Kf) :: S     ! Overlap matrix
        REAL*8, dimension(Kf,Kf) :: X     ! Transformation matrix

        REAL*8, dimension(Kf,Kf) :: Hc    ! Core Hamiltonian

        REAL*8, dimension(Kf,Kf) :: Pold  ! Old density matrix
        REAL*8, dimension(Kf,Kf) :: Pnew  ! New density matrix

        REAl*8, dimension(Kf,Kf,Kf,Kf) :: ee ! List of electron-electron integrals

        REAL*8, dimension(Kf,Kf) :: F     ! Fock matrix
        REAL*8, dimension(Kf) :: E           ! Orbital energies

        LOGICAL :: converged = .FALSE.      ! Convergence parameter
        INTEGER, PARAMETER :: maxiter = 100 ! Maximal number of iterations
        INTEGER :: step = 0                 ! SCF steps counter

    END SUBROUTINE

END MODULE FORCE
