MODULE RHF

    USE UTILS, only: EIGS, print_real_matrix
    USE DENSITY, only: P_density
    USE ELECTRONIC, only: G_ee

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE RHF_step(Kf,Ne,H,X,ee,Pold,Pnew,F,orbitalE,verbose)

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                           ! Number of basis functions
            !INTEGER, intent(in) :: Nn                           ! Number of nuclei
            INTEGER, intent(in) :: Ne                           ! Number of electrons
            !REAL*8, dimension(Kf,3), intent(in) :: basis_D      ! Contraction coefficients
            !REAL*8, dimension(Kf,3), intent(in) :: basis_A      ! Contraction exponential coefficients
            !INTEGER, dimension(Kf,3), intent(in) :: basis_L     ! Basis function angular momenta
            !REAL*8, dimension(Kf,3), intent(in) :: basis_R      ! Basis set centers
            !REAL*8, dimension(Nn,3), intent(in) :: Rn           ! Nuclear coordinates
            !INTEGER, dimension(Nn), intent(in) :: Zn            ! Nuclear charges
            REAL*8, dimension(Kf,Kf), intent(in) :: H           ! Core Hamiltonian
            REAL*8, dimension(Kf,Kf), intent(in) :: X           ! Transformation maxrix
            REAL*8, dimension(Kf,Kf), intent(in) :: Pold        ! Old density matrix
            REAL*8, dimension(Kf,Kf,Kf,Kf), intent(in) :: ee    ! Electron-electron list
            LOGICAL, intent(in) :: verbose                      ! Flag to print matrices

            ! INTERMEDIATE VARIABLES
            REAL*8, dimension(Kf,Kf) :: G                       ! Electron-electron repulsion matrix
            REAL*8, dimension(Kf,Kf) :: Fx                      ! Fock matrix in the orthogonal basis set
            REAL*8, dimension(Kf,Kf) :: Cx                      ! Coefficient matrix in the orthogonal basis set
            REAL*8, dimension(Kf,Kf) :: C                       ! Coefficient matrix in the original basis set

            ! INPUT / OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: Pnew       ! New density matrix
            REAL*8, dimension(Kf,Kf), intent(out) :: F          ! Fock operator
            REAL*8, dimension(Kf), intent(out) :: orbitalE      ! Orbital energies

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Density matrix P:"
                CALL print_real_matrix(Kf,Kf,Pold)
            END IF

            CALL G_ee(Kf,ee,Pold,G) ! Compute new electron-electron repulsion matrix G

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Electron-electron repulsion matrix G:"
                CALL print_real_matrix(Kf,Kf,G)
            END IF

            F = H + G ! Compute new Fock operator

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix F:"
                CALL print_real_matrix(Kf,Kf,F)
            END IF

            Fx = MATMUL(TRANSPOSE(X),MATMUL(F,X))

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Fock matrix in orthogonal orbital basis Fx:"
                CALL print_real_matrix(Kf,Kf,Fx)
            END IF

            ! ------------------------------------
            ! Solve orthogonalized Roothan eqution
            ! ------------------------------------
            CALL EIGS(Kf,Fx,Cx,orbitalE) ! Compute coefficients (in orthonormal basis) and orbital energies

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients in orthogonal orbital basis Cx:"
                CALL print_real_matrix(Kf,Kf,Cx)
            END IF

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Orbital energies:"
                CALL print_real_matrix(Kf,1,orbitalE)
            END IF

            C = MATMUL(X,Cx) ! Compute coefficients in the original basis

            IF (verbose) THEN
                WRITE(*,*)
                WRITE(*,*) "Coefficients:"
                CALL print_real_matrix(Kf,Kf,C)
            END IF

            CALL P_density(Kf,Ne,C,Pnew) ! Compute new density

        END SUBROUTINE RHF_step


END MODULE RHF
