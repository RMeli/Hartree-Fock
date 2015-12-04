PROGRAM HF_N2

    USE INPUT
    USE RHF
    USE FORCES

    IMPLICIT NONE

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO)

    ! System and basis set informations
    REAL*8, allocatable, dimension(:,:) :: Rn           ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Zn            ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
    INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D      ! Conttaction linear coefficients

    REAL*8, allocatable, dimension(:,:) :: F

    ! ------------------------------------------------
    ! LOAD SYSTEM AND BASIS SET INFORMATIONS FROM FILE
    ! ------------------------------------------------

    CALL load("tests/H2_f.in",Ne,Nn,K,c,Rn,Zn,basis_R,basis_L,basis_A,basis_D)

    ! -----
    ! FORCE
    ! -----

    ALLOCATE(F(Nn,3))

    ! Force on atom 1
    CALL force_fd(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,F,1D-4)

    WRITE(*,*) "######"
    WRITE(*,*) "FORCES"
    WRITE(*,*) "######"
    CALL print_real_matrix(Nn,3,F)

    ! ---
    ! END
    ! ---

    DEALLOCATE(Rn,Zn,basis_R,basis_L,basis_A,basis_D,F) ! Deallocate allocated memory

END PROGRAM HF_N2
