PROGRAM HF

    USE INPUT
    USE RHF

    IMPLICIT NONE

    CHARACTER (len=20) :: fname ! Input filename
    INTEGER :: nargs            ! Number of input arguments

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

    REAL*8 :: final_E ! Total converged energy

    ! --------------------------
    ! READ COMMAND LINE ARGUMENT
    ! --------------------------

    nargs = IARGC()

    IF (nargs .LT. 1) THEN
        WRITE(*,*) "NO INPUT FILE!"
        CALL EXIT(-1)
    END IF

    CALL GETARG(1, fname)

    WRITE(*,*) fname

    ! ------------------------------------------------
    ! LOAD SYSTEM AND BASIS SET INFORMATIONS FROM FILE
    ! ------------------------------------------------

    CALL load(fname,Ne,Nn,K,c,Rn,Zn,basis_R,basis_L,basis_A,basis_D)

    ! ------------------------
    ! TOTAL ENERGY CALCULATION
    ! ------------------------

    CALL RHF_SCF(K,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,final_E,.TRUE.)

    DEALLOCATE(Rn,Zn,basis_R,basis_L,basis_A,basis_D) ! Deallocate allocated memory

END PROGRAM HF
