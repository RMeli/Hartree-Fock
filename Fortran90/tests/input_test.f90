PROGRAM input_test

    USE INPUT
    USE UTILS

    IMPLICIT NONE

    ! System an basis set sizes
    INTEGER :: Ne   ! Number of electrons
    INTEGER :: Nn   ! Number of nuclei
    INTEGER :: K    ! Basis set size
    INTEGER :: c    ! Contractions (TODO: split valence basis set)

    ! System and basis set informations
    REAL*8, allocatable, dimension(:,:) :: R            ! Atomic potisions
    INTEGER, allocatable, dimension(:) :: Z             ! Atomic charges
    REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
    INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
    REAL*8, allocatable, dimension(:,:) :: basis_D      ! Conttaction linear coefficients

    CALL load("tests/N2_f.in",Ne,Nn,K,c,R,Z,basis_R,basis_L,basis_A,basis_D)

    CALL print_real_matrix(Nn,3,R)
    CALL print_integer_matrix(Nn,1,Z)
    CALL print_real_matrix(K,3,basis_R)
    CALL print_integer_matrix(K,3,basis_L)
    CALL print_real_matrix(K,c,basis_A)
    CALL print_real_matrix(K,c,basis_D)

    DEALLOCATE(R,Z,basis_R,basis_L,basis_A,basis_D)


END PROGRAM input_test
