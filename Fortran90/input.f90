MODULE INPUT

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE load(fname,Ne,Nn,K,c,R,basis_R,basis_L,basis_A,basis_D)

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(out) :: Ne              ! Number of electrons
            INTEGER, intent(out) :: Nn              ! Number of nuclei
            INTEGER, intent(out) :: K               ! Basis set size
            INTEGER, intent(out) :: c               ! Contractions (TODO: split valence basis set)
            CHARACTER(len=*), intent(in) :: fname   ! Input file name

            ! INTERMEDIATE VARIABLES
            INTEGER :: i ! Loop index

            ! ALLOCATABLE
            REAL*8, allocatable, dimension(:,:) :: R            ! Atomic potisions
            REAL*8, allocatable, dimension(:,:) :: basis_R      ! Basis functions' centers
            INTEGER, allocatable, dimension(:,:) :: basis_L     ! Basis functions' angular momenta
            REAL*8, allocatable, dimension(:,:) :: basis_A      ! Contraction exponential coefficients
            REAL*8, allocatable, dimension(:,:) :: basis_D      ! Conttaction linear coefficients

            ! Open file containing system and basis set informations
            OPEN(unit=100,file=fname,form="formatted",status="old",action="read")

            READ(100,*) Ne
            READ(100,*) Nn
            READ(100,*) K
            READ(100,*) C

            ! Allocate arrays once we know system's and basis set's specifications (Ne,Nn,K,c)
            ALLOCATE(R(Nn,3))
            ALLOCATE(basis_R(K,3))
            ALLOCATE(basis_L(K,3))
            ALLOCATE(basis_A(K,c))
            ALLOCATE(basis_D(K,c))

            ! Read atomic positions
            DO i = 1, Nn
                READ(100,*) R(i,:)
            END DO

            ! Read basis set informations
            DO i = 1, K
                READ(100,*) basis_R(i,:), basis_L(i,:), basis_A(i,:), basis_D(i,:)
            END DO

            CLOSE(unit=100) ! Close the file

        END SUBROUTINE load

END MODULE INPUT
