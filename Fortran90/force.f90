MODULE FORCES

    USE RHF

    IMPLICIT NONE

    CONTAINS

        FUNCTION finite_difference(f1,f2,x1,x2) result(df)

            IMPLICIT NONE

            ! INPUT
            REAL*8, intent(in) :: f1, f2, x1, x2

            ! OUTPUT
            REAL*8 :: df

            df = (f2 - f1) / (x2 - x1)

        END FUNCTION finite_difference

        SUBROUTINE force_fd(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rn,idx,force)
            ! ------------------------------------------------------------------------
            ! Compute forces on atom IDX by finite difference (computationally costly)
            ! ------------------------------------------------------------------------

            ! INPUT
            INTEGER, intent(in) :: Kf                                   ! Basis set size
            INTEGER, intent(in) :: Ne                                   ! Number of electrons
            INTEGER, intent(in) :: Nn                                   ! Number of nuclei
            INTEGER, intent(in) :: c                                    ! Number of contractions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L             ! Angular momenta of basis set Gaussians
            REAL*8, dimension(Kf,3), intent(in) :: basis_R              ! Centers of basis set Gaussians
            REAL*8, dimension(Kf,c), intent(in) :: basis_D, basis_A     ! Basis set coefficients
            REAL*8, dimension(Nn,3), intent(in) :: Rn                   ! Nuclear positions
            INTEGER, dimension(Nn), intent(in) :: Zn                    ! Nuclear charges
            INTEGER, intent(in) :: idx                                  ! Atom index (compute forces for atom idx)

            ! OUTPUT
            REAL*8, dimension(Nn,3), intent(out) :: force               ! Vector of forces (on every atom)

            ! INTERMEDIATE VARIABLES
            REAL*8, PARAMETER :: dd = 1e-9
            REAL*8, dimension(3) :: dx = (/dd,0.0D0,0.0D0/)
            REAL*8, dimension(3) :: dy = (/0.0D0,dd,0.0D0/)
            REAL*8, dimension(3) :: dz = (/0.0D0,0.0D0,dd/)
            REAL*8 :: dEp, dEm
            REAL*8, dimension(Nn,3) :: Rnp, Rnm

            ! -------------
            ! FORCE ALONG X
            ! -------------

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dx

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnp,dEp,.FALSE.)

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dx

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnm,dEm,.FALSE.)

            force(idx,1) = - finite_difference(dEp,dEm,Rnp(idx,1),Rnm(idx,1))

            ! -------------
            ! FORCE ALONG Y
            ! -------------

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dy

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnp,dEp,.FALSE.)

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dy

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnm,dEm,.FALSE.)

            force(idx,2) = - finite_difference(dEp,dEm,Rnp(idx,2),Rnm(idx,2))

            ! -------------
            ! FORCE ALONG Z
            ! -------------

            Rnp = Rn
            Rnp(idx,:) = Rn(idx,:) + dz

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnp,dEp,.FALSE.)

            Rnm = Rn
            Rnm(idx,:) = Rn(idx,:) - dz

            CALL RHF_SCF(Kf,c,Ne,Nn,basis_D,basis_A,basis_L,basis_R,Zn,Rnm,dEm,.FALSE.)

            force(idx,3) = - finite_difference(dEp,dEm,Rnp(idx,3),Rnm(idx,3))

        END SUBROUTINE



END MODULE FORCES
