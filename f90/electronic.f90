MODULE ELECTRONIC

    USE FACT
    USE NUCLEAR, only: f

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
            t = t / (factorial(r) * factorial(l - 2*r))

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
            B = B * (-1)**(ll) * theta(l,l1,l2,Rp-Ra,Rp-Rb,r,g1)
            B = B * theta(ll,l3,l4,Rq-Rc,Rq-Rd,rr,g2) * (-1)**i
            B = B * (2*delta)**(2*(r+rr)) * factorial(l + ll - 2 * r - 2 * rr)
            B = B * delta**i * (Rp-Rq)**(l+ll - 2*(r+rr+i))
            B = B / ( (4*delta)**(l+ll) * factorial(i) * factorial(l+ll - 2*(r+rr+i)) )

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

            ! INPUT


            ! INTERMEDIATE VARIABLES

            ! OUTPUT

        END FUNCTION electronic_coeff

END MODULE ELECTRONIC
