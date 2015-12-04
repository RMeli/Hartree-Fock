MODULE DYNAMICS

    USE FORCES, only: force_fd

    IMPLICIT NONE

    CONTAINS

    ! --------------------
    ! POSITION PROPAGATION
    ! --------------------
    SUBROUTINE position_prop(Nn,mass,pos,vel,force_t,dt)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Nn                           ! Number of nuclei (number of moving particles)
        REAL*8, dimension(Nn) :: mass                       ! Nuclear masses
        REAL*8, intent(in) :: dt                            ! Integration time step
        REAL*8, dimension(Nn,3), intent(in) :: force_t      ! Force acting on nuclei (PES gradient) at time t
        REAL*8, dimension(Nn,3), intent(in) :: vel       ! Nuclear velocities at time t

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: pos       ! Nuclear positions

        ! INTERMEDIATE VARIABLES
        INTEGER :: i                                    ! Loop index

        ! Update positions of all atoms
        DO i = 1, Nn
            pos(i,:) = pos(i,:) + vel(i,:) * dt + 0.5D0 * force_t(i,:) / mass(i) * dt**2.0D0
        END DO

    END SUBROUTINE

    ! --------------------
    ! VELOCITY PROPAGATION
    ! --------------------
    SUBROUTINE velocity_prop(Nn,mass,vel,force_t,force_tdt,dt)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Nn                           ! Number of nuclei (number of moving particles)
        REAL*8, dimension(Nn) :: mass                       ! Nuclear masses
        REAL*8, intent(in) :: dt                            ! Integration time step
        REAL*8, dimension(Nn,3), intent(in) :: force_t      ! Force acting on nuclei (PES gradient) at time t
        REAL*8, dimension(Nn,3), intent(in) :: force_tdt    ! Force acting on nuclei (PES gradient) at time t+dt

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: vel       ! Nuclear velocities

        ! INTERMEDIATE VARIABLES
        INTEGER :: i                                        ! Loop index

        ! Update velocities of all atoms
        DO i = 1, Nn
            vel(i,:) = vel(i,:) + 0.5D0 * (force_t(i,:) + force_tdt(i,:)) / mass(i) * dt
        END DO

    END SUBROUTINE


    ! ----------------------------
    ! CLASSICAL MOLECULAR DYNAMICS
    ! ----------------------------
    SUBROUTINE classical_step(Nn,mass,pos,vel,force_t,dU,dt)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Nn                           ! Number of nuclei (number of moving particles)
        REAL*8, dimension(Nn) :: mass                       ! Nuclear masses
        REAL*8, intent(in) :: dt                            ! Integration time step

        ! INPUT/OUTPUT
        REAL*8, dimension(Nn,3), intent(inout) :: pos       ! Nuclear positions
        REAL*8, dimension(Nn,3), intent(inout) :: vel       ! Nuclear velocities
        REAL*8, dimension(Nn,3), intent(inout) :: force_t   ! Force acting on nuclei (PES gradient) at time t

        ! INTERMEDIATE VARIABLES
        REAL*8, dimension(Nn,3) :: force_tdt                ! Force acting on nuclei (PES gradient) at time t+dt

        ! For classical MD we need an analytical description of the PES
        INTERFACE
            FUNCTION dU(Nn,pos)
                INTEGER, intent(in) :: Nn
                REAL*8, dimension(Nn,3), intent(in) :: pos
                REAL*8, dimension(Nn,3) :: dU
            END FUNCTION
        END INTERFACE

        ! Update positions
        CALL position_prop(Nn,mass,pos,vel,force_t,dt) ! Use forces at time t

        ! Compute forces at time t+dt
        force_tdt = (-1) * dU(Nn,pos) ! POS has been updated by POSITION_PROP

        ! Update velocities
        CALL velocity_prop(Nn,mass,vel,force_t,force_tdt,dt)

        force_t = force_tdt ! Store the force at t+dt for the next step

    END SUBROUTINE classical_step

    ! -----------------------------------
    ! BORN-OPPENHEIMER MOLECULAR DYNAMICS
    ! -----------------------------------
    SUBROUTINE BO_step()

        IMPLICIT NONE

    END SUBROUTINE BO_step

END MODULE DYNAMICS
