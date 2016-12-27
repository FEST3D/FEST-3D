module state
    !-------------------------------------------------------------------
    ! The state module contains the state variables and the methods that
    ! act on them. 
    !
    ! The state of the system is defined using the density, velocity and
    ! pressure (primitive variables qp) at the grid points. This current
    ! version assumes the grid is in atmost two dimensions. 
    !-------------------------------------------------------------------
    
    use global, only: FILE_NAME_LENGTH, STATE_FILE_UNIT, OUT_FILE_UNIT, &
            DESCRIPTION_STRING_LENGTH, STRING_BUFFER_LENGTH
    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
    use global_vars, only : grid_x
    use global_vars, only : grid_y
    use global_vars, only : grid_z

    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z

    use global_vars, only : n_var
    use global_vars, only : qp
    use global_vars, only : qp_inf
    use global_vars, only : density
    use global_vars, only : x_speed
    use global_vars, only : y_speed
    use global_vars, only : z_speed
    use global_vars, only : pressure
    use global_vars, only : density_inf
    use global_vars, only : x_speed_inf
    use global_vars, only : y_speed_inf
    use global_vars, only : z_speed_inf
    use global_vars, only : pressure_inf
    use global_vars, only : tk
    use global_vars, only : tw
    use global_vars, only : tk_inf
    use global_vars, only : tw_inf
    use global_vars, only : gm
    use global_vars, only : mu_ref
    use global_vars, only : supersonic_flag
    use global_vars, only : turbulence
    
    use utils, only: alloc, dealloc, dmsg
    use layout, only: process_id
    use string
!    use grid, only: imx, jmx, kmx, grid_x, grid_y, grid_z
!    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz

    implicit none
    private

!    ! State variables
!    real, public, dimension(:, :, :, :), allocatable, target :: qp
!    ! Infinity variables (free stream conditions)
!    real, public, dimension(:), allocatable, target :: qp_inf
!    ! State variable component aliases
!    integer, public :: n_var
!    real, public, dimension(:, :, :), pointer :: density
!    real, public, dimension(:, :, :), pointer :: x_speed
!    real, public, dimension(:, :, :), pointer :: y_speed
!    real, public, dimension(:, :, :), pointer :: z_speed
!    real, public, dimension(:, :, :), pointer :: pressure
!    real, public, pointer :: density_inf
!    real, public, pointer :: x_speed_inf
!    real, public, pointer :: y_speed_inf
!    real, public, pointer :: z_speed_inf
!    real, public, pointer :: pressure_inf
!    ! Supersonic flag
!    logical :: supersonic_flag
!    ! Ratio of specific heats (gamma)
!    real, public :: gm
!    ! Specific gas constant
!    real, public :: R_gas
!    ! Constants related to viscosity
!    real, public :: mu_ref, T_ref, Sutherland_temp, Pr
!    integer, public :: ilimiter_switch
!
    ! Including Turbulence variables
!    include "turbulence_models/include/state/variables_deceleration.inc"

    ! Public methods
    public :: setup_state
    public :: destroy_state
!   public :: set_ghost_cell_data
    public :: writestate_vtk

    contains

        subroutine link_aliases()
            implicit none
            call dmsg(1, 'state', 'link_aliases')
            density(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 1)
            x_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 2)
            y_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 3)
            z_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 4)
            pressure(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 5)
            density_inf => qp_inf(1)
            x_speed_inf => qp_inf(2)
            y_speed_inf => qp_inf(3)
            z_speed_inf => qp_inf(4)
            pressure_inf => qp_inf(5)
            include "turbulence_models/include/state/link_aliases.inc"
        end subroutine link_aliases

        subroutine unlink_aliases()
            implicit none
            call dmsg(1, 'state', 'unlink_aliases')
            nullify(density)
            nullify(x_speed)
            nullify(y_speed)
            nullify(z_speed)
            nullify(pressure)
            nullify(density_inf)
            nullify(x_speed_inf)
            nullify(y_speed_inf)
            nullify(z_speed_inf)
            nullify(pressure_inf)
            include "turbulence_models/include/state/unlink_aliases.inc"
        end subroutine unlink_aliases

        subroutine allocate_memory()
            !-----------------------------------------------------------
            ! Allocate memory for the state variables
            !
            ! This assumes that imx and jmx (the grid size) has been set
            ! within the state module.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'state', 'allocate_memory')

            ! The state of the system is defined by the primitive 
            ! variables (density, velocity and pressure) at the grid
            ! cell centers. 
            ! There are (imx - 1) x (jmx - 1) grid cells within the 
            ! domain. We require a row of ghost cells on each boundary.
            ! This current implementation is for a 2D/1D case. 
            call alloc(qp, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for state ' // &
                        'variable qp.')
            call alloc(qp_inf, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for state ' // &
                        'variable qp_inf.')
            !!!!include only when turblent variables are differetn than qp(6:7)
            !include "turbulence_models/sst/state/allocate_memory.inc"
        end subroutine allocate_memory

        subroutine deallocate_memory()

            implicit none

            call dmsg(1, 'state', 'deallocate_memory')

            call dealloc(qp)

        end subroutine deallocate_memory

        subroutine setup_state(free_stream_density, free_stream_x_speed, &
                free_stream_y_speed, free_stream_z_speed, &
                free_stream_pressure, state_file_level)
            !-----------------------------------------------------------
            ! Setup the state module.
            !
            ! This subroutine should be run before the state variables
            ! are initilized. This subroutine allocates the memory for 
            ! state variables and sets up the aliases to refer to the 
            ! components of the state.
            !-----------------------------------------------------------

            implicit none
            real, intent(in) :: free_stream_density
            real, intent(in) :: free_stream_x_speed, free_stream_y_speed, &
                                free_stream_z_speed
            real, intent(in) :: free_stream_pressure
            integer          :: state_file_level

            call dmsg(1, 'state', 'setup_state')

            call allocate_memory()
            call link_aliases()
            call init_infinity_values(free_stream_density, &
                    free_stream_x_speed, free_stream_y_speed, &
                    free_stream_z_speed, free_stream_pressure)
            call set_supersonic_flag()
            call initstate(state_file_level)
            !!!!include only when turblent variables are differetn than qp(6:7)
            !include "turbulence_models/sst/state/deallocate_memory.inc"

        end subroutine setup_state

        subroutine destroy_state()
            !-----------------------------------------------------------
            ! Destroy the state module.
            !
            ! This subroutine destroys the state module which includes
            ! unlinking the aliases for the state components and 
            ! deallocating the memory held by the state variables. 
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'state', 'destroy_state')

            call unlink_aliases()
            call deallocate_memory()

        end subroutine destroy_state

!       subroutine set_inlet_and_exit_state_variables()
!           !-----------------------------------------------------------
!           ! Set / extrapolate inlet and exit state variables
!           !
!           ! The inlet and exit variables are set based on the 
!           ! prescribed infinity (free stream) values and extrapolation
!           ! of the neighbouring cells. The decision to impose or 
!           ! extrapolate is based on the wave speeds and whether or not
!           ! the flow is supersonic.
!           !-----------------------------------------------------------

!           implicit none
!           
!           call dmsg(1, 'state', 'set_inlet_and_exit_state_variables')

!           ! Impose the density, x_speed and y_speed at the inlet
!           density(0, :, :) = density_inf
!           x_speed(0, :, :) = x_speed_inf
!           y_speed(0, :, :) = y_speed_inf
!           z_speed(0, :, :) = z_speed_inf
!           ! Extrapolate these quantities at the exit
!           density(imx, :, :) = density(imx - 1, :, :)
!           x_speed(imx, :, :) = x_speed(imx - 1, :, :)
!           y_speed(imx, :, :) = y_speed(imx - 1, :, :)
!           z_speed(imx, :, :) = z_speed(imx - 1, :, :)
!           ! If the flow is subsonic, impose the back pressure
!           ! Else, impose the inlet pressure
!           ! Extrapolate at the other end
!           if (supersonic_flag .eqv. .TRUE.) then
!               pressure(0, :, :) = pressure_inf
!               pressure(imx, :, :) = pressure(imx - 1, :, :)
!           else
!               pressure(imx, :, :) = pressure_inf
!               pressure(0, :, :) = pressure(1, :, :)
!           end if
!           
!       end subroutine set_inlet_and_exit_state_variables

!       subroutine set_front_and_back_ghost_cell_data()
!           !-----------------------------------------------------------
!           ! Set the state variables for the front and back ghosh cells
!           !
!           ! The pressure and density for the front and back ghost 
!           ! cells is extrapolated from the neighboring interior cells.
!           ! The velocity components are computed by applying the 
!           ! flow tangency conditions.
!           !-----------------------------------------------------------

!           implicit none

!           call dmsg(1, 'state', 'set_front_and_back_ghost_cell_data')

!           pressure(:, 0, :) = pressure(:, jmx-1, :)
!           pressure(:, jmx, :) = pressure(:, 1, :)
!           density(:, 0, :) = density(:, jmx-1, :)
!           density(:, jmx, :) = density(:, 1, :)
!        
!           call apply_eta_flow_tangency_conditions()

!       end subroutine set_front_and_back_ghost_cell_data

!       subroutine set_top_and_bottom_ghost_cell_data()
!           !-----------------------------------------------------------
!           ! Set the state variables for the top and bottom ghosh cells
!           !
!           ! The pressure and density for the top and bottom ghost 
!           ! cells is extrapolated from the neighboring interior cells.
!           ! The velocity components are computed by applying the 
!           ! flow tangency conditions.
!           !-----------------------------------------------------------

!           implicit none

!           call dmsg(1, 'state', 'set_top_and_bottom_ghost_cell_data')

!           pressure(:, :,  0) = pressure(:, :, kmx-1)
!           pressure(:, :, kmx) = pressure(:, :, 1)
!           density(:, :, 0) = density(:, :,  kmx-1)
!           density(:, :, kmx) = density(:, :, 1)
!           call apply_zeta_flow_tangency_conditions()

!       end subroutine set_top_and_bottom_ghost_cell_data
!       
!       subroutine apply_eta_flow_tangency_conditions()
!           !-----------------------------------------------------------
!           ! Apply the flow tangency conditions for the eta face
!           !
!           ! The flow tangency conditions ensure that there is no flow
!           ! across the boundaries. This is done by ensuring that the
!           ! flow is parallel to the boundary.
!           !-----------------------------------------------------------

!           implicit none
!           
!           call dmsg(1, 'state', 'apply_eta_flow_tangency_conditions')

!           if (mu_ref .eq. 0.0) then
!             ! For the back cells
!             x_speed(1:imx-1, jmx, 1:kmx-1) = x_speed(1:imx-1, jmx-1, 1:kmx-1) - &
!                     (2. * &
!                         ((x_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             ynx(1:imx-1, jmx, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) &
!                         + (y_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             yny(1:imx-1, jmx, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) &
!                         + (z_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             ynz(1:imx-1, jmx, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) &
!                         ) &
!                     )
!             y_speed(1:imx-1, jmx, 1:kmx-1) = y_speed(1:imx-1, jmx-1, 1:kmx-1) - &
!                     (2. * &
!                         ((x_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             ynx(1:imx-1, jmx, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) &
!                         + (y_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             yny(1:imx-1, jmx, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) &
!                         + (z_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             ynz(1:imx-1, jmx, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) &
!                         ) &
!                     )
!             z_speed(1:imx-1, jmx, 1:kmx-1) = z_speed(1:imx-1, jmx-1, 1:kmx-1) - &
!                     (2. * &
!                         ((x_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             ynx(1:imx-1, jmx, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
!                         + (y_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             yny(1:imx-1, jmx, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
!                         + (z_speed(1:imx-1, jmx-1, 1:kmx-1) * &
!                             ynz(1:imx-1, jmx, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
!                         ) &
!                     )

!             ! For the front cells
!             x_speed(1:imx-1, 0, 1:kmx-1) = x_speed(1:imx-1, 1, 1:kmx-1) - &
!                  (2. * &
!                      ((x_speed(1:imx-1, 1, 1:kmx-1) * &
!                          ynx(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) &
!                      + (y_speed(1:imx-1, 1, 1:kmx-1) * &
!                          yny(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) &
!                      + (z_speed(1:imx-1, 1, 1:kmx-1) * &
!                          ynz(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) &
!                      ) &
!                  )
!             y_speed(1:imx-1, 0, 1:kmx-1) = y_speed(1:imx-1, 1, 1:kmx-1) - &
!                  (2. * &
!                      ((x_speed(1:imx-1, 1, 1:kmx-1) * &
!                          ynx(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) &
!                      + (y_speed(1:imx-1, 1, 1:kmx-1) * &
!                          yny(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) &
!                      + (z_speed(1:imx-1, 1, 1:kmx-1) * &
!                          ynz(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) &
!                      ) &
!                  )
!             z_speed(1:imx-1, 0, 1:kmx-1) = z_speed(1:imx-1, 1, 1:kmx-1) - &
!                  (2. * &
!                      ((x_speed(1:imx-1, 1, 1:kmx-1) * &
!                          ynx(1:imx-1, 1, 1:kmx-1) * znx(1:imx-1, 1, 1:kmx-1)) &
!                      + (y_speed(1:imx-1, 1, 1:kmx-1) * &
!                          yny(1:imx-1, 1, 1:kmx-1) * znx(1:imx-1, 1, 1:kmx-1)) &
!                      + (z_speed(1:imx-1, 1, 1:kmx-1) * &
!                          ynz(1:imx-1, 1, 1:kmx-1) * znx(1:imx-1, 1, 1:kmx-1)) &
!                      ) &
!                  )
!           else
!             x_speed(1:imx-1, jmx, 1:kmx-1) = - x_speed(1:imx-1, jmx-1, 1:kmx-1)
!             y_speed(1:imx-1, jmx, 1:kmx-1) = - y_speed(1:imx-1, jmx-1, 1:kmx-1)
!             z_speed(1:imx-1, jmx, 1:kmx-1) = - z_speed(1:imx-1, jmx-1, 1:kmx-1)
!        
!             x_speed(1:imx-1, 0, 1:kmx-1) = - x_speed(1:imx-1, 1, 1:kmx-1)
!             y_speed(1:imx-1, 0, 1:kmx-1) = - y_speed(1:imx-1, 1, 1:kmx-1)
!             z_speed(1:imx-1, 0, 1:kmx-1) = - z_speed(1:imx-1, 1, 1:kmx-1)
!           end if
!           
!           x_speed(1:imx-1, 0, 1:kmx-1) = x_speed(1:imx-1, jmx-1, 1:kmx-1)
!           y_speed(1:imx-1, 0, 1:kmx-1) = y_speed(1:imx-1, jmx-1, 1:kmx-1)
!           z_speed(1:imx-1, 0, 1:kmx-1) = z_speed(1:imx-1, jmx-1, 1:kmx-1)
!           x_speed(1:imx-1, jmx, 1:kmx-1) = x_speed(1:imx-1, 1, 1:kmx-1)
!           y_speed(1:imx-1, jmx, 1:kmx-1) = y_speed(1:imx-1, 1, 1:kmx-1)
!           z_speed(1:imx-1, jmx, 1:kmx-1) = z_speed(1:imx-1, 1, 1:kmx-1)

!       end subroutine apply_eta_flow_tangency_conditions

!       subroutine apply_zeta_flow_tangency_conditions()
!           !-----------------------------------------------------------
!           ! Apply the flow tangency conditions
!           !
!           ! The flow tangency conditions ensure that there is no flow
!           ! across the boundaries. This is done by ensuring that the
!           ! flow is parallel to the boundary.
!           !-----------------------------------------------------------

!           implicit none
!           
!           call dmsg(1, 'state', 'apply_zeta_flow_tangency_conditions')

!           ! For the top cells
!           x_speed(1:imx-1, 1:jmx-1, kmx) = x_speed(1:imx-1, 1:jmx-1, kmx-1) - &
!                   (2. * &
!                       ((x_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           znx(1:imx-1, 1:jmx-1, kmx) * znx(1:imx-1, 1:jmx-1, kmx)) &
!                       + (y_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           zny(1:imx-1, 1:jmx-1, kmx) * znx(1:imx-1, 1:jmx-1, kmx)) &
!                       + (z_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           znz(1:imx-1, 1:jmx-1, kmx) * znx(1:imx-1, 1:jmx-1, kmx)) &
!                       ) &
!                   )
!           y_speed(1:imx-1, 1:jmx-1, kmx) = y_speed(1:imx-1, 1:jmx-1, kmx-1) - &
!                   (2. * &
!                       ((x_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           znx(1:imx-1, 1:jmx-1, kmx) * zny(1:imx-1, 1:jmx-1, kmx)) &
!                       + (y_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           zny(1:imx-1, 1:jmx-1, kmx) * zny(1:imx-1, 1:jmx-1, kmx)) &
!                       + (z_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           znz(1:imx-1, 1:jmx-1, kmx) * zny(1:imx-1, 1:jmx-1, kmx)) &
!                       ) &
!                   )
!           z_speed(1:imx-1, 1:jmx-1, kmx) = z_speed(1:imx-1, 1:jmx-1, kmx-1) - &
!                   (2. * &
!                       ((x_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           znx(1:imx-1, 1:jmx-1, kmx) * znz(1:imx-1, 1:jmx-1, kmx)) &
!                       + (y_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           zny(1:imx-1, 1:jmx-1, kmx) * znz(1:imx-1, 1:jmx-1, kmx)) &
!                       + (z_speed(1:imx-1, 1:jmx-1, kmx-1) * &
!                           znz(1:imx-1, 1:jmx-1, kmx) * znz(1:imx-1, 1:jmx-1, kmx)) &
!                       ) &
!                   )
!           
!           x_speed(1:imx-1, 1:jmx-1, kmx) = - x_speed(1:imx-1, 1:jmx-1, kmx-1)
!           y_speed(1:imx-1, 1:jmx-1, kmx) = - y_speed(1:imx-1, 1:jmx-1, kmx-1)
!           z_speed(1:imx-1, 1:jmx-1, kmx) = - z_speed(1:imx-1, 1:jmx-1, kmx-1)

!           ! For the bottom cells
!           x_speed(1:imx-1, 1:jmx-1, 0) = x_speed(1:imx-1, 1:jmx-1, 1) - &
!                   (2. * &
!                       ((x_speed(1:imx-1, 1:jmx-1, 1) * &
!                           znx(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) &
!                       + (y_speed(1:imx-1, 1:jmx-1, 1) * &
!                           zny(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) &
!                       + (z_speed(1:imx-1, 1:jmx-1, 1) * &
!                           znz(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) &
!                       ) &
!                   )
!           y_speed(1:imx-1, 1:jmx-1, 0) = y_speed(1:imx-1, 1:jmx-1, 1) - &
!                   (2. * &
!                       ((x_speed(1:imx-1, 1:jmx-1, 1) * &
!                           znx(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) &
!                       + (y_speed(1:imx-1, 1:jmx-1, 1) * &
!                           zny(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) &
!                       + (z_speed(1:imx-1, 1:jmx-1, 1) * &
!                           znz(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) &
!                       ) &
!                   )
!           z_speed(1:imx-1, 1:jmx-1, 0) = z_speed(1:imx-1, 1:jmx-1, 1) - &
!                   (2. * &
!                       ((x_speed(1:imx-1, 1:jmx-1, 1) * &
!                           znx(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
!                       + (y_speed(1:imx-1, 1:jmx-1, 1) * &
!                           zny(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
!                       + (z_speed(1:imx-1, 1:jmx-1, 1) * &
!                           znz(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
!                       ) &
!                   )

!           x_speed(1:imx-1, 1:jmx-1, 0) = - x_speed(1:imx-1, 1:jmx-1, 1)
!           y_speed(1:imx-1, 1:jmx-1, 0) = - y_speed(1:imx-1, 1:jmx-1, 1)
!           z_speed(1:imx-1, 1:jmx-1, 0) = - z_speed(1:imx-1, 1:jmx-1, 1)

!           x_speed(1:imx-1, 1:jmx-1, 0) = x_speed(1:imx-1, 1:jmx-1, kmx-1)
!           y_speed(1:imx-1, 1:jmx-1, 0) = y_speed(1:imx-1, 1:jmx-1, kmx-1)
!           z_speed(1:imx-1, 1:jmx-1, 0) = z_speed(1:imx-1, 1:jmx-1, kmx-1)
!           x_speed(1:imx-1, 1:jmx-1, kmx) = x_speed(1:imx-1, 1:jmx-1, 1)
!           y_speed(1:imx-1, 1:jmx-1, kmx) = y_speed(1:imx-1, 1:jmx-1, 1)
!           z_speed(1:imx-1, 1:jmx-1, kmx) = z_speed(1:imx-1, 1:jmx-1, 1)

!       end subroutine apply_zeta_flow_tangency_conditions
!       
!       subroutine set_ghost_cell_data()
!           !-----------------------------------------------------------
!           ! Set the data in the ghost cells
!           !
!           ! The ghost cell data is either imposed or extrapolated from
!           ! the neighboring cells inside the domain. The decision to
!           ! impose or extrapolate is taken based on certain 
!           ! conditions.
!           !-----------------------------------------------------------

!           implicit none

!           call dmsg(1, 'state', 'set_ghost_cell_data')

!           call set_inlet_and_exit_state_variables()
!           call set_front_and_back_ghost_cell_data()
!           call set_top_and_bottom_ghost_cell_data()

!       end subroutine set_ghost_cell_data

        subroutine init_infinity_values(free_stream_density, &
                free_stream_x_speed, free_stream_y_speed, &
                free_stream_z_speed, free_stream_pressure)
            !-----------------------------------------------------------
            ! Set the values of the infinity variables
            !-----------------------------------------------------------

            implicit none
            real, intent(in) :: free_stream_density
            real, intent(in) :: free_stream_x_speed
            real, intent(in) :: free_stream_y_speed
            real, intent(in) :: free_stream_z_speed
            real, intent(in) :: free_stream_pressure
            
            call dmsg(1, 'state', 'init_infinity_values')

            density_inf = free_stream_density
            x_speed_inf = free_stream_x_speed
            y_speed_inf = free_stream_y_speed
            z_speed_inf = free_stream_z_speed
            pressure_inf = free_stream_pressure
            include "turbulence_models/include/state/init_infinity_values.inc"

        end subroutine init_infinity_values
        
        function sound_speed_inf() result(a)
            !-----------------------------------------------------------
            ! Return the free stream speed of sound.
            !-----------------------------------------------------------

            implicit none
            real :: a

            a = sqrt(gm * pressure_inf / density_inf)

        end function sound_speed_inf

        subroutine set_supersonic_flag()
            !-----------------------------------------------------------
            ! Set the supersonic flag based on the infinity conditions.
            !
            ! In the ghost cells, the values of the primitive variables
            ! are set either based on the neighbouring cells' values or 
            ! the infinity values. This splitting is based on the wave
            ! speeds. The supersonic flag is used as an indication as 
            ! to how these cells should get their values. 
            !-----------------------------------------------------------

            implicit none
            real :: avg_inlet_mach

            call dmsg(1, 'state', 'set_supersonic_flag')

            avg_inlet_mach = sqrt(x_speed_inf ** 2. + y_speed_inf ** 2. + &
                                  z_speed_inf ** 2.) / sound_speed_inf()

            if (avg_inlet_mach >= 1) then
                supersonic_flag = .TRUE.
            else
                supersonic_flag = .FALSE.
            end if

            call dmsg(5, 'state', 'set_supersonic_flag', &
                    'Supersonic flag set to ' + supersonic_flag)

        end subroutine set_supersonic_flag

        subroutine initstate(state_file_level)
            !-----------------------------------------------------------
            ! Initialize the state
            !
            ! If state_file is a tilde (~), then the state should be 
            ! set to the infinity values. Otherwise, read the state_file
            ! to get the state values.
            !-----------------------------------------------------------

            implicit none
            integer                         :: state_file_level
            character(len=FILE_NAME_LENGTH) :: state_file
            
            call dmsg(1, 'state', 'initstate')

            if (state_file_level .eq. 0) then
                ! Set the state to the infinity values
                call init_state_with_infinity_values()
            else
                write(state_file,'(a,i2.2,a,i5.5,a)') &
                  "results/process_",process_id,"/output",state_file_level,".vtk"
                call readstate_vtk(state_file)
            end if

        end subroutine initstate

        subroutine init_state_with_infinity_values()
            !-----------------------------------------------------------
            ! Initialize the state based on the infinity values.
            !-----------------------------------------------------------
            
            implicit none
            integer :: i
            
            call dmsg(1, 'state', 'init_state_with_infinity_values')
            
            do i = 1,n_var
                qp(:, :, :, i) = qp_inf(i)
            end do 
            !!!!include only when turblent variables are differetn than qp(6:7)
            !include "turbulence_models/sst/state/init_state_with_infinity_values.inc"
            
        end subroutine init_state_with_infinity_values

        subroutine writestate_vtk(outfile, comment, extravar, extravarname)
            !-----------------------------------------------------------
            ! Write the state of the system to a file
            !-----------------------------------------------------------

            implicit none
            character(len=FILE_NAME_LENGTH), intent(in) :: outfile
            character(len=DESCRIPTION_STRING_LENGTH), optional, intent(in) :: &
                    comment
            real, dimension(:, :, :), optional, intent(in) :: extravar
            character(len=DESCRIPTION_STRING_LENGTH), optional, intent(in) :: extravarname
            integer :: i, j, k
            
            call dmsg(1, 'state', 'writestate_vtk')

            open(OUT_FILE_UNIT, file=outfile + '.part')

            write(OUT_FILE_UNIT, fmt='(a)') '# vtk DataFile Version 3.1'

            if (present(comment)) then
                write(OUT_FILE_UNIT, fmt='(a)') & 
                        trim('cfd-iitm output (' + comment + ')')
            else
                write(OUT_FILE_UNIT, '(a)') 'cfd-iitm output'
            end if

            write(OUT_FILE_UNIT, '(a)') 'ASCII'
            write(OUT_FILE_UNIT, '(a)') 'DATASET STRUCTURED_GRID'
            write(OUT_FILE_UNIT, *) 

            write(OUT_FILE_UNIT, fmt='(a, i0, a, i0, a, i0)') &
                'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
            write(OUT_FILE_UNIT, fmt='(a, i0, a)') &
                'POINTS ', imx*jmx*kmx, ' DOUBLE'

            do k = 1, kmx
             do j = 1, jmx
              do i = 1, imx
                write(OUT_FILE_UNIT, fmt='(f0.16, a, f0.16, a, f0.16)') &
                    grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k)
              end do
             end do
            end do
            write(OUT_FILE_UNIT, *) 

            ! Cell data
            write(OUT_FILE_UNIT, fmt='(a, i0)') &
                'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)

            ! Writing Velocity
            write(OUT_FILE_UNIT, '(a)') 'VECTORS Velocity FLOAT'
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                write(OUT_FILE_UNIT, fmt='(f0.16, a, f0.16, a, f0.16)') &
                    x_speed(i, j, k), ' ', y_speed(i, j, k), ' ', z_speed(i, j, k)
              end do
             end do
            end do
            write(OUT_FILE_UNIT, *) 

            ! Writing Density
            write(OUT_FILE_UNIT, '(a)') 'SCALARS Density FLOAT'
            write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                write(OUT_FILE_UNIT, fmt='(f0.16)') density(i, j, k)
              end do
             end do
            end do
            write(OUT_FILE_UNIT, *) 

            ! Writing Pressure
            write(OUT_FILE_UNIT, '(a)') 'SCALARS Pressure FLOAT'
            write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                write(OUT_FILE_UNIT, fmt='(f0.16)') pressure(i, j, k)
              end do
             end do
            end do
            write(OUT_FILE_UNIT, *) 

            ! Signed distance if IB
            if (present(extravar)) then
                write(OUT_FILE_UNIT, '(a)') &
                    trim('SCALARS '), trim(extravarname), ' FLOAT'
                write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
                do k = 1, kmx - 1
                 do j = 1, jmx - 1
                  do i = 1, imx - 1
                    write(OUT_FILE_UNIT, fmt='(f0.16)') extravar(i, j, k)
                  end do
                 end do
                end do 
                write(OUT_FILE_UNIT, *) 
            end if

            include "turbulence_models/include/state/writestate_vtk.inc"
            
            close(OUT_FILE_UNIT)

            call rename(outfile + '.part', outfile)

        end subroutine writestate_vtk

        subroutine readstate_vtk(state_file)
            !-----------------------------------------------------------
            ! Read the state of the system from a file
            !-----------------------------------------------------------

            implicit none
            character(len=FILE_NAME_LENGTH), intent(in) :: state_file
            integer :: i, j, k
            
            call dmsg(1, 'state', 'readstate_vtk')

            open(OUT_FILE_UNIT, file=state_file)

            read(OUT_FILE_UNIT, *) ! Skip first line
            read(OUT_FILE_UNIT, *) ! Skip comment
            read(OUT_FILE_UNIT, *) ! Skip ASCII
            read(OUT_FILE_UNIT, *) ! Skip DATASET
            read(OUT_FILE_UNIT, *) ! Skip Extra line

            read(OUT_FILE_UNIT, *) ! Skip DIMENSIONS
            read(OUT_FILE_UNIT, *) ! Skip POINTS
            do k = 1, kmx
             do j = 1, jmx
              do i = 1, imx
                read(OUT_FILE_UNIT, *) ! Skip grid points
              end do
             end do
            end do
            read(OUT_FILE_UNIT, *) ! Skip blank space

            ! Cell data
            read(OUT_FILE_UNIT, *) ! Skip CELL_DATA
            read(OUT_FILE_UNIT, *) ! Skip VECTORS Velocity
 
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                read(OUT_FILE_UNIT, *) x_speed(i, j, k), y_speed(i, j, k), z_speed(i, j, k)
              end do
             end do
            end do

            read(OUT_FILE_UNIT, *) ! Skip Blank line
            read(OUT_FILE_UNIT, *) ! Skip SCALARS DENSITY
            read(OUT_FILE_UNIT, *) ! Skip LOOKUP_TABLE
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                read(OUT_FILE_UNIT, *) density(i, j, k)
              end do
             end do
            end do

            read(OUT_FILE_UNIT, *) ! Skip Blank line
            read(OUT_FILE_UNIT, *) ! Skip SCALARS Pressure
            read(OUT_FILE_UNIT, *) ! Skip LOOKUP_TABLE
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                read(OUT_FILE_UNIT, *) pressure(i, j, k)
              end do
             end do
            end do
            ! Extra var not needed for state
            include "turbulence_models/include/state/readstate_vtk.inc"

            close(OUT_FILE_UNIT)

        end subroutine readstate_vtk

end module state
