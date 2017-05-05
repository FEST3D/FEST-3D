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
    use global_vars, only : start_from
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
    use global_vars, only : infile
    
    use global_vars, only  : free_stream_density
    use global_vars, only  : free_stream_x_speed
    use global_vars, only  : free_stream_y_speed
    use global_vars, only  : free_stream_z_speed
    use global_vars, only  : free_stream_pressure
    use global_vars, only  : free_stream_tk
    use global_vars, only  : free_stream_tw

  use global_vars, only: mass_residue
  use global_vars, only: x_mom_residue
  use global_vars, only: y_mom_residue
  use global_vars, only: z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only:    vis_resnorm
  use global_vars, only:   cont_resnorm
  use global_vars, only:  x_mom_resnorm
  use global_vars, only:  y_mom_resnorm
  use global_vars, only:  z_mom_resnorm
  use global_vars, only: energy_resnorm

    use utils, only: alloc, dealloc, dmsg
    use layout, only: process_id
    use string
  use read_output, only: read_file

    implicit none
    private


    real :: speed_inf
    ! Public methods
    public :: setup_state
    public :: destroy_state
!   public :: set_ghost_cell_data

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
        end subroutine allocate_memory

        subroutine deallocate_memory()

            implicit none

            call dmsg(1, 'state', 'deallocate_memory')

            call dealloc(qp)

        end subroutine deallocate_memory

        subroutine setup_state()
            !-----------------------------------------------------------
            ! Setup the state module.
            !
            ! This subroutine should be run before the state variables
            ! are initilized. This subroutine allocates the memory for 
            ! state variables and sets up the aliases to refer to the 
            ! components of the state.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'state', 'setup_state')

            call allocate_memory()
            call link_aliases()
            call init_infinity_values()
            call set_supersonic_flag()
            call initstate()

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

        subroutine init_infinity_values()
            !-----------------------------------------------------------
            ! Set the values of the infinity variables
            !-----------------------------------------------------------

            implicit none
            
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

        subroutine initstate()
            !-----------------------------------------------------------
            ! Initialize the state
            !
            ! If state_file is a tilde (~), then the state should be 
            ! set to the infinity values. Otherwise, read the state_file
            ! to get the state values.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'state', 'initstate')

            if (start_from .eq. 0) then
                ! Set the state to the infinity values
                call init_state_with_infinity_values()
            else
                write(infile,'(a,i4.4,a,i2.2)') &
                  "time_directories/",start_from,"/process_",process_id
                !call readstate_vtk(state_file)
                call read_file()

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

!        subroutine readstate_vtk(state_file)
!            !-----------------------------------------------------------
!            ! Read the state of the system from a file
!            !-----------------------------------------------------------
!
!            implicit none
!            character(len=FILE_NAME_LENGTH), intent(in) :: state_file
!            integer :: i, j, k
!            
!            call dmsg(1, 'state', 'readstate_vtk')
!
!            open(OUT_FILE_UNIT, file=state_file)
!
!            read(OUT_FILE_UNIT, *) ! Skip first line
!            read(OUT_FILE_UNIT, *) ! Skip comment
!            read(OUT_FILE_UNIT, *) ! Skip ASCII
!            read(OUT_FILE_UNIT, *) ! Skip DATASET
!            read(OUT_FILE_UNIT, *) ! Skip Extra line
!
!            read(OUT_FILE_UNIT, *) ! Skip DIMENSIONS
!            read(OUT_FILE_UNIT, *) ! Skip POINTS
!            do k = 1, kmx
!             do j = 1, jmx
!              do i = 1, imx
!                read(OUT_FILE_UNIT, *) ! Skip grid points
!              end do
!             end do
!            end do
!            read(OUT_FILE_UNIT, *) ! Skip blank space
!
!            ! Cell data
!            read(OUT_FILE_UNIT, *) ! Skip CELL_DATA
!            read(OUT_FILE_UNIT, *) ! Skip VECTORS Velocity
! 
!            do k = 1, kmx - 1
!             do j = 1, jmx - 1
!              do i = 1, imx - 1
!                read(OUT_FILE_UNIT, *) x_speed(i, j, k), y_speed(i, j, k), z_speed(i, j, k)
!              end do
!             end do
!            end do
!
!            read(OUT_FILE_UNIT, *) ! Skip Blank line
!            read(OUT_FILE_UNIT, *) ! Skip SCALARS DENSITY
!            read(OUT_FILE_UNIT, *) ! Skip LOOKUP_TABLE
!            do k = 1, kmx - 1
!             do j = 1, jmx - 1
!              do i = 1, imx - 1
!                read(OUT_FILE_UNIT, *) density(i, j, k)
!              end do
!             end do
!            end do
!
!            read(OUT_FILE_UNIT, *) ! Skip Blank line
!            read(OUT_FILE_UNIT, *) ! Skip SCALARS Pressure
!            read(OUT_FILE_UNIT, *) ! Skip LOOKUP_TABLE
!            do k = 1, kmx - 1
!             do j = 1, jmx - 1
!              do i = 1, imx - 1
!                read(OUT_FILE_UNIT, *) pressure(i, j, k)
!              end do
!             end do
!            end do
!            ! Extra var not needed for state
!            include "turbulence_models/include/state/readstate_vtk.inc"
!
!            close(OUT_FILE_UNIT)
!
!        end subroutine readstate_vtk

end module state
