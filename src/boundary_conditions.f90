module boundary_conditions
    !-------------------------------------------------------------------
    ! The boundary_conditions module contains methods that will enforce
    ! the required boundary conditions on the state variables. 
    !
    ! The boundary conditions can be grouped into the following major
    ! categories:
    ! a. Fix
    !    This involves fixing a parameter at the given boundary.
    ! b. Copy
    !    This involves copying over the value of the parameter from the 
    !    neighbouring cell.
    ! c. Misc
    !    The following are the other types of conditions implemented.
    !    i. Using equations to compute parameters:
    !       After fixing and copying some of the variables, equations like
    !       Bernoulli's equation or total enthalpy conservation can be used
    !       to compute the remaining parameters.
    !    ii. Flow tangency
    !       The flow tangency condition is a type of vector-symmetry 
    !       condition. It is applied to the flow velocity to make the 
    !       velocity at the wall parallel to it. 
    !
    ! This module provides subroutines that can be used to effect these
    ! boundary conditions. The region over which the conditions are 
    ! enforced needs to be specified. Currently, this can be done using
    ! the face argument. Valid options for the face argument are:
    !   "imin", "imax", "jmin", "jmax"
    ! where the first character (i / j) denotes the direction (xi / eta)
    ! and the final three (min / max) denote which side of the domain.
    !-------------------------------------------------------------------
    
    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx

    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z

    use global_vars, only : n_var
    use global_vars, only : qp
    use global_vars, only : density
    use global_vars, only : x_speed
    use global_vars, only : y_speed
    use global_vars, only : z_speed
    use global_vars, only : pressure
    use global_vars, only : tk
    use global_vars, only : tw
    use global_vars, only : density_inf
    use global_vars, only : x_speed_inf
    use global_vars, only : y_speed_inf
    use global_vars, only : z_speed_inf
    use global_vars, only : pressure_inf
    use global_vars, only : tk_inf
    use global_vars, only : tw_inf
    use global_vars, only : mu_ref
    use global_vars, only : mu_t
    use global_vars, only : R_gas
    use global_vars, only : T_ref
    use global_vars, only : Sutherland_temp
    use global_vars, only : turbulence
    use global_vars, only : dist
    use global_vars, only : process_id

    use global_sst , only : beta1

    use utils, only: alloc, dealloc, dmsg
    use bitwise
    use global, only: BOUNDARY_CONDITIONS_FILE_UNIT, STRING_BUFFER_LENGTH, &
            FILE_NAME_LENGTH
!    use grid, only: imx, jmx, kmx
!    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz
!    use state, only: density, x_speed, y_speed, z_speed, pressure, &
!            density_inf, x_speed_inf, y_speed_inf, z_speed_inf, &
!            pressure_inf, qp, n_var
    use face_interpolant, only: x_x_speed_left, x_x_speed_right, &
            x_y_speed_left, x_y_speed_right, x_z_speed_left, x_z_speed_right, &
            y_x_speed_left, y_x_speed_right, y_y_speed_left, y_y_speed_right, &
            y_z_speed_left, y_z_speed_right, z_x_speed_left, z_x_speed_right, &
            z_y_speed_left, z_y_speed_right, z_z_speed_left, z_z_speed_right

    use global_sst, only : sst_F1

    include "turbulence_models/include/bc/import_module.inc" 
    implicit none
    private
    ! Boundary condition (bc) descriptor variables
    ! If a common boundary condition it to be applied to all the faces
    ! in a set, a single element integer array will be used. 
    ! Otherwise, there will be as many integers as faces in the set, 
    ! each denoting the boundary conditions applicable to that face.
    integer, public, dimension(:, :), allocatable :: &
            bc_imn, bc_imx, &
            bc_jmn, bc_jmx, &
            bc_kmn, bc_kmx
    ! Variables to store extra parameters for the FIX type conditions
    ! (e.g. pressure to fix the value to)
    real, dimension(:, :, :), allocatable :: &
            bc_imn_fix_values, bc_imx_fix_values, &
            bc_jmn_fix_values, bc_jmx_fix_values, &
            bc_kmn_fix_values, bc_kmx_fix_values
    ! Boundary conditions list
    include "bclist_declaration.inc"

    ! Public methods
    public :: setup_boundary_conditions
    public :: destroy_boundary_conditions
    public :: apply_boundary_conditions
    public :: set_wall_bc_at_faces

    contains

        subroutine setup_boundary_conditions(bc_filename)
            implicit none
            character(len=*), intent(in) :: bc_filename
            call dmsg(1, 'boundary_conditions', 'setup_boundary_conditions')
            include "bclist_definition.inc"
            open(BOUNDARY_CONDITIONS_FILE_UNIT, file=bc_filename)
            ! Ignore the file header
            read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
            read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
            read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
            ! Continue setting up sides until there are sides to setup
            do while (setup_next_side() .neqv. .FALSE.)
            end do
            call initialize_remaining()
            close(BOUNDARY_CONDITIONS_FILE_UNIT)
            call dmsg(1, 'boundary_conditions', 'setup_boundary_conditions', 'done')
        end subroutine setup_boundary_conditions

        subroutine destroy_boundary_conditions()

            implicit none
            
            call dmsg(1, 'boundary_conditions', 'destroy_boundary_conditions')

            call dealloc(bc_imn)
            call dealloc(bc_imx)
            call dealloc(bc_jmn)
            call dealloc(bc_jmx)
            call dealloc(bc_kmn)
            call dealloc(bc_kmx)
            call dealloc(bc_imn_fix_values)
            call dealloc(bc_imx_fix_values)
            call dealloc(bc_jmn_fix_values)
            call dealloc(bc_jmx_fix_values)
            call dealloc(bc_kmn_fix_values)
            call dealloc(bc_kmx_fix_values)

        end subroutine destroy_boundary_conditions

        subroutine initialize_remaining()
            implicit none
            ! Allocate memory to store just a single integer
            if (allocated(bc_imn) .eqv. .FALSE.) then
                call alloc(bc_imn, 1, 1, 1, 1, &
                        errmsg='Error: Unable to allocate memory for side_bc.')
                bc_imn = 0
            end if
            if (allocated(bc_imx) .eqv. .FALSE.) then
                call alloc(bc_imx, 1, 1, 1, 1, &
                        errmsg='Error: Unable to allocate memory for side_bc.')
                bc_imx = 0
            end if
            if (allocated(bc_jmn) .eqv. .FALSE.) then
                call alloc(bc_jmn, 1, 1, 1, 1, &
                        errmsg='Error: Unable to allocate memory for side_bc.')
                bc_jmn = 0
            end if
            if (allocated(bc_jmx) .eqv. .FALSE.) then
                call alloc(bc_jmx, 1, 1, 1, 1, &
                        errmsg='Error: Unable to allocate memory for side_bc.')
                bc_jmx = 0
            end if
            if (allocated(bc_kmn) .eqv. .FALSE.) then
                call alloc(bc_kmn, 1, 1, 1, 1, &
                        errmsg='Error: Unable to allocate memory for side_bc.')
                bc_kmn = 0
            end if
            if (allocated(bc_kmx) .eqv. .FALSE.) then
                call alloc(bc_kmx, 1, 1, 1, 1, &
                        errmsg='Error: Unable to allocate memory for side_bc.')
                bc_kmx = 0
            end if
        end subroutine initialize_remaining

        function get_bc_token_value(token) result(val)
            implicit none
            character(len=*), intent(in) :: token
            integer :: val
            !TODO: Change this if-ladder implementation to a switch case
            include "bc_token_value_map_implementation.inc"
        end function get_bc_token_value

        subroutine get_next_boundary_condition_value(bc_val, fix_val_type, fix_val)
            implicit none
            character(len=STRING_BUFFER_LENGTH) :: buf
            integer, intent(out) :: bc_val
            integer, intent(out) :: fix_val_type
            real, intent(out) :: fix_val
            integer :: ios
            !TODO: Update this to allow reading a value with the condition also (done) -> Update this to make the statements more readable. 
            ! Idea: The condition string is of the form: 
            !     - condition value
            ! Find the where the second space occurs to separate the value from the condition
            ! Extract the value and convert the condition string into its numerical value
            read(BOUNDARY_CONDITIONS_FILE_UNIT, '(A)', iostat=ios) buf
!            print *, "Read boundary condition: ", trim(buf)
            if (ios /= 0) then
                print *, 'Error while reading a boundary condition name' // &
                        ' in the boundary conditions file.'
                print *, 'Current buffer length is set to: ', &
                        STRING_BUFFER_LENGTH
                stop
            end if
            if (buf(1:2) == '- ') then
                read(buf(index(buf(3:), ' ')+3:), *, iostat=ios) fix_val
                if (ios /= 0) then
                    fix_val = -1
                end if
                write(6,"(A,f0.4)") "  >> Read bc: "//trim(buf)//&
                                 ", Read value:- ", fix_val
                if (.false.) then
                    bc_val = 0
                end if
                bc_val = get_bc_token_value(buf(3:index(buf(3:), ' ')+1))
                if (bc_val == BC_FIX_DENSITY) then
                    fix_val_type = BC_FIX_DENSITY
                else if (bc_val == BC_FIX_X_SPEED) then
                    fix_val_type = BC_FIX_X_SPEED
                else if (bc_val == BC_FIX_Y_SPEED) then
                    fix_val_type = BC_FIX_Y_SPEED
                else if (bc_val == BC_FIX_Z_SPEED) then
                    fix_val_type = BC_FIX_Z_SPEED
                else if (bc_val == BC_FIX_PRESSURE) then
                    fix_val_type = BC_FIX_PRESSURE
                else
                    fix_val_type = -1
                end if
            else
                bc_val = 0
                fix_val_type = -1
            end if
            if (fix_val_type > 0) then
                ! The variable number in the array is the logarithm of fix_val_type 
                ! (to the base 2) because fix_val_type holds the integer corresponding 
                ! to the boundary condition we are looking at. 
                fix_val_type = nint(log(real(fix_val_type)) / log(2.)) + 1
                if (fix_val == -1) then
!                    print *, 'Fix value not specified. Using initial condition'
                    if (bc_val == BC_FIX_DENSITY) then
                        fix_val = density_inf
                    else if (bc_val == BC_FIX_X_SPEED) then
                        fix_val = x_speed_inf 
                    else if (bc_val == BC_FIX_Y_SPEED) then
                        fix_val = y_speed_inf 
                    else if (bc_val == BC_FIX_Z_SPEED) then
                        fix_val = z_speed_inf 
                    else if (bc_val == BC_FIX_PRESSURE) then
                        fix_val = pressure_inf
                    end if
                end if
            end if
        end subroutine get_next_boundary_condition_value

        subroutine setup_common_boundary_conditions(side_bc, side_bc_fix_vals)
            implicit none
            integer, dimension(:, :), allocatable, intent(out) :: side_bc
            real, dimension(:, :, :), allocatable, intent(out) :: side_bc_fix_vals
            integer :: bc_val
            integer :: fix_val_type
            real :: bc_fix_val
            ! Allocate memory to store just a single integer
            call alloc(side_bc, 1, 1, 1, 1, &
                    errmsg='Error: Unable to allocate memory for side_bc.')
            ! Allocate memory to store the fix type conditions values
            call alloc(side_bc_fix_vals, 1, 1, 1, 1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for side_bc_fix_vals')
            ! Initialize the boundary condition variables
            side_bc = 0
            side_bc_fix_vals = 0
            ! bc_val represents the integer corresponding to the 
            ! boundary condition just read from the file.
            ! Initialize it to -1 so that entering the loop is 
            ! possible.
            bc_val = -1
            do while (bc_val /= 0)
                call get_next_boundary_condition_value(bc_val, fix_val_type, bc_fix_val)
                side_bc(1, 1) = side_bc(1, 1) .or. bc_val
                if (fix_val_type > 0) then
                    side_bc_fix_vals(1, 1, fix_val_type) = bc_fix_val
                end if
            end do
        end subroutine setup_common_boundary_conditions
        
 !      subroutine get_next_range(startidx, stopidx, dim_max)
 !          !-----------------------------------------------------------
 !          ! Get next range from the boundary conditions config file
 !          !
 !          ! In the case of range based boundary conditions, the config
 !          ! file will have different boundary conditions for different
 !          ! ranges. This subroutine will extract the range of face 
 !          ! numbers that the following conditions will be applicable 
 !          ! to.
 !          !-----------------------------------------------------------
 !          implicit none
 !          integer, intent(inout) :: startidx, stopidx
 !          integer, intent(in) :: dim_max
 !          character(len=STRING_BUFFER_LENGTH) :: buf
 !          character(len=3) :: temp
 !          integer :: ios

 !          ! Read the range line
 !          read(BOUNDARY_CONDITIONS_FILE_UNIT, '(A)', iostat=ios) buf
 !          print *, "Range line: ", buf
 !          if (ios /= 0) then
 !              print *, 'Error while reading a boundary condition name' // &
 !                      ' in the boundary conditions file.'
 !              print *, 'Current buffer length is set to: ', &
 !                      STRING_BUFFER_LENGTH
 !              stop
 !          end if

 !          startidx = 0
 !          stopidx  0

 !          ! Try tokenizing the line assuming all components are there
 !          read(buf, *, iostat=ios) temp, startidx, temp, stopidx
 !          if (ios /= 0) then
 !              ! At least one component is not present
 !              ! Assume no start idx is present
 !              read(buf, *, iostat=ios) temp, temp, stopidx
 !              if (ios /= 0) then
 !                  ! Something went wrong
 !                  ! Assume stopidx is not present
 !                  read(buf, *, iostat=ios) temp, startidx, temp
 !                  if (ios /= 0) then
 !                      ! Something went wrong again
 !                      ! Check if this line is a valid range line
 !                      if (buf == '##R :') then
 !                          ! This is for the entire domain
 !                          ! Why didn't this guy just use the common 
 !                          ! type boundary condition?! Would have saved
 !                          ! him some memory... But, hey, if he wants 
 !                          ! it this way, sure.
 !                          startidx = 1
 !                          stopidx = dim_max
 !                          print *, "For entire domain"
 !                      else
 !                          ! No more ranges to read
 !                          stopidx = 0
 !                          startidx = 0
 !                          print *, "No more ranges to read"
 !                      end if
 !                  else
 !                      stopidx = dim_max
 !                      print *, "stopidx was not present"
 !                  end if
 !              else
 !                  startidx = 1
 !                  print *, "startidx was not present"
 !              endif
 !          end if
 !          print *, "startidx = ", startidx
 !          print *, "stopidx = ", stopidx
 !          
 !      end subroutine get_next_range

 !      subroutine setup_range_based_boundary_conditions(side_bc, side_bc_fix_vals, dim_max)
 !          implicit none
 !          integer, dimension(:), allocatable, intent(out) :: side_bc
 !          real, dimension(:, :), allocatable, intent(out) :: side_bc_fix_vals
 !          integer, intent(in) :: dim_max
 !          integer :: startidx, stopidx
 !          integer :: bc_val
 !          integer :: fix_val_type
 !          real :: bc_fix_val
 !          ! There are multiple ranges. Allocate enough memory to 
 !          ! store boundary conditions individually for faces. 
 !          call alloc(side_bc, 1, dim_max, &
 !                  errmsg='Error: Unable to allocate memory for side_bc.')
 !          call alloc(side_bc_fix_vals, 1, dim_max, 1, 4, &
 !                  errmsg='Error: Unable to allocate memory for side_bc_fix_vals')
 !          ! Initialize each faces' boundary conditions
 !          side_bc = 0
 !          side_bc_fix_vals = 0
 !          startidx = -1
 !          stopidx = -1
 !          call get_next_range(startidx, stopidx, dim_max)
 !          do while (stopidx /= 0)
 !              bc_val = -1
 !              do while (bc_val /= 0)
 !                  call get_next_boundary_condition_value(bc_val, fix_val_type, bc_fix_val)
 !                  side_bc(startidx:stopidx) = &
 !                          side_bc(startidx:stopidx) .or. bc_val
 !                  if (fix_val_type > 0) then
 !                      side_bc_fix_vals(startidx:stopidx, fix_val_type) = bc_fix_val
 !                  end if
 !              end do
 !              call get_next_range(startidx, stopidx, dim_max)
 !          end do
 !      end subroutine setup_range_based_boundary_conditions

        function setup_next_side() result(found_new_boundary)
            implicit none
            character(len=STRING_BUFFER_LENGTH) :: buf
            integer :: ios
            logical :: found_new_boundary
            ! Let the boundary directive line remain in the buffer so 
            ! that the setup_side subroutine can read it.
            read(BOUNDARY_CONDITIONS_FILE_UNIT, "(A)", iostat=ios) buf
            write(6,'(A)') ">>> Boundary directive line: "//trim(buf)
            select case (buf)
                case ("# imn")
                    call setup_common_boundary_conditions(bc_imn, bc_imn_fix_values)
                    found_new_boundary = .TRUE.
                case ("# imx")
                    call setup_common_boundary_conditions(bc_imx, bc_imx_fix_values)
                    found_new_boundary = .TRUE.
                case ("# jmn")
                    call setup_common_boundary_conditions(bc_jmn, bc_jmn_fix_values)
                    found_new_boundary = .TRUE.
                case ("# jmx")
                    call setup_common_boundary_conditions(bc_jmx, bc_jmx_fix_values)
                    found_new_boundary = .TRUE.
                case ("# kmn")
                    call setup_common_boundary_conditions(bc_kmn, bc_kmn_fix_values)
                    found_new_boundary = .TRUE.
                case ("# kmx")
                    call setup_common_boundary_conditions(bc_kmx, bc_kmx_fix_values)
                    found_new_boundary = .TRUE.
 !              case ("# imn +")
 !                  call setup_range_based_boundary_conditions(bc_imn, bc_imn_fix_values, jmx)
 !                  found_new_boundary = .TRUE.
 !              case ("# imx +")
 !                  call setup_range_based_boundary_conditions(bc_imx, bc_imx_fix_values, jmx)
 !                  found_new_boundary = .TRUE.
 !              case ("# jmn +")
 !                  call setup_range_based_boundary_conditions(bc_jmn, bc_jmn_fix_values, imx)
 !                  found_new_boundary = .TRUE.
 !              case ("# jmx +")
 !                  call setup_range_based_boundary_conditions(bc_jmx, bc_jmx_fix_values, imx)
 !                  found_new_boundary = .TRUE.
                case default
                    found_new_boundary = .FALSE.
            end select
        end function setup_next_side

        subroutine apply_boundary_conditions()
            implicit none
            integer :: current_condition, i, j, cell_ind, m, n
            integer, dimension(1:2) :: n_shape
            integer :: imin_flow_tangency_flag=0
            integer :: imax_flow_tangency_flag=0
            integer :: jmin_flow_tangency_flag=0
            integer :: jmax_flow_tangency_flag=0
            integer :: kmin_flow_tangency_flag=0
            integer :: kmax_flow_tangency_flag=0

            call dmsg(1, 'boundary_conditions', 'apply_boundary conditions')
            current_condition = 1
            do while (current_condition /= 0)
                ! Check if current_condition needs to be applied
                ! Apply it
                ! Go to the next condition (bit shift to left)
                n_shape = shape(bc_imn)
                m = n_shape(1)
                n = n_shape(2)
                if ((m == 1) .and. (n == 1)) then
                    cell_ind = -1
                else
                    cell_ind = 1
                end if
                do j = 1, n
                 do i = 1, m
                    if ((bc_imn(i, j) .and. current_condition) .eq. BC_FIX_DENSITY) then
                        call fix_density("imin", i * cell_ind, bc_imn_fix_values(i, j, 1))
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_FIX_X_SPEED) then
                        call fix_x_speed("imin", i * cell_ind, bc_imn_fix_values(i, j, 2))
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_FIX_Y_SPEED) then
                        call fix_y_speed("imin", i * cell_ind, bc_imn_fix_values(i, j, 3))
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_FIX_Z_SPEED) then
                        call fix_z_speed("imin", i * cell_ind, bc_imn_fix_values(i, j, 4))
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_FIX_PRESSURE) then
                        call fix_pressure("imin", i * cell_ind, bc_imn_fix_values(i, j, 5))
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_COPY_DENSITY) then
                        call copy_density("imin", i * cell_ind)
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_COPY_X_SPEED) then
                        call copy_x_speed("imin", i * cell_ind)
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_COPY_Y_SPEED) then
                        call copy_y_speed("imin", i * cell_ind)
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_COPY_Z_SPEED) then
                        call copy_z_speed("imin", i * cell_ind)
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_COPY_PRESSURE) then
                        call copy_pressure("imin", i * cell_ind)
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_FLOW_TANGENCY) then
                        call flow_tangency("imin", i * cell_ind)
                        imin_flow_tangency_flag=1
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_NO_SLIP) then
                        call no_slip("imin", i * cell_ind)
                    else if ((bc_imn(i, j) .and. current_condition) .eq. BC_PERIODIC) then
                        call periodic("imin", i * cell_ind)
                    end if

                    include "turbulence_models/include/bc/apply_boundary_condition_imin.inc"

                    if ((bc_imn(i, j) .and. current_condition) .ne. BC_INTERFACE) then
                      if (imin_flow_tangency_flag==0)then
                        !call extra_ghost_cells("imin")
                      end if
                    end if
                 end do
                end do
                
                n_shape = shape(bc_imx)
                m = n_shape(1)
                n = n_shape(2)
                if ((m == 1) .and. (n == 1)) then
                    cell_ind = -1
                else
                    cell_ind = 1
                end if
                do j = 1, n
                 do i = 1, m
                    if ((bc_imx(i, j) .and. current_condition) .eq. BC_FIX_DENSITY) then
                        call fix_density("imax", i * cell_ind, bc_imx_fix_values(i, j, 1))
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_FIX_X_SPEED) then
                        call fix_x_speed("imax", i * cell_ind, bc_imx_fix_values(i, j, 2))
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_FIX_Y_SPEED) then
                        call fix_y_speed("imax", i * cell_ind, bc_imx_fix_values(i, j, 3))
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_FIX_Z_SPEED) then
                        call fix_z_speed("imax", i * cell_ind, bc_imx_fix_values(i, j, 4))
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_FIX_PRESSURE) then
                        call fix_pressure("imax", i * cell_ind, bc_imx_fix_values(i, j, 5))
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_COPY_DENSITY) then
                        call copy_density("imax", i * cell_ind)
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_COPY_X_SPEED) then
                        call copy_x_speed("imax", i * cell_ind)
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_COPY_Y_SPEED) then
                        call copy_y_speed("imax", i * cell_ind)
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_COPY_Z_SPEED) then
                        call copy_z_speed("imax", i * cell_ind)
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_COPY_PRESSURE) then
                        call copy_pressure("imax", i * cell_ind)
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_FLOW_TANGENCY) then
                        call flow_tangency("imax", i * cell_ind)
                        imax_flow_tangency_flag=1
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_NO_SLIP) then
                        call no_slip("imax", i * cell_ind)
                    else if ((bc_imx(i, j) .and. current_condition) .eq. BC_PERIODIC) then
                        call periodic("imax", i * cell_ind)
                    end if
                    include "turbulence_models/include/bc/apply_boundary_condition_imax.inc"

                    if ((bc_imx(i, j) .and. current_condition) .ne. BC_INTERFACE) then
                      if (imax_flow_tangency_flag==0)then
                        !call extra_ghost_cells("imax")
                      end if
                    end if
                 end do
                end do
                
                n_shape = shape(bc_jmn)
                m = n_shape(1)
                n = n_shape(2)
                if ((m == 1) .and. (n == 1)) then
                    cell_ind = -1
                else
                    cell_ind = 1
                end if
                do j = 1, n
                 do i = 1, m
                    if ((bc_jmn(i, j) .and. current_condition) .eq. BC_FIX_DENSITY) then
                        call fix_density("jmin", i * cell_ind, bc_jmn_fix_values(i, j, 1))
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_FIX_X_SPEED) then
                        call fix_x_speed("jmin", i * cell_ind, bc_jmn_fix_values(i, j, 2))
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_FIX_Y_SPEED) then
                        call fix_y_speed("jmin", i * cell_ind, bc_jmn_fix_values(i, j, 3))
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_FIX_Z_SPEED) then
                        call fix_z_speed("jmin", i * cell_ind, bc_jmn_fix_values(i, j, 4))
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_FIX_PRESSURE) then
                        call fix_pressure("jmin", i * cell_ind, bc_jmn_fix_values(i, j, 5))
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_COPY_DENSITY) then
                        call copy_density("jmin", i * cell_ind)
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_COPY_X_SPEED) then
                        call copy_x_speed("jmin", i * cell_ind)
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_COPY_Y_SPEED) then
                        call copy_y_speed("jmin", i * cell_ind)
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_COPY_Z_SPEED) then
                        call copy_z_speed("jmin", i * cell_ind)
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_COPY_PRESSURE) then
                        call copy_pressure("jmin", i * cell_ind)
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_FLOW_TANGENCY) then
                        call flow_tangency("jmin", i * cell_ind)
                        jmin_flow_tangency_flag=1
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_NO_SLIP) then
                        call no_slip("jmin", i * cell_ind)
                    else if ((bc_jmn(i, j) .and. current_condition) .eq. BC_PERIODIC) then
                        call periodic("jmin", i * cell_ind)
                    end if
                    include "turbulence_models/include/bc/apply_boundary_condition_jmin.inc"

                    if ((bc_jmn(i, j) .and. current_condition) .ne. BC_INTERFACE) then
                      if (jmin_flow_tangency_flag==0)then
                        !call extra_ghost_cells("jmin")
                      end if
                    end if
                 end do
                end do

                n_shape = shape(bc_jmx)
                m = n_shape(1)
                n = n_shape(2)
                if ((m == 1) .and. (n == 1)) then
                    cell_ind = -1
                else
                    cell_ind = 1
                end if
                do j = 1, n
                 do i = 1, m
                    if ((bc_jmx(i, j) .and. current_condition) .eq. BC_FIX_DENSITY) then
                        call fix_density("jmax", i * cell_ind, bc_jmx_fix_values(i, j, 1))
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_FIX_X_SPEED) then
                        call fix_x_speed("jmax", i * cell_ind, bc_jmx_fix_values(i, j, 2))
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_FIX_Y_SPEED) then
                        call fix_y_speed("jmax", i * cell_ind, bc_jmx_fix_values(i, j, 3))
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_FIX_Z_SPEED) then
                        call fix_z_speed("jmax", i * cell_ind, bc_jmx_fix_values(i, j, 4))
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_FIX_PRESSURE) then
                        call fix_pressure("jmax", i * cell_ind, bc_jmx_fix_values(i, j, 5))
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_COPY_DENSITY) then
                        call copy_density("jmax", i * cell_ind)
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_COPY_X_SPEED) then
                        call copy_x_speed("jmax", i * cell_ind)
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_COPY_Y_SPEED) then
                        call copy_y_speed("jmax", i * cell_ind)
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_COPY_Z_SPEED) then
                        call copy_z_speed("jmax", i * cell_ind)
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_COPY_PRESSURE) then
                        call copy_pressure("jmax", i * cell_ind)
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_FLOW_TANGENCY) then
                        call flow_tangency("jmax", i * cell_ind)
                        jmax_flow_tangency_flag=1
                    elseif ((bc_jmx(i, j) .and. current_condition) .eq. BC_NO_SLIP) then
                        call no_slip("jmax", i * cell_ind)
                    else if ((bc_jmx(i, j) .and. current_condition) .eq. BC_PERIODIC) then
                        call periodic("jmax", i * cell_ind)
                    end if
                    include "turbulence_models/include/bc/apply_boundary_condition_jmax.inc"

                    if ((bc_jmx(i, j) .and. current_condition) .ne. BC_INTERFACE) then
                      if (jmax_flow_tangency_flag==0)then
                        call extra_ghost_cells("jmax")
                      end if
                    end if
                 end do
                end do

                n_shape = shape(bc_kmn)
                m = n_shape(1)
                n = n_shape(2)
                if ((m == 1) .and. (n == 1)) then
                    cell_ind = -1
                else
                    cell_ind = 1
                end if
                do j = 1, n
                 do i = 1, m
                    if ((bc_kmn(i, j) .and. current_condition) .eq. BC_FIX_DENSITY) then
                        call fix_density("kmin", i * cell_ind, bc_kmn_fix_values(i, j, 1))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_FIX_X_SPEED) then
                        call fix_x_speed("kmin", i * cell_ind, bc_kmn_fix_values(i, j, 2))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_FIX_Y_SPEED) then
                        call fix_y_speed("kmin", i * cell_ind, bc_kmn_fix_values(i, j, 3))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_FIX_Z_SPEED) then
                        call fix_z_speed("kmin", i * cell_ind, bc_kmn_fix_values(i, j, 4))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_FIX_PRESSURE) then
                        call fix_pressure("kmin", i * cell_ind, bc_kmn_fix_values(i, j, 5))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_COPY_DENSITY) then
                        call copy_density("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_COPY_X_SPEED) then
                        call copy_x_speed("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_COPY_Y_SPEED) then
                        call copy_y_speed("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_COPY_Z_SPEED) then
                        call copy_z_speed("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_COPY_PRESSURE) then
                        call copy_pressure("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_FLOW_TANGENCY) then
                        call flow_tangency("kmin", i * cell_ind)
                        kmin_flow_tangency_flag=1
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_NO_SLIP) then
                        call no_slip("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. BC_PERIODIC) then
                        call periodic("kmin", i * cell_ind)
                    end if
                   include "turbulence_models/include/bc/apply_boundary_condition_kmin.inc"

                    if ((bc_kmn(i, j) .and. current_condition) .ne. BC_INTERFACE) then
                      if (kmin_flow_tangency_flag==0)then
                        !call extra_ghost_cells("kmin")
                      end if
                    end if
                 end do
                end do

                n_shape = shape(bc_kmx)
                m = n_shape(1)
                n = n_shape(2)
                if ((m == 1) .and. (n == 1)) then
                    cell_ind = -1
                else
                    cell_ind = 1
                end if
                do j = 1, n
                 do i = 1, m
                    if ((bc_kmx(i, j) .and. current_condition) .eq. BC_FIX_DENSITY) then
                        call fix_density("kmax", i * cell_ind, bc_kmx_fix_values(i, j, 1))
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_FIX_X_SPEED) then
                        call fix_x_speed("kmax", i * cell_ind, bc_kmx_fix_values(i, j, 2))
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_FIX_Y_SPEED) then
                        call fix_y_speed("kmax", i * cell_ind, bc_kmx_fix_values(i, j, 3))
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_FIX_Z_SPEED) then
                        call fix_z_speed("kmax", i * cell_ind, bc_kmx_fix_values(i, j, 4))
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_FIX_PRESSURE) then
                        call fix_pressure("kmax", i * cell_ind, bc_kmx_fix_values(i, j, 5))
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_COPY_DENSITY) then
                        call copy_density("kmax", i * cell_ind)
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_COPY_X_SPEED) then
                        call copy_x_speed("kmax", i * cell_ind)
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_COPY_Y_SPEED) then
                        call copy_y_speed("kmax", i * cell_ind)
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_COPY_Z_SPEED) then
                        call copy_z_speed("kmax", i * cell_ind)
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_COPY_PRESSURE) then
                        call copy_pressure("kmax", i * cell_ind)
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_FLOW_TANGENCY) then
                        call flow_tangency("kmax", i * cell_ind)
                        kmax_flow_tangency_flag=1
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_NO_SLIP) then
                        call no_slip("kmax", i * cell_ind)
                    else if ((bc_kmx(i, j) .and. current_condition) .eq. BC_PERIODIC) then
                        call periodic("kmax", i * cell_ind)
                    end if
                    include "turbulence_models/include/bc/apply_boundary_condition_kmax.inc"

                    if ((bc_kmx(i, j) .and. current_condition) .ne. BC_INTERFACE) then
                      if (kmax_flow_tangency_flag==0)then
                        !call extra_ghost_cells("kmax")
                      end if
                    end if
                 end do
                end do

                if (current_condition > 2**13) then !TODO change it to 17 for turbulence integeraitokn
                    current_condition = 0
                end if
                current_condition = current_condition * 2

            end do

        end subroutine apply_boundary_conditions

        subroutine set_wall_bc_at_faces

            implicit none

            integer :: i, j, k
            real :: ul_temp, vl_temp, wl_temp, ur_temp, vr_temp, &
                wr_temp, nx, ny, nz, dot_prod_l, dot_prod_r

            call dmsg(1, 'boundary_conditions', 'set_wall_bc_at_faces')

            !TODO: Assumes common type boundary condition

            if ((bc_imn(1,1) .and. BC_FLOW_TANGENCY) .eq. BC_FLOW_TANGENCY) then
                do k = 1, kmx-1
                 do j = 1, jmx-1
                    ul_temp = x_x_speed_left(1, j, k)
                    ur_temp = x_x_speed_right(1, j, k)
                    vl_temp = x_y_speed_left(1, j, k)
                    vr_temp = x_y_speed_right(1, j, k)
                    wl_temp = x_z_speed_left(1, j, k)
                    wr_temp = x_z_speed_right(1, j, k)
                    nx = xnx(1, j, k)
                    ny = xny(1, j, k)
                    nz = xnz(1, j, k)
                    dot_prod_l = (ul_temp * nx) + (vl_temp * ny) + (wl_temp * nz)
                    dot_prod_r = (ur_temp * nx) + (vr_temp * ny) + (wr_temp * nz)

                    x_x_speed_left(1, j, k) = ul_temp - (dot_prod_l * nx)
                    x_y_speed_left(1, j, k) = vl_temp - (dot_prod_l * ny)
                    x_z_speed_left(1, j, k) = wl_temp - (dot_prod_l * nz)
                    x_x_speed_right(1, j, k) = ur_temp - (dot_prod_r * nx)
                    x_y_speed_right(1, j, k) = vr_temp - (dot_prod_r * ny)
                    x_z_speed_right(1, j, k) = wr_temp - (dot_prod_r * nz)
                 end do
                end do
            else if ((bc_imx(1,1) .and. BC_FLOW_TANGENCY) .eq. BC_FLOW_TANGENCY) then
                do k = 1, kmx-1
                 do j = 1, jmx-1
                    ul_temp = x_x_speed_left(imx, j, k)
                    ur_temp = x_x_speed_right(imx, j, k)
                    vl_temp = x_y_speed_left(imx, j, k)
                    vr_temp = x_y_speed_right(imx, j, k)
                    wl_temp = x_z_speed_left(imx, j, k)
                    wr_temp = x_z_speed_right(imx, j, k)
                    nx = xnx(imx, j, k)
                    ny = xny(imx, j, k)
                    nz = xnz(imx, j, k)
                    dot_prod_l = (ul_temp * nx) + (vl_temp * ny) + (wl_temp * nz)
                    dot_prod_r = (ur_temp * nx) + (vr_temp * ny) + (wr_temp * nz)

                    x_x_speed_left(imx, j, k) = ul_temp - (dot_prod_l * nx)
                    x_y_speed_left(imx, j, k) = vl_temp - (dot_prod_l * ny)
                    x_z_speed_left(imx, j, k) = wl_temp - (dot_prod_l * nz)
                    x_x_speed_right(imx, j, k) = ur_temp - (dot_prod_r * nx)
                    x_y_speed_right(imx, j, k) = vr_temp - (dot_prod_r * ny)
                    x_z_speed_right(imx, j, k) = wr_temp - (dot_prod_r * nz)
                 end do
                end do
            else if ((bc_jmn(1,1) .and. BC_FLOW_TANGENCY) .eq. BC_FLOW_TANGENCY) then
                do k = 1, kmx-1
                 do i = 1, imx-1
                    ul_temp = y_x_speed_left(i, 1, k)
                    ur_temp = y_x_speed_right(i, 1, k)
                    vl_temp = y_y_speed_left(i, 1, k)
                    vr_temp = y_y_speed_right(i, 1, k)
                    wl_temp = y_z_speed_left(i, 1, k)
                    wr_temp = y_z_speed_right(i, 1, k)
                    nx = ynx(i, 1, k)
                    ny = yny(i, 1, k)
                    nz = ynz(i, 1, k)
                    dot_prod_l = (ul_temp * nx) + (vl_temp * ny) + (wl_temp * nz)
                    dot_prod_r = (ur_temp * nx) + (vr_temp * ny) + (wr_temp * nz)

                    y_x_speed_left(i, 1, k) = ul_temp - (dot_prod_l * nx)
                    y_y_speed_left(i, 1, k) = vl_temp - (dot_prod_l * ny)
                    y_z_speed_left(i, 1, k) = wl_temp - (dot_prod_l * nz)
                    y_x_speed_right(i, 1, k) = ur_temp - (dot_prod_r * nx)
                    y_y_speed_right(i, 1, k) = vr_temp - (dot_prod_r * ny)
                    y_z_speed_right(i, 1, k) = wr_temp - (dot_prod_r * nz)
                 end do
                end do
            else if ((bc_jmx(1,1) .and. BC_FLOW_TANGENCY) .eq. BC_FLOW_TANGENCY) then
                do k = 1, kmx-1
                 do i = 1, imx-1
                    ul_temp = y_x_speed_left(i, jmx, k)
                    ur_temp = y_x_speed_right(i, jmx, k)
                    vl_temp = y_y_speed_left(i, jmx, k)
                    vr_temp = y_y_speed_right(i, jmx, k)
                    wl_temp = y_z_speed_left(i, jmx, k)
                    wr_temp = y_z_speed_right(i, jmx, k)
                    nx = ynx(i, jmx, k)
                    ny = yny(i, jmx, k)
                    nz = ynz(i, jmx, k)
                    dot_prod_l = (ul_temp * nx) + (vl_temp * ny) + (wl_temp * nz)
                    dot_prod_r = (ur_temp * nx) + (vr_temp * ny) + (wr_temp * nz)

                    y_x_speed_left(i, jmx, k) = ul_temp - (dot_prod_l * nx)
                    y_y_speed_left(i, jmx, k) = vl_temp - (dot_prod_l * ny)
                    y_z_speed_left(i, jmx, k) = wl_temp - (dot_prod_l * nz)
                    y_x_speed_right(i, jmx, k) = ur_temp - (dot_prod_r * nx)
                    y_y_speed_right(i, jmx, k) = vr_temp - (dot_prod_r * ny)
                    y_z_speed_right(i, jmx, k) = wr_temp - (dot_prod_r * nz)
                 end do
                end do
            else if ((bc_kmn(1,1) .and. BC_FLOW_TANGENCY) .eq. BC_FLOW_TANGENCY) then
                do j = 1, jmx-1
                 do i = 1, imx-1
                    ul_temp = z_x_speed_left(i, j, 1)
                    ur_temp = z_x_speed_right(i, j, 1)
                    vl_temp = z_y_speed_left(i, j, 1)
                    vr_temp = z_y_speed_right(i, j, 1)
                    wl_temp = z_z_speed_left(i, j, 1)
                    wr_temp = z_z_speed_right(i, j, 1)
                    nx = znx(i, j, 1)
                    ny = zny(i, j, 1)
                    nz = znz(i, j, 1)
                    dot_prod_l = (ul_temp * nx) + (vl_temp * ny) + (wl_temp * nz)
                    dot_prod_r = (ur_temp * nx) + (vr_temp * ny) + (wr_temp * nz)

                    z_x_speed_left(i, j, 1) = ul_temp - (dot_prod_l * nx)
                    z_y_speed_left(i, j, 1) = vl_temp - (dot_prod_l * ny)
                    z_z_speed_left(i, j, 1) = wl_temp - (dot_prod_l * nz)
                    z_x_speed_right(i, j, 1) = ur_temp - (dot_prod_r * nx)
                    z_y_speed_right(i, j, 1) = vr_temp - (dot_prod_r * ny)
                    z_z_speed_right(i, j, 1) = wr_temp - (dot_prod_r * nz)
                 end do
                end do
            else if ((bc_kmx(1,1) .and. BC_FLOW_TANGENCY) .eq. BC_FLOW_TANGENCY) then
                do j = 1, jmx-1
                 do i = 1, imx-1
                    ul_temp = z_x_speed_left(i, j, kmx)
                    ur_temp = z_x_speed_right(i, j, kmx)
                    vl_temp = z_y_speed_left(i, j, kmx)
                    vr_temp = z_y_speed_right(i, j, kmx)
                    wl_temp = z_z_speed_left(i, j, kmx)
                    wr_temp = z_z_speed_right(i, j, kmx)
                    nx = znx(i, j, kmx)
                    ny = zny(i, j, kmx)
                    nz = znz(i, j, kmx)
                    dot_prod_l = (ul_temp * nx) + (vl_temp * ny) + (wl_temp * nz)
                    dot_prod_r = (ur_temp * nx) + (vr_temp * ny) + (wr_temp * nz)

                    z_x_speed_left(i, j, kmx) = ul_temp - (dot_prod_l * nx)
                    z_y_speed_left(i, j, kmx) = vl_temp - (dot_prod_l * ny)
                    z_z_speed_left(i, j, kmx) = wl_temp - (dot_prod_l * nz)
                    z_x_speed_right(i, j, kmx) = ur_temp - (dot_prod_r * nx)
                    z_y_speed_right(i, j, kmx) = vr_temp - (dot_prod_r * ny)
                    z_z_speed_right(i, j, kmx) = wr_temp - (dot_prod_r * nz)
                 end do
                end do
            end if

            if ((bc_imn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                x_x_speed_left(1, :, :) = 0.
                x_y_speed_left(1, :, :) = 0.
                x_z_speed_left(1, :, :) = 0.
                x_x_speed_right(1, :, :) = 0.
                x_y_speed_right(1, :, :) = 0.
                x_z_speed_right(1, :, :) = 0.
            else if ((bc_imx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                x_x_speed_left(imx, :, :) = 0.
                x_y_speed_left(imx, :, :) = 0.
                x_z_speed_left(imx, :, :) = 0.
                x_x_speed_right(imx, :, :) = 0.
                x_y_speed_right(imx, :, :) = 0.
                x_z_speed_right(imx, :, :) = 0.
            else if ((bc_jmn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                y_x_speed_left(:, 1, :) = 0.
                y_y_speed_left(:, 1, :) = 0.
                y_z_speed_left(:, 1, :) = 0.
                y_x_speed_right(:, 1, :) = 0.
                y_y_speed_right(:, 1, :) = 0.
                y_z_speed_right(:, 1, :) = 0.
            else if ((bc_jmx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                y_x_speed_left(:, jmx, :) = 0.
                y_y_speed_left(:, jmx, :) = 0.
                y_z_speed_left(:, jmx, :) = 0.
                y_x_speed_right(:, jmx, :) = 0.
                y_y_speed_right(:, jmx, :) = 0.
                y_z_speed_right(:, jmx, :) = 0.
            else if ((bc_kmn(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                z_x_speed_left(:, :, 1) = 0.
                z_y_speed_left(:, :, 1) = 0.
                z_z_speed_left(:, :, 1) = 0.
                z_x_speed_right(:, :, 1) = 0.
                z_y_speed_right(:, :, 1) = 0.
                z_z_speed_right(:, :, 1) = 0.
            else if ((bc_kmx(1,1) .and. BC_NO_SLIP) .eq. BC_NO_SLIP) then
                z_x_speed_left(:, :, kmx) = 0.
                z_y_speed_left(:, :, kmx) = 0.
                z_z_speed_left(:, :, kmx) = 0.
                z_x_speed_right(:, :, kmx) = 0.
                z_y_speed_right(:, :, kmx) = 0.
                z_z_speed_right(:, :, kmx) = 0.
            end if

            include "turbulence_models/include/bc/set_wall_bc_at_faces.inc"

        end subroutine set_wall_bc_at_faces

        !---------------------------------------------------------------
        ! "Fix" type conditions
        !
        ! These subroutines are named "fix_<parameter>" where 
        ! <parameter> can be density, x_speed, y_speed, pressure. Each 
        ! of the subroutines accept an argument face as described above. 
        ! The subroutine will set the value of the ghost cell at the 
        ! face specified. 
        !
        ! The subroutine also accepts an optional val. If val is 
        ! specified, the parameter is fixed to that value, else the 
        ! infinity values are used. 
        !---------------------------------------------------------------

        subroutine fix_density(face, cell_ind, optional_fixed_density)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_density
            real :: fixed_density

            if (present(optional_fixed_density)) then
                fixed_density = optional_fixed_density
            else
                fixed_density = density_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    density(0, 1:jmx-1, 1:kmx-1) = fixed_density
                    density(-1, 1:jmx-1, 1:kmx-1) = fixed_density
                    density(-2, 1:jmx-1, 1:kmx-1) = fixed_density
                else if (face == "imax") then
                    density(imx, 1:jmx-1, 1:kmx-1) = fixed_density
                    density(imx+1, 1:jmx-1, 1:kmx-1) = fixed_density
                    density(imx+2, 1:jmx-1, 1:kmx-1) = fixed_density
                else if (face == "jmin") then
                    density(1:imx-1, 0, 1:kmx-1) = fixed_density
                    density(1:imx-1, -1, 1:kmx-1) = fixed_density
                    density(1:imx-1, -2, 1:kmx-1) = fixed_density
                else if (face == "jmax") then
                    density(1:imx-1, jmx, 1:kmx-1) = fixed_density
                    density(1:imx-1, jmx+1, 1:kmx-1) = fixed_density
                    density(1:imx-1, jmx+2, 1:kmx-1) = fixed_density
                else if (face == "kmin") then
                    density(1:imx-1, 1:jmx-1, 0) = fixed_density
                    density(1:imx-1, 1:jmx-1, -1) = fixed_density
                    density(1:imx-1, 1:jmx-1, -2) = fixed_density
                else if (face == "kmax") then
                    density(1:imx-1, 1:jmx-1, kmx) = fixed_density
                    density(1:imx-1, 1:jmx-1, kmx+1) = fixed_density
                    density(1:imx-1, 1:jmx-1, kmx+2) = fixed_density
                end if
       !    else
       !        if (face == "imin") then
       !            density(0, :, cell_ind) = fixed_density
       !        else if (face == "imax") then
       !            density(imx, :, cell_ind) = fixed_density
       !        else if (face == "jmin") then
       !            density(cell_ind, :, 0) = fixed_density
       !        else if (face == "jmax") then
       !            density(cell_ind, :, jmx) = fixed_density
       !        end if
            end if

        end subroutine fix_density

        subroutine fix_x_speed(face, cell_ind, optional_fixed_x_speed)
            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_x_speed
            real :: fixed_x_speed

            if (present(optional_fixed_x_speed)) then
                fixed_x_speed = optional_fixed_x_speed
            else
                fixed_x_speed = x_speed_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    x_speed(0, 1:jmx-1, 1:kmx-1) = fixed_x_speed
                else if (face == "imax") then
                    x_speed(imx, 1:jmx-1, 1:kmx-1) = fixed_x_speed
                else if (face == "jmin") then
                    x_speed(1:imx-1, 0, 1:kmx-1) = fixed_x_speed
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx, 1:kmx-1) = fixed_x_speed
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, 0) = fixed_x_speed
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx) = fixed_x_speed
                end if
       !    else
       !        if (face == "imin") then
       !            x_speed(0, cell_ind) = fixed_x_speed
       !        else if (face == "imax") then
       !            x_speed(imx, cell_ind) = fixed_x_speed
       !        else if (face == "jmin") then
       !            x_speed(cell_ind, 0) = fixed_x_speed
       !        else if (face == "jmax") then
       !            x_speed(cell_ind, jmx) = fixed_x_speed
       !        end if
            end if

        end subroutine fix_x_speed

        subroutine fix_y_speed(face, cell_ind, optional_fixed_y_speed)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_y_speed
            real :: fixed_y_speed

            if (present(optional_fixed_y_speed)) then
                fixed_y_speed = optional_fixed_y_speed
            else
                fixed_y_speed = y_speed_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    y_speed(0, 1:jmx-1, 1:kmx-1) = fixed_y_speed
                else if (face == "imax") then
                    y_speed(imx, 1:jmx-1, 1:kmx-1) = fixed_y_speed
                else if (face == "jmin") then
                    y_speed(1:imx-1, 0, 1:kmx-1) = fixed_y_speed
                else if (face == "jmax") then
                    y_speed(1:imx-1, jmx, 1:kmx-1) = fixed_y_speed
                else if (face == "kmin") then
                    y_speed(1:imx-1, 1:jmx-1, 0) = fixed_y_speed
                else if (face == "kmax") then
                    y_speed(1:imx-1, 1:jmx-1, kmx) = fixed_y_speed
                end if
        !   else
        !       if (face == "imin") then
        !           y_speed(0, cell_ind) = fixed_y_speed
        !       else if (face == "imax") then
        !           y_speed(imx, cell_ind) = fixed_y_speed
        !       else if (face == "jmin") then
        !           y_speed(cell_ind, 0) = fixed_y_speed
        !       else if (face == "jmax") then
        !           y_speed(cell_ind, jmx) = fixed_y_speed
        !       end if
            end if

        end subroutine fix_y_speed

        subroutine fix_z_speed(face, cell_ind, optional_fixed_z_speed)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_z_speed
            real :: fixed_z_speed

            if (present(optional_fixed_z_speed)) then
                fixed_z_speed = optional_fixed_z_speed
            else
                fixed_z_speed = z_speed_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    z_speed(0, 1:jmx-1, 1:kmx-1) = fixed_z_speed
                else if (face == "imax") then
                    z_speed(imx, 1:jmx-1, 1:kmx-1) = fixed_z_speed
                else if (face == "jmin") then
                    z_speed(1:imx-1, 0, 1:kmx-1) = fixed_z_speed
                else if (face == "jmax") then
                    z_speed(1:imx-1, jmx, 1:kmx-1) = fixed_z_speed
                else if (face == "kmin") then
                    z_speed(1:imx-1, 1:jmx-1, 0) = fixed_z_speed
                else if (face == "kmax") then
                    z_speed(1:imx-1, 1:jmx-1, kmx) = fixed_z_speed
                end if
        !   else
        !       if (face == "imin") then
        !           y_speed(0, cell_ind) = fixed_y_speed
        !       else if (face == "imax") then
        !           y_speed(imx, cell_ind) = fixed_y_speed
        !       else if (face == "jmin") then
        !           y_speed(cell_ind, 0) = fixed_y_speed
        !       else if (face == "jmax") then
        !           y_speed(cell_ind, jmx) = fixed_y_speed
        !       end if
            end if

        end subroutine fix_z_speed

        subroutine fix_pressure(face, cell_ind, optional_fixed_pressure)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind
            real, intent(in), optional :: optional_fixed_pressure
            real :: fixed_pressure

            if (present(optional_fixed_pressure)) then
                fixed_pressure = optional_fixed_pressure
            else
                fixed_pressure = pressure_inf
            end if
            
            if (cell_ind == -1) then
                if (face == "imin") then
                    pressure( 0, 1:jmx-1, 1:kmx-1)    = fixed_pressure
                    pressure(-1, 1:jmx-1, 1:kmx-1)    = fixed_pressure
                    pressure(-2, 1:jmx-1, 1:kmx-1)    = fixed_pressure
                else if (face == "imax") then
                    pressure(imx  , 1:jmx-1, 1:kmx-1) = fixed_pressure
                    pressure(imx+1, 1:jmx-1, 1:kmx-1) = fixed_pressure
                    pressure(imx+2, 1:jmx-1, 1:kmx-1) = fixed_pressure
                else if (face == "jmin") then
                    pressure(1:imx-1,  0, 1:kmx-1)    = fixed_pressure
                    pressure(1:imx-1, -1, 1:kmx-1)    = fixed_pressure
                    pressure(1:imx-1, -2, 1:kmx-1)    = fixed_pressure
                else if (face == "jmax") then
                    pressure(1:imx-1, jmx  , 1:kmx-1) = fixed_pressure
                    pressure(1:imx-1, jmx+1, 1:kmx-1) = fixed_pressure
                    pressure(1:imx-1, jmx+2, 1:kmx-1) = fixed_pressure
                else if (face == "kmin") then
                    pressure(1:imx-1, 1:jmx-1,  0)    = fixed_pressure
                    pressure(1:imx-1, 1:jmx-1, -1)    = fixed_pressure
                    pressure(1:imx-1, 1:jmx-1, -2)    = fixed_pressure
                else if (face == "kmax") then
                    pressure(1:imx-1, 1:jmx-1, kmx  ) = fixed_pressure
                    pressure(1:imx-1, 1:jmx-1, kmx+1) = fixed_pressure
                    pressure(1:imx-1, 1:jmx-1, kmx+2) = fixed_pressure
                end if
        !   else
        !       if (face == "imin") then
        !           pressure(0, cell_ind) = fixed_pressure
        !       else if (face == "imax") then
        !           pressure(imx, cell_ind) = fixed_pressure
        !       else if (face == "jmin") then
        !           pressure(cell_ind, 0) = fixed_pressure
        !       else if (face == "jmax") then
        !           pressure(cell_ind, jmx) = fixed_pressure
        !       end if
            end if

        end subroutine fix_pressure

        !---------------------------------------------------------------
        ! "Copy" type conditions
        !
        ! These subroutines are named "copy_<parameter>" where 
        ! <parameter> can be density, x_speed, y_speed, pressure. Each 
        ! of the subroutines accept an argument face as described above.
        ! The subroutine will set the value of the ghost cell at the 
        ! face specified. 
        !---------------------------------------------------------------

        subroutine copy_density(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    density( 0, 1:jmx-1, 1:kmx-1)    = density(1, 1:jmx-1, 1:kmx-1)
                    density(-1, 1:jmx-1, 1:kmx-1)    = density(2, 1:jmx-1, 1:kmx-1)
                    density(-2, 1:jmx-1, 1:kmx-1)    = density(3, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    density(imx  , 1:jmx-1, 1:kmx-1) = density(imx-1, 1:jmx-1, 1:kmx-1)
                    density(imx+1, 1:jmx-1, 1:kmx-1) = density(imx-2, 1:jmx-1, 1:kmx-1)
                    density(imx+2, 1:jmx-1, 1:kmx-1) = density(imx-3, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    density(1:imx-1,  0, 1:kmx-1)    = density(1:imx-1, 1, 1:kmx-1)
                    density(1:imx-1, -1, 1:kmx-1)    = density(1:imx-1, 2, 1:kmx-1)
                    density(1:imx-1, -2, 1:kmx-1)    = density(1:imx-1, 3, 1:kmx-1)
                else if (face == "jmax") then
                    density(1:imx-1, jmx  , 1:kmx-1) = density(1:imx-1, jmx-1, 1:kmx-1)
                    density(1:imx-1, jmx+1, 1:kmx-1) = density(1:imx-1, jmx-2, 1:kmx-1)
                    density(1:imx-1, jmx+2, 1:kmx-1) = density(1:imx-1, jmx-3, 1:kmx-1)
                else if (face == "kmin") then
                    density(1:imx-1, 1:jmx-1,  0)    = density(1:imx-1, 1:jmx-1, 1)
                    density(1:imx-1, 1:jmx-1, -1)    = density(1:imx-1, 1:jmx-1, 2)
                    density(1:imx-1, 1:jmx-1, -2)    = density(1:imx-1, 1:jmx-1, 3)
                else if (face == "kmax") then
                    density(1:imx-1, 1:jmx-1, kmx  ) = density(1:imx-1, 1:jmx-1, kmx-1)
                    density(1:imx-1, 1:jmx-1, kmx+1) = density(1:imx-1, 1:jmx-1, kmx-2)
                    density(1:imx-1, 1:jmx-1, kmx+2) = density(1:imx-1, 1:jmx-1, kmx-3)
                end if
    !       else
    !           if (face == "imin") then
    !               density(0, cell_ind) = density(1, cell_ind)
    !           else if (face == "imax") then
    !               density(imx, cell_ind) = density(imx-1, cell_ind)
    !           else if (face == "jmin") then
    !               density(cell_ind, 0) = density(cell_ind, 1)
    !           else if (face == "jmax") then
    !               density(cell_ind, jmx) = density(cell_ind, jmx-1)
    !           end if
            end if

        end subroutine copy_density

        subroutine copy_x_speed(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    x_speed( 0, 1:jmx-1, 1:kmx-1)    = x_speed( 1, 1:jmx-1, 1:kmx-1)
                    x_speed(-1, 1:jmx-1, 1:kmx-1)    = x_speed( 0, 1:jmx-1, 1:kmx-1)
                    x_speed(-2, 1:jmx-1, 1:kmx-1)    = x_speed(-1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    x_speed(imx  , 1:jmx-1, 1:kmx-1) = x_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    x_speed(imx+1, 1:jmx-1, 1:kmx-1) = x_speed(imx  , 1:jmx-1, 1:kmx-1)
                    x_speed(imx+2, 1:jmx-1, 1:kmx-1) = x_speed(imx+1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    x_speed(1:imx-1, 0, 1:kmx-1)     = x_speed(1:imx-1, 1, 1:kmx-1)
                    x_speed(1:imx-1,-1, 1:kmx-1)     = x_speed(1:imx-1, 0, 1:kmx-1)
                    x_speed(1:imx-1,-2, 1:kmx-1)     = x_speed(1:imx-1,-1, 1:kmx-1)
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx  , 1:kmx-1) = x_speed(1:imx-1, jmx-1, 1:kmx-1)
                    x_speed(1:imx-1, jmx+1, 1:kmx-1) = x_speed(1:imx-1, jmx  , 1:kmx-1)
                    x_speed(1:imx-1, jmx+2, 1:kmx-1) = x_speed(1:imx-1, jmx+1, 1:kmx-1)
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, 0)     = x_speed(1:imx-1, 1:jmx-1,  1)
                    x_speed(1:imx-1, 1:jmx-1,-1)     = x_speed(1:imx-1, 1:jmx-1,  0)
                    x_speed(1:imx-1, 1:jmx-1,-2)     = x_speed(1:imx-1, 1:jmx-1, -1)
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx  ) = x_speed(1:imx-1, 1:jmx-1, kmx-1)
                    x_speed(1:imx-1, 1:jmx-1, kmx+1) = x_speed(1:imx-1, 1:jmx-1, kmx  )
                    x_speed(1:imx-1, 1:jmx-1, kmx+2) = x_speed(1:imx-1, 1:jmx-1, kmx+1)
                end if
      !     else
      !         if (face == "imin") then
      !             x_speed(0, cell_ind) = x_speed(1, cell_ind)
      !         else if (face == "imax") then
      !             x_speed(imx, cell_ind) = x_speed(imx-1, cell_ind)
      !         else if (face == "jmin") then
      !             x_speed(cell_ind, 0) = x_speed(cell_ind, 1)
      !         else if (face == "jmax") then
      !             x_speed(cell_ind, jmx) = x_speed(cell_ind, jmx-1)
      !         end if
            end if

        end subroutine copy_x_speed

        subroutine copy_y_speed(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    y_speed( 0, 1:jmx-1, 1:kmx-1)    = y_speed( 1, 1:jmx-1, 1:kmx-1)
                    y_speed(-1, 1:jmx-1, 1:kmx-1)    = y_speed( 0, 1:jmx-1, 1:kmx-1)
                    y_speed(-2, 1:jmx-1, 1:kmx-1)    = y_speed(-1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    y_speed(imx  , 1:jmx-1, 1:kmx-1) = y_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    y_speed(imx+1, 1:jmx-1, 1:kmx-1) = y_speed(imx  , 1:jmx-1  , 1:kmx-1)
                    y_speed(imx+2, 1:jmx-1, 1:kmx-1) = y_speed(imx+1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    y_speed(1:imx-1,  0, 1:kmx-1)    = y_speed(1:imx-1, 1, 1:kmx-1)
                    y_speed(1:imx-1, -1, 1:kmx-1)    = y_speed(1:imx-1, 0, 1:kmx-1)
                    y_speed(1:imx-1, -2, 1:kmx-1)    = y_speed(1:imx-1,-1, 1:kmx-1)
                else if (face == "jmax") then
                    y_speed(1:imx-1, jmx  , 1:kmx-1) = y_speed(1:imx-1, jmx-1, 1:kmx-1)
                    y_speed(1:imx-1, jmx+1, 1:kmx-1) = y_speed(1:imx-1, jmx  , 1:kmx-1)
                    y_speed(1:imx-1, jmx+2, 1:kmx-1) = y_speed(1:imx-1, jmx+1, 1:kmx-1)
                else if (face == "kmin") then
                    y_speed(1:imx-1, 1:jmx-1, 0)     = y_speed(1:imx-1, 1:jmx-1,  1)
                    y_speed(1:imx-1, 1:jmx-1,-1)     = y_speed(1:imx-1, 1:jmx-1,  0)
                    y_speed(1:imx-1, 1:jmx-1,-2)     = y_speed(1:imx-1, 1:jmx-1, -1)
                else if (face == "kmax") then
                    y_speed(1:imx-1, 1:jmx-1, kmx  ) = y_speed(1:imx-1, 1:jmx-1, kmx-1)
                    y_speed(1:imx-1, 1:jmx-1, kmx+1) = y_speed(1:imx-1, 1:jmx-1, kmx  )
                    y_speed(1:imx-1, 1:jmx-1, kmx+2) = y_speed(1:imx-1, 1:jmx-1, kmx+1)
                end if
      !     else
      !         if (face == "imin") then
      !             x_speed(0, cell_ind) = x_speed(1, cell_ind)
      !         else if (face == "imax") then
      !             x_speed(imx, cell_ind) = x_speed(imx-1, cell_ind)
      !         else if (face == "jmin") then
      !             x_speed(cell_ind, 0) = x_speed(cell_ind, 1)
      !         else if (face == "jmax") then
      !             x_speed(cell_ind, jmx) = x_speed(cell_ind, jmx-1)
      !         end if
            end if

        end subroutine copy_y_speed

        subroutine copy_z_speed(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    z_speed( 0, 1:jmx-1, 1:kmx-1)    = z_speed( 1, 1:jmx-1, 1:kmx-1)
                    z_speed(-1, 1:jmx-1, 1:kmx-1)    = z_speed( 0, 1:jmx-1, 1:kmx-1)
                    z_speed(-2, 1:jmx-1, 1:kmx-1)    = z_speed(-1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    z_speed(imx  , 1:jmx-1, 1:kmx-1) = z_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    z_speed(imx+1, 1:jmx-1, 1:kmx-1) = z_speed(imx  , 1:jmx-1, 1:kmx-1)
                    z_speed(imx+2, 1:jmx-1, 1:kmx-1) = z_speed(imx+1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    z_speed(1:imx-1,  0, 1:kmx-1)    = z_speed(1:imx-1,  1, 1:kmx-1)
                    z_speed(1:imx-1, -1, 1:kmx-1)    = z_speed(1:imx-1,  0, 1:kmx-1)
                    z_speed(1:imx-1, -2, 1:kmx-1)    = z_speed(1:imx-1, -1, 1:kmx-1)
                else if (face == "jmax") then
                    z_speed(1:imx-1, jmx  , 1:kmx-1) = z_speed(1:imx-1, jmx-1, 1:kmx-1)
                    z_speed(1:imx-1, jmx+1, 1:kmx-1) = z_speed(1:imx-1, jmx  , 1:kmx-1)
                    z_speed(1:imx-1, jmx+2, 1:kmx-1) = z_speed(1:imx-1, jmx+1, 1:kmx-1)
                else if (face == "kmin") then
                    z_speed(1:imx-1, 1:jmx-1,  0)    = z_speed(1:imx-1, 1:jmx-1,  1)
                    z_speed(1:imx-1, 1:jmx-1, -1)    = z_speed(1:imx-1, 1:jmx-1,  0)
                    z_speed(1:imx-1, 1:jmx-1, -2)    = z_speed(1:imx-1, 1:jmx-1, -1)
                else if (face == "kmax") then
                    z_speed(1:imx-1, 1:jmx-1, kmx  ) = z_speed(1:imx-1, 1:jmx-1, kmx-1)
                    z_speed(1:imx-1, 1:jmx-1, kmx+1) = z_speed(1:imx-1, 1:jmx-1, kmx  )
                    z_speed(1:imx-1, 1:jmx-1, kmx+2) = z_speed(1:imx-1, 1:jmx-1, kmx+1)
                end if
      !     else
      !         if (face == "imin") then
      !             x_speed(0, cell_ind) = x_speed(1, cell_ind)
      !         else if (face == "imax") then
      !             x_speed(imx, cell_ind) = x_speed(imx-1, cell_ind)
      !         else if (face == "jmin") then
      !             x_speed(cell_ind, 0) = x_speed(cell_ind, 1)
      !         else if (face == "jmax") then
      !             x_speed(cell_ind, jmx) = x_speed(cell_ind, jmx-1)
      !         end if
            end if

        end subroutine copy_z_speed

        subroutine copy_pressure(face, cell_ind)

            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    pressure(0, 1:jmx-1, 1:kmx-1) = pressure(1, 1:jmx-1, 1:kmx-1)
                    pressure(-1, 1:jmx-1, 1:kmx-1) = pressure(2, 1:jmx-1, 1:kmx-1)
                    pressure(-2, 1:jmx-1, 1:kmx-1) = pressure(3, 1:jmx-1, 1:kmx-1)
                 else if (face == "imax") then
                    pressure(imx, 1:jmx-1, 1:kmx-1) = pressure(imx-1, 1:jmx-1, 1:kmx-1)
                    pressure(imx+1, 1:jmx-1, 1:kmx-1) = pressure(imx-2, 1:jmx-1, 1:kmx-1)
                    pressure(imx+2, 1:jmx-1, 1:kmx-1) = pressure(imx-3, 1:jmx-1, 1:kmx-1)
                 else if (face == "jmin") then
                    pressure(1:imx-1, 0, 1:kmx-1) = pressure(1:imx-1, 1, 1:kmx-1)
                    pressure(1:imx-1, -1, 1:kmx-1) = pressure(1:imx-1, 2, 1:kmx-1)
                    pressure(1:imx-1, -2, 1:kmx-1) = pressure(1:imx-1, 3, 1:kmx-1)
                else if (face == "jmax") then
                    pressure(1:imx-1, jmx, 1:kmx-1) = pressure(1:imx-1, jmx-1, 1:kmx-1)
                    pressure(1:imx-1, jmx+1, 1:kmx-1) = pressure(1:imx-1, jmx-2, 1:kmx-1)
                    pressure(1:imx-1, jmx+2, 1:kmx-1) = pressure(1:imx-1, jmx-3, 1:kmx-1)
                else if (face == "kmin") then
                    pressure(1:imx-1, 1:jmx-1, 0) = pressure(1:imx-1, 1:jmx-1, 1)
                    pressure(1:imx-1, 1:jmx-1, -1) = pressure(1:imx-1, 1:jmx-1, 2)
                    pressure(1:imx-1, 1:jmx-1, -2) = pressure(1:imx-1, 1:jmx-1, 3)
                else if (face == "kmax") then
                    pressure(1:imx-1, 1:jmx-1, kmx) = pressure(1:imx-1, 1:jmx-1, kmx-1)
                    pressure(1:imx-1, 1:jmx-1, kmx+1) = pressure(1:imx-1, 1:jmx-1, kmx-2)
                    pressure(1:imx-1, 1:jmx-1, kmx+2) = pressure(1:imx-1, 1:jmx-1, kmx-3)
                end if
        !   else
        !       if (face == "imin") then
        !           pressure(0, cell_ind) = pressure(1, cell_ind)
        !       else if (face == "imax") then
        !           pressure(imx, cell_ind) = pressure(imx-1, cell_ind)
        !       else if (face == "jmin") then
        !           pressure(cell_ind, 0) = pressure(cell_ind, 1)
        !       else if (face == "jmax") then
        !           pressure(cell_ind, jmx) = pressure(cell_ind, jmx-1)
        !       end if
            end if

        end subroutine copy_pressure

        !---------------------------------------------------------------
        ! Misc conditions
        !
        ! All other condition are grouped under "Misc". These have been 
        ! described above.
        !---------------------------------------------------------------

        subroutine flow_tangency(face, cell_ind)
            !-----------------------------------------------------------
            ! The flow tangency condition is an enforcement of the slip
            ! wall condition, i.e., the velocity at the wall should be 
            ! parallel to it. This is achieved by setting the velocity 
            ! of the ghost cell to be the mirror image (about the 
            ! boundary) of its neighbouring cell's velocity.
            !-----------------------------------------------------------
        
            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    x_speed(0, 1:jmx-1, 1:kmx-1) = x_speed(1, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1, 1:jmx-1, 1:kmx-1) * xnx(1, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(1, 1:jmx-1, 1:kmx-1) * xny(1, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(1, 1:jmx-1, 1:kmx-1) * xnz(1, 1:jmx-1, 1:kmx-1)) &
                                ) * xnx(1, 1:jmx-1, 1:kmx-1) &
                            )
                    y_speed(0, 1:jmx-1, 1:kmx-1) = y_speed(1, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1, 1:jmx-1, 1:kmx-1) * xnx(1, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(1, 1:jmx-1, 1:kmx-1) * xny(1, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(1, 1:jmx-1, 1:kmx-1) * xnz(1, 1:jmx-1, 1:kmx-1)) &
                                ) * xny(1, 1:jmx-1, 1:kmx-1) &
                            )
                    z_speed(0, 1:jmx-1, 1:kmx-1) = z_speed(1, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1, 1:jmx-1, 1:kmx-1) * xnx(1, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(1, 1:jmx-1, 1:kmx-1) * xny(1, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(1, 1:jmx-1, 1:kmx-1) * xnz(1, 1:jmx-1, 1:kmx-1)) &
                                ) * xnz(1, 1:jmx-1, 1:kmx-1) &
                            )
                else if (face == "imax") then
                    x_speed(imx, 1:jmx-1, 1:kmx-1) = x_speed(imx-1, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-1, 1:jmx-1, 1:kmx-1) * xnx(imx, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-1, 1:jmx-1, 1:kmx-1) * xny(imx, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-1, 1:jmx-1, 1:kmx-1) * xnz(imx, 1:jmx-1, 1:kmx-1)) &
                                ) * xnx(imx, 1:jmx-1, 1:kmx-1) &
                            )
                    y_speed(imx, 1:jmx-1, 1:kmx-1) = y_speed(imx-1, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-1, 1:jmx-1, 1:kmx-1) * xnx(imx, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-1, 1:jmx-1, 1:kmx-1) * xny(imx, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-1, 1:jmx-1, 1:kmx-1) * xnz(imx, 1:jmx-1, 1:kmx-1)) &
                                ) * xny(imx, 1:jmx-1, 1:kmx-1) &
                            )
                    z_speed(imx, 1:jmx-1, 1:kmx-1) = z_speed(imx-1, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-1, 1:jmx-1, 1:kmx-1) * xnx(imx, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-1, 1:jmx-1, 1:kmx-1) * xny(imx, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-1, 1:jmx-1, 1:kmx-1) * xnz(imx, 1:jmx-1, 1:kmx-1)) &
                                ) * xnz(imx, 1:jmx-1, 1:kmx-1) &
                            )
                else if (face == "jmin") then
                    x_speed(1:imx-1, 0, 1:kmx-1) = x_speed(1:imx-1, 1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 1, 1:kmx-1) * ynz(1:imx-1, 1, 1:kmx-1)) &
                                ) * ynx(1:imx-1, 1, 1:kmx-1) &
                            )
                    y_speed(1:imx-1, 0, 1:kmx-1) = y_speed(1:imx-1, 1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 1, 1:kmx-1) * ynz(1:imx-1, 1, 1:kmx-1)) &
                                ) * yny(1:imx-1, 1, 1:kmx-1) &
                            )
                    z_speed(1:imx-1, 0, 1:kmx-1) = z_speed(1:imx-1, 1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 1, 1:kmx-1) * ynz(1:imx-1, 1, 1:kmx-1)) &
                                ) * ynz(1:imx-1, 1, 1:kmx-1) &
                            )
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx, 1:kmx-1) = x_speed(1:imx-1, jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-1, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-1, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-1, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
                                ) * ynx(1:imx-1, jmx, 1:kmx-1) &
                            )
                    y_speed(1:imx-1, jmx, 1:kmx-1) = y_speed(1:imx-1, jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-1, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-1, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-1, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
                                ) * yny(1:imx-1, jmx, 1:kmx-1) &
                            )
                    z_speed(1:imx-1, jmx, 1:kmx-1) = z_speed(1:imx-1, jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-1, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-1, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-1, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
                                ) * ynz(1:imx-1, jmx, 1:kmx-1) &
                            )
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, 0) = x_speed(1:imx-1, 1:jmx-1, 1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
                                ) * znx(1:imx-1, 1:jmx-1, 1) &
                            )
                    y_speed(1:imx-1, 1:jmx-1, 0) = y_speed(1:imx-1, 1:jmx-1, 1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
                                ) * zny(1:imx-1, 1:jmx-1, 1) &
                            )
                    z_speed(1:imx-1, 1:jmx-1, 0) = z_speed(1:imx-1, 1:jmx-1, 1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
                                ) * znz(1:imx-1, 1:jmx-1, 1) &
                            )
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx) = x_speed(1:imx-1, 1:jmx-1, kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-1) * znx(1:imx-1, 1:jmx-1, kmx)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-1) * zny(1:imx-1, 1:jmx-1, kmx)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-1) * znz(1:imx-1, 1:jmx-1, kmx)) &
                                ) * znx(1:imx-1, 1:jmx-1, kmx) &
                            )
                    y_speed(1:imx-1, 1:jmx-1, kmx) = y_speed(1:imx-1, 1:jmx-1, kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-1) * znx(1:imx-1, 1:jmx-1, kmx)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-1) * zny(1:imx-1, 1:jmx-1, kmx)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-1) * znz(1:imx-1, 1:jmx-1, kmx)) &
                                ) * zny(1:imx-1, 1:jmx-1, kmx) &
                            )
                    z_speed(1:imx-1, 1:jmx-1, kmx) = z_speed(1:imx-1, 1:jmx-1, kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-1) * znx(1:imx-1, 1:jmx-1, kmx)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-1) * zny(1:imx-1, 1:jmx-1, kmx)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-1) * znz(1:imx-1, 1:jmx-1, kmx)) &
                                ) * znz(1:imx-1, 1:jmx-1, kmx) &
                            )
                end if
                !!!!-----------------------------------------------------------------------
                !!~!!  GHOST Cell Changes !!~!!
                !!!!-----------------------------------------------------------------------
                if (face == "imin") then
                    x_speed(-1, 1:jmx-1, 1:kmx-1) = x_speed(2, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(2, 1:jmx-1, 1:kmx-1) * xnx(2, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(2, 1:jmx-1, 1:kmx-1) * xny(2, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(2, 1:jmx-1, 1:kmx-1) * xnz(2, 1:jmx-1, 1:kmx-1)) &
                                ) * xnx(2, 1:jmx-1, 1:kmx-1) &
                            )
                    y_speed(-1, 1:jmx-1, 1:kmx-1) = y_speed(2, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(2, 1:jmx-1, 1:kmx-1) * xnx(2, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(2, 1:jmx-1, 1:kmx-1) * xny(2, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(2, 1:jmx-1, 1:kmx-1) * xnz(2, 1:jmx-1, 1:kmx-1)) &
                                ) * xny(2, 1:jmx-1, 1:kmx-1) &
                            )
                    z_speed(-1, 1:jmx-1, 1:kmx-1) = z_speed(2, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(2, 1:jmx-1, 1:kmx-1) * xnx(2, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(2, 1:jmx-1, 1:kmx-1) * xny(2, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(2, 1:jmx-1, 1:kmx-1) * xnz(2, 1:jmx-1, 1:kmx-1)) &
                                ) * xnz(2, 1:jmx-1, 1:kmx-1) &
                            )
                else if (face == "imax") then
                    x_speed(imx+1, 1:jmx-1, 1:kmx-1) = x_speed(imx-2, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-2, 1:jmx-1, 1:kmx-1) * xnx(imx-1, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-2, 1:jmx-1, 1:kmx-1) * xny(imx-1, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-2, 1:jmx-1, 1:kmx-1) * xnz(imx-1, 1:jmx-1, 1:kmx-1)) &
                                ) * xnx(imx-1, 1:jmx-1, 1:kmx-1) &
                            )
                    y_speed(imx+1, 1:jmx-1, 1:kmx-1) = y_speed(imx-2, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-2, 1:jmx-1, 1:kmx-1) * xnx(imx-1, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-2, 1:jmx-1, 1:kmx-1) * xny(imx-1, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-2, 1:jmx-1, 1:kmx-1) * xnz(imx-1, 1:jmx-1, 1:kmx-1)) &
                                ) * xny(imx-1, 1:jmx-1, 1:kmx-1) &
                            )
                    z_speed(imx+1, 1:jmx-1, 1:kmx-1) = z_speed(imx-2, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-2, 1:jmx-1, 1:kmx-1) * xnx(imx-1, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-2, 1:jmx-1, 1:kmx-1) * xny(imx-1, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-2, 1:jmx-1, 1:kmx-1) * xnz(imx-1, 1:jmx-1, 1:kmx-1)) &
                                ) * xnz(imx-1, 1:jmx-1, 1:kmx-1) &
                            )
                else if (face == "jmin") then
                    x_speed(1:imx-1, -1, 1:kmx-1) = x_speed(1:imx-1, 2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 2, 1:kmx-1) * ynx(1:imx-1, 2, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 2, 1:kmx-1) * yny(1:imx-1, 2, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 2, 1:kmx-1) * ynz(1:imx-1, 2, 1:kmx-1)) &
                                ) * ynx(1:imx-1, 2, 1:kmx-1) &
                            )
                    y_speed(1:imx-1, -1, 1:kmx-1) = y_speed(1:imx-1, 2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 2, 1:kmx-1) * ynx(1:imx-1, 2, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 2, 1:kmx-1) * yny(1:imx-1, 2, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 2, 1:kmx-1) * ynz(1:imx-1, 2, 1:kmx-1)) &
                                ) * yny(1:imx-1, 2, 1:kmx-1) &
                            )
                    z_speed(1:imx-1, -1, 1:kmx-1) = z_speed(1:imx-1, 2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 2, 1:kmx-1) * ynx(1:imx-1, 2, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 2, 1:kmx-1) * yny(1:imx-1, 2, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 2, 1:kmx-1) * ynz(1:imx-1, 2, 1:kmx-1)) &
                                ) * ynz(1:imx-1, 2, 1:kmx-1) &
                            )
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx+1, 1:kmx-1) = x_speed(1:imx-1, jmx-2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-2, 1:kmx-1) * ynx(1:imx-1, jmx-1, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-2, 1:kmx-1) * yny(1:imx-1, jmx-1, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-2, 1:kmx-1) * ynz(1:imx-1, jmx-1, 1:kmx-1)) &
                                ) * ynx(1:imx-1, jmx-1, 1:kmx-1) &
                            )
                    y_speed(1:imx-1, jmx+1, 1:kmx-1) = y_speed(1:imx-1, jmx-2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-2, 1:kmx-1) * ynx(1:imx-1, jmx-1, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-2, 1:kmx-1) * yny(1:imx-1, jmx-1, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-2, 1:kmx-1) * ynz(1:imx-1, jmx-1, 1:kmx-1)) &
                                ) * yny(1:imx-1, jmx-1, 1:kmx-1) &
                            )
                    z_speed(1:imx-1, jmx+1, 1:kmx-1) = z_speed(1:imx-1, jmx-2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-2, 1:kmx-1) * ynx(1:imx-1, jmx-1, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-2, 1:kmx-1) * yny(1:imx-1, jmx-1, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-2, 1:kmx-1) * ynz(1:imx-1, jmx-1, 1:kmx-1)) &
                                ) * ynz(1:imx-1, jmx-1, 1:kmx-1) &
                            )
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, -1) = x_speed(1:imx-1, 1:jmx-1, 2) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 2) * znx(1:imx-1, 1:jmx-1, 2)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 2) * zny(1:imx-1, 1:jmx-1, 2)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 2) * znz(1:imx-1, 1:jmx-1, 2)) &
                                ) * znx(1:imx-1, 1:jmx-1, 2) &
                            )
                    y_speed(1:imx-1, 1:jmx-1, -1) = y_speed(1:imx-1, 1:jmx-1, 2) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 2) * znx(1:imx-1, 1:jmx-1, 2)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 2) * zny(1:imx-1, 1:jmx-1, 2)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 2) * znz(1:imx-1, 1:jmx-1, 2)) &
                                ) * zny(1:imx-1, 1:jmx-1, 2) &
                            )
                    z_speed(1:imx-1, 1:jmx-1, -1) = z_speed(1:imx-1, 1:jmx-1, 2) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 2) * znx(1:imx-1, 1:jmx-1, 2)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 2) * zny(1:imx-1, 1:jmx-1, 2)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 2) * znz(1:imx-1, 1:jmx-1, 2)) &
                                ) * znz(1:imx-1, 1:jmx-1, 2) &
                            )
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx+1) = x_speed(1:imx-1, 1:jmx-1, kmx-2) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-2) * znx(1:imx-1, 1:jmx-1, kmx-1)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-2) * zny(1:imx-1, 1:jmx-1, kmx-1)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-2) * znz(1:imx-1, 1:jmx-1, kmx-1)) &
                                ) * znx(1:imx-1, 1:jmx-1, kmx-1) &
                            )
                    y_speed(1:imx-1, 1:jmx-1, kmx+1) = y_speed(1:imx-1, 1:jmx-1, kmx-2) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-2) * znx(1:imx-1, 1:jmx-1, kmx-1)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-2) * zny(1:imx-1, 1:jmx-1, kmx-1)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-2) * znz(1:imx-1, 1:jmx-1, kmx-1)) &
                                ) * zny(1:imx-1, 1:jmx-1, kmx-1) &
                            )
                    z_speed(1:imx-1, 1:jmx-1, kmx+1) = z_speed(1:imx-1, 1:jmx-1, kmx-2) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-2) * znx(1:imx-1, 1:jmx-1, kmx-1)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-2) * zny(1:imx-1, 1:jmx-1, kmx-1)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-2) * znz(1:imx-1, 1:jmx-1, kmx-1)) &
                                ) * znz(1:imx-1, 1:jmx-1, kmx-1) &
                            )
                end if
                if (face == "imin") then
                    x_speed(-2, 1:jmx-1, 1:kmx-1) = x_speed(3, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(3, 1:jmx-1, 1:kmx-1) * xnx(3, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(3, 1:jmx-1, 1:kmx-1) * xny(3, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(3, 1:jmx-1, 1:kmx-1) * xnz(3, 1:jmx-1, 1:kmx-1)) &
                                ) * xnx(3, 1:jmx-1, 1:kmx-1) &
                            )
                    y_speed(-2, 1:jmx-1, 1:kmx-1) = y_speed(3, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(3, 1:jmx-1, 1:kmx-1) * xnx(3, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(3, 1:jmx-1, 1:kmx-1) * xny(3, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(3, 1:jmx-1, 1:kmx-1) * xnz(3, 1:jmx-1, 1:kmx-1)) &
                                ) * xny(3, 1:jmx-1, 1:kmx-1) &
                            )
                    z_speed(-2, 1:jmx-1, 1:kmx-1) = z_speed(3, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(3, 1:jmx-1, 1:kmx-1) * xnx(3, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(3, 1:jmx-1, 1:kmx-1) * xny(3, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(3, 1:jmx-1, 1:kmx-1) * xnz(3, 1:jmx-1, 1:kmx-1)) &
                                ) * xnz(3, 1:jmx-1, 1:kmx-1) &
                            )
                else if (face == "imax") then
                    x_speed(imx+2, 1:jmx-1, 1:kmx-1) = x_speed(imx-3, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-3, 1:jmx-1, 1:kmx-1) * xnx(imx-2, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-3, 1:jmx-1, 1:kmx-1) * xny(imx-2, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-3, 1:jmx-1, 1:kmx-1) * xnz(imx-2, 1:jmx-1, 1:kmx-1)) &
                                ) * xnx(imx-2, 1:jmx-1, 1:kmx-1) &
                            )
                    y_speed(imx+2, 1:jmx-1, 1:kmx-1) = y_speed(imx-3, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-3, 1:jmx-1, 1:kmx-1) * xnx(imx-2, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-3, 1:jmx-1, 1:kmx-1) * xny(imx-2, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-3, 1:jmx-1, 1:kmx-1) * xnz(imx-2, 1:jmx-1, 1:kmx-1)) &
                                ) * xny(imx-2, 1:jmx-1, 1:kmx-1) &
                            )
                    z_speed(imx+2, 1:jmx-1, 1:kmx-1) = z_speed(imx-3, 1:jmx-1, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(imx-3, 1:jmx-1, 1:kmx-1) * xnx(imx-2, 1:jmx-1, 1:kmx-1)) + &
                                (y_speed(imx-3, 1:jmx-1, 1:kmx-1) * xny(imx-2, 1:jmx-1, 1:kmx-1)) + &
                                (z_speed(imx-3, 1:jmx-1, 1:kmx-1) * xnz(imx-2, 1:jmx-1, 1:kmx-1)) &
                                ) * xnz(imx-2, 1:jmx-1, 1:kmx-1) &
                            )
                else if (face == "jmin") then
                    x_speed(1:imx-1, -2, 1:kmx-1) = x_speed(1:imx-1, 3, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 3, 1:kmx-1) * ynx(1:imx-1, 3, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 3, 1:kmx-1) * yny(1:imx-1, 3, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 3, 1:kmx-1) * ynz(1:imx-1, 3, 1:kmx-1)) &
                                ) * ynx(1:imx-1, 3, 1:kmx-1) &
                            )
                    y_speed(1:imx-1, -2, 1:kmx-1) = y_speed(1:imx-1, 3, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 3, 1:kmx-1) * ynx(1:imx-1, 3, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 3, 1:kmx-1) * yny(1:imx-1, 3, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 3, 1:kmx-1) * ynz(1:imx-1, 3, 1:kmx-1)) &
                                ) * yny(1:imx-1, 3, 1:kmx-1) &
                            )
                    z_speed(1:imx-1, -2, 1:kmx-1) = z_speed(1:imx-1, 3, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 3, 1:kmx-1) * ynx(1:imx-1, 3, 1:kmx-1)) + &
                                (y_speed(1:imx-1, 3, 1:kmx-1) * yny(1:imx-1, 3, 1:kmx-1)) + &
                                (z_speed(1:imx-1, 3, 1:kmx-1) * ynz(1:imx-1, 3, 1:kmx-1)) &
                                ) * ynz(1:imx-1, 3, 1:kmx-1) &
                            )
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx+2, 1:kmx-1) = x_speed(1:imx-1, jmx-3, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-3, 1:kmx-1) * ynx(1:imx-1, jmx-2, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-3, 1:kmx-1) * yny(1:imx-1, jmx-2, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-3, 1:kmx-1) * ynz(1:imx-1, jmx-2, 1:kmx-1)) &
                                ) * ynx(1:imx-1, jmx-2, 1:kmx-1) &
                            )
                    y_speed(1:imx-1, jmx+2, 1:kmx-1) = y_speed(1:imx-1, jmx-2, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-3, 1:kmx-1) * ynx(1:imx-1, jmx-2, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-3, 1:kmx-1) * yny(1:imx-1, jmx-2, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-3, 1:kmx-1) * ynz(1:imx-1, jmx-2, 1:kmx-1)) &
                                ) * yny(1:imx-1, jmx-2, 1:kmx-1) &
                            )
                    z_speed(1:imx-1, jmx+2, 1:kmx-1) = z_speed(1:imx-1, jmx-3, 1:kmx-1) - &
                            (2. * ( &
                                (x_speed(1:imx-1, jmx-3, 1:kmx-1) * ynx(1:imx-1, jmx-2, 1:kmx-1)) + &
                                (y_speed(1:imx-1, jmx-3, 1:kmx-1) * yny(1:imx-1, jmx-2, 1:kmx-1)) + &
                                (z_speed(1:imx-1, jmx-3, 1:kmx-1) * ynz(1:imx-1, jmx-2, 1:kmx-1)) &
                                ) * ynz(1:imx-1, jmx-2, 1:kmx-1) &
                            )
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, -2) = x_speed(1:imx-1, 1:jmx-1, 3) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 3) * znx(1:imx-1, 1:jmx-1, 3)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 3) * zny(1:imx-1, 1:jmx-1, 3)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 3) * znz(1:imx-1, 1:jmx-1, 3)) &
                                ) * znx(1:imx-1, 1:jmx-1, 2) &
                            )
                    y_speed(1:imx-1, 1:jmx-1, -2) = y_speed(1:imx-1, 1:jmx-1, 3) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 3) * znx(1:imx-1, 1:jmx-1, 3)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 3) * zny(1:imx-1, 1:jmx-1, 3)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 3) * znz(1:imx-1, 1:jmx-1, 3)) &
                                ) * zny(1:imx-1, 1:jmx-1, 3) &
                            )
                    z_speed(1:imx-1, 1:jmx-1, -2) = z_speed(1:imx-1, 1:jmx-1, 3) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, 3) * znx(1:imx-1, 1:jmx-1, 3)) + &
                                (y_speed(1:imx-1, 1:jmx-1, 3) * zny(1:imx-1, 1:jmx-1, 3)) + &
                                (z_speed(1:imx-1, 1:jmx-1, 3) * znz(1:imx-1, 1:jmx-1, 3)) &
                                ) * znz(1:imx-1, 1:jmx-1, 3) &
                            )
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx+2) = x_speed(1:imx-1, 1:jmx-1, kmx-3) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-3) * znx(1:imx-1, 1:jmx-1, kmx-2)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-3) * zny(1:imx-1, 1:jmx-1, kmx-2)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-3) * znz(1:imx-1, 1:jmx-1, kmx-2)) &
                                ) * znx(1:imx-1, 1:jmx-1, kmx-2) &
                            )
                    y_speed(1:imx-1, 1:jmx-1, kmx+2) = y_speed(1:imx-1, 1:jmx-1, kmx-3) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-3) * znx(1:imx-1, 1:jmx-1, kmx-2)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-3) * zny(1:imx-1, 1:jmx-1, kmx-2)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-3) * znz(1:imx-1, 1:jmx-1, kmx-2)) &
                                ) * zny(1:imx-1, 1:jmx-1, kmx-2) &
                            )
                    z_speed(1:imx-1, 1:jmx-1, kmx+2) = z_speed(1:imx-1, 1:jmx-1, kmx-3) - &
                            (2. * ( &
                                (x_speed(1:imx-1, 1:jmx-1, kmx-3) * znx(1:imx-1, 1:jmx-1, kmx-2)) + &
                                (y_speed(1:imx-1, 1:jmx-1, kmx-3) * zny(1:imx-1, 1:jmx-1, kmx-2)) + &
                                (z_speed(1:imx-1, 1:jmx-1, kmx-3) * znz(1:imx-1, 1:jmx-1, kmx-2)) &
                                ) * znz(1:imx-1, 1:jmx-1, kmx-2) &
                            )
                end if
     !      else
     !          if (face == "imin") then
     !              x_speed(0, cell_ind) = x_speed(1, cell_ind) - &
     !                      (2. * ( &
     !                          (x_speed(1, cell_ind) * xnx(1, cell_ind)) + &
     !                          (y_speed(1, cell_ind) * xny(1, cell_ind)) &
     !                          ) * xnx(1, cell_ind) &
     !                      )
     !              y_speed(0, cell_ind) = y_speed(1, cell_ind) - &
     !                      (2. * ( &
     !                          (x_speed(1, cell_ind) * xnx(1, cell_ind)) + &
     !                          (y_speed(1, cell_ind) * xny(1, cell_ind)) &
     !                          ) * xny(1, cell_ind) &
     !                      )
     !          else if (face == "imax") then
     !              x_speed(imx, cell_ind) = x_speed(imx-1, cell_ind) - &
     !                      (2. * ( &
     !                          (x_speed(imx-1, cell_ind) * xnx(imx, cell_ind)) + &
     !                          (y_speed(imx-1, cell_ind) * xny(imx, cell_ind)) &
     !                          ) * xnx(imx, cell_ind) &
     !                      )
     !              y_speed(imx, cell_ind) = y_speed(imx-1, cell_ind) - &
     !                      (2. * ( &
     !                          (x_speed(imx-1, cell_ind) * xnx(imx, cell_ind)) + &
     !                          (y_speed(imx-1, cell_ind) * xny(imx, cell_ind)) &
     !                          ) * xny(imx, cell_ind) &
     !                      )
     !          else if (face == "jmin") then
     !              x_speed(cell_ind, 0) = x_speed(cell_ind, 1) - &
     !                      (2. * ( &
     !                          (x_speed(cell_ind, 1) * ynx(cell_ind, 1)) + &
     !                          (y_speed(cell_ind, 1) * yny(cell_ind, 1)) &
     !                          ) * ynx(cell_ind, 1) &
     !                      )
     !              y_speed(cell_ind, 0) = y_speed(cell_ind, 1) - &
     !                      (2. * ( &
     !                          (x_speed(cell_ind, 1) * ynx(cell_ind, 1)) + &
     !                          (y_speed(cell_ind, 1) * yny(cell_ind, 1)) &
     !                          ) * yny(cell_ind, 1) &
     !                      )
     !          else if (face == "jmax") then
     !              x_speed(cell_ind, jmx) = x_speed(cell_ind, jmx-1) - &
     !                      (2. * ( &
     !                          (x_speed(cell_ind, jmx-1) * ynx(cell_ind, jmx)) + &
     !                          (y_speed(cell_ind, jmx-1) * yny(cell_ind, jmx)) &
     !                          ) * ynx(cell_ind, jmx) &
     !                      )
     !              y_speed(cell_ind, jmx) = y_speed(cell_ind, jmx-1) - &
     !                      (2. * ( &
     !                          (x_speed(cell_ind, jmx-1) * ynx(cell_ind, jmx)) + &
     !                          (y_speed(cell_ind, jmx-1) * yny(cell_ind, jmx)) &
     !                          ) * yny(cell_ind, jmx) &
     !                      )
     !          end if
            end if

        end subroutine flow_tangency

        subroutine no_slip(face, cell_ind)
            !-----------------------------------------------------------
            ! The no slip condition is an enforcement of the no slip
            ! wall condition, i.e., the velocity at the wall should be 
            ! zero. This is achieved by setting the velocity 
            ! of the ghost cell to be the negative of its neighbouring 
            ! cell's velocity.
            !-----------------------------------------------------------
        
            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    x_speed(0, 1:jmx-1, 1:kmx-1) = - x_speed(1, 1:jmx-1, 1:kmx-1)
                    y_speed(0, 1:jmx-1, 1:kmx-1) = - y_speed(1, 1:jmx-1, 1:kmx-1)
                    z_speed(0, 1:jmx-1, 1:kmx-1) = - z_speed(1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    x_speed(imx, 1:jmx-1, 1:kmx-1) = - x_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    y_speed(imx, 1:jmx-1, 1:kmx-1) = - y_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    z_speed(imx, 1:jmx-1, 1:kmx-1) = - z_speed(imx-1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    x_speed(1:imx-1, 0, 1:kmx-1) = - x_speed(1:imx-1, 1, 1:kmx-1)
                    y_speed(1:imx-1, 0, 1:kmx-1) = - y_speed(1:imx-1, 1, 1:kmx-1)
                    z_speed(1:imx-1, 0, 1:kmx-1) = - z_speed(1:imx-1, 1, 1:kmx-1)
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx, 1:kmx-1) = - x_speed(1:imx-1, jmx-1, 1:kmx-1)
                    y_speed(1:imx-1, jmx, 1:kmx-1) = - y_speed(1:imx-1, jmx-1, 1:kmx-1)
                    z_speed(1:imx-1, jmx, 1:kmx-1) = - z_speed(1:imx-1, jmx-1, 1:kmx-1)
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, 0) = - x_speed(1:imx-1, 1:jmx-1, 1)
                    y_speed(1:imx-1, 1:jmx-1, 0) = - y_speed(1:imx-1, 1:jmx-1, 1)
                    z_speed(1:imx-1, 1:jmx-1, 0) = - z_speed(1:imx-1, 1:jmx-1, 1)
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx) = - x_speed(1:imx-1, 1:jmx-1, kmx-1)
                    y_speed(1:imx-1, 1:jmx-1, kmx) = - y_speed(1:imx-1, 1:jmx-1, kmx-1)
                    z_speed(1:imx-1, 1:jmx-1, kmx) = - z_speed(1:imx-1, 1:jmx-1, kmx-1)
                end if
                !----------------------------------------------------------------
                ! Extra Ghost cells
                !---------------------------------------------------------------
                
                if (face == "imin") then
                    x_speed(-1, 1:jmx-1, 1:kmx-1) = - x_speed(2, 1:jmx-1, 1:kmx-1)
                    y_speed(-1, 1:jmx-1, 1:kmx-1) = - y_speed(2, 1:jmx-1, 1:kmx-1)
                    z_speed(-1, 1:jmx-1, 1:kmx-1) = - z_speed(2, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    x_speed(imx+1, 1:jmx-1, 1:kmx-1) = - x_speed(imx-2, 1:jmx-1, 1:kmx-1)
                    y_speed(imx+1, 1:jmx-1, 1:kmx-1) = - y_speed(imx-2, 1:jmx-1, 1:kmx-1)
                    z_speed(imx+1, 1:jmx-1, 1:kmx-1) = - z_speed(imx-2, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    x_speed(1:imx-1, -1, 1:kmx-1) = - x_speed(1:imx-1, 2, 1:kmx-1)
                    y_speed(1:imx-1, -1, 1:kmx-1) = - y_speed(1:imx-1, 2, 1:kmx-1)
                    z_speed(1:imx-1, -1, 1:kmx-1) = - z_speed(1:imx-1, 2, 1:kmx-1)
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx+1, 1:kmx-1) = - x_speed(1:imx-1, jmx-2, 1:kmx-1)
                    y_speed(1:imx-1, jmx+1, 1:kmx-1) = - y_speed(1:imx-1, jmx-2, 1:kmx-1)
                    z_speed(1:imx-1, jmx+1, 1:kmx-1) = - z_speed(1:imx-1, jmx-2, 1:kmx-1)
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, -1) = - x_speed(1:imx-1, 1:jmx-1, 2)
                    y_speed(1:imx-1, 1:jmx-1, -1) = - y_speed(1:imx-1, 1:jmx-1, 2)
                    z_speed(1:imx-1, 1:jmx-1, -1) = - z_speed(1:imx-1, 1:jmx-1, 2)
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx+1) = - x_speed(1:imx-1, 1:jmx-1, kmx-2)
                    y_speed(1:imx-1, 1:jmx-1, kmx+1) = - y_speed(1:imx-1, 1:jmx-1, kmx-2)
                    z_speed(1:imx-1, 1:jmx-1, kmx+1) = - z_speed(1:imx-1, 1:jmx-1, kmx-2)
                end if

                !------------------------------------------------------------------
                ! 3rd ghost cell
                !-----------------------------------------------------------------

                if (face == "imin") then
                    x_speed(-2, 1:jmx-1, 1:kmx-1) = - x_speed(3, 1:jmx-1, 1:kmx-1)
                    y_speed(-2, 1:jmx-1, 1:kmx-1) = - y_speed(3, 1:jmx-1, 1:kmx-1)
                    z_speed(-2, 1:jmx-1, 1:kmx-1) = - z_speed(3, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    x_speed(imx+2, 1:jmx-1, 1:kmx-1) = - x_speed(imx-3, 1:jmx-1, 1:kmx-1)
                    y_speed(imx+2, 1:jmx-1, 1:kmx-1) = - y_speed(imx-3, 1:jmx-1, 1:kmx-1)
                    z_speed(imx+2, 1:jmx-1, 1:kmx-1) = - z_speed(imx-3, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    x_speed(1:imx-1, -2, 1:kmx-1) = - x_speed(1:imx-1, 3, 1:kmx-1)
                    y_speed(1:imx-1, -2, 1:kmx-1) = - y_speed(1:imx-1, 3, 1:kmx-1)
                    z_speed(1:imx-1, -2, 1:kmx-1) = - z_speed(1:imx-1, 3, 1:kmx-1)
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx+2, 1:kmx-1) = - x_speed(1:imx-1, jmx-3, 1:kmx-1)
                    y_speed(1:imx-1, jmx+2, 1:kmx-1) = - y_speed(1:imx-1, jmx-3, 1:kmx-1)
                    z_speed(1:imx-1, jmx+2, 1:kmx-1) = - z_speed(1:imx-1, jmx-3, 1:kmx-1)
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, -2) = - x_speed(1:imx-1, 1:jmx-1, 3)
                    y_speed(1:imx-1, 1:jmx-1, -2) = - y_speed(1:imx-1, 1:jmx-1, 3)
                    z_speed(1:imx-1, 1:jmx-1, -2) = - z_speed(1:imx-1, 1:jmx-1, 3)
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx+2) = - x_speed(1:imx-1, 1:jmx-1, kmx-3)
                    y_speed(1:imx-1, 1:jmx-1, kmx+2) = - y_speed(1:imx-1, 1:jmx-1, kmx-3)
                    z_speed(1:imx-1, 1:jmx-1, kmx+2) = - z_speed(1:imx-1, 1:jmx-1, kmx-3)
                end if


            end if

            include "turbulence_models/include/bc/no_slip.inc"

        end subroutine no_slip

        subroutine periodic(face, cell_ind)
            !-----------------------------------------------------------
            ! The periodic boundary condition is used when the domain
            ! wraps around itself. In those cases, the ghost cells of the
            ! extremes are the same as the interior cells in the opposite
            ! face.
            !-----------------------------------------------------------
        
            implicit none
            character(len=4), intent(in) :: face
            integer, intent(in) :: cell_ind

            if (cell_ind == -1) then
                if (face == "imin") then
                    x_speed(0, 1:jmx-1, 1:kmx-1) = x_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    y_speed(0, 1:jmx-1, 1:kmx-1) = y_speed(imx-1, 1:jmx-1, 1:kmx-1)
                    z_speed(0, 1:jmx-1, 1:kmx-1) = z_speed(imx-1, 1:jmx-1, 1:kmx-1)
                else if (face == "imax") then
                    x_speed(imx, 1:jmx-1, 1:kmx-1) = x_speed(1, 1:jmx-1, 1:kmx-1)
                    y_speed(imx, 1:jmx-1, 1:kmx-1) = y_speed(1, 1:jmx-1, 1:kmx-1)
                    z_speed(imx, 1:jmx-1, 1:kmx-1) = z_speed(1, 1:jmx-1, 1:kmx-1)
                else if (face == "jmin") then
                    x_speed(1:imx-1, 0, 1:kmx-1) = x_speed(1:imx-1, jmx-1, 1:kmx-1)
                    y_speed(1:imx-1, 0, 1:kmx-1) = y_speed(1:imx-1, jmx-1, 1:kmx-1)
                    z_speed(1:imx-1, 0, 1:kmx-1) = z_speed(1:imx-1, jmx-1, 1:kmx-1)
                else if (face == "jmax") then
                    x_speed(1:imx-1, jmx, 1:kmx-1) = x_speed(1:imx-1, 1, 1:kmx-1)
                    y_speed(1:imx-1, jmx, 1:kmx-1) = y_speed(1:imx-1, 1, 1:kmx-1)
                    z_speed(1:imx-1, jmx, 1:kmx-1) = z_speed(1:imx-1, 1, 1:kmx-1)
                else if (face == "kmin") then
                    x_speed(1:imx-1, 1:jmx-1, 0) = x_speed(1:imx-1, 1:jmx-1, kmx-1)
                    y_speed(1:imx-1, 1:jmx-1, 0) = y_speed(1:imx-1, 1:jmx-1, kmx-1)
                    z_speed(1:imx-1, 1:jmx-1, 0) = z_speed(1:imx-1, 1:jmx-1, kmx-1)
                else if (face == "kmax") then
                    x_speed(1:imx-1, 1:jmx-1, kmx) = x_speed(1:imx-1, 1:jmx-1, 1)
                    y_speed(1:imx-1, 1:jmx-1, kmx) = y_speed(1:imx-1, 1:jmx-1, 1)
                    z_speed(1:imx-1, 1:jmx-1, kmx) = z_speed(1:imx-1, 1:jmx-1, 1)
                end if
            end if

        end subroutine periodic

        subroutine extra_ghost_cells(face)
            !-----------------------------------------------------------
            ! Since higher order methods are used, the interface boundary
            ! condition should require higher order accurate reconstructions
            ! even at the physical boundary. To do that, it is estimated
            ! that totally 3 layers of ghost cells are required.
            ! The above boundary conditions are for the first layer. The
            ! boundary conditions for the other layers, as described in this
            ! subroutine, assume a linear extrapolation using interior and
            ! first ghost cell.
            !-----------------------------------------------------------
        
            implicit none
            character(len=4), intent(in) :: face

            if (face == "imin") then
                qp(-1, 1:jmx-1, 1:kmx-1, :) = - qp(1, 1:jmx-1, 1:kmx-1, :) + &
                                            2 * qp(0, 1:jmx-1, 1:kmx-1, :)
                qp(-2, 1:jmx-1, 1:kmx-1, :) = - 2 * qp(1, 1:jmx-1, 1:kmx-1, :) + &
                                                3 * qp(0, 1:jmx-1, 1:kmx-1, :)
            else if (face == "imax") then
                qp(imx+1, 1:jmx-1, 1:kmx-1, :) = - qp(imx-1, 1:jmx-1, 1:kmx-1, :) + &
                                               2 * qp(imx, 1:jmx-1, 1:kmx-1, :)
                qp(imx+2, 1:jmx-1, 1:kmx-1, :) = - 2 * qp(imx-1, 1:jmx-1, 1:kmx-1, :) + &
                                                   3 * qp(imx, 1:jmx-1, 1:kmx-1, :)
            else if (face == "jmin") then
                qp(1:imx-1, -1, 1:kmx-1, :) = - qp(1:imx-1, 1, 1:kmx-1, :) + &
                                            2 * qp(1:imx-1, 0, 1:kmx-1, :)
                qp(1:imx-1, -2, 1:kmx-1, :) = - 2 * qp(1:imx-1, 1, 1:kmx-1, :) + &
                                                3 * qp(1:imx-1, 0, 1:kmx-1, :)
            else if (face == "jmax") then
                qp(1:imx-1, jmx+1, 1:kmx-1, :) = - qp(1:imx-1, jmx-1, 1:kmx-1, :) + &
                                               2 * qp(1:imx-1, jmx, 1:kmx-1, :)
                qp(1:imx-1, jmx+2, 1:kmx-1, :) = - 2 * qp(1:imx-1, jmx-1, 1:kmx-1, :) + &
                                                   3 * qp(1:imx-1, jmx, 1:kmx-1, :)
            else if (face == "kmin") then
                qp(1:imx-1, 1:jmx-1, -1, :) = - qp(1:imx-1, 1:jmx-1, 1, :) + &
                                            2 * qp(1:imx-1, 1:jmx-1, 0, :)
                qp(1:imx-1, 1:jmx-1, -2, :) = - 2 * qp(1:imx-1, 1:jmx-1, 1, :) + &
                                                3 * qp(1:imx-1, 1:jmx-1, 0, :)
            else if (face == "kmax") then
                qp(1:imx-1, 1:jmx-1, kmx+1, :) = - qp(1:imx-1, 1:jmx-1, kmx-1, :) + &
                                               2 * qp(1:imx-1, 1:jmx-1, kmx, :)
                qp(1:imx-1, 1:jmx-1, kmx+2, :) = - 2 * qp(1:imx-1, 1:jmx-1, kmx-1, :) + &
                                                   3 * qp(1:imx-1, 1:jmx-1, kmx, :)
            end if

        end subroutine extra_ghost_cells

        include "turbulence_models/include/bc/fix_and_copy_turb_var.inc"


end module boundary_conditions
