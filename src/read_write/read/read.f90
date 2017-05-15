module read
  !-------------------------------------------------------
  ! This module read input control files which include:
  !   1. control.md
  !   2. fvscheme.md
  !   3. flow.md
  !   4. res_control.md ! todo
  !   5. state_read_write_control.md ! todo
  !------------------------------------------------------


  use global, only: CONTROL_FILE_UNIT
  use global, only:  SCHEME_FILE_UNIT
  use global, only:    FLOW_FILE_UNIT
  use global, only: control_file
  use global, only:  scheme_file
  use global, only:    flow_file
  use global, only: STRING_BUFFER_LENGTH
  use global, only : OUTIN_FILE_UNIT
  use global, only : outin_file

  use global_vars, only: CFL
  use global_vars, only: max_iters
  use global_vars, only: start_from
  use global_vars, only: checkpoint_iter
  use global_vars, only: res_write_interval
  use global_vars, only: write_file_format
  use global_vars, only: write_data_format
  use global_vars, only: write_percision
  use global_vars, only: purge_write
  use global_vars, only: tolerance

  use global_vars, only: time_stepping_method
  use global_vars, only: time_step_accuracy
  use global_vars, only: global_time_step

  use global_vars, only: n_var
  use global_vars, only: free_stream_density
  use global_vars, only: free_stream_x_speed
  use global_vars, only: free_stream_y_speed
  use global_vars, only: free_stream_z_speed
  use global_vars, only: free_stream_pressure
  use global_vars, only: free_stream_tk
  use global_vars, only: free_stream_tw
  use global_vars, only: gm    !gamma
  use global_vars, only: R_gas !univarsal gas constant
  use global_vars, only: mu_ref !viscoity
  use global_vars, only: mu_variation !viscoity variation type
  use global_vars, only: T_ref
  use global_vars, only: Sutherland_temp
  use global_vars, only: Pr !prandtl number
  use global_vars, only: ilimiter_switch
  use global_vars, only: PB_switch
  use global_vars, only: accur
  
  use global_vars, only: interpolant
  use global_vars, only: scheme_name
  use global_vars, only: turbulence
  use global_vars, only: rw_list
  use global_vars, only: n_write
  use utils      , only: DEBUG_LEVEL
  use utils      , only: dmsg
  use string

  implicit none
  private

  public :: read_input_and_controls

    contains

      subroutine read_input_and_controls()
        implicit none
        call read_controls()
        call read_scheme()
        call read_flow()
        call read_output_control()
      end subroutine read_input_and_controls


      subroutine get_next_token(token_file_unit, buf)
        !-----------------------------------------------------------
        ! Extract the next token from the config file
        !
        ! Each token is on a separate line.
        ! There may be multiple comments (lines beginning with #) 
        ! and blank lines in between.
        ! The purpose of this subroutine is to ignore all these 
        ! lines and return the next "useful" line.
        !-----------------------------------------------------------

        implicit none
        integer                            , intent(in)  :: token_file_unit
        character(len=STRING_BUFFER_LENGTH), intent(out) :: buf
        integer :: ios

        do
            read(token_file_unit, '(A)', iostat=ios) buf
            if (ios /= 0) then
                print *, 'Error while reading config file.'
                print *, 'Current buffer length is set to: ', &
                        STRING_BUFFER_LENGTH
                stop
            end if
            if (index(buf, '#') == 1) then
                ! The current line begins with a hash
                ! Ignore it
                continue
            else if (len_trim(buf) == 0) then
                ! The current line is empty
                ! Ignore it
                continue
            else
                ! A new token has been found
                ! Break out
                exit
            end if
        end do
        call dmsg(0, 'read', 'get_next_token', 'Returning: ' // trim(buf))

      end subroutine get_next_token



      subroutine read_controls()
        !---------------------------------------------
        ! read control.md file
        !---------------------------------------------
        implicit none
        character(len=STRING_BUFFER_LENGTH) :: buf

        call dmsg(1, 'read', 'read_controls')

        open(CONTROL_FILE_UNIT, file=control_file, status='old', action='read')

        !ignoring file header
        read(CONTROL_FILE_UNIT,*)
        read(CONTROL_FILE_UNIT,*)
        read(CONTROL_FILE_UNIT,*)

        ! READ CFL
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) CFL
        call dmsg(5, 'read', 'read_controls', &
                msg='CFL = ' + CFL)

        ! READ start_from
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) start_from
        call dmsg(5, 'read', 'read_controls', &
                msg='Simlulation  start from  level = ' + start_from)

        ! READ max_iters
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) max_iters
        call dmsg(5, 'read', 'read_controls', &
                msg=' Stop at iteration = ' + max_iters)

        ! READ checkpoint_iter
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) checkpoint_iter
        call dmsg(5, 'read', 'read_controls', &
                msg=' Solution write interval = ' + checkpoint_iter)

        ! READ write_file_format
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) write_file_format
        call dmsg(5, 'read', 'read_controls', &
                msg='Solution file format  = ' + write_file_format)

        ! READ write_data_format
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) write_data_format
        call dmsg(5, 'read', 'read_controls', &
                msg='solution file data format = ' + write_data_format)

        ! READ write_percision
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) write_percision
        call dmsg(5, 'read', 'read_controls', &
                msg='File write percision = ' + write_percision)

        ! READ purge_write
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) purge_write
        call dmsg(5, 'read', 'read_controls', &
                msg='Purge folder more then  = ' + purge_write)

        ! READ res_write_interval
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) res_write_interval
        call dmsg(5, 'read', 'read_controls', &
                msg='resnorm write interval  = ' + res_write_interval)

        ! READ tolerance
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) tolerance
        call dmsg(5, 'read', 'read_controls', &
                msg='Tolerance  = ' + tolerance)

        ! READ DEBUG_LEVEL
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) DEBUG_LEVEL
        call dmsg(5, 'read', 'read_controls', &
                msg='DEBUG_LEVEL = ' + DEBUG_LEVEL)

        close(CONTROL_FILE_UNIT)

      end subroutine read_controls


      subroutine read_scheme()
        !---------------------------------------------
        ! read fvscheme.md control file
        !--------------------------------------------
        implicit none
        character(len=STRING_BUFFER_LENGTH) :: buf
        integer                             :: ios

        call dmsg(1, 'read', 'read_scheme')

        open(SCHEME_FILE_UNIT, file=scheme_file, status='old', action='read')

        ! ignoring file header
        read(SCHEME_FILE_UNIT,*)
        read(SCHEME_FILE_UNIT,*)
        read(SCHEME_FILE_UNIT,*)
       
        ! read scheme name
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme_name
        call dmsg(5, 'read', 'read_scheme', &
                msg='scheme_name = ' + scheme_name)

        ! read interpolant
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) interpolant
        interpolant = trim(interpolant)
        call dmsg(5, 'read', 'read_scheme', &
                msg='interpolant = ' + interpolant)

        ! read ilimiter and PB switch
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) ilimiter_switch, PB_switch
        call dmsg(5, 'read', 'read_scheme', &
                msg='limiter switch = ' + ilimiter_switch )
        call dmsg(5, 'read', 'read_scheme', &
                  msg='PB switch = ' + PB_switch )

          ! read turbulence model
          call get_next_token(SCHEME_FILE_UNIT, buf)
          read(buf, *) turbulence
          call dmsg(5, 'read', 'read_scheme', &
                  msg='Turbulence Model = ' + turbulence)

        ! read time stepping method
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *, iostat=ios) time_stepping_method, global_time_step
        if (ios /= 0) then
            read(buf, *) time_stepping_method
            global_time_step = -1
        end if
        call dmsg(5, 'read', 'read_scheme', &
                msg='time_stepping_method = ' + time_stepping_method)
        call dmsg(5, 'read', 'read_scheme', &
                msg='global_time_step = ' + global_time_step)

        ! read time integration method
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) time_step_accuracy
        call dmsg(5, 'read', 'read_scheme', &
                msg='time_step_accuracy  = ' + time_step_accuracy)

        ! read higher order boundary
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) accur
        call dmsg(5, 'read', 'read_scheme', &
                msg='higher order boundary  = ' + accur)


        close(SCHEME_FILE_UNIT)

      end subroutine read_scheme

      subroutine read_flow()
        !---------------------------------------------
        ! read flow.md control file
        !--------------------------------------------
        implicit none

        character(len=STRING_BUFFER_LENGTH) :: buf

        call dmsg(1, 'read', 'read_flow')

        open(FLOW_FILE_UNIT, file=flow_file, status='old', action='read')

        ! ignoring file header
        read(FLOW_FILE_UNIT,*)
        read(FLOW_FILE_UNIT,*)
        read(FLOW_FILE_UNIT,*)
       
        ! read number of variable
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) n_var
        call dmsg(5, 'read', 'read_flow', &
                msg='Number of variables = ' + n_var)

        ! read rho_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_density
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_density = ' + free_stream_density)

        ! read u_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_x_speed
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_x_speed = ' + free_stream_x_speed)

        ! read v_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_y_speed
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_y_speed = ' + free_stream_y_speed)

        ! read w_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_z_speed
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_z_speed = ' + free_stream_z_speed)

        ! read P_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_pressure
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_pressure = ' + free_stream_pressure)

        ! read TKE_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_tk
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_TKE = ' + free_stream_tk)

        ! read omega_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_tw
        call dmsg(5, 'read', 'read_flow', &
                msg='free_stream_omega = ' + free_stream_tw)

        ! read reference viscosity
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) mu_ref
        call dmsg(5, 'read', 'read_flow', &
                msg='mu_reference = ' + mu_ref)

        ! Type of variation for viscosity
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) mu_variation
        call dmsg(5, 'read', 'read_flow', &
                msg='mu_variation = ' + mu_variation)

        ! read T_red
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) T_ref
        call dmsg(5, 'read', 'read_flow', &
                msg='T_reference = ' + T_ref)

        ! read Sutherland temp
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) Sutherland_temp
        call dmsg(5, 'read', 'read_flow', &
                msg='Sutherland temperature = ' + Sutherland_temp)

        ! read prandtl number
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) Pr
        call dmsg(5, 'read', 'read_flow', &
                msg='Prandtl Number = ' + Pr)

        ! read gamma
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) gm
        call dmsg(5, 'read', 'read_flow', &
                msg='gamma = ' + gm)

        ! read universal gas constant
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) R_gas
        call dmsg(5, 'read', 'read_flow', &
                msg='R_gas = ' + R_gas)
          

        close(FLOW_FILE_UNIT)

      end subroutine read_flow

        

      subroutine read_output_control()
        implicit none
        integer           :: i
        character(len=64) :: buf

        open(OUTIN_FILE_UNIT, file=outin_file, status='old', action='read')

        ! reading only counter first for dimension
        read(OUTIN_FILE_UNIT, *)
        read(OUTIN_FILE_UNIT, *)
        buf=" not } "
        n_write = 0
        do while (.true.)
          read(OUTIN_FILE_UNIT, *) buf
          if (trim(buf)=='}') EXIT
          n_write = n_write + 1
        end do

        allocate(rw_list(1:n_write))

        ! reading data types only, dumping n_write
        rewind(OUTIN_FILE_UNIT)
        read(OUTIN_FILE_UNIT, *)
        read(OUTIN_FILE_UNIT, *)
        buf=" "
        do i = 1,n_write
          read(OUTIN_FILE_UNIT, *) buf
          read(buf,*) rw_list(i)
        end do
        print*, rw_list

        close(OUTIN_FILE_UNIT)

      end subroutine read_output_control

end module read
