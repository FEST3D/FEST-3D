  !< This module read input control files which include:
  !<   1. control.md
  !<   2. fvscheme.md
  !<   3. flow.md
  !<   4. res_control.md
  !<   5. state_read_write_control.md
module read
  !< This module read input control files which include:
  !<   1. control.md
  !<   2. fvscheme.md
  !<   3. flow.md
  !<   4. res_control.md
  !<   5. state_read_write_control.md
  !------------------------------------------------------

#include "../../debug.h"
  use vartypes
  implicit none
  private

  public :: read_input_and_controls

    contains

      subroutine read_input_and_controls(files, control, scheme, flow)
        !< Read all the input control files
        implicit none
        type(filetype), intent(in) :: files
        type(controltype), intent(inout) :: control
        type(schemetype), intent(inout) :: scheme
        type(flowtype), intent(inout) :: flow
        call read_controls(files, control)
        call read_scheme(files, scheme)
        call read_flow(files, control, flow)
        call read_output_control(files, control)
        call read_Res_list(files, control)
      end subroutine read_input_and_controls


      subroutine get_next_token(token_file_unit, buf)
        !< Extract the next token from the config file
        !<
        !< Each token is on a separate line.
        !< There may be multiple comments (lines beginning with #) 
        !< and blank lines in between.
        !< The purpose of this subroutine is to ignore all these 
        !< lines and return the next "useful" line.
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

      end subroutine get_next_token



      subroutine read_controls(files, control)
        !< Read control.md file
        !---------------------------------------------
        implicit none
        type(filetype), intent(in) :: files
        type(controltype), intent(inout) :: control
        character(len=STRING_BUFFER_LENGTH) :: buf

        DebugCall('read_controls')

        open(files%CONTROL_FILE_UNIT, file=files%control_file, status='old', action='read')

        !ignoring file header
        read(files%CONTROL_FILE_UNIT,*)
        read(files%CONTROL_FILE_UNIT,*)
        read(files%CONTROL_FILE_UNIT,*)

        ! READ CFL
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%CFL
        DebugInfo("CFL = "//trim(buf))

        ! READ start_from
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%start_from
        DebugInfo('Start from  level = '//trim(buf))

        ! READ max_iters
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%max_iters
        DebugInfo('Stop at iteration = '//trim(buf))

        ! READ checkpoint_iter
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%checkpoint_iter
        DebugInfo(' Solution write interval = '//trim(buf))

        ! READ write_file_format
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%write_file_format
        DebugInfo('Solution file format  = '//trim(buf))

        ! READ write_data_format
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%write_data_format
        DebugInfo('solution file data format = '//trim(buf))

        ! READ read_file_format
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%read_file_format
        DebugInfo('Restart file format  = '//trim(buf))

        ! READ_read data_format
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%read_data_format
        DebugInfo('Restart file data format = '//trim(buf))

        ! READ write_percision
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%write_percision
        DebugInfo('File write percision = '//trim(buf))

        ! READ purge_write
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%purge_write
        DebugInfo('Purge folder more then  = '//trim(buf))

        ! READ res_write_interval
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%res_write_interval
        DebugInfo('resnorm write interval  = '//trim(buf))

        ! READ tolerance
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%tolerance, control%tolerance_type
        DebugInfo(trim(control%tolerance_type)//' Tolerance  = '//trim(buf))

        ! READ DEBUG_LEVEL
        call get_next_token(files%CONTROL_FILE_UNIT, buf)
        read(buf, *) control%DEBUG_LEVEL
        DebugInfo('DEBUG_LEVEL = '//trim(buf))

        close(files%CONTROL_FILE_UNIT)

      end subroutine read_controls


      subroutine read_scheme(files, scheme)
        !< Read fvscheme.md control file
        !--------------------------------------------
        implicit none
        type(filetype), intent(in) ::  files
        type(schemetype), intent(inout) :: scheme
        character(len=STRING_BUFFER_LENGTH) :: buf
        integer                             :: ios

        DebugCall('read_scheme')

        open(files%SCHEME_FILE_UNIT, file=files%scheme_file, status='old', action='read')

        ! ignoring file header
        read(files%SCHEME_FILE_UNIT,*)
        read(files%SCHEME_FILE_UNIT,*)
        read(files%SCHEME_FILE_UNIT,*)
       
        ! read scheme name
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%scheme_name
        DebugInfo('scheme_name = '//trim(buf))

        ! read interpolant
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%interpolant
        scheme%interpolant = trim(scheme%interpolant)
        DebugInfo('interpolant = '//trim(buf))

        ! read ilimiter and PB switch
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%ilimiter_switch, scheme%jlimiter_switch, scheme%klimiter_switch, &
                     scheme%iPB_switch, scheme%jPB_switch, scheme%kPB_switch
        DebugInfo('ilimiter switch = '//trim(buf) )
        DebugInfo('jlimiter switch = '//trim(buf) )
        DebugInfo('klimiter switch = '//trim(buf) )
          DebugInfo('PB switch = '//trim(buf) )

        ! read turbulent limiter
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%itlimiter_switch,scheme%jtlimiter_switch,scheme%ktlimiter_switch 
        DebugInfo('ilimiter switch = '//trim(buf) )
        DebugInfo('jlimiter switch = '//trim(buf) )
        DebugInfo('klimiter switch = '//trim(buf) )

        ! read turbulence model
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%turbulence
        DebugInfo('Turbulence Model = '//trim(buf))

        ! read transition model
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%transition
        DebugInfo('Transition Model = '//trim(buf))

        ! read time stepping method
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *, iostat=ios) scheme%time_stepping_method, scheme%global_time_step
        if (ios /= 0) then
            read(buf, *) scheme%time_stepping_method
            scheme%global_time_step = -1
        end if
        DebugInfo('time_stepping_method = '//trim(buf))
        DebugInfo('global_time_step = '//trim(buf))

        ! read time integration method
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%time_step_accuracy
        DebugInfo('time_step_accuracy  = '//trim(buf))

        ! read higher order boundary
        call get_next_token(files%SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme%accur
        DebugInfo('higher order boundary  = '//trim(buf))


        close(files%SCHEME_FILE_UNIT)

      end subroutine read_scheme

      subroutine read_flow(files, control, flow)
        !< Read flow.md control file
        !--------------------------------------------
        implicit none
        type(filetype), intent(in) :: files
        type(controltype), intent(inout) :: control
        type(flowtype), intent(inout) :: flow

        character(len=STRING_BUFFER_LENGTH) :: buf

        DebugCall('read_flow')

        open(files%FLOW_FILE_UNIT, file=files%flow_file, status='old', action='read')

        ! ignoring file header
        read(files%FLOW_FILE_UNIT,*)
        read(files%FLOW_FILE_UNIT,*)
        read(files%FLOW_FILE_UNIT,*)
       
        ! read number of variable
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) control%n_var
        DebugInfo('Number of variables = '//trim(buf))

        ! read rho_inf
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%density_inf
        DebugInfo('free_stream_density = '//trim(buf))

        ! read u_inf
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%x_speed_inf
        DebugInfo('free_stream_x_speed = '//trim(buf))

        ! read v_inf
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%y_speed_inf
        DebugInfo('free_stream_y_speed = '//trim(buf))

        ! read w_inf
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%z_speed_inf
        DebugInfo('free_stream_z_speed = '//trim(buf))

        ! read P_inf
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%pressure_inf
        DebugInfo('free_stream_pressure = '//trim(buf))

        ! read turbulence intensity in percentage
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%tu_inf
        DebugInfo('free_stream_Turb_intensity = '//trim(buf))

        ! read viscosity ratio
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%mu_ratio_inf
        DebugInfo('free_stream_mu_ratio = '//trim(buf))

        ! read intermittency
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%tgm_inf
        DebugInfo('free_stream_Intermittency = '//trim(buf))

        ! read reference viscosity
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%mu_ref
        DebugInfo('mu_reference = '//trim(buf))

        ! Type of variation for viscosity
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%mu_variation
        DebugInfo('mu_variation = '//trim(buf))

        ! read T_red
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%T_ref
        DebugInfo('T_reference = '//trim(buf))

        ! read Sutherland temp
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%Sutherland_temp
        DebugInfo('Sutherland temperature = '//trim(buf))

        ! read prandtl number
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%Pr, flow%tPr
        DebugInfo('Prandtl Number = '//trim(buf))

        ! read gamma
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%gm
        DebugInfo('gamma = '//trim(buf))

        ! read universal gas constant
        call get_next_token(files%FLOW_FILE_UNIT, buf)
        read(buf, *) flow%R_gas
        DebugInfo('R_gas = '//trim(buf))
          

        close(files%FLOW_FILE_UNIT)

      end subroutine read_flow

        

      subroutine read_output_control(files, control)
        !< Read output_contorl.md file
        implicit none
        type(filetype), intent(in) :: files
        type(controltype), intent(inout) :: control
        integer           :: i
        character(len=64) :: buf
        integer :: ios
        logical :: ok
        
        call get_rw_count(files, control)
        inquire(files%OUTIN_FILE_UNIT, opened=ok)
        if(ok)  close(files%OUTIN_FILE_UNIT)
        !call close_file(files%OUTIN_FILE_UNIT)
        open(files%OUTIN_FILE_UNIT, file=files%outin_file, status='old', action='read')

        ! variables to write
        do while(.true.)
          read(files%OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        do i = 1,control%w_count
          read(files%OUTIN_FILE_UNIT, *) buf
          read(buf,*) control%w_list(i)
        end do

        ! restart variables to read
        do while(.true.)
          read(files%OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        do i = 1,control%r_count
          read(files%OUTIN_FILE_UNIT, *) buf
          read(buf,*) control%r_list(i)
        end do
        if(control%r_count==0) control%r_list=control%w_list

        close(files%OUTIN_FILE_UNIT)

      end subroutine read_output_control

      subroutine get_rw_count(files, control)
        !< Get read/write count
        implicit none
        type(filetype), intent(in) :: files
        type(controltype), intent(inout) :: control
        integer :: ios
        character(len=64) :: buf
        logical :: ok

        control%r_count=0
        control%w_count=0
        inquire(files%OUTIN_FILE_UNIT, opened=ok)
        if(ok)  close(files%OUTIN_FILE_UNIT)
        !call close_file(files%OUTIN_FILE_UNIT)
        open(files%OUTIN_FILE_UNIT, file=files%outin_file, status='old', action='read')

        ! write list dimension
        do while(.true.)
          read(files%OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        control%w_count = 0
        do while (.true.)
          read(files%OUTIN_FILE_UNIT, *, iostat=ios) buf
          if (trim(buf)=='}') EXIT
          if(is_iostat_end(ios)) EXIT
          control%w_count = control%w_count + 1
        end do

        if(control%w_count>0) then
          allocate(control%w_list(1:control%w_count))
        else
          control%w_count=3
          allocate(control%w_list(1:control%w_count))
          control%w_list(1) = "Velocity"
          control%w_list(2) = "Density"
          control%w_list(3) = "Pressure"
        end if

        ! read list dimesnion 
        do while(.true.)
          read(files%OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        control%r_count = 0
        do while (.true.)
          read(files%OUTIN_FILE_UNIT, *, iostat=ios) buf
          if (trim(buf)=='}') EXIT
          if(is_iostat_end(ios)) EXIT
          control%r_count = control%r_count + 1
        end do
        if(control%r_count==0) then
          allocate(control%r_list(1:control%w_count))
        else
          allocate(control%r_list(1:control%r_count))
        end if

        close(files%OUTIN_FILE_UNIT)

      end subroutine get_rw_count


      subroutine get_count_within_braces(handler, count)
        !< Get number of variables between two curly braces
        implicit none
        integer, intent(in) :: handler
        !< File handler from which list number is extracted
        integer, intent(out) :: count
        !< Extracted count
        integer :: skip

        ! skipping lines outside braces
        skip  = get_number_of_line('{', handler)
        ! finding actual count if any
        count = get_number_of_line('}', handler)

      end subroutine get_count_within_braces


      function get_number_of_line(till, infile) result(number)
        !< Get number of lines till some character like "#"
        implicit none
        integer          ,intent(in) :: infile
        character(len= 1),intent(in) :: till
        character(len=64)   :: buf
        integer             :: ios
        integer             :: number
        number=0
        do while(.true.)
          read(infile, *, iostat=ios) buf
          if(trim(buf)==till) EXIT
          if(is_iostat_end(ios)) EXIT
          number = number + 1
        end do
      end function get_number_of_line


      subroutine read_Res_list(files, control)
        !< Read Residual file: res_control.md
        implicit none
        type(filetype), intent(in) :: files
        type(controltype), intent(inout) :: control
        integer           :: i
        integer           :: skip
        logical :: ok

        open(files%RES_CONTROL_FILE_UNIT, file=files%res_control_file, status='old', action='read')
        call get_count_within_braces(files%RES_CONTROL_FILE_UNIT, control%Res_count)
        !call close_file(files%RES_CONTROL_FILE_UNIT)
        inquire(files%RES_CONTROL_FILE_UNIT, opened=ok)
        if(ok)  close(files%RES_CONTROL_FILE_UNIT)

        open(files%RES_CONTROL_FILE_UNIT, file=files%res_control_file, status='old', action='read')
        ! skipping line
        skip  = get_number_of_line('{',files%RES_CONTROL_FILE_UNIT)

        !reading vaules
        if(control%Res_count==0)then
          allocate(control%Res_list(1:2))
          control%Res_count=2
          control%Res_list(1)="Mass_abs"
          control%Res_list(2)="Resnorm_abs"
        else
          allocate(control%Res_list(1:control%Res_count))
        end if
        do i = 1,control%Res_count
          read(files%RES_CONTROL_FILE_UNIT, *) control%Res_list(i)
        end do

        close(files%RES_CONTROL_FILE_UNIT)

      end subroutine read_Res_list

end module read
