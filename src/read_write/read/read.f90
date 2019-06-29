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

  use global, only: CONTROL_FILE_UNIT
  use global, only:  SCHEME_FILE_UNIT
  use global, only:    FLOW_FILE_UNIT
  use global, only: control_file
  use global, only:  scheme_file
  use global, only:    flow_file
  use global, only: STRING_BUFFER_LENGTH
  use global, only: OUTIN_FILE_UNIT
  use global, only: outin_file
  use global, only: RES_CONTROL_FILE_UNIT
  use global, only: res_control_file

  use global_vars, only: CFL
  use global_vars, only: max_iters
  use global_vars, only: start_from
  use global_vars, only: checkpoint_iter
  use global_vars, only: res_write_interval
  use global_vars, only: write_file_format
  use global_vars, only: write_data_format
  use global_vars, only: read_file_format
  use global_vars, only: read_data_format
  use global_vars, only: write_percision
  use global_vars, only: purge_write
  use global_vars, only: tolerance
  use global_vars, only: tolerance_type
  use global_vars, only: process_id

  use global_vars, only: time_stepping_method
  use global_vars, only: time_step_accuracy
  use global_vars, only: global_time_step

  use global_vars, only: n_var
  use global_vars, only: free_stream_density
  use global_vars, only: free_stream_x_speed
  use global_vars, only: free_stream_y_speed
  use global_vars, only: free_stream_z_speed
  use global_vars, only: free_stream_pressure
  use global_vars, only: free_stream_tu
  use global_vars, only: free_stream_tgm
  use global_vars, only: mu_ratio_inf
  use global_vars, only: gm    !gamma
  use global_vars, only: R_gas !univarsal gas constant
  use global_vars, only: mu_ref !viscoity
  use global_vars, only: mu_variation !viscoity variation type
  use global_vars, only: T_ref
  use global_vars, only: Sutherland_temp
  use global_vars, only: Pr !prandtl number
  use global_vars, only: tPr !Turbulent prandtl number
  use global_vars, only: ilimiter_switch
  use global_vars, only: jlimiter_switch
  use global_vars, only: klimiter_switch
  use global_vars, only: itlimiter_switch
  use global_vars, only: jtlimiter_switch
  use global_vars, only: ktlimiter_switch
  use global_vars, only: iPB_switch
  use global_vars, only: jPB_switch
  use global_vars, only: kPB_switch
  use global_vars, only: accur
  
  use global_vars, only: interpolant
  use global_vars, only: scheme_name
  use global_vars, only: turbulence
  use global_vars, only: transition
  use global_vars, only: r_list
  use global_vars, only: w_list
  use global_vars, only: r_count
  use global_vars, only: w_count
  use global_vars, only: Res_list
  use global_vars, only: Res_count
  use utils      , only: DEBUG_LEVEL
  use utils      , only: dmsg
  use string
  use fclose     , only: close_file


  implicit none
  private

  public :: read_input_and_controls

    contains

      subroutine read_input_and_controls()
        !< Read all the input control files
        implicit none
        call read_controls()
        call read_scheme()
        call read_flow()
        call read_output_control()
        call read_Res_list()
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
        call dmsg(0, 'read', 'get_next_token', 'Returning: ' // trim(buf))

      end subroutine get_next_token



      subroutine read_controls()
        !< Read control.md file
        !---------------------------------------------
        implicit none
        character(len=STRING_BUFFER_LENGTH) :: buf

        DebugCall('read_controls')

        open(CONTROL_FILE_UNIT, file=control_file, status='old', action='read')

        !ignoring file header
        read(CONTROL_FILE_UNIT,*)
        read(CONTROL_FILE_UNIT,*)
        read(CONTROL_FILE_UNIT,*)

        ! READ CFL
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) CFL
        DebugInfo("CFL = "//trim(buf))

        ! READ start_from
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) start_from
        DebugInfo('Start from  level = '//trim(buf))

        ! READ max_iters
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) max_iters
        DebugInfo('Stop at iteration = '//trim(buf))

        ! READ checkpoint_iter
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) checkpoint_iter
        DebugInfo(' Solution write interval = '//trim(buf))

        ! READ write_file_format
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) write_file_format
        DebugInfo('Solution file format  = '//trim(buf))

        ! READ write_data_format
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) write_data_format
        DebugInfo('solution file data format = '//trim(buf))

        ! READ read_file_format
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) read_file_format
        DebugInfo('Restart file format  = '//trim(buf))

        ! READ_read data_format
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) read_data_format
        DebugInfo('Restart file data format = '//trim(buf))

        ! READ write_percision
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) write_percision
        DebugInfo('File write percision = '//trim(buf))

        ! READ purge_write
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) purge_write
        DebugInfo('Purge folder more then  = '//trim(buf))

        ! READ res_write_interval
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) res_write_interval
        DebugInfo('resnorm write interval  = '//trim(buf))

        ! READ tolerance
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) tolerance, tolerance_type
        DebugInfo(trim(tolerance_type)//' Tolerance  = '//trim(buf))

        ! READ DEBUG_LEVEL
        call get_next_token(CONTROL_FILE_UNIT, buf)
        read(buf, *) DEBUG_LEVEL
        DebugInfo('DEBUG_LEVEL = '//trim(buf))

        close(CONTROL_FILE_UNIT)

      end subroutine read_controls


      subroutine read_scheme()
        !< Read fvscheme.md control file
        !--------------------------------------------
        implicit none
        character(len=STRING_BUFFER_LENGTH) :: buf
        integer                             :: ios

        DebugCall('read_scheme')

        open(SCHEME_FILE_UNIT, file=scheme_file, status='old', action='read')

        ! ignoring file header
        read(SCHEME_FILE_UNIT,*)
        read(SCHEME_FILE_UNIT,*)
        read(SCHEME_FILE_UNIT,*)
       
        ! read scheme name
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) scheme_name
        DebugInfo('scheme_name = '//trim(buf))

        ! read interpolant
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) interpolant
        interpolant = trim(interpolant)
        DebugInfo('interpolant = '//trim(buf))

        ! read ilimiter and PB switch
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) ilimiter_switch,jlimiter_switch,klimiter_switch, &
                     iPB_switch, jPB_switch, kPB_switch
        DebugInfo('ilimiter switch = '//trim(buf) )
        DebugInfo('jlimiter switch = '//trim(buf) )
        DebugInfo('klimiter switch = '//trim(buf) )
          DebugInfo('PB switch = '//trim(buf) )

        ! read turbulent limiter
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) itlimiter_switch,jtlimiter_switch,ktlimiter_switch 
        DebugInfo('ilimiter switch = '//trim(buf) )
        DebugInfo('jlimiter switch = '//trim(buf) )
        DebugInfo('klimiter switch = '//trim(buf) )

        ! read turbulence model
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) turbulence
        DebugInfo('Turbulence Model = '//trim(buf))

        ! read transition model
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) transition
        DebugInfo('Transition Model = '//trim(buf))

        ! read time stepping method
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *, iostat=ios) time_stepping_method, global_time_step
        if (ios /= 0) then
            read(buf, *) time_stepping_method
            global_time_step = -1
        end if
        DebugInfo('time_stepping_method = '//trim(buf))
        DebugInfo('global_time_step = '//trim(buf))

        ! read time integration method
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) time_step_accuracy
        DebugInfo('time_step_accuracy  = '//trim(buf))

        ! read higher order boundary
        call get_next_token(SCHEME_FILE_UNIT, buf)
        read(buf, *) accur
        DebugInfo('higher order boundary  = '//trim(buf))


        close(SCHEME_FILE_UNIT)

      end subroutine read_scheme

      subroutine read_flow()
        !< Read flow.md control file
        !--------------------------------------------
        implicit none

        character(len=STRING_BUFFER_LENGTH) :: buf

        DebugCall('read_flow')

        open(FLOW_FILE_UNIT, file=flow_file, status='old', action='read')

        ! ignoring file header
        read(FLOW_FILE_UNIT,*)
        read(FLOW_FILE_UNIT,*)
        read(FLOW_FILE_UNIT,*)
       
        ! read number of variable
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) n_var
        DebugInfo('Number of variables = '//trim(buf))

        ! read rho_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_density
        DebugInfo('free_stream_density = '//trim(buf))

        ! read u_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_x_speed
        DebugInfo('free_stream_x_speed = '//trim(buf))

        ! read v_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_y_speed
        DebugInfo('free_stream_y_speed = '//trim(buf))

        ! read w_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_z_speed
        DebugInfo('free_stream_z_speed = '//trim(buf))

        ! read P_inf
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_pressure
        DebugInfo('free_stream_pressure = '//trim(buf))

        ! read turbulence intensity in percentage
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_tu
        DebugInfo('free_stream_Turb_intensity = '//trim(buf))

        ! read viscosity ratio
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) mu_ratio_inf
        DebugInfo('free_stream_mu_ratio = '//trim(buf))

        ! read intermittency
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) free_stream_tgm
        DebugInfo('free_stream_Intermittency = '//trim(buf))

        ! read reference viscosity
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) mu_ref
        DebugInfo('mu_reference = '//trim(buf))

        ! Type of variation for viscosity
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) mu_variation
        DebugInfo('mu_variation = '//trim(buf))

        ! read T_red
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) T_ref
        DebugInfo('T_reference = '//trim(buf))

        ! read Sutherland temp
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) Sutherland_temp
        DebugInfo('Sutherland temperature = '//trim(buf))

        ! read prandtl number
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) Pr, tPr
        DebugInfo('Prandtl Number = '//trim(buf))

        ! read gamma
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) gm
        DebugInfo('gamma = '//trim(buf))

        ! read universal gas constant
        call get_next_token(FLOW_FILE_UNIT, buf)
        read(buf, *) R_gas
        DebugInfo('R_gas = '//trim(buf))
          

        close(FLOW_FILE_UNIT)

      end subroutine read_flow

        

      subroutine read_output_control()
        !< Read output_contorl.md file
        implicit none
        integer           :: i
        character(len=64) :: buf
        integer :: ios
        
        call get_rw_count()
        call close_file(OUTIN_FILE_UNIT)
        open(OUTIN_FILE_UNIT, file=outin_file, status='old', action='read')

        ! variables to write
        do while(.true.)
          read(OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        do i = 1,w_count
          read(OUTIN_FILE_UNIT, *) buf
          read(buf,*) w_list(i)
        end do

        ! restart variables to read
        do while(.true.)
          read(OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        do i = 1,r_count
          read(OUTIN_FILE_UNIT, *) buf
          read(buf,*) r_list(i)
        end do
        if(r_count==0) r_list=w_list

        close(OUTIN_FILE_UNIT)

      end subroutine read_output_control

      subroutine get_rw_count()
        !< Get read/write count
        implicit none
        integer :: ios
        character(len=64) :: buf

        r_count=0
        w_count=0
        call close_file(OUTIN_FILE_UNIT)
        open(OUTIN_FILE_UNIT, file=outin_file, status='old', action='read')

        ! write list dimension
        do while(.true.)
          read(OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        w_count = 0
        do while (.true.)
          read(OUTIN_FILE_UNIT, *, iostat=ios) buf
          if (trim(buf)=='}') EXIT
          if(is_iostat_end(ios)) EXIT
          w_count = w_count + 1
        end do

        if(w_count>0) then
          allocate(w_list(1:w_count))
        else
          w_count=3
          allocate(w_list(1:w_count))
          w_list(1) = "Velocity"
          w_list(2) = "Density"
          w_list(3) = "Pressure"
        end if

        ! read list dimesnion 
        do while(.true.)
          read(OUTIN_FILE_UNIT, *, iostat=ios) buf
          if(trim(buf)=='{') EXIT
          if(is_iostat_end(ios)) EXIT
        end do
        r_count = 0
        do while (.true.)
          read(OUTIN_FILE_UNIT, *, iostat=ios) buf
          if (trim(buf)=='}') EXIT
          if(is_iostat_end(ios)) EXIT
          r_count = r_count + 1
        end do
        if(r_count==0) then
          allocate(r_list(1:w_count))
        else
          allocate(r_list(1:r_count))
        end if

        close(OUTIN_FILE_UNIT)

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


      subroutine read_Res_list()
        !< Read Residual file: res_control.md
        implicit none
        integer           :: i
        integer           :: skip

        open(RES_CONTROL_FILE_UNIT, file=res_control_file, status='old', action='read')
        call get_count_within_braces(RES_CONTROL_FILE_UNIT, Res_count)
        call close_file(RES_CONTROL_FILE_UNIT)

        open(RES_CONTROL_FILE_UNIT, file=res_control_file, status='old', action='read')
        ! skipping line
        skip  = get_number_of_line('{', RES_CONTROL_FILE_UNIT)

        !reading vaules
        if(Res_count==0)then
          allocate(Res_list(1:2))
          Res_count=2
          Res_list(1)="Mass_abs"
          Res_list(2)="Resnorm_abs"
        else
          allocate(Res_list(1:Res_count))
        end if
        do i = 1,Res_count
          read(RES_CONTROL_FILE_UNIT, *) Res_list(i)
        end do


        call close_file(RES_CONTROL_FILE_UNIT)

      end subroutine read_Res_list

end module read
