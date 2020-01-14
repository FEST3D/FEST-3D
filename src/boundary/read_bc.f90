  !< Get all the fixed values from the bc_**.md file 
module read_bc
  !< Get all the fixed values from the bc_**.md file 
  !-----------------------------------------------------
#include "../error.h"
  use vartypes
  implicit none
  private

  character(len=STRING_BUFFER_LENGTH) :: buf
  !< String to extract single line from the file
  public :: read_fixed_values


  contains

    subroutine read_fixed_values(files, scheme, flow, bc)
      !< Read fixed values for each block face
      implicit none
      type(filetype), intent(in) :: files
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) ::flow
      type(boundarytype), intent(inout) :: bc
      integer :: count=0

      call fill_fixed_values(scheme, flow, bc)

      open(unit=files%BOUNDARY_CONDITIONS_FILE_UNIT, file=files%bcfile)
            read(files%BOUNDARY_CONDITIONS_FILE_UNIT, *)
            read(files%BOUNDARY_CONDITIONS_FILE_UNIT, *)
            read(files%BOUNDARY_CONDITIONS_FILE_UNIT, *)
      do while(count<6)
        read(files%BOUNDARY_CONDITIONS_FILE_UNIT, "(A)") buf
        if(buf(1:1)=='#')then
          count=count+1
          call get_fixed_values(files, scheme,flow, bc, count)
        end if
      end do
      close(files%BOUNDARY_CONDITIONS_FILE_UNIT)

    end subroutine read_fixed_values

    subroutine get_fixed_values(files, scheme, flow, bc, count)
      !< Extract fixed value from the bc_**.md file
      implicit none
      type(filetype), intent(in) :: files
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(boundarytype), intent(inout) :: bc
      integer, intent(in) :: count
      real(wp) :: fix_val
      integer :: ios
      do while(.true.)
        read(files%BOUNDARY_CONDITIONS_FILE_UNIT,"(A)") buf
        if(buf(1:2) == '- ') then 
          read(buf(index(buf(3:), ' ')+3:), *, iostat=ios) fix_val
          select case(buf(3:index(buf(3:), " ")+1))
            case ("FIX_DENSITY")
              call set_value(bc%fixed_density , fix_val, flow%density_inf , count, ios)

            case ("FIX_X_SPEED")
              call set_value(bc%fixed_x_speed , fix_val, flow%x_speed_inf , count, ios)

            case ("FIX_Y_SPEED")
              call set_value(bc%fixed_y_speed , fix_val, flow%y_speed_inf , count, ios)

            case ("FIX_Z_SPEED")
              call set_value(bc%fixed_z_speed , fix_val, flow%z_speed_inf , count, ios)

            case ("FIX_PRESSURE")
              call set_value(bc%fixed_pressure, fix_val, flow%pressure_inf, count, ios)

            case ("WALL_TEMPERATURE")
              call set_value(bc%fixed_wall_temperature, fix_val, 0.0, count, ios)

            case ("TOTAL_TEMPERATURE")
              call set_value(bc%fixed_Ttemperature, fix_val, 0.0, count, ios)

            case ("TOTAL_PRESSURE")
              call set_value(bc%fixed_Tpressure, fix_val, 0.0, count, ios)

          end select

          select case (scheme%turbulence)

            case ("none")
              !do nothing
              continue

            case ("sst", 'tw', 'sst2003')
              select case(buf(3:index(buf(3:), " ")+1))
                case ("FIX_tk")
                  call set_value(bc%fixed_tk      , fix_val, flow%tk_inf      , count, ios)
                case ("FIX_tw")
                  call set_value(bc%fixed_tw      , fix_val, flow%tw_inf      , count, ios)
                case DEFAULT
                  ! no a value to fix
                  continue
              end select
              
            case ("kkl")
              select case(buf(3:index(buf(3:), " ")+1))
                case ("FIX_tk")
                  call set_value(bc%fixed_tk      , fix_val, flow%tk_inf      , count, ios)
                case ("FIX_tkl")
                  call set_value(bc%fixed_tkl     , fix_val, flow%tkl_inf     , count, ios)
                case DEFAULT
                  ! no a value to fix
                  continue
              end select

            case ("sa", "saBC")
              select case(buf(3:index(buf(3:), " ")+1))
                case ("FIX_tv")
                  call set_value(bc%fixed_tk      , fix_val, flow%tv_inf      , count, ios)
                case DEFAULT
                  ! no a value to fix
                  continue
              end select

          end select

        else
          exit
        end if
      end do

    end subroutine get_fixed_values

    
    subroutine fill_fixed_values(scheme, flow, bc)
      !< Fill the Fixed_var array with with free-stream value
      !< or default values.
      implicit none
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(boundarytype), intent(inout) :: bc
      integer :: count
      integer :: ios=-1

      do count = 1,6
        !case ("FIX_DENSITY")
          call set_value(bc%fixed_density , flow%density_inf, flow%density_inf , count, ios)

        !case ("FIX_X_SPEED")
          call set_value(bc%fixed_x_speed , flow%x_speed_inf, flow%x_speed_inf , count, ios)

        !case ("FIX_Y_SPEED")
          call set_value(bc%fixed_y_speed , flow%y_speed_inf, flow%y_speed_inf , count, ios)

        !case ("FIX_Z_SPEED")
          call set_value(bc%fixed_z_speed , flow%z_speed_inf, flow%z_speed_inf , count, ios)

        !case ("FIX_PRESSURE")
          call set_value(bc%fixed_pressure, flow%pressure_inf, flow%pressure_inf, count, ios)

        !case ("WALL_TEMPERATURE")
          call set_value(bc%fixed_wall_temperature, 0.0, 0.0, count, ios)

        !case ("TOTAL_TEMPERATURE")
          call set_value(bc%fixed_Ttemperature, 0.0, 0.0, count, ios)

        !case ("TOTAL_PRESSURE")
          call set_value(bc%fixed_Tpressure, 0.0, 0.0, count, ios)


        select case (scheme%turbulence)

          case ("none")
            !do nothing
            continue

          case ("sst", 'tw', 'sst2003')
            !case ("FIX_tk")
              call set_value(bc%fixed_tk      , flow%tk_inf, flow%tk_inf      , count, ios)
            !case ("FIX_tw")
              call set_value(bc%fixed_tw      , flow%tw_inf, flow%tw_inf      , count, ios)
            
          case ("kkl")
            !case ("FIX_tk")
              call set_value(bc%fixed_tk      , flow%tk_inf, flow%tk_inf      , count, ios)
            !case ("FIX_tkl")
              call set_value(bc%fixed_tkl     , flow%tkl_inf, flow%tkl_inf     , count, ios)

          case ("sa", "saBC")
            !case ("FIX_tv")
              call set_value(bc%fixed_tk      , flow%tv_inf, flow%tv_inf      , count, ios)

          case DEFAULT
            Fatal_error

        end select
      end do

    end subroutine fill_fixed_values



    subroutine set_value(fixed_var, fix_val, inf_val, count, ios)
      !< Set particular value to the Fixed_var variable
      implicit none
      integer, intent(in) :: ios
      integer, intent(in) :: count
      real(wp)   , intent(in) :: fix_val
      real(wp)   , intent(in) :: inf_val
      real(wp)   , intent(out), dimension(:) :: fixed_var
      if(ios==0)then
        fixed_var(count) = fix_val
      else 
        fixed_var(count) = inf_val
      end if
    end subroutine set_value
end module read_bc
