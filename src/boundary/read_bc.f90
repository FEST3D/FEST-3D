  !< Get all the fixed values from the bc_**.md file 
module read_bc
  !< Get all the fixed values from the bc_**.md file 
  !-----------------------------------------------------
  ! 170516  Jatinder Pal Singh Sandhu
  ! Aim : get all the fixed valed from bc_**.md file
  !-----------------------------------------------------
#include "../error.inc"
  use global     , only: BOUNDARY_CONDITIONS_FILE_UNIT
  use global     , only: STRING_BUFFER_LENGTH
  use global_vars, only: density_inf
  use global_vars, only: x_speed_inf
  use global_vars, only: y_speed_inf
  use global_vars, only: z_speed_inf
  use global_vars, only: pressure_inf
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only: tv_inf
  use global_vars, only: te_inf
  use global_vars, only: tkl_inf
  use global_vars, only: fixed_density
  use global_vars, only: fixed_x_speed
  use global_vars, only: fixed_y_speed
  use global_vars, only: fixed_z_speed
  use global_vars, only: fixed_pressure
  use global_vars, only: fixed_tk
  use global_vars, only: fixed_tw
  use global_vars, only: fixed_te
  use global_vars, only: fixed_tv
  use global_vars, only: fixed_tkl
  use global_vars, only: fixed_Tpressure
  use global_vars, only: fixed_Ttemperature
  use global_vars, only: fixed_wall_temperature
  use global_vars, only: process_id
  use global_vars, only: turbulence
  use layout     , only: bc_file

  implicit none
  private

  character(len=STRING_BUFFER_LENGTH) :: buf
  !< String to extract single line from the file
  public :: read_fixed_values


  contains

    subroutine read_fixed_values()
      !< Read fixed values for each block face
      implicit none
      integer :: count=0

      call fill_fixed_values()

      open(unit=BOUNDARY_CONDITIONS_FILE_UNIT, file=bc_file)
            read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
            read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
            read(BOUNDARY_CONDITIONS_FILE_UNIT, *)
      do while(count<6)
        read(BOUNDARY_CONDITIONS_FILE_UNIT, "(A)") buf
        if(buf(1:1)=='#')then
          count=count+1
          call get_fixed_values(count)
        end if
      end do
      close(BOUNDARY_CONDITIONS_FILE_UNIT)

    end subroutine read_fixed_values

    subroutine get_fixed_values(count)
      !< Extract fixed value from the bc_**.md file
      implicit none
      integer, intent(in) :: count
      real :: fix_val
      integer :: ios
      do while(.true.)
        read(BOUNDARY_CONDITIONS_FILE_UNIT,"(A)") buf
        if(buf(1:2) == '- ') then 
          read(buf(index(buf(3:), ' ')+3:), *, iostat=ios) fix_val
          select case(buf(3:index(buf(3:), " ")+1))
            case ("FIX_DENSITY")
              call set_value(fixed_density , fix_val, density_inf , count, ios)

            case ("FIX_X_SPEED")
              call set_value(fixed_x_speed , fix_val, x_speed_inf , count, ios)

            case ("FIX_Y_SPEED")
              call set_value(fixed_y_speed , fix_val, y_speed_inf , count, ios)

            case ("FIX_Z_SPEED")
              call set_value(fixed_z_speed , fix_val, z_speed_inf , count, ios)

            case ("FIX_PRESSURE")
              call set_value(fixed_pressure, fix_val, pressure_inf, count, ios)

            case ("WALL_TEMPERATURE")
              call set_value(fixed_wall_temperature, fix_val, 0.0, count, ios)

            case ("TOTAL_TEMPERATURE")
              call set_value(fixed_Ttemperature, fix_val, 0.0, count, ios)

            case ("TOTAL_PRESSURE")
              call set_value(fixed_Tpressure, fix_val, 0.0, count, ios)

          end select

          select case (turbulence)

            case ("none")
              !do nothing
              continue

            case ("sst", 'tw', 'sst2003')
              select case(buf(3:index(buf(3:), " ")+1))
                case ("FIX_tk")
                  call set_value(fixed_tk      , fix_val, tk_inf      , count, ios)
                case ("FIX_tw")
                  call set_value(fixed_tw      , fix_val, tw_inf      , count, ios)
                case DEFAULT
                  ! no a value to fix
                  continue
              end select
              
            case ("kkl")
              select case(buf(3:index(buf(3:), " ")+1))
                case ("FIX_tk")
                  call set_value(fixed_tk      , fix_val, tk_inf      , count, ios)
                case ("FIX_tkl")
                  call set_value(fixed_tkl     , fix_val, tkl_inf     , count, ios)
                case DEFAULT
                  ! no a value to fix
                  continue
              end select

            case ("sa", "saBC")
              select case(buf(3:index(buf(3:), " ")+1))
                case ("FIX_tv")
                  call set_value(fixed_tk      , fix_val, tv_inf      , count, ios)
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

    
    subroutine fill_fixed_values()
      !< fill the Fixed_var array with with free-stream value
      !< or default values.
      implicit none
      integer :: count
      integer :: ios=-1

      do count = 1,6
        !case ("FIX_DENSITY")
          call set_value(fixed_density , density_inf, density_inf , count, ios)

        !case ("FIX_X_SPEED")
          call set_value(fixed_x_speed , x_speed_inf, x_speed_inf , count, ios)

        !case ("FIX_Y_SPEED")
          call set_value(fixed_y_speed , y_speed_inf, y_speed_inf , count, ios)

        !case ("FIX_Z_SPEED")
          call set_value(fixed_z_speed , z_speed_inf, z_speed_inf , count, ios)

        !case ("FIX_PRESSURE")
          call set_value(fixed_pressure, pressure_inf, pressure_inf, count, ios)

        !case ("WALL_TEMPERATURE")
          call set_value(fixed_wall_temperature, 0.0, 0.0, count, ios)

        !case ("TOTAL_TEMPERATURE")
          call set_value(fixed_Ttemperature, 0.0, 0.0, count, ios)

        !case ("TOTAL_PRESSURE")
          call set_value(fixed_Tpressure, 0.0, 0.0, count, ios)


        select case (turbulence)

          case ("none")
            !do nothing
            continue

          case ("sst", 'tw', 'sst2003')
            !case ("FIX_tk")
              call set_value(fixed_tk      , tk_inf, tk_inf      , count, ios)
            !case ("FIX_tw")
              call set_value(fixed_tw      , tw_inf, tw_inf      , count, ios)
            
          case ("kkl")
            !case ("FIX_tk")
              call set_value(fixed_tk      , tk_inf, tk_inf      , count, ios)
            !case ("FIX_tkl")
              call set_value(fixed_tkl     , tkl_inf, tkl_inf     , count, ios)

          case ("sa", "saBC")
            !case ("FIX_tv")
              call set_value(fixed_tk      , tv_inf, tv_inf      , count, ios)

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
      real   , intent(in) :: fix_val
      real   , intent(in) :: inf_val
      real   , intent(out), dimension(:) :: fixed_var
      if(ios==0)then
        fixed_var(count) = fix_val
      else 
        fixed_var(count) = inf_val
      end if
    end subroutine set_value
end module read_bc
