module read_bc
  !-----------------------------------------------------
  ! 170516  Jatinder Pal Singh Sandhu
  ! Aim : get all the fixed valed from bc_**.md file
  !-----------------------------------------------------
  use global     , only: BOUNDARY_CONDITIONS_FILE_UNIT
  use global     , only: STRING_BUFFER_LENGTH
  use global_vars, only: density_inf
  use global_vars, only: x_speed_inf
  use global_vars, only: y_speed_inf
  use global_vars, only: z_speed_inf
  use global_vars, only: pressure_inf
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only: fixed_density
  use global_vars, only: fixed_x_speed
  use global_vars, only: fixed_y_speed
  use global_vars, only: fixed_z_speed
  use global_vars, only: fixed_pressure
  use global_vars, only: fixed_tk
  use global_vars, only: fixed_tw
  use global_vars, only: process_id
  use layout     , only: bc_file

  implicit none
  private

  character(len=STRING_BUFFER_LENGTH) :: buf
  public :: read_fixed_values


  contains

    subroutine read_fixed_values
      implicit none
      integer :: count=0

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

            case ("FIX_tk")
              call set_value(fixed_tk      , fix_val, tk_inf      , count, ios)

            case ("FIX_tw")
              call set_value(fixed_tw      , fix_val, tw_inf      , count, ios)

            case DEFAULT
              ! no a value to fix
              continue

          end select

        else
          exit
        end if
      end do
    end subroutine get_fixed_values

    subroutine set_value(fixed_var, fix_val, inf_val, count, ios)
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
