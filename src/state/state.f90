  !< Allocate memory to the state variables and initialize them
module state
  !< Allocate memory to the state variables and initialize them
  !< The state of the system is defined using the density, velocity and
  !< pressure (primitive variables qp), and trubulent and transition
  !< variables at the cell-center points.
    !-------------------------------------------------------------------
    ! The state module contains the state variables and the methods that
    ! act on them. 
    !-------------------------------------------------------------------
    
#include "../debug.h"
#include "../error.h"

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
    use global_vars, only : sst_n_var
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
    use global_vars, only : te
    use global_vars, only : tv
    use global_vars, only : tkl
    use global_vars, only : tgm
    use global_vars, only : tk_inf
    use global_vars, only : tw_inf
    use global_vars, only : te_inf
    use global_vars, only : tv_inf
    use global_vars, only : tkl_inf
    use global_vars, only : tgm_inf
    use global_vars, only : gm
    use global_vars, only : mu_ref
    use global_vars, only : turbulence
    use global_vars, only : transition
    use global_vars, only : infile
    use global_vars, only : intermittency
    use global_vars, only : ExtraVar1
    use global_vars, only : ExtraVar2
    use global_vars, only : ExtraVar3
    use global_vars, only : ExtraVar4
    use global_vars, only : ExtraVar5
    
    use global_vars, only  : free_stream_density
    use global_vars, only  : free_stream_x_speed
    use global_vars, only  : free_stream_y_speed
    use global_vars, only  : free_stream_z_speed
    use global_vars, only  : free_stream_pressure
    use global_vars, only  : free_stream_tk
    use global_vars, only  : free_stream_tw
    use global_vars, only  : free_stream_tu
    use global_vars, only  : free_stream_tgm
    use global_vars, only  : vel_mag
    use global_vars, only  : MInf
    use global_vars, only  : Reynolds_number
    use global_vars, only  : mu_ratio_inf
    use global_vars, only  : Turb_intensity_inf


    use utils,       only: alloc, dealloc, dmsg
    use layout,      only: process_id
    use string
    use read_output, only: read_file

    use check_output_control, only : verify_write_control

    implicit none
    private


    ! Public methods
    public :: setup_state
    public :: destroy_state

    contains


        subroutine link_aliases()
          !< Setup state variable pointers

            implicit none

            DebugCall("link_aliases")

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

            select case (trim(turbulence))

                case ("none")
                    !include nothing
                    continue
                
                case ("sst", "sst2003", "bsl", "des-sst", "kw")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    tw(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
                    tk_inf => qp_inf(6)
                    tw_inf => qp_inf(7)

                case ("kkl")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    tkl(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
                    tk_inf => qp_inf(6)
                    tkl_inf => qp_inf(7)

                case ("sa", "saBC")
                    tv(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    tv_inf => qp_inf(6)

                case ("ke")
                    tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
                    te(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
                    tk_inf => qp_inf(6)
                    te_inf => qp_inf(7)

                case ("les")
                  continue
                  ! todo

                case DEFAULT
                  Fatal_error

            end select


            ! Transition modeling
            select case(trim(transition))

              case('lctm2015')
                tgm(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, n_var)
                tgm_inf => qp_inf(n_var)

              case('bc', 'none')
                !do nothing
                continue

              case DEFAULT
                Fatal_error

            end Select

        end subroutine link_aliases



        subroutine unlink_aliases()
          !< Nullify the pointer link

            implicit none

            DebugCall("unlink_aliases")

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

            select case (trim(turbulence))

                case ("none")
                    continue

                  case ("sst", "sst2003", "bsl", "kw", "des-sst")
                    nullify(tk)
                    nullify(tw)
                    nullify(tk_inf)
                    nullify(tw_inf)

                case ("kkl")
                    nullify(tk)
                    nullify(tkl)
                    nullify(tk_inf)
                    nullify(tkl_inf)

                case ("sa", "saBC")
                    nullify(tv)
                    nullify(tv_inf)

                case ("ke")
                    nullify(tk)
                    nullify(te)
                    nullify(tk_inf)
                    nullify(te_inf)

                case ("les")
                    continue
                    ! todo

                case DEFAULT
                  Fatal_error

            end select


            !Transition modeling
            select case(trim(transition))

              case('lctm2015')
                nullify(tgm)
                nullify(tgm_inf)

              case('bc', 'none')
                !do nothing
                continue

              case DEFAULT
                Fatal_error

            end Select

        end subroutine unlink_aliases




        subroutine allocate_memory()
            !< Allocate memory to the state variables
            !-----------------------------------------------------------
            implicit none

            DebugCall("allocate_memory")

            ! The state of the system is defined by the primitive 
            ! variables (density, velocity and pressure) at the grid
            ! cell centers. 
            call alloc(qp, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, n_var, AErrMsg("qp"))
            call alloc(qp_inf, 1, n_var, AErrMsg("qp_inf"))

        end subroutine allocate_memory



        subroutine deallocate_memory()
          !< Deallocate memory from the state variable

            implicit none

            DebugCall("allocate_memory")

            call dealloc(qp)

        end subroutine deallocate_memory



        subroutine setup_state()
            !< Setup the state module.
            !< This subroutine should be run before the state variables
            !< are initilized. This subroutine allocates the memory for 
            !< state variables and sets up the aliases to refer to the 
            !< components of the state
            !-----------------------------------------------------------

            implicit none

            DebugCall("setup_state")

            call set_n_var_value()
            call allocate_memory()
            call link_aliases()
            call init_infinity_values()
            call initstate()

        end subroutine setup_state



        subroutine destroy_state()
            !< Destroy the state module.
            !< This subroutine destroys the state module which includes
            !< unlinking the aliases for the state components and 
            !< deallocating the memory held by the state variables
            !-----------------------------------------------------------

            implicit none
            
            DebugCall("destroy_state")

            call unlink_aliases()
            call deallocate_memory()

        end subroutine destroy_state



        subroutine init_infinity_values()
            !< Set the values of the infinity variables "qp_inf"
            !-----------------------------------------------------------

            implicit none
            
            DebugCall("init_infinity_values")

            density_inf = free_stream_density
            x_speed_inf = free_stream_x_speed
            y_speed_inf = free_stream_y_speed
            z_speed_inf = free_stream_z_speed
            pressure_inf = free_stream_pressure
            vel_mag = sqrt(x_speed_inf**2 + y_speed_inf**2 + z_speed_inf**2)
            MInf    = Vel_mag/sqrt(gm*pressure_inf/density_inf)
            Reynolds_number = density_inf*vel_mag*1.0/mu_ref
            Turb_intensity_inf = free_stream_tu/100

            select case (trim(turbulence))
                
                case ("none")
                    continue

                case ("sst", "sst2003", "bsl")
                    tk_inf = 1.5*((Vel_mag*Turb_Intensity_inf)**2)
                    tw_inf = density_inf*tk_inf/(mu_ref*mu_ratio_inf)

                case ("kkl")
                    tk_inf = 9*(1e-9)*(sound_speed_inf()**2)
                    tkl_inf = 1.5589*(1e-6)*(mu_ref*sound_speed_inf())/density_inf

                case ("sa")
                     tv_inf = mu_ratio_inf*mu_ref/density_inf

                case ("saBC")
                    tv_inf = 0.005*mu_ratio_inf*mu_ref/density_inf

                case ("kw")
                    tk_inf = 1.5*((Vel_mag*Turb_Intensity_inf)**2)
                    tw_inf = density_inf*tk_inf/(mu_ref*mu_ratio_inf)

                case ("ke")
                    tk_inf = 1.5*((Vel_mag*Turb_Intensity_inf)**2)
                    tw_inf = 0.09*density_inf*tk_inf*tk_inf/(mu_ref*mu_ratio_inf)

                case ("des-sst")
                    tk_inf = 1.5*((Vel_mag*Turb_Intensity_inf)**2)
                    tw_inf = density_inf*tk_inf/(mu_ref*mu_ratio_inf)

                case ("les")
                  continue
                  ! todo

                case DEFAULT
                  Fatal_error

            end select


            !Transition modeling
            select case(trim(transition))

              case('lctm2015')
                tgm_inf = free_stream_tgm

              case('bc', 'none')
                !do nothing
                continue

              case DEFAULT
                Fatal_error

            end Select

        end subroutine init_infinity_values


        
        function sound_speed_inf() result(a)
            !< Return the free stream speed of sound.
            !-----------------------------------------------------------

            implicit none
            real :: a

            a = sqrt(gm * pressure_inf / density_inf)

        end function sound_speed_inf



        subroutine initstate()
            !< Initialize the state.
            !< If load file(start_from) is 0, then the state should be 
            !< set to the infinity values. Otherwise, read the state_file
            !< to get the state values
            !-----------------------------------------------------------

            implicit none
            
            DebugCall("initstate")

            call  verify_write_control()

            if (start_from .eq. 0) then
                ! Set the state to the infinity values
                call init_state_with_infinity_values()
                !!----------------------------------------
                !!following are added spefically for 
                !! shock tube test case
                !!---------------------------------------
                !if(process_id<2) then
                !  pressure = 1.0
                !  density = 1.0
                !  x_speed = 0.0
                !  y_speed = 0.0
                !  z_speed = 0.0
                !else
                !  pressure = 0.1
                !  density = 0.125
                !  x_speed = 0.0
                !  y_speed = 0.0
                !  z_speed = 0.0
                !end if
                !X_speed = 0.0
                !Y_speed = 0.0
                !Z_speed = 0.0
            else
                write(infile,'(a,i4.4,a,i2.2)') &
                  "time_directories/",start_from,"/process_",process_id
                ! Set the state to the infinity values so if some
                ! variable are not restart variable they get free_stream value
                call init_state_with_infinity_values()
                call read_file()

            end if

        end subroutine initstate



        subroutine init_state_with_infinity_values()
            !< Initialize the state based on the infinity values
            !-----------------------------------------------------------
            
            implicit none
            integer :: i
            
            DebugCall("init_state_with_infinity_values")
            
            do i = 1,n_var
                qp(:, :, :, i) = qp_inf(i)
            end do 
            
        end subroutine init_state_with_infinity_values



        subroutine set_n_var_value()
          !< Set number of variable to solver for based on
          !< the tubulence and transition model being used
          implicit none

          DebugCall("set_n_var_value")

          select case (trim(turbulence))
            case('none')
              n_var=5

            case('sa', 'saBC')
              n_var=6
              
            case('sst', "sst2003", 'bsl', 'kw', 'ke', 'kkl', 'Des-kw')
              n_var=7

            case DEFAULT
              n_var=5

          end select


          !Transition modeling
          select case(trim(transition))
            case('lctm2015')
              n_var = n_var + 1

            case('bc', 'none')
              n_var = n_var + 0

            case DEFAULT
              Fatal_error

          end Select

        end subroutine set_n_var_value

end module state
