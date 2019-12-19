module face_interpolant

#include "../../debug.h"
#include "../../error.h"
    use vartypes
!    use global, only: INTERPOLANT_NAME_LENGTH
    use global_vars, only : qp
!    use global_vars, only : gm
!    use global_vars, only : turbulence

!    use global_vars, only : interpolant

    use utils, only: alloc, dealloc
    use muscl, only: setup_scheme_muscl => setup_scheme, &
            destroy_scheme_muscl => destroy_scheme, &
            compute_muscl_states, &
            x_qp_left_muscl => x_qp_left, &
            x_qp_right_muscl => x_qp_right, &
            y_qp_left_muscl => y_qp_left, &
            y_qp_right_muscl => y_qp_right, &
            z_qp_left_muscl => z_qp_left, &
            z_qp_right_muscl => z_qp_right
    use ppm, only: setup_scheme_ppm => setup_scheme, &
            destroy_scheme_ppm => destroy_scheme, &
            compute_ppm_states, &
            x_qp_left_ppm => x_qp_left, &
            x_qp_right_ppm => x_qp_right, &
            y_qp_left_ppm => y_qp_left, &
            y_qp_right_ppm => y_qp_right, &
            z_qp_left_ppm => z_qp_left, &
            z_qp_right_ppm => z_qp_right
    use weno, only: setup_scheme_weno => setup_scheme, &
            destroy_scheme_weno => destroy_scheme, &
            compute_weno_states, &
             x_qp_left_weno => x_qp_left, &
            x_qp_right_weno => x_qp_right, &
             y_qp_left_weno => y_qp_left, &
            y_qp_right_weno => y_qp_right, &
             z_qp_left_weno => z_qp_left, &
            z_qp_right_weno => z_qp_right
    use weno_NM, only: setup_scheme_weno_NM => setup_scheme, &
            destroy_scheme_weno_NM => destroy_scheme, &
            compute_weno_NM_states, &
             x_qp_left_weno_NM => x_qp_left, &
            x_qp_right_weno_NM => x_qp_right, &
             y_qp_left_weno_NM => y_qp_left, &
            y_qp_right_weno_NM => y_qp_right, &
             z_qp_left_weno_NM => z_qp_left, &
            z_qp_right_weno_NM => z_qp_right
    !include "turbulence_models/include/face_interpolant/import_module.inc"

    implicit none
    private


    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :, :), allocatable, target :: y_qp_left, y_qp_right
    real, dimension(:, :, :, :), allocatable, target :: z_qp_left, z_qp_right
    real, dimension(:, :, :), pointer :: x_density_left, x_density_right
    real, dimension(:, :, :), pointer :: x_x_speed_left, x_x_speed_right
    real, dimension(:, :, :), pointer :: x_y_speed_left, x_y_speed_right
    real, dimension(:, :, :), pointer :: x_z_speed_left, x_z_speed_right
    real, dimension(:, :, :), pointer :: x_pressure_left, x_pressure_right
    real, dimension(:, :, :), pointer :: y_density_left, y_density_right
    real, dimension(:, :, :), pointer :: y_x_speed_left, y_x_speed_right
    real, dimension(:, :, :), pointer :: y_y_speed_left, y_y_speed_right
    real, dimension(:, :, :), pointer :: y_z_speed_left, y_z_speed_right
    real, dimension(:, :, :), pointer :: y_pressure_left, y_pressure_right
    real, dimension(:, :, :), pointer :: z_density_left, z_density_right
    real, dimension(:, :, :), pointer :: z_x_speed_left, z_x_speed_right
    real, dimension(:, :, :), pointer :: z_y_speed_left, z_y_speed_right
    real, dimension(:, :, :), pointer :: z_z_speed_left, z_z_speed_right
    real, dimension(:, :, :), pointer :: z_pressure_left, z_pressure_right

    integer :: imx, jmx, kmx, n_var

    !turbulent variable left and right with public deceleration
    !include "turbulence_models/include/face_interpolant/variables_deceleration.inc"
!    real, dimension(:, :, :), pointer :: x_tk_left, x_tk_right
!    real, dimension(:, :, :), pointer :: y_tk_left, y_tk_right
!    real, dimension(:, :, :), pointer :: z_tk_left, z_tk_right
!    real, dimension(:, :, :), pointer :: x_tw_left, x_tw_right
!    real, dimension(:, :, :), pointer :: y_tw_left, y_tw_right
!    real, dimension(:, :, :), pointer :: z_tw_left, z_tw_right
!    real, dimension(:, :, :), pointer :: x_tkl_left, x_tkl_right
!    real, dimension(:, :, :), pointer :: y_tkl_left, y_tkl_right
!    real, dimension(:, :, :), pointer :: z_tkl_left, z_tkl_right
!    real, dimension(:, :, :), pointer :: x_tv_left, x_tv_right
!    real, dimension(:, :, :), pointer :: y_tv_left, y_tv_right
!    real, dimension(:, :, :), pointer :: z_tv_left, z_tv_right
!    real, dimension(:, :, :), pointer :: x_te_left, x_te_right
!    real, dimension(:, :, :), pointer :: y_te_left, y_te_right
!    real, dimension(:, :, :), pointer :: z_te_left, z_te_right
!    real, dimension(:, :, :), pointer :: x_tgm_left, x_tgm_right
!    real, dimension(:, :, :), pointer :: y_tgm_left, y_tgm_right
!    real, dimension(:, :, :), pointer :: z_tgm_left, z_tgm_right
!
!    !public member
!    public :: x_tk_left, x_tk_right
!    public :: y_tk_left, y_tk_right
!    public :: z_tk_left, z_tk_right
!    public :: x_tw_left, x_tw_right
!    public :: y_tw_left, y_tw_right
!    public :: z_tw_left, z_tw_right
!    public :: x_tkl_left, x_tkl_right
!    public :: y_tkl_left, y_tkl_right
!    public :: z_tkl_left, z_tkl_right
!    public :: x_tv_left, x_tv_right
!    public :: y_tv_left, y_tv_right
!    public :: z_tv_left, z_tv_right
!    public :: x_te_left, x_te_right
!    public :: y_te_left, y_te_right
!    public :: z_te_left, z_te_right
!    public :: x_tgm_left, x_tgm_right
!    public :: y_tgm_left, y_tgm_right
!    public :: z_tgm_left, z_tgm_right

    ! Public members
!    public :: interpolant
    public :: setup_interpolant_scheme
    public :: extrapolate_cell_averages_to_faces
    public :: destroy_interpolant_scheme
    public :: compute_face_interpolant
!    public :: x_sound_speed_left
!    public :: x_sound_speed_right
!    public :: y_sound_speed_left
!    public :: y_sound_speed_right
!    public :: z_sound_speed_left
!    public :: z_sound_speed_right
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right
    public :: x_density_left, x_density_right
    public :: x_x_speed_left, x_x_speed_right
    public :: x_y_speed_left, x_y_speed_right
    public :: x_z_speed_left, x_z_speed_right
    public :: x_pressure_left, x_pressure_right
    public :: y_density_left, y_density_right
    public :: y_x_speed_left, y_x_speed_right
    public :: y_y_speed_left, y_y_speed_right
    public :: y_z_speed_left, y_z_speed_right
    public :: y_pressure_left, y_pressure_right
    public :: z_density_left, z_density_right
    public :: z_x_speed_left, z_x_speed_right
    public :: z_y_speed_left, z_y_speed_right
    public :: z_z_speed_left, z_z_speed_right
    public :: z_pressure_left, z_pressure_right

    contains

        subroutine allocate_memory()
            implicit none
            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')
            call alloc(z_qp_left, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_left.')
            call alloc(z_qp_right, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_right.')
        end subroutine allocate_memory

        subroutine link_aliases()
            implicit none

            ! Link xi faces left pointers
            x_density_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 1)
            x_x_speed_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 2)
            x_y_speed_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 3)
            x_z_speed_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 4)
            x_pressure_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 5)

            ! Link xi faces right pointers
            x_density_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 1)
            x_x_speed_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 2)
            x_y_speed_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 3)
            x_z_speed_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 4)
            x_pressure_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 5)

            ! Link eta faces left pointers
            y_density_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 1)
            y_x_speed_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 2)
            y_y_speed_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 3)
            y_z_speed_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 4)
            y_pressure_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 5)

            ! Link eta faces right pointers
            y_density_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 1)
            y_x_speed_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 2)
            y_y_speed_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 3)
            y_z_speed_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 4)
            y_pressure_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 5)
            
            ! Link zeta faces left pointers
            z_density_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 1)
            z_x_speed_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 2)
            z_y_speed_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 3)
            z_z_speed_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 4)
            z_pressure_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 5)

            ! Link zeta faces right pointers
            z_density_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 1)
            z_x_speed_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 2)
            z_y_speed_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 3)
            z_z_speed_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 4)
            z_pressure_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 5)

            !turbulent variable tk and tw linking _qp_ (6:7)
            !include "turbulence_models/include/face_interpolant/link_aliases.inc"
  
            !select case (turbulence)
       
            !    case ("none")
            !        !include nothing
            !        continue
            !   
            !    case ("sst", "sst2003")
            !        !include "turbulence_models/sst/face_interpolant/link_aliases.inc"
            !        ! Link xi faces left pointers
            !        x_tk_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)
            !        x_tw_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 7)

            !        ! Link xi faces right pointers
            !        x_tk_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)
            !        x_tw_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 7)

            !        ! Link eta faces left pointers
            !        y_tk_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)
            !        y_tw_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 7)

            !        ! Link eta faces right pointers
            !        y_tk_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            !        y_tw_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 7)
            !        
            !        ! Link zeta faces left pointers
            !        z_tk_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)
            !        z_tw_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 7)

            !        ! Link zeta faces right pointers
            !        z_tk_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
            !        z_tw_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 7)

  
            !    case ("kkl")
            !        !include "turbulence_models/kkl/face_interpolant/link_aliases.inc"
            !        ! Link xi faces left pointers
            !        x_tk_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)
            !        x_tkl_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 7)

            !        ! Link xi faces right pointers
            !        x_tk_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)
            !        x_tkl_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 7)

            !        ! Link eta faces left pointers
            !        y_tk_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)
            !        y_tkl_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 7)

            !        ! Link eta faces right pointers
            !        y_tk_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            !        y_tkl_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 7)
            !        
            !        ! Link zeta faces left pointers
            !        z_tk_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)
            !        z_tkl_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 7)

            !        ! Link zeta faces right pointers
            !        z_tk_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
            !        z_tkl_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 7)

  
            !    case ("sa", "saBC")
            !        !include "turbulence_models/sa/face_interpolant/link_aliases.inc"
            !        ! Link xi faces left pointers
            !        x_tv_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)

            !        ! Link xi faces right pointers
            !        x_tv_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)

            !        ! Link eta faces left pointers
            !        y_tv_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)

            !        ! Link eta faces right pointers
            !        y_tv_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            !        
            !        ! Link zeta faces left pointers
            !        z_tv_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)

            !        ! Link zeta faces right pointers
            !        z_tv_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
  
            !    case ("ke")
            !        !include "turbulence_models/ke/face_interpolant/link_aliases.inc"
            !        ! Link xi faces left pointers
            !        x_tk_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)
            !        x_te_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 7)

            !        ! Link xi faces right pointers
            !        x_tk_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)
            !        x_te_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 7)

            !        ! Link eta faces left pointers
            !        y_tk_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)
            !        y_te_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 7)

            !        ! Link eta faces right pointers
            !        y_tk_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            !        y_te_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 7)
            !        
            !        ! Link zeta faces left pointers
            !        z_tk_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)
            !        z_te_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 7)

            !        ! Link zeta faces right pointers
            !        z_tk_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
            !        z_te_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 7)
  
            !    case ("kw")
            !        !include "turbulence_models/kw/face_interpolant/link_aliases.inc"
            !        ! Link xi faces left pointers
            !        x_tk_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)
            !        x_tw_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 7)

            !        ! Link xi faces right pointers
            !        x_tk_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)
            !        x_tw_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 7)

            !        ! Link eta faces left pointers
            !        y_tk_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)
            !        y_tw_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 7)

            !        ! Link eta faces right pointers
            !        y_tk_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            !        y_tw_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 7)
            !        
            !        ! Link zeta faces left pointers
            !        z_tk_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)
            !        z_tw_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 7)

            !        ! Link zeta faces right pointers
            !        z_tk_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
            !        z_tw_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 7)
  
            !    case ("des-sst")
            !        !include "turbulence_models/des-sst/face_interpolant/link_aliases.inc"
            !        ! Link xi faces left pointers
            !        x_tk_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)
            !        x_tw_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 7)

            !        ! Link xi faces right pointers
            !        x_tk_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)
            !        x_tw_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 7)

            !        ! Link eta faces left pointers
            !        y_tk_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)
            !        y_tw_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 7)

            !        ! Link eta faces right pointers
            !        y_tk_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            !        y_tw_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 7)
            !        
            !        ! Link zeta faces left pointers
            !        z_tk_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)
            !        z_tw_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 7)

            !        ! Link zeta faces right pointers
            !        z_tk_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
            !        z_tw_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 7)

  
            !    case DEFAULT
            !       stop
  
            !end select
        end subroutine link_aliases

        subroutine setup_interpolant_scheme(control, scheme, dims)
            implicit none
            type(controltype), intent(in) :: control
            type(schemetype), intent(in) :: scheme
            type(extent), intent(in) :: dims

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            n_var = control%n_var

            select case (scheme%interpolant)
                case ("none")
                    ! Do nothing
                    continue
                case ("ppm")
                    call setup_scheme_ppm(control, dims)
                case ("muscl")
                    call setup_scheme_muscl(control, dims)
                case ("weno")
                    call setup_scheme_weno(control, dims)
                case ("weno_NM")
                    call setup_scheme_weno_NM(control, dims)
                case default
                    Fatal_error
            end select
            call allocate_memory()
            call link_aliases()
        end subroutine setup_interpolant_scheme

        subroutine deallocate_memory()
            implicit none
            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)
        end subroutine deallocate_memory

!        subroutine unlink_aliases()
!            implicit none
!
!            ! Unlink xi faces left pointers
!            nullify(x_density_left)
!            nullify(x_x_speed_left)
!            nullify(x_y_speed_left)
!            nullify(x_z_speed_left)
!            nullify(x_pressure_left)
!
!            ! Unlink xi faces right pointers
!            nullify(x_density_right)
!            nullify(x_x_speed_right)
!            nullify(x_y_speed_right)
!            nullify(x_z_speed_right)
!            nullify(x_pressure_right)
!
!            ! Unlink eta faces left pointers
!            nullify(y_density_left)
!            nullify(y_x_speed_left)
!            nullify(y_y_speed_left)
!            nullify(y_z_speed_left)
!            nullify(y_pressure_left)
!
!            ! Unlink eta faces right pointers
!            nullify(y_density_right)
!            nullify(y_x_speed_right)
!            nullify(y_y_speed_right)
!            nullify(y_z_speed_right)
!            nullify(y_pressure_right)
!            
!            ! Unlink tau faces left pointers
!            nullify(z_density_left)
!            nullify(z_x_speed_left)
!            nullify(z_y_speed_left)
!            nullify(z_z_speed_left)
!            nullify(z_pressure_left)
!
!            ! Unlink tau faces right pointers
!            nullify(z_density_right)
!            nullify(z_x_speed_right)
!            nullify(z_y_speed_right)
!            nullify(z_z_speed_right)
!            nullify(z_pressure_right)
!
!            include "turbulence_models/include/face_interpolant/unlink_aliases.inc" 
!        end subroutine unlink_aliases

        subroutine destroy_interpolant_scheme(scheme)
            implicit none
            type(schemetype), intent(in) :: scheme
            !call unlink_aliases()
            call deallocate_memory()
            select case (scheme%interpolant)
                case ("none")
                    ! Do nothing
                    continue
                case ("ppm")
                    call destroy_scheme_ppm()
                case ("muscl")
                    call destroy_scheme_muscl()
                case ("weno")
                    call destroy_scheme_weno()
                case ("weno_NM")
                    call destroy_scheme_weno_NM()
                case default
                    Fatal_error
            end select
        end subroutine destroy_interpolant_scheme

        subroutine extrapolate_cell_averages_to_faces()
            implicit none

            DebugCall('extrapolate_cell_averages_to_faces')

            x_qp_left(:, :, :, :) = qp(-1:imx, 1:jmx-1, 1:kmx-1, 1:n_var)
            x_qp_right(:, :, :, :) = qp(0:imx+1, 1:jmx-1, 1:kmx-1, 1:n_var)
            y_qp_left(:, :, :, :) = qp(1:imx-1, -1:jmx, 1:kmx-1, 1:n_var)
            y_qp_right(:, :, :, :) = qp(1:imx-1, 0:jmx+1, 1:kmx-1, 1:n_var)
            z_qp_left(:, :, :, :) = qp(1:imx-1, 1:jmx-1, -1:kmx, 1:n_var)
            z_qp_right(:, :, :, :) = qp(1:imx-1, 1:jmx-1, 0:kmx+1, 1:n_var)
        end subroutine extrapolate_cell_averages_to_faces

        subroutine compute_face_interpolant(scheme, flow)
            implicit none
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(in) :: flow
            select case (scheme%interpolant)
                case ("none")
                    call extrapolate_cell_averages_to_faces()
                case ("ppm")
                    call compute_ppm_states(scheme, flow)
                    x_qp_left(:, :, :, :) = x_qp_left_ppm(:, :, :, :)
                    x_qp_right(:, :, :, :) = x_qp_right_ppm(:, :, :, :)
                    y_qp_left(:, :, :, :) = y_qp_left_ppm(:, :, :, :)
                    y_qp_right(:, :, :, :) = y_qp_right_ppm(:, :, :, :)
                    z_qp_left(:, :, :, :) = z_qp_left_ppm(:, :, :, :)
                    z_qp_right(:, :, :, :) = z_qp_right_ppm(:, :, :, :)
                case ("muscl")
                    call compute_muscl_states(scheme, flow)
                    x_qp_left = x_qp_left_muscl
                    x_qp_right = x_qp_right_muscl
                    y_qp_left = y_qp_left_muscl
                    y_qp_right = y_qp_right_muscl
                    z_qp_left = z_qp_left_muscl
                    z_qp_right = z_qp_right_muscl
                    !x_qp_left( :, :, :, 6:) = qp(-1:imx, 1:jmx-1, 1:kmx-1, 6:n_var)
                    !x_qp_right(:, :, :, 6:) = qp(0:imx+1, 1:jmx-1, 1:kmx-1, 6:n_var)
                    !y_qp_left( :, :, :, 6:) = qp(1:imx-1, -1:jmx, 1:kmx-1, 6:n_var)
                    !y_qp_right(:, :, :, 6:) = qp(1:imx-1, 0:jmx+1, 1:kmx-1, 6:n_var)
                    !z_qp_left( :, :, :, 6:) = qp(1:imx-1, 1:jmx-1, -1:kmx, 6:n_var)
                    !z_qp_right(:, :, :, 6:) = qp(1:imx-1, 1:jmx-1, 0:kmx+1, 6:n_var)
                case ("weno")
                    call compute_weno_states()
                    x_qp_left  =  x_qp_left_weno
                    x_qp_right = x_qp_right_weno
                    y_qp_left  =  y_qp_left_weno
                    y_qp_right = y_qp_right_weno
                    z_qp_left  =  z_qp_left_weno
                    z_qp_right = z_qp_right_weno
                    x_qp_left( :, :, :, 6:) = qp(-1:imx, 1:jmx-1, 1:kmx-1, 6:n_var)
                    x_qp_right(:, :, :, 6:) = qp(0:imx+1, 1:jmx-1, 1:kmx-1, 6:n_var)
                    y_qp_left( :, :, :, 6:) = qp(1:imx-1, -1:jmx, 1:kmx-1, 6:n_var)
                    y_qp_right(:, :, :, 6:) = qp(1:imx-1, 0:jmx+1, 1:kmx-1, 6:n_var)
                    z_qp_left( :, :, :, 6:) = qp(1:imx-1, 1:jmx-1, -1:kmx, 6:n_var)
                    z_qp_right(:, :, :, 6:) = qp(1:imx-1, 1:jmx-1, 0:kmx+1, 6:n_var)
                case ("weno_NM")
                    call compute_weno_NM_states()
                    x_qp_left  =  x_qp_left_weno_NM
                    x_qp_right = x_qp_right_weno_NM
                    y_qp_left  =  y_qp_left_weno_NM
                    y_qp_right = y_qp_right_weno_NM
                    z_qp_left  =  z_qp_left_weno_NM
                    z_qp_right = z_qp_right_weno_NM
                case default
                    Fatal_error
            end select
        end subroutine compute_face_interpolant

end module face_interpolant
