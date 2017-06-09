module source
  !----------------------------------------------------------------------------
  !compute residue created by the source term
  !loop over cell centers.
  !----------------------------------------------------------------------------

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx

    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z
    use global_vars, only : xA, yA, zA    !face area
    use global_vars, only : vol => volume

    use global_vars, only : gm
    use global_vars, only : n_var
    use global_vars, only : R_gas
    use global_vars, only : mu_ref
    use global_vars, only : T_ref
    use global_vars, only : Pr
    use global_vars, only : Sutherland_temp
    use global_vars, only : density
    use global_vars, only : x_speed
    use global_vars, only : y_speed
    use global_vars, only : z_speed
    use global_vars, only : pressure
    use global_vars, only : tk
    use global_vars, only : tw
    use global_vars, only : turbulence
    use global_vars, only : TKE_residue
    use global_vars, only : omega_residue
    use global_vars, only : dist
    use global_vars  ,only :   gradu_x
    use global_vars  ,only :   gradu_y
    use global_vars  ,only :   gradu_z
    use global_vars  ,only :   gradv_x
    use global_vars  ,only :   gradv_y
    use global_vars  ,only :   gradv_z
    use global_vars  ,only :   gradw_x
    use global_vars  ,only :   gradw_y
    use global_vars  ,only :   gradw_z
    use global_vars  ,only :   gradT_x
    use global_vars  ,only :   gradT_y
    use global_vars  ,only :   gradT_z
    use global_vars  ,only :   gradtk_x
    use global_vars  ,only :   gradtk_y
    use global_vars  ,only :   gradtk_z
    use global_vars  ,only :   gradtw_x
    use global_vars  ,only :   gradtw_y
    use global_vars  ,only :   gradtw_z
    use global_vars  ,only :   gradqp_z
    use global_vars  ,only :   process_id
    use utils, only: alloc, dealloc, dmsg
    use utils, only: turbulence_read_error
    use string
    use face_interpolant, only: x_density_left, x_density_right, &
        y_density_left, y_density_right, z_density_left, z_density_right, &
        x_x_speed_left, x_x_speed_right, x_y_speed_left, x_y_speed_right, &
        x_z_speed_left, x_z_speed_right, y_x_speed_left, y_x_speed_right, &
        y_y_speed_left, y_y_speed_right, y_z_speed_left, y_z_speed_right, &
        z_x_speed_left, z_x_speed_right, z_y_speed_left, z_y_speed_right, &
        z_z_speed_left, z_z_speed_right, x_pressure_left, x_pressure_right, &
        y_pressure_left, y_pressure_right, z_pressure_left, z_pressure_right
    use sst_source  ,only : add_sst_source
    !use turbulent_source , only : setup_turbulent_grad
    !use turbulent_source , only : destroy_turbulent_grad
    use turbulent_source , only : compute_turbulent_grad
!  include "turbulence_models/include/source/import_module.inc" 

  implicit none
  private

!    real, public, dimension(:, :, :), allocatable :: gradu_x, gradu_y, &
!    gradu_z, gradv_x, gradv_y, gradv_z, gradw_x, gradw_y, gradw_z, &
!    gradT_x, gradT_y, gradT_z
!    include "turbulence_models/include/source/variable_deceleration.inc"  

!    public :: setup_source
!    public :: destroy_source
    public :: add_source_term_residue
    public :: compute_gradients_cell_centre

    contains

!    subroutine setup_source()
!
!        implicit none
!
!        call dmsg(1, 'source', 'setup_source')
!
!      ! call alloc(gradu_y_f, 1, imx, 1, jmx-1, 1, kmx-1, &
!      !         errmsg='Error: Unable to allocate memory for ' // &
!      !             'Gradu_y_f - source')
!        call alloc(gradu_x, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradu_x - source')
!        call alloc(gradu_y, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradu_y - source')
!        call alloc(gradu_z, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradu_z - source')
!
!        call alloc(gradv_x, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradv_x - source')
!        call alloc(gradv_y, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradv_y - source')
!        call alloc(gradv_z, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradv_z - source')
!
!        call alloc(gradw_x, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradw_x - source')
!        call alloc(gradw_y, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradw_y - source')
!        call alloc(gradw_z, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'Gradw_z - source')
!
!        call alloc(gradT_x, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'GradT_x - source')
!        call alloc(gradT_y, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'GradT_y - source')
!        call alloc(gradT_z, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for ' // &
!                    'GradT_z - source')
!
!
!        call setup_turbulent_grad()
!!        include "turbulence_models/include/source/setup_source.inc"
!
!    end subroutine setup_source
!
!    subroutine destroy_source()
!
!        implicit none
!
!        call dmsg(1, 'source', 'destroy_source')
!
!        call dealloc(gradu_x)
!        call dealloc(gradu_y)
!        call dealloc(gradu_z)
!
!        call dealloc(gradv_x)
!        call dealloc(gradv_y)
!        call dealloc(gradv_z)
!        
!        call dealloc(gradw_x)
!        call dealloc(gradw_y)
!        call dealloc(gradw_z)
!        
!        call dealloc(gradT_x)
!        call dealloc(gradT_y)
!        call dealloc(gradT_z)
!
! 
!        call destroy_turbulent_grad
!!        include "turbulence_models/include/source/destroy_source.inc"
!
!    end subroutine destroy_source
    
    subroutine compute_gradients_cell_centre()
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the cell
    ! centre
    !-----------------------------------------------------------------

        implicit none

        integer :: i, j, k
        real :: T_r, T_l, T_face

        !include "turbulence_models/include/source/gradient_init.inc"
        
        call dmsg(1, 'source', 'compute_gradient_cell_center')

        gradu_x = 0.0
        gradu_y = 0.0
        gradu_z = 0.0
        gradv_x = 0.0
        gradv_y = 0.0
        gradv_z = 0.0
        gradw_x = 0.0
        gradw_y = 0.0
        gradw_z = 0.0
        gradT_x = 0.0
        gradT_y = 0.0
        gradT_z = 0.0



        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            
            ! Solving for gradu
            gradu_x(i,j,k) = ( &
                             - ((x_x_speed_left(i, j, k) + x_x_speed_right(i, j, k)) * &
                                  xnx(i,j,k) * xA(i,j,k)) &
                             + ((x_x_speed_left(i+1, j, k) + x_x_speed_right(i+1, j, k)) * &
                                  xnx(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_x_speed_left(i, j, k) + y_x_speed_right(i, j, k)) * &
                                  ynx(i, j, k) * yA(i, j, k)) &
                             + ((y_x_speed_left(i, j+1, k) + y_x_speed_right(i, j+1, k)) * &
                                  ynx(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_x_speed_left(i, j, k) + z_x_speed_right(i, j, k)) * &
                                  znx(i, j, k) * zA(i, j, k)) &
                             + ((z_x_speed_left(i, j, k+1) + z_x_speed_right(i, j, k+1)) * &
                                  znx(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            gradu_y(i,j,k) = ( &
                             - ((x_x_speed_left(i, j, k) + x_x_speed_right(i, j, k)) * &
                                  xny(i,j,k) * xA(i,j,k)) &
                             + ((x_x_speed_left(i+1, j, k) + x_x_speed_right(i+1, j, k)) * &
                                  xny(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_x_speed_left(i, j, k) + y_x_speed_right(i, j, k)) * &
                                  yny(i, j, k) * yA(i, j, k)) &
                             + ((y_x_speed_left(i, j+1, k) + y_x_speed_right(i, j+1, k)) * &
                                  yny(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_x_speed_left(i, j, k) + z_x_speed_right(i, j, k)) * &
                                  zny(i, j, k) * zA(i, j, k)) &
                             + ((z_x_speed_left(i, j, k+1) + z_x_speed_right(i, j, k+1)) * &
                                  zny(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            gradu_z(i,j,k) = ( &
                             - ((x_x_speed_left(i, j, k) + x_x_speed_right(i, j, k)) * &
                                  xnz(i,j,k) * xA(i,j,k)) &
                             + ((x_x_speed_left(i+1, j, k) + x_x_speed_right(i+1, j, k)) * &
                                  xnz(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_x_speed_left(i, j, k) + y_x_speed_right(i, j, k)) * &
                                  ynz(i, j, k) * yA(i, j, k)) &
                             + ((y_x_speed_left(i, j+1, k) + y_x_speed_right(i, j+1, k)) * &
                                  ynz(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_x_speed_left(i, j, k) + z_x_speed_right(i, j, k)) * &
                                  znz(i, j, k) * zA(i, j, k)) &
                             + ((z_x_speed_left(i, j, k+1) + z_x_speed_right(i, j, k+1)) * &
                                  znz(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)


            ! Solving for gradv
            gradv_x(i,j,k) = ( &
                             - ((x_y_speed_left(i, j, k) + x_y_speed_right(i, j, k)) * &
                                  xnx(i,j,k) * xA(i,j,k)) &
                             + ((x_y_speed_left(i+1, j, k) + x_y_speed_right(i+1, j, k)) * &
                                  xnx(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_y_speed_left(i, j, k) + y_y_speed_right(i, j, k)) * &
                                  ynx(i, j, k) * yA(i, j, k)) &
                             + ((y_y_speed_left(i, j+1, k) + y_y_speed_right(i, j+1, k)) * &
                                  ynx(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_y_speed_left(i, j, k) + z_y_speed_right(i, j, k)) * &
                                  znx(i, j, k) * zA(i, j, k)) &
                             + ((z_y_speed_left(i, j, k+1) + z_y_speed_right(i, j, k+1)) * &
                                  znx(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            gradv_y(i,j,k) = ( &
                             - ((x_y_speed_left(i, j, k) + x_y_speed_right(i, j, k)) * &
                                  xny(i,j,k) * xA(i,j,k)) &
                             + ((x_y_speed_left(i+1, j, k) + x_y_speed_right(i+1, j, k)) * &
                                  xny(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_y_speed_left(i, j, k) + y_y_speed_right(i, j, k)) * &
                                  yny(i, j, k) * yA(i, j, k)) &
                             + ((y_y_speed_left(i, j+1, k) + y_y_speed_right(i, j+1, k)) * &
                                  yny(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_y_speed_left(i, j, k) + z_y_speed_right(i, j, k)) * &
                                  zny(i, j, k) * zA(i, j, k)) &
                             + ((z_y_speed_left(i, j, k+1) + z_y_speed_right(i, j, k+1)) * &
                                  zny(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            gradv_z(i,j,k) = ( &
                             - ((x_y_speed_left(i, j, k) + x_y_speed_right(i, j, k)) * &
                                  xnz(i,j,k) * xA(i,j,k)) &
                             + ((x_y_speed_left(i+1, j, k) + x_y_speed_right(i+1, j, k)) * &
                                  xnz(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_y_speed_left(i, j, k) + y_y_speed_right(i, j, k)) * &
                                  ynz(i, j, k) * yA(i, j, k)) &
                             + ((y_y_speed_left(i, j+1, k) + y_y_speed_right(i, j+1, k)) * &
                                  ynz(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_y_speed_left(i, j, k) + z_y_speed_right(i, j, k)) * &
                                  znz(i, j, k) * zA(i, j, k)) &
                             + ((z_y_speed_left(i, j, k+1) + z_y_speed_right(i, j, k+1)) * &
                                  znz(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            
            ! Solving for gradw
            gradw_x(i,j,k) = ( &
                             - ((x_z_speed_left(i, j, k) + x_z_speed_right(i, j, k)) * &
                                  xnx(i,j,k) * xA(i,j,k)) &
                             + ((x_z_speed_left(i+1, j, k) + x_z_speed_right(i+1, j, k)) * &
                                  xnx(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_z_speed_left(i, j, k) + y_z_speed_right(i, j, k)) * &
                                  ynx(i, j, k) * yA(i, j, k)) &
                             + ((y_z_speed_left(i, j+1, k) + y_z_speed_right(i, j+1, k)) * &
                                  ynx(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_z_speed_left(i, j, k) + z_z_speed_right(i, j, k)) * &
                                  znx(i, j, k) * zA(i, j, k)) &
                             + ((z_z_speed_left(i, j, k+1) + z_z_speed_right(i, j, k+1)) * &
                                  znx(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            gradw_y(i,j,k) = ( &
                             - ((x_z_speed_left(i, j, k) + x_z_speed_right(i, j, k)) * &
                                  xny(i,j,k) * xA(i,j,k)) &
                             + ((x_z_speed_left(i+1, j, k) + x_z_speed_right(i+1, j, k)) * &
                                  xny(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_z_speed_left(i, j, k) + y_z_speed_right(i, j, k)) * &
                                  yny(i, j, k) * yA(i, j, k)) &
                             + ((y_z_speed_left(i, j+1, k) + y_z_speed_right(i, j+1, k)) * &
                                  yny(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_z_speed_left(i, j, k) + z_z_speed_right(i, j, k)) * &
                                  zny(i, j, k) * zA(i, j, k)) &
                             + ((z_z_speed_left(i, j, k+1) + z_z_speed_right(i, j, k+1)) * &
                                  zny(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            gradw_z(i,j,k) = ( &
                             - ((x_z_speed_left(i, j, k) + x_z_speed_right(i, j, k)) * &
                                  xnz(i,j,k) * xA(i,j,k)) &
                             + ((x_z_speed_left(i+1, j, k) + x_z_speed_right(i+1, j, k)) * &
                                  xnz(i+1, j, k) * xA(i+1, j, k)) &
                             - ((y_z_speed_left(i, j, k) + y_z_speed_right(i, j, k)) * &
                                  ynz(i, j, k) * yA(i, j, k)) &
                             + ((y_z_speed_left(i, j+1, k) + y_z_speed_right(i, j+1, k)) * &
                                  ynz(i, j+1, k) * yA(i, j+1, k)) &
                             - ((z_z_speed_left(i, j, k) + z_z_speed_right(i, j, k)) * &
                                  znz(i, j, k) * zA(i, j, k)) &
                             + ((z_z_speed_left(i, j, k+1) + z_z_speed_right(i, j, k+1)) * &
                                  znz(i, j, k+1) * zA(i, j, k+1)) &
                             ) * 0.5 / vol(i, j, k)

            ! Finding grad T
            ! Since T is not stored, each face value is calculated in-situ
            ! Note that grad T is in the direction of the face normal.
            ! Hence each term in the gradT_x, gradT_y and gradT_z are
            ! same and are different only by which component of the face
            ! normal they are multiplied by.
            ! In the below formulation, Tface varies with different indices
            ! and face directions used. Hence, the six terms that make up 
            ! gradT_x, gradT_y, gradT_z are of indices: 
            ! x_face(i,j,k)
            ! x_face(i+1, j, k)
            ! y_face(i, j, k)
            ! y_face(i, j+1, k)
            ! z_face(i, j, k)
            ! z_face(i, j, k+1)
            !
            ! NOTE: The factor of 0.5 will be multiplied in the end

            ! x_face(i, j, k)
            T_l = x_pressure_left(i, j, k) / (x_density_left(i, j, k) * R_gas)
            T_r = x_pressure_right(i, j, k) / (x_density_right(i, j, k) * R_gas)
            T_face = (T_l + T_r)
            gradT_x(i, j, k) = gradT_x(i, j, k) - &
                               (T_face * xnx(i, j, k) * xA(i, j, k))
            gradT_y(i, j, k) = gradT_y(i, j, k) - &
                               (T_face * xny(i, j, k) * xA(i, j, k))
            gradT_z(i, j, k) = gradT_z(i, j, k) - &
                               (T_face * xnz(i, j, k) * xA(i, j, k))

            ! x_face(i+1, j, k)
            T_l = x_pressure_left(i+1, j, k) / (x_density_left(i+1, j, k) * R_gas)
            T_r = x_pressure_right(i+1, j, k) / (x_density_right(i+1, j, k) * R_gas)
            T_face = (T_l + T_r)
            gradT_x(i, j, k) = gradT_x(i, j, k) + &
                               (T_face * xnx(i+1, j, k) * xA(i+1, j, k))
            gradT_y(i, j, k) = gradT_y(i, j, k) + &
                               (T_face * xny(i+1, j, k) * xA(i+1, j, k))
            gradT_z(i, j, k) = gradT_z(i, j, k) + &
                               (T_face * xnz(i+1, j, k) * xA(i+1, j, k))

            ! y_face(i, j, k)
            T_l = y_pressure_left(i, j, k) / (y_density_left(i, j, k) * R_gas)
            T_r = y_pressure_right(i, j, k) / (y_density_right(i, j, k) * R_gas)
            T_face = (T_l + T_r)
            gradT_x(i, j, k) = gradT_x(i, j, k) - &
                               (T_face * ynx(i, j, k) * yA(i, j, k))
            gradT_y(i, j, k) = gradT_y(i, j, k) - &
                               (T_face * yny(i, j, k) * yA(i, j, k))
            gradT_z(i, j, k) = gradT_z(i, j, k) - &
                               (T_face * ynz(i, j, k) * yA(i, j, k))

            ! y_face(i, j+1, k)
            T_l = y_pressure_left(i, j+1, k) / (y_density_left(i, j+1, k) * R_gas)
            T_r = y_pressure_right(i, j+1, k) / (y_density_right(i, j+1, k) * R_gas)
            T_face = (T_l + T_r)
            gradT_x(i, j, k) = gradT_x(i, j, k) + &
                               (T_face * ynx(i, j+1, k) * yA(i, j+1, k))
            gradT_y(i, j, k) = gradT_y(i, j, k) + &
                               (T_face * yny(i, j+1, k) * yA(i, j+1, k))
            gradT_z(i, j, k) = gradT_z(i, j, k) + &
                               (T_face * ynz(i, j+1, k) * yA(i, j+1, k))

            ! z_face(i, j, k)
            T_l = z_pressure_left(i, j, k) / (z_density_left(i, j, k) * R_gas)
            T_r = z_pressure_right(i, j, k) / (z_density_right(i, j, k) * R_gas)
            T_face = (T_l + T_r)
            gradT_x(i, j, k) = gradT_x(i, j, k) - &
                               (T_face * znx(i, j, k) * zA(i, j, k))
            gradT_y(i, j, k) = gradT_y(i, j, k) - &
                               (T_face * zny(i, j, k) * zA(i, j, k))
            gradT_z(i, j, k) = gradT_z(i, j, k) - &
                               (T_face * znz(i, j, k) * zA(i, j, k))

            ! z_face(i, j, k+1)
            T_l = z_pressure_left(i, j, k+1) / (z_density_left(i, j, k+1) * R_gas)
            T_r = z_pressure_right(i, j, k+1) / (z_density_right(i, j, k+1) * R_gas)
            T_face = (T_l + T_r)
            gradT_x(i, j, k) = gradT_x(i, j, k) + &
                               (T_face * znx(i, j, k+1) * zA(i, j, k+1))
            gradT_y(i, j, k) = gradT_y(i, j, k) + &
                               (T_face * zny(i, j, k+1) * zA(i, j, k+1))
            gradT_z(i, j, k) = gradT_z(i, j, k) + &
                               (T_face * znz(i, j, k+1) * zA(i, j, k+1))

            ! Factor of volume
            gradT_x(i, j, k) = gradT_x(i, j, k) * 0.5 / vol(i, j, k)
            gradT_y(i, j, k) = gradT_y(i, j, k) * 0.5 / vol(i, j, k)
            gradT_z(i, j, k) = gradT_z(i, j, k) * 0.5 / vol(i, j, k)

!            include "turbulence_models/include/source/gradient_find.inc"
          if(gradu_y(i,j,k)<0.) print*, process_id, i,j,k, gradu_y(i,j,k)
          end do
         end do
        end do

    end subroutine compute_gradients_cell_centre

    
    subroutine add_source_term_residue()

      implicit none

      call dmsg(1, 'source', 'compute_source_term_residue')

!      call compute_gradients_cell_centre()
!      call compute_gradient(gradu_x, x_speed, 'x')
!      call compute_gradient(gradu_y, x_speed, 'y')
!      call compute_gradient(gradu_z, x_speed, 'z')
!      call compute_gradient(gradv_x, y_speed, 'x')
!      call compute_gradient(gradv_y, y_speed, 'y')
!      call compute_gradient(gradv_z, y_speed, 'z')
!      call compute_gradient(gradw_x, z_speed, 'x')
!      call compute_gradient(gradw_y, z_speed, 'y')
!      call compute_gradient(gradw_z, z_speed, 'z')
!      call compute_gradient_temp(gradT_x    , 'x')
!      call compute_gradient_temp(gradT_y    , 'y')
!      call compute_gradient_temp(gradT_z    , 'z')
!      call compute_turbulent_grad()
!      if(kmx==2) gradqp_z=0.0

!      do k = 1, kmx-1
!        do j = 1, jmx-1
!          do i = 1, imx-1
!            include "turbulence_models/include/source/compute.inc"
!          end do
!        end do
!      end do
      select case (turbulence)

        case ('none')
          !do nothing
          continue

        case ('sst')
          call add_sst_source()

        case DEFAULT
          call turbulence_read_error()

      end select

    end subroutine add_source_term_residue


    subroutine compute_gradient(grad, var, dir)
      implicit none
      real, dimension( 0:imx  , 0:jmx  , 0:kmx  ), intent(out) :: grad
      real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(in) :: var
      character(len=*)                           , intent(in) :: dir
      
      real, dimension(:,:,:), pointer  :: nx
      real, dimension(:,:,:), pointer  :: ny
      real, dimension(:,:,:), pointer  :: nz

      integer :: i
      integer :: j
      integer :: k

      select case(dir)
        case('x')
          nx(1:imx  ,1:jmx-1,1:kmx-1) => xnx
          ny(1:imx-1,1:jmx  ,1:kmx-1) => ynx
          nz(1:imx-1,1:jmx-1,1:kmx  ) => znx
        case('y')
          nx(1:imx  ,1:jmx-1,1:kmx-1) => xny
          ny(1:imx-1,1:jmx  ,1:kmx-1) => yny
          nz(1:imx-1,1:jmx-1,1:kmx  ) => zny
        case('z')
          nx(1:imx  ,1:jmx-1,1:kmx-1) => xnz
          ny(1:imx-1,1:jmx  ,1:kmx-1) => ynz
          nz(1:imx-1,1:jmx-1,1:kmx  ) => znz
        case DEFAULT
          print*, "ERROR: gradient direction error"
      end select
      grad = 0.0

      do k=1,kmx-1
        do j=1,jmx-1
          do i=1,imx-1
            grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*nx(i,j,k)*xA(i,j,k) &
                          -(var(i  ,j-1,k  )+var(i,j,k))*ny(i,j,k)*yA(i,j,k) &
                          -(var(i  ,j  ,k-1)+var(i,j,k))*nz(i,j,k)*zA(i,j,k) &
                          +(var(i+1,j  ,k  )+var(i,j,k))*nx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                          +(var(i  ,j+1,k  )+var(i,j,k))*ny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                          +(var(i  ,j  ,k+1)+var(i,j,k))*nz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                         )/(2*vol(i,j,k))
          end do
        end do
      end do

    end subroutine compute_gradient

    subroutine compute_gradient_temp(grad, dir)

      implicit none
      real, dimension( 0:imx  , 0:jmx  , 0:kmx  ), intent(out) :: grad
      character(len=*)                           , intent(in) :: dir
      
      real, dimension(6)               :: T
      real                             :: cell_T
      real, dimension(:,:,:), pointer  :: nx
      real, dimension(:,:,:), pointer  :: ny
      real, dimension(:,:,:), pointer  :: nz

      integer :: i
      integer :: j
      integer :: k

      select case(dir)
        case('x')
          nx(1:imx  ,1:jmx-1,1:kmx-1) => xnx
          ny(1:imx-1,1:jmx  ,1:kmx-1) => ynx
          nz(1:imx-1,1:jmx-1,1:kmx  ) => znx
        case('y')
          nx(1:imx  ,1:jmx-1,1:kmx-1) => xny
          ny(1:imx-1,1:jmx  ,1:kmx-1) => yny
          nz(1:imx-1,1:jmx-1,1:kmx  ) => zny
        case('z')
          nx(1:imx  ,1:jmx-1,1:kmx-1) => xnz
          ny(1:imx-1,1:jmx  ,1:kmx-1) => ynz
          nz(1:imx-1,1:jmx-1,1:kmx  ) => znz
        case DEFAULT
          print*, "ERROR: gradient direction error"
      end select
      grad = 0.0

      do k=1,kmx-1
        do j=1,jmx-1
          do i=1,imx-1

            cell_T = (pressure(i,j,k)/density(i,j,k))/R_gas

            T(1)   = (pressure(i-1,j,k)/density(i-1,j,k))/R_gas + cell_T
            T(2)   = (pressure(i,j-1,k)/density(i,j-1,k))/R_gas + cell_T
            T(3)   = (pressure(i,j,k-1)/density(i,j,k-1))/R_gas + cell_T
            T(4)   = (pressure(i+1,j,k)/density(i+1,j,k))/R_gas + cell_T
            T(5)   = (pressure(i,j+1,k)/density(i,j+1,k))/R_gas + cell_T
            T(6)   = (pressure(i,j,k+1)/density(i,j,k+1))/R_gas + cell_T

            grad(i,j,k) =(-T(1)*nx(i,j,k)*xA(i,j,k) &
                          -T(2)*ny(i,j,k)*yA(i,j,k) &
                          -T(3)*nz(i,j,k)*zA(i,j,k) &
                          +T(4)*nx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                          +T(5)*ny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                          +T(6)*nz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                         )/(2*vol(i,j,k))
          end do
        end do
      end do

    end subroutine compute_gradient_temp


end module source


