module source
  !----------------------------------------------------------------------------
  !compute residue created by the source term
  !loop over cell centers.
  !----------------------------------------------------------------------------
  use utils, only: alloc, dealloc, dmsg
  use string
    use grid, only: imx, jmx, kmx
    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz, &
                        xA, yA, zA
    use geometry, only: vol => volume
    use state, only: gm, n_var, R_gas, mu_ref, T_ref, Pr, Sutherland_temp, &
                     density, x_speed, y_speed, z_speed, pressure
    use face_interpolant, only: x_density_left, x_density_right, &
        y_density_left, y_density_right, z_density_left, z_density_right, &
        x_x_speed_left, x_x_speed_right, x_y_speed_left, x_y_speed_right, &
        x_z_speed_left, x_z_speed_right, y_x_speed_left, y_x_speed_right, &
        y_y_speed_left, y_y_speed_right, y_z_speed_left, y_z_speed_right, &
        z_x_speed_left, z_x_speed_right, z_y_speed_left, z_y_speed_right, &
        z_z_speed_left, z_z_speed_right, x_pressure_left, x_pressure_right, &
        y_pressure_left, y_pressure_right, z_pressure_left, z_pressure_right
  include "turbulence_models/include/source/import_module.inc" 

  implicit none
  private

    real, public, dimension(:, :, :), allocatable :: gradu_x, gradu_y, &
    gradu_z, gradv_x, gradv_y, gradv_z, gradw_x, gradw_y, gradw_z, &
    gradT_x, gradT_y, gradT_z
    include "turbulence_models/include/source/variable_deceleration.inc"  

    public :: setup_source
    public :: destroy_source
    public :: add_source_term_residue
    public :: compute_gradients_cell_centre

    contains

    subroutine setup_source()

        implicit none

        call dmsg(1, 'source', 'setup_source')

      ! call alloc(gradu_y_f, 1, imx, 1, jmx-1, 1, kmx-1, &
      !         errmsg='Error: Unable to allocate memory for ' // &
      !             'Gradu_y_f - source')
        call alloc(gradu_x, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradu_x - source')
        call alloc(gradu_y, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradu_y - source')
        call alloc(gradu_z, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradu_z - source')

        call alloc(gradv_x, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradv_x - source')
        call alloc(gradv_y, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradv_y - source')
        call alloc(gradv_z, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradv_z - source')

        call alloc(gradw_x, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradw_x - source')
        call alloc(gradw_y, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradw_y - source')
        call alloc(gradw_z, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradw_z - source')

        call alloc(gradT_x, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'GradT_x - source')
        call alloc(gradT_y, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'GradT_y - source')
        call alloc(gradT_z, 0, imx, 0, jmx, 0, kmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'GradT_z - source')


        include "turbulence_models/include/source/setup_source.inc"

    end subroutine setup_source

    subroutine destroy_source()

        implicit none

        call dmsg(1, 'source', 'destroy_source')

        call dealloc(gradu_x)
        call dealloc(gradu_y)
        call dealloc(gradu_z)

        call dealloc(gradv_x)
        call dealloc(gradv_y)
        call dealloc(gradv_z)
        
        call dealloc(gradw_x)
        call dealloc(gradw_y)
        call dealloc(gradw_z)
        
        call dealloc(gradT_x)
        call dealloc(gradT_y)
        call dealloc(gradT_z)

        
        include "turbulence_models/include/source/destroy_source.inc"

    end subroutine destroy_source
    
    subroutine compute_gradients_cell_centre()
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the cell
    ! centre
    !-----------------------------------------------------------------

        implicit none

        integer :: i, j, k
        real :: T_r, T_l, T_face

        include "turbulence_models/include/source/gradient_init.inc"
        
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

            include "turbulence_models/include/source/gradient_find.inc"
          end do
         end do
        end do

    end subroutine compute_gradients_cell_centre

    
    subroutine add_source_term_residue()

      implicit none
      integer :: i, j, k

      call dmsg(1, 'source', 'compute_source_term_residue')

      call compute_gradients_cell_centre()

      do k = 1, kmx-1
        do j = 1, jmx-1
          do i = 1, imx-1
            include "turbulence_models/include/source/compute.inc"
          end do
        end do
      end do

    end subroutine add_source_term_residue


end module source


