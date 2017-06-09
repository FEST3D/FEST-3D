module summon_grad_evaluation
  !----------------------------------------------------------
  !170608  -Jatinder Pal Singh Sandhu
  ! Aim : call is made to all the required gradients
  !                         based on input conditions
  !----------------------------------------------------------

  use global_vars, only : kmx
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : tk
  use global_vars, only : tw
  use global_vars, only : turbulence
  use global_vars, only : gradu_x
  use global_vars, only : gradu_y
  use global_vars, only : gradu_z
  use global_vars, only : gradv_x
  use global_vars, only : gradv_y
  use global_vars, only : gradv_z
  use global_vars, only : gradw_x
  use global_vars, only : gradw_y
  use global_vars, only : gradw_z
  use global_vars, only : gradT_x
  use global_vars, only : gradT_y
  use global_vars, only : gradT_z
  use global_vars, only : gradtk_x
  use global_vars, only : gradtk_y
  use global_vars, only : gradtk_z
  use global_vars, only : gradtw_x
  use global_vars, only : gradtw_y
  use global_vars, only : gradtw_z
  use global_vars, only : gradqp_z
  use global_vars, only : process_id
  use utils      , only : alloc
  use utils      , only : dealloc
  use utils      , only : dmsg
  use utils      , only : turbulence_read_error
  use string
  use evaluate_grad  , only : compute_gradient_G
  use evaluate_grad  , only : compute_gradient_T
  use ghost_gradients, only : apply_gradient_bc 

  implicit none
  private

  public :: evaluate_all_gradients

  contains

    subroutine evaluate_all_gradients()

      implicit none

      call dmsg(1, 'summon_grad_evaluation', 'evaluate_all_gradients')

      call compute_gradient_G(gradu_x, x_speed, 'x')
      call compute_gradient_G(gradv_x, y_speed, 'x')
      call compute_gradient_G(gradw_x, z_speed, 'x')
      call compute_gradient_T(gradT_x         , 'x')
      call compute_gradient_G(gradu_y, x_speed, 'y')
      call compute_gradient_G(gradv_y, y_speed, 'y')
      call compute_gradient_G(gradw_y, z_speed, 'y')
      call compute_gradient_T(gradT_y         , 'y')
      if(kmx>2) then
      call compute_gradient_G(gradu_z, x_speed, 'z')
      call compute_gradient_G(gradv_z, y_speed, 'z')
      call compute_gradient_G(gradw_z, z_speed, 'z')
      call compute_gradient_T(gradT_z         , 'z')
      else
       gradqp_z=0.0
      end if

      select case (turbulence)

        case ('none')
          !do nothing
          continue

        case ('sst')
          call compute_gradient_G(gradtk_x, tk, 'x')
          call compute_gradient_G(gradtw_x, tw, 'x')
          call compute_gradient_G(gradtk_y, tk, 'y')
          call compute_gradient_G(gradtw_y, tw, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtk_z, tk, 'z')
          call compute_gradient_G(gradtw_z, tw, 'z')
          end if

        case DEFAULT
          call turbulence_read_error()

      end select

      !applying boundary condition to gradients
      call apply_gradient_bc()

    end subroutine evaluate_all_gradients


end module summon_grad_evaluation
