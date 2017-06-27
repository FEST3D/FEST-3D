module gradients
  !------------------------------------------------------------------
  ! 170509  Jatinder Pal Singh Sandhu
  !         - first build
  !aim : 1) allocate memory to laminar gradients if flow is viscous
  !      2) allocate memory to tubulence gradients base on model used
  !-------------------------------------------------------------------

  use global_vars,  only : imx
  use global_vars,  only : jmx
  use global_vars,  only : kmx
  use global_vars,  only : mu_ref
  use global_vars,  only : turbulence 
  use global_vars,  only : n_grad
  use global_vars,  only : sst_n_grad
  use global_vars,  only : gradqp_x
  use global_vars,  only : gradqp_y
  use global_vars,  only : gradqp_z

  use utils,        only : dmsg
  use utils,        only : dealloc
  use utils,        only : alloc
  use utils,        only : turbulence_read_error

  use laminar_gradients,  only : setup_laminar_grad
  use laminar_gradients,  only : destroy_laminar_grad

  use sst_gradients,      only : setup_sst_grad
  use sst_gradients,      only : destroy_sst_grad
  use kkl_gradients,      only : setup_kkl_grad
  use kkl_gradients,      only : destroy_kkl_grad

  implicit none
  private

  public :: setup_gradients
  public :: destroy_gradients

  contains

    subroutine setup_gradients
      implicit none

      if(mu_ref/=0)then

        call get_n_grad()
        call allocate_memory()

        ! linking pointer to laminar gradients
        call setup_laminar_grad()

        ! linking pointer to turbulent gradients
        select case (trim(turbulence))
          
          case('none')
            !do nothing
            continue

          case('sst')
            call setup_sst_grad()

          case('kkl')
            call setup_kkl_grad()

          case DEFAULT
            call turbulence_read_error()

        end select

      end if
    end subroutine setup_gradients

    subroutine destroy_gradients
      implicit none

      if(mu_ref/=0)then

        call destroy_memory()

        ! unlink laminar grad pointer
        call destroy_laminar_grad()

        !unlink turublent grad pointer
        select case (trim(turbulence))
          
          case('none')
            !do nothing
            continue

          case('sst')
            call destroy_sst_grad()

          case('kkl')
            call destroy_kkl_grad()

          case DEFAULT
            call turbulence_read_error()

        end select

      end if

    end subroutine destroy_gradients


    subroutine get_n_grad()
      implicit none

      select case (trim(turbulence))
        
        case('none')
          !do nothing
          continue

        case('sst')
          n_grad = n_grad + sst_n_grad

        case('kkl')
          n_grad = n_grad + 2

        case DEFAULT
          call turbulence_read_error()

      end select

    end subroutine get_n_grad


    subroutine allocate_memory()
      implicit none

      call dmsg(1, 'gradients', 'allocate_memory')

      call alloc(gradqp_x, 0, imx, 0, jmx, 0, kmx, 1, n_grad,&
              errmsg='Error: Unable to allocate memory for ' // &
                  'Gradu_x - gradients')
      call alloc(gradqp_y, 0, imx, 0, jmx, 0, kmx, 1, n_grad, &
              errmsg='Error: Unable to allocate memory for ' // &
                  'Gradu_y - gradients')
      call alloc(gradqp_z, 0, imx, 0, jmx, 0, kmx, 1, n_grad, &
              errmsg='Error: Unable to allocate memory for ' // &
                  'Gradu_z - gradients')

    end subroutine allocate_memory


    subroutine destroy_memory()
      implicit none

      call dmsg(1, 'gradients', 'destroy_memory')

      call dealloc(gradqp_x)
      call dealloc(gradqp_y)
      call dealloc(gradqp_z)

    end subroutine destroy_memory
    

end module gradients
