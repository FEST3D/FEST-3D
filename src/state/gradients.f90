  !< allocate memory to laminar gradients if flow is viscous and
  !< allocate memory to tubulence gradients base upon the model being used
module gradients
  !< allocate memory to laminar gradients if flow is viscous and
  !< allocate memory to tubulence gradients base upon the model being used
  !------------------------------------------------------------------
  ! 170509  Jatinder Pal Singh Sandhu
  !         - first build
  !-------------------------------------------------------------------
#include "../error.h"
#include "../debug.h"

  use global_vars,  only : imx
  use global_vars,  only : jmx
  use global_vars,  only : kmx
  use global_vars,  only : mu_ref
  use global_vars,  only : turbulence 
  use global_vars,  only : transition 
  use global_vars,  only : n_grad
  use global_vars,  only : sst_n_grad
  use global_vars,  only : gradqp_x
  use global_vars,  only : gradqp_y
  use global_vars,  only : gradqp_z
  use global_vars,  only : process_id

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
  use sa_gradients ,      only : setup_sa_grad
  use sa_gradients ,      only : destroy_sa_grad
  use lctm2015_gradients ,      only : setup_lctm2015_grad
  use lctm2015_gradients ,      only : destroy_lctm2015_grad

  implicit none
  private

  public :: setup_gradients
  public :: destroy_gradients

  contains


    subroutine setup_gradients
      !< Memoery allocation to the gradient variables and 
      !< setup pointer to the slice to the main gradient variable
      !< based on the various models being used.

      implicit none

      DebugCall("setup_gradients")

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

          case('sa', 'saBC')
            call setup_sa_grad()

          case('sst', 'sst2003')
            call setup_sst_grad()

          case('kkl')
            call setup_kkl_grad()

          case DEFAULT
            !call turbulence_read_error()
            Fatal_error

        end select

        !Transition modeling
        select case(trim(transition))

          case('lctm2015')
            call setup_lctm2015_grad()

          case('none','bc')
            !do nothing
            continue

          case DEFAULT
            Fatal_error

        end Select

      end if
    end subroutine setup_gradients



    subroutine destroy_gradients
      !< Deallocate memoery and nullify pointers
      !< to the gradient variables.
      
      implicit none

      DebugCall("destroy_gradients")

      if(mu_ref/=0)then

        call destroy_memory()

        ! unlink laminar grad pointer
        call destroy_laminar_grad()

        !unlink turublent grad pointer
        select case (trim(turbulence))
          
          case('none')
            !do nothing
            continue

          case('sa', 'saBC')
            call destroy_sa_grad()

          case('sst', 'sst2003')
            call destroy_sst_grad()

          case('kkl')
            call destroy_kkl_grad()

          case DEFAULT
           ! call turbulence_read_error()
           Fatal_error

        end select

        !Transition modeling
        select case(trim(transition))

          case('lctm2015')
            call destroy_lctm2015_grad()

          case('none','bc')
            !do nothing
            continue

          case DEFAULT
            Fatal_error

        end Select

      end if

    end subroutine destroy_gradients



    subroutine get_n_grad()
      !< Set number of variables for which
      !< gradient is required based on the
      !< being used

      implicit none

      DebugCall("get_n_grad")

      select case (trim(turbulence))
        
        case('none')
          !do nothing
          continue

        case ('sa', 'saBC')
          n_grad = 5

        case('sst', 'sst2003')
          n_grad = 6

        case('kkl')
          n_grad = 6

        case DEFAULT
          !call turbulence_read_error()
          Fatal_error

      end select


      !Transition modeling
      select case(trim(transition))

        case('lctm2015')
          n_grad = n_grad + 1

        case('none','bc')
          n_grad = n_grad + 0

        case DEFAULT
          Fatal_error

      end Select

    end subroutine get_n_grad



    subroutine allocate_memory()
      !< allocating memory to the gradient variable being used

      implicit none

      DebugCall("allocate_memory")

      call alloc(gradqp_x, 0, imx, 0, jmx, 0, kmx, 1, n_grad, AErrMsg("gradqp_x"))
      call alloc(gradqp_y, 0, imx, 0, jmx, 0, kmx, 1, n_grad, AErrMsg("gradqp_y"))
      call alloc(gradqp_z, 0, imx, 0, jmx, 0, kmx, 1, n_grad, AErrMsg("gradqp_z"))

    end subroutine allocate_memory



    subroutine destroy_memory()
      !< deallocate memeory from the gradient variables

      implicit none

      DebugCall("deallocate_memory")

      call dealloc(gradqp_x)
      call dealloc(gradqp_y)
      call dealloc(gradqp_z)

    end subroutine destroy_memory
    


end module gradients
