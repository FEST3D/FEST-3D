  !< Setup and nullify pointers for LCTM2015 transition 
  !<model to the main array which stores gradient of all variables
module lctm2015_gradients
  !< Setup and nullify pointers for LCTM2015 transition 
  !<model to the main array which stores gradient of all variables
  !---------------------------------------------------------------
  ! 1705009  Jatinder Pal Singh Sandhu
  !          - first build
  ! aim : link LCTM2015 pointer to allocated memory for gradients
  !---------------------------------------------------------------

#include "../debug.h"
#include "../error.h"

  use global_vars,  only : process_id
  use global_vars,  only : transition
  use global_vars,  only : imx
  use global_vars,  only : jmx
  use global_vars,  only : kmx
  use global_vars,  only : n_grad
  use global_vars,  only : gradqp_x
  use global_vars,  only : gradqp_y
  use global_vars,  only : gradqp_z

  use global_vars,  only : gradtgm_x
  use global_vars,  only : gradtgm_y
  use global_vars,  only : gradtgm_z

  use utils,        only : dmsg
  implicit none
  private

  public :: setup_lctm2015_grad
  public :: destroy_lctm2015_grad

  contains

    subroutine setup_lctm2015_grad()
      !< setup Pointer to the main array which stores gradient 
      !< all variables with x, y, z

      implicit none

      DebugCall('setup_sst_grad')

      select case(trim(transition))

        case('lctm2015')
          gradtgm_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, n_grad)
          gradtgm_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, n_grad)
          gradtgm_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, n_grad)

        case( 'bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error

      end Select

    end subroutine setup_lctm2015_grad


    subroutine destroy_lctm2015_grad()
      !< nullify all the gradient pointer setup for k-kL model
      implicit none

      DebugCall('destroy_sst_grad')

      select case(trim(transition))

        case('lctm2015')
          nullify(gradtgm_x)
          nullify(gradtgm_y)
          nullify(gradtgm_z)

        case('bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error

      end Select

    end subroutine destroy_lctm2015_grad


end module lctm2015_gradients
