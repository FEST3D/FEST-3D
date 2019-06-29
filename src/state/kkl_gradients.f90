  !< Setup and nullify pointers for k-kL model to the main array which stores gradient of all variables
module kkl_gradients
  !< Setup and nullify pointers for k-kL model to the main array which stores gradient of all variables
  !---------------------------------------------------------------
  ! 1705009  Jatinder Pal Singh Sandhu
  !          - first build
  ! aim : link kkl pointer to allocated memory for gradients
  !---------------------------------------------------------------
#include "../debug.h"

  use global_vars,  only : process_id
  use global_vars,  only : imx
  use global_vars,  only : jmx
  use global_vars,  only : kmx
  use global_vars,  only : gradqp_x
  use global_vars,  only : gradqp_y
  use global_vars,  only : gradqp_z

  use global_vars,  only : gradtk_x
  use global_vars,  only : gradtk_y
  use global_vars,  only : gradtk_z 
  use global_vars,  only : gradtkl_x
  use global_vars,  only : gradtkl_y
  use global_vars,  only : gradtkl_z

  use utils,        only : dmsg
  implicit none
  private

  public :: setup_kkl_grad
  public :: destroy_kkl_grad

  contains

    subroutine setup_kkl_grad()
      !< Setup Pointer to the main array which stores gradient 
      !< all variables with x, y, z

      implicit none

      DebugCall('setup_kkl_grad')

      gradtk_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 5)
      gradtkl_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 6)

      gradtk_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 5)
      gradtkl_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 6)

      gradtk_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 5)
      gradtkl_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 6)

    end subroutine setup_kkl_grad


    subroutine destroy_kkl_grad()
      !< Nullify all the gradient pointer setup for k-kL model

      implicit none

      DebugCall('destroy_kkl_grad')

      nullify(gradtk_x)
      nullify(gradtkl_x)

      nullify(gradtk_y)
      nullify(gradtkl_y)

      nullify(gradtk_z)
      nullify(gradtkl_z)

    end subroutine destroy_kkl_grad


end module kkl_gradients
