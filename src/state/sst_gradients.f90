  !< Setup and nullify pointers for SST model to the main array which stores gradient of all variables
module sst_gradients
  !< Setup and nullify pointers for SST model to the main array which stores gradient of all variables
  !---------------------------------------------------------------
  ! 1705009  Jatinder Pal Singh Sandhu
  !          - first build
  ! aim : link sst pointer to allocated memory for gradients
  !---------------------------------------------------------------

#include "../debug.h"
#include "../error.h"

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
  use global_vars,  only : gradv_x 
  use global_vars,  only : gradv_y
  use global_vars,  only : gradv_z
  use global_vars,  only : gradtw_x
  use global_vars,  only : gradtw_y
  use global_vars,  only : gradtw_z
  use global_vars,  only : gradT_x
  use global_vars,  only : gradT_y
  use global_vars,  only : gradT_z
  use global_vars,  only : gradtgm_x
  use global_vars,  only : gradtgm_y
  use global_vars,  only : gradtgm_z

  use utils,        only : dmsg
  implicit none
  private

  public :: setup_sst_grad
  public :: destroy_sst_grad

  contains

    subroutine setup_sst_grad()
      !< Setup Pointer to the main array which stores gradient 
      !< all variables with x, y, z

      implicit none

      DebugCall('setup_sst_grad')

      gradtk_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 5)
      gradtw_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 6)

      gradtk_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 5)
      gradtw_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 6)

      gradtk_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 5)
      gradtw_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 6)

    end subroutine setup_sst_grad


    subroutine destroy_sst_grad()
      !< Nullify all the gradient pointer setup for k-kL model
      implicit none

      DebugCall('destroy_sst_grad')

      nullify(gradtk_x)
      nullify(gradtw_x)

      nullify(gradtk_y)
      nullify(gradtw_y)

      nullify(gradtk_z)
      nullify(gradtw_z)


    end subroutine destroy_sst_grad


end module sst_gradients
