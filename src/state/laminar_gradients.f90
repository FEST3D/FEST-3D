  !< Link/Pointer to the gradient of u,v,w, and Temperature
module laminar_gradients
  !< Link/Pointer to the gradient of u,v,w, and Temperature
  !---------------------------------------------------------------
  ! 1705009  Jatinder Pal Singh Sandhu
  !          - first build
  ! aim : link laminar pointer to allocated memory for gradients
  !---------------------------------------------------------------

#include "../debug.h"
  
  use global_vars,  only : process_id
  use global_vars,  only : imx
  use global_vars,  only : jmx
  use global_vars,  only : kmx
  use global_vars,  only : gradqp_x
  use global_vars,  only : gradqp_y
  use global_vars,  only : gradqp_z

  use global_vars,  only : n_grad
  use global_vars,  only : gradu_x
  use global_vars,  only : gradu_y
  use global_vars,  only : gradu_z 
  use global_vars,  only : gradv_x 
  use global_vars,  only : gradv_y
  use global_vars,  only : gradv_z
  use global_vars,  only : gradw_x
  use global_vars,  only : gradw_y
  use global_vars,  only : gradw_z
  use global_vars,  only : gradT_x
  use global_vars,  only : gradT_y
  use global_vars,  only : gradT_z

  use utils,        only : dmsg

  implicit none

  contains

    subroutine setup_laminar_grad()
      !< Setup pointer to the gradient of U, v, w, Temperature
      !< with respect to x, y, and z

      implicit none

      DebugCall('setup_laminar_grad')

      gradu_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 1)
      gradv_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 2)
      gradw_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 3)
      gradT_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 4)

      gradu_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 1)
      gradv_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 2)
      gradw_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 3)
      gradT_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 4)

      gradu_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 1)
      gradv_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 2)
      gradw_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 3)
      gradT_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 4)

    end subroutine setup_laminar_grad


    subroutine destroy_laminar_grad()
      !< Unlink the laminar gradient pointers

      implicit none

      DebugCall('destroy_laminar_grad')

      nullify(gradu_x)
      nullify(gradv_x)
      nullify(gradw_x)
      nullify(gradT_x)

      nullify(gradu_y)
      nullify(gradv_y)
      nullify(gradw_y)
      nullify(gradT_y)

      nullify(gradu_z)
      nullify(gradv_z)
      nullify(gradw_z)
      nullify(gradT_z)

    end subroutine destroy_laminar_grad


end module laminar_gradients
