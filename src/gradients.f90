  !< Allocate memory to laminar gradients if flow is viscous and
  !< allocate memory to tubulence gradients base upon the model being used
module gradients
  !< Allocate memory to laminar gradients if flow is viscous and
  !< allocate memory to tubulence gradients base upon the model being used
  !------------------------------------------------------------------
  ! 170509  Jatinder Pal Singh Sandhu
  !         - first build
  !-------------------------------------------------------------------
#include "error.h"
#include "debug.h"
  use vartypes
  use utils,        only : alloc


  implicit none
  !private
  ! gradients
  integer                                           :: n_grad=4 
  !< Number of variable to store gradient for
  real(wp), dimension(:, :, :, :), allocatable, target  :: gradqp_x 
  !< Store gradient of n_grad variables with respect to direction x
  real(wp), dimension(:, :, :, :), allocatable, target  :: gradqp_y 
  !< Store gradient of n_grad variables with respect to direction y
  real(wp), dimension(:, :, :, :), allocatable, target  :: gradqp_z 
  !< Store gradient of n_grad variables with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradu_x  
  !< Gradient of variable U with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradu_y  
  !< Gradient of variable U with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradu_z  
  !< Gradient of variable U with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradv_x  
  !< Gradient of variable V with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradv_y  
  !< Gradient of variable V with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradv_z  
  !< Gradient of variable V with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradw_x  
  !< Gradient of variable W with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradw_y  
  !< Gradient of variable W with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradw_z  
  !< Gradient of variable W with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradT_x  
  !< Gradient of variable Temperature with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradT_y  
  !< Gradient of variable Temperature with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradT_z  
  !< Gradient of variable Temperature with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradtk_x 
  !< Gradient of variable turbulent kinetic energy with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradtk_y 
  !< Gradient of variable turbulent kinetic energy with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradtk_z 
  !< Gradient of variable turbulent kinetic energy with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradtw_x 
  !< Gradient of variable dissipation rate with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradtw_y 
  !< Gradient of variable dissipation rate with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradtw_z 
  !< Gradient of variable dissipation rate with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradtkl_x
  !< Gradient of variable kL with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradtkl_y
  !< Gradient of variable kL with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradtkl_z
  !< Gradient of variable kL with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradte_x 
  !< Gradient of variable turbulent energy dissiaption with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradte_y 
  !< Gradient of variable turbulent energy dissiaption with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradte_z 
  !< Gradient of variable turbulent energy dissiaption with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradtv_x 
  !< Gradient of variable turbulenct visocity(SA mode) with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradtv_y 
  !< Gradient of variable turbulenct visocity(SA mode) with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradtv_z 
  !< Gradient of variable turbulenct visocity(SA mode) with respect to direction z
  real(wp), dimension(:, :, :),                 pointer :: gradtgm_x
  !< Gradient of variable intermittency with respect to direction x
  real(wp), dimension(:, :, :),                 pointer :: gradtgm_y
  !< Gradient of variable intermittency with respect to direction y
  real(wp), dimension(:, :, :),                 pointer :: gradtgm_z
  !< Gradient of variable intermittency with respect to direction z
  real(wp) :: R_gas

  integer :: imx, jmx, kmx, n_var

  type :: singlesub
    integer :: imin, imax, il, iu
    integer :: jmin, jmax, jl, ju
    integer :: kmin, kmax, kl, ku
    integer :: sig=1
  end type singlesub

  public :: setup_gradients
  public :: evaluate_all_gradients
  !public :: destroy_gradients

  contains


    subroutine setup_gradients(control, scheme, flow, dims)
      !< Memoery allocation to the gradient variables and 
      !< setup pointer to the slice to the main gradient variable
      !< based on the various models being used.

      implicit none
      type(controltype), intent(in) :: control
      !< Control parameters
      type(schemetype) , intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype)   , intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx

      DebugCall("setup_gradients")

        imx = dims%imx
        jmx = dims%jmx
        kmx = dims%kmx

        n_var = control%n_var
        R_gas = flow%R_gas

      if(flow%mu_ref/=0)then

        call get_n_grad(scheme)
        !call allocate_memory()
        call alloc(gradqp_x, 0, imx, 0, jmx, 0, kmx, 1, n_grad, AErrMsg("gradqp_x"))
        call alloc(gradqp_y, 0, imx, 0, jmx, 0, kmx, 1, n_grad, AErrMsg("gradqp_y"))
        call alloc(gradqp_z, 0, imx, 0, jmx, 0, kmx, 1, n_grad, AErrMsg("gradqp_z"))

        ! Linking pointer to laminar gradients
        !call setup_laminar_grad()
        DebugCall('setup_laminar_grad')
        !< Setup Pointer to the main array which stores gradient 
        !< all variables with x, y, z

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

        ! Linking pointer to turbulent gradients
        select case (trim(scheme%turbulence))
          
          case('none')
            !do nothing
            continue

          case('sa', 'saBC')
            !call setup_sa_grad()
            DebugCall("setup_sa_grad")
            gradtv_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 5)
            gradtv_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 5)
            gradtv_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 5)

          case('sst', 'sst2003')
            !< Setup Pointer to the main array which stores gradient 
            !< all variables with x, y, z
            DebugCall('setup_sst_grad')

            gradtk_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 5)
            gradtw_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 6)

            gradtk_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 5)
            gradtw_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 6)

            gradtk_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 5)
            gradtw_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 6)

          case('kkl')
            !call setup_kkl_grad()
            DebugCall('setup_kkl_grad')

            gradtk_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 5)
            gradtkl_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, 6)

            gradtk_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 5)
            gradtkl_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, 6)

            gradtk_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 5)
            gradtkl_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, 6)

          case DEFAULT
            Fatal_error

        end select

        !Transition modeling
        select case(trim(scheme%transition))

          case('lctm2015')
            !call setup_lctm2015_grad()
            gradtgm_x(0:imx, 0:jmx, 0:kmx) => gradqp_x(:, :, :, n_grad)
            gradtgm_y(0:imx, 0:jmx, 0:kmx) => gradqp_y(:, :, :, n_grad)
            gradtgm_z(0:imx, 0:jmx, 0:kmx) => gradqp_z(:, :, :, n_grad)

          case('none','bc')
            !do nothing
            continue

          case DEFAULT
            Fatal_error

        end Select

      end if
    end subroutine setup_gradients



    subroutine get_n_grad(scheme)
      !< Set number of variables for which
      !< gradient is required based on the
      !< being used

      implicit none
      type(schemetype) , intent(in) :: scheme
      !< finite-volume Schemes

      DebugCall("get_n_grad")

      select case (trim(scheme%turbulence))
        
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
          Fatal_error

      end select


      !Transition modeling
      select case(trim(scheme%transition))

        case('lctm2015')
          n_grad = n_grad + 1

        case('none','bc')
          n_grad = n_grad + 0

        case DEFAULT
          Fatal_error

      end Select

    end subroutine get_n_grad


    subroutine evaluate_all_gradients(qp, Temp, cells, Ifaces, Jfaces, Kfaces, scheme, bc, dims)
      !< Call to all the required gradients and 
      !< apply boundary condition for ghost cell
      !< gradients

      implicit none
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), Target :: qp
      !< Store primitive variable at cell center
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2), intent(in) :: Temp
      !< Store Temperature variable at cell center
      type(schemetype) , intent(in) :: scheme
      !< finite-volume Schemes
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      type(boundarytype), intent(in) :: bc
      !< boundary conditions and fixed values
      real(wp), dimension(:, :, :), pointer :: x_speed      
       !< U pointer, point to slice of qp (:,:,:,2) 
      real(wp), dimension(:, :, :), pointer :: y_speed      
       !< V pointer, point to slice of qp (:,:,:,3) 
      real(wp), dimension(:, :, :), pointer :: z_speed      
       !< W pointer, point to slice of qp (:,:,:,4)
      real(wp), dimension(:, :, :), pointer :: tk   !< TKE/mass
      real(wp), dimension(:, :, :), pointer :: tw   !< Omega
      real(wp), dimension(:, :, :), pointer :: tv   !< SA visocity
      real(wp), dimension(:, :, :), pointer :: tkl  !< KL K-KL method
      real(wp), dimension(:, :, :), pointer :: tgm  !< Intermittency of LCTM2015

      DebugCall('evaluate_all_gradients')

      x_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 2)
      y_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 3)
      z_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 4)

      !call compute_gradient_G(gradu_x, qp(:,:,:,2), cells, Ifaces, Jfaces, Kfaces, dims, 'x')
      call compute_gradient_G(gradu_x, x_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
      call compute_gradient_G(gradv_x, y_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
      call compute_gradient_G(gradw_x, z_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
      call compute_gradient_G(gradT_x, Temp,    cells, Ifaces, Jfaces, Kfaces, dims, 'x')
      !call compute_gradient_G(gradu_y, qp(:,:,:,2), cells, Ifaces, Jfaces, Kfaces, dims, 'y')
      call compute_gradient_G(gradu_y, x_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
      call compute_gradient_G(gradv_y, y_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
      call compute_gradient_G(gradw_y, z_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
      call compute_gradient_G(gradT_y, Temp   , cells, Ifaces, Jfaces, Kfaces, dims, 'y')
      if(dims%kmx>2) then
        call compute_gradient_G(gradu_z, x_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
        call compute_gradient_G(gradv_z, y_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
        call compute_gradient_G(gradw_z, z_speed, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
        !call compute_gradient_T(gradT_z         , cells, Ifaces, Jfaces, Kfaces, dims, 'z')
        call compute_gradient_G(gradT_z, Temp   , cells, Ifaces, Jfaces, Kfaces, dims, 'z')
      else
       gradqp_z=0.0
      end if

      select case (trim(scheme%turbulence))

        case ('none')
          !do nothing
          continue

        case ('sa', 'saBC')
          tv(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 6)
          call compute_gradient_G(gradtv_x, tv, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
          call compute_gradient_G(gradtv_y, tv, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtv_z, tv, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
          end if

        case ('sst', 'sst2003')
          tk(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 6)
          tw(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 7)
          call compute_gradient_G(gradtk_x, tk, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
          call compute_gradient_G(gradtw_x, tw, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
          call compute_gradient_G(gradtk_y, tk, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
          call compute_gradient_G(gradtw_y, tw, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtk_z, tk, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
          call compute_gradient_G(gradtw_z, tw, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
          end if

        case ('kkl')
          tk(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 6)
          tkl(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 7)
          call compute_gradient_G(gradtk_x , tk , cells, Ifaces, Jfaces, Kfaces, dims, 'x')
          call compute_gradient_G(gradtkl_x, tkl, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
          call compute_gradient_G(gradtk_y , tk , cells, Ifaces, Jfaces, Kfaces, dims, 'y')
          call compute_gradient_G(gradtkl_y, tkl, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtk_z , tk , cells, Ifaces, Jfaces, Kfaces, dims, 'z')
          call compute_gradient_G(gradtkl_z, tkl, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
          end if

        case DEFAULT
          Fatal_error

      end select


      select case(trim(scheme%transition))
        case('lctm2015')
          tgm(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 8)
          call compute_gradient_G(gradtgm_x, tgm, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
          call compute_gradient_G(gradtgm_y, tgm, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
          if(kmx>2)then
          call compute_gradient_G(gradtgm_z, tgm, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
          end if

        case('bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error

      end Select

      call apply_gradient_bc(qp, temp, cells, Ifaces, Jfaces, Kfaces, bc, dims)

    end subroutine evaluate_all_gradients


    subroutine compute_gradient_G(grad, var, cells, Ifaces, Jfaces, Kfaces, dims, dir)
      !<  Compute gradient of any input scalar
      implicit none
      type(extent), intent(in) :: dims
      real(wp), dimension( 0:dims%imx  , 0:dims%jmx  , 0:dims%kmx  ), intent(out) :: grad
      !< Output variable storing the graident of var
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: var
      !< Input variable of which graident is required
      character(len=*)                           , intent(in) :: dir
      !< Direction with respect to which gradients are calculated
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal

      integer :: i
      integer :: j
      integer :: k

      DebugCall('compute_gradient_G')
      grad(:,:,:) = 0.0
      select case(dir)
        case('x')
          do k=0,dims%kmx
            do j=0,dims%jmx
              do i=0,dims%imx
                grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*Ifaces(i,j,k)%nx*Ifaces(i,j,k)%A &
                              -(var(i  ,j-1,k  )+var(i,j,k))*Jfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                              -(var(i  ,j  ,k-1)+var(i,j,k))*Kfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                              +(var(i+1,j  ,k  )+var(i,j,k))*Ifaces(i+1,j  ,k  )%nx*Ifaces(i+1,j  ,k  )%A &
                              +(var(i  ,j+1,k  )+var(i,j,k))*Jfaces(i  ,j+1,k  )%nx*Jfaces(i  ,j+1,k  )%A &
                              +(var(i  ,j  ,k+1)+var(i,j,k))*Kfaces(i  ,j  ,k+1)%nx*Kfaces(i  ,j  ,k+1)%A &
                             )/(2*cells(i,j,k)%volume)
              end do
            end do
          end do
        case('y')
          do k=0,dims%kmx
            do j=0,dims%jmx
              do i=0,dims%imx
                grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*Ifaces(i,j,k)%ny*Ifaces(i,j,k)%A &
                              -(var(i  ,j-1,k  )+var(i,j,k))*Jfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                              -(var(i  ,j  ,k-1)+var(i,j,k))*Kfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                              +(var(i+1,j  ,k  )+var(i,j,k))*Ifaces(i+1,j  ,k  )%ny*Ifaces(i+1,j  ,k  )%A &
                              +(var(i  ,j+1,k  )+var(i,j,k))*Jfaces(i  ,j+1,k  )%ny*Jfaces(i  ,j+1,k  )%A &
                              +(var(i  ,j  ,k+1)+var(i,j,k))*Kfaces(i  ,j  ,k+1)%ny*Kfaces(i  ,j  ,k+1)%A &
                             )/(2*cells(i,j,k)%volume)
              end do
            end do
          end do
        case('z')
          do k=0,dims%kmx
            do j=0,dims%jmx
              do i=0,dims%imx
                grad(i,j,k) =(-(var(i-1,j  ,k  )+var(i,j,k))*Ifaces(i,j,k)%nz*Ifaces(i,j,k)%A &
                              -(var(i  ,j-1,k  )+var(i,j,k))*Jfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                              -(var(i  ,j  ,k-1)+var(i,j,k))*Kfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                              +(var(i+1,j  ,k  )+var(i,j,k))*Ifaces(i+1,j  ,k  )%nz*Ifaces(i+1,j  ,k  )%A &
                              +(var(i  ,j+1,k  )+var(i,j,k))*Jfaces(i  ,j+1,k  )%nz*Jfaces(i  ,j+1,k  )%A &
                              +(var(i  ,j  ,k+1)+var(i,j,k))*Kfaces(i  ,j  ,k+1)%nz*Kfaces(i  ,j  ,k+1)%A &
                             )/(2*cells(i,j,k)%volume)
              end do
            end do
          end do
        case DEFAULT
          print*, "ERROR: gradient direction error"
          Fatal_error
      end select

      if(any(isnan(grad)))then
        Fatal_error
      end if

    end subroutine compute_gradient_G



    subroutine apply_gradient_bc(qp, temp, cells, Ifaces, Jfaces, Kfaces, bc, dims)
      !< Call same subroutine for all the face
      !< Apply/set value of all gradient in the ghost cells
      !< gradqp_G = (qp_I - qp_G)*Area_W*unit_normal_G/(volume_G)
      !< volume_G = volume_I
      !-----------------------------------------------------------
      implicit none
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Input variable of which graident is required
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2), intent(in) :: Temp
      !< Intput Temperature variable
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      type(boundarytype), intent(in) :: bc
      type(singlesub) :: domain

      DebugCall('apply_gradient_bc')


      domain%imin = 1
      domain%imax = 1
      domain%jmin = 1
      domain%jmax = dims%jmx-1
      domain%kmin = 1
      domain%kmax = dims%kmx-1
      domain%il   = 1; domain%jl = 0; domain%kl = 0
      domain%iu   = 0; domain%ju = 0; domain%ku = 0
      domain%sig  = 1
      call apply_gradient_bc_face(qp, temp, cells, Ifaces, dims, domain, bc%imin_id, bc%fixed_wall_temperature(1))

      domain%imin = dims%imx
      domain%imax = dims%imx
      domain%jmin = 1
      domain%jmax = dims%jmx-1
      domain%kmin = 1
      domain%kmax = dims%kmx-1
      domain%il   = 0; domain%jl = 0; domain%kl = 0
      domain%iu   = 1; domain%ju = 0; domain%ku = 0
      domain%sig  = -1
      call apply_gradient_bc_face(qp, temp, cells, Ifaces, dims, domain, bc%imax_id, bc%fixed_wall_temperature(2))

      domain%imin = 1
      domain%imax = dims%imx-1
      domain%jmin = 1
      domain%jmax = 1
      domain%kmin = 1
      domain%kmax = dims%kmx-1
      domain%il   = 0; domain%jl = 1; domain%kl = 0
      domain%iu   = 0; domain%ju = 0; domain%ku = 0
      domain%sig  = 1
      call apply_gradient_bc_face(qp, temp, cells, Jfaces, dims, domain, bc%jmin_id, bc%fixed_wall_temperature(3))

      domain%imin = 1
      domain%imax = dims%imx-1
      domain%jmin = dims%jmx
      domain%jmax = dims%jmx
      domain%kmin = 1
      domain%kmax = dims%kmx-1
      domain%il   = 0; domain%jl = 0; domain%kl = 0
      domain%iu   = 0; domain%ju = 1; domain%ku = 0
      domain%sig  = -1
      call apply_gradient_bc_face(qp, temp, cells, Jfaces, dims, domain, bc%jmax_id, bc%fixed_wall_temperature(4))

      domain%imin = 1
      domain%imax = dims%imx-1
      domain%jmin = 1
      domain%jmax = dims%jmx-1
      domain%kmin = 1
      domain%kmax = 1
      domain%il   = 0; domain%jl = 0; domain%kl = 1
      domain%iu   = 0; domain%ju = 0; domain%ku = 0
      domain%sig  = 1
      call apply_gradient_bc_face(qp, temp, cells, Kfaces, dims, domain, bc%kmin_id, bc%fixed_wall_temperature(5))

      domain%imin = 1
      domain%imax = dims%imx-1
      domain%jmin = 1
      domain%jmax = dims%jmx-1
      domain%kmin = dims%kmx
      domain%kmax = dims%kmx
      domain%il   = 0; domain%jl = 0; domain%kl = 0
      domain%iu   = 0; domain%ju = 0; domain%ku = 1
      domain%sig  = -1
      call apply_gradient_bc_face(qp, temp, cells, Kfaces, dims, domain, bc%kmax_id, bc%fixed_wall_temperature(6))
          
    end subroutine apply_gradient_bc

    

    subroutine apply_gradient_bc_face(qp, temp, cells, faces, dims, domain, bc_id, fixed_temp)
      !< Call same subroutine for all the face
      !< Apply/set value of all gradient in the ghost cells
      !< gradqp_G = (qp_I - qp_G)*Area_W*unit_normal_G/(volume_G)
      !< volume_G = volume_I
      !-----------------------------------------------------------
      implicit none
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(singlesub), intent(in) :: domain
      !< flags for direction 
      real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Input variable of which graident is required
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2), intent(in) :: Temp
      !< Intput Temperature variable
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: faces
      !< Input varaible which stores any(I,J,K) faces' area and unit normal
      integer, intent(in) :: bc_id
      real(wp), intent(in) :: fixed_temp
      real(wp), dimension(n_grad) :: qp_I
      real(wp), dimension(n_grad) :: qp_G
      real(wp)    :: T_I
      real(wp)    :: T_G 
      real(wp)    :: c_x
      real(wp)    :: c_y
      real(wp)    :: c_z
      integer :: i,j,k,l
      real(wp)    :: nx
      real(wp)    :: ny
      real(wp)    :: nz
      real(wp)    :: dot
      integer  :: il, jl, kl
      integer  :: iu, ju, ku

      il = domain%il
      jl = domain%jl
      kl = domain%kl
      iu = domain%iu
      ju = domain%ju
      ku = domain%ku

      do k = domain%kmin, domain%kmax
        do j = domain%jmin, domain%jmax
          do i = domain%imin, domain%imax
            nx   = faces(i,j,k)%nx
            ny   = faces(i,j,k)%ny
            nz   = faces(i,j,k)%nz
            c_x  = faces(i,j,k)%A*nx/cells(i-iu,j-ju,k-ku)%volume
            c_y  = faces(i,j,k)%A*ny/cells(i-iu,j-ju,k-ku)%volume
            c_z  = faces(i,j,k)%A*nz/cells(i-iu,j-ju,k-ku)%volume
            T_I  = Temp(i-iu,j-ju,k-ku)
            T_G  = Temp(i-il,j-jl,k-kl)
            qp_I = qp(i-iu,j-ju,k-ku,2:dims%n_var)
            qp_G = qp(i-il,j-jl,k-kl,2:dims%n_var)

            ! normal component of gradient
            gradqp_x(i-il,j-jl,k-kl,:) = domain%sig*(qp_I - qp_G)*c_x 
            gradqp_y(i-il,j-jl,k-kl,:) = domain%sig*(qp_I - qp_G)*c_y
            gradqp_z(i-il,j-jl,k-kl,:) = domain%sig*(qp_I - qp_G)*c_z
            gradqp_x(i-il,j-jl,k-kl,4) = domain%sig*( T_I -  T_G)*c_x
            gradqp_y(i-il,j-jl,k-kl,4) = domain%sig*( T_I -  T_G)*c_y
            gradqp_z(i-il,j-jl,k-kl,4) = domain%sig*( T_I -  T_G)*c_z
            if(bc_id==-5 .and. (fixed_temp<1. .and. fixed_temp>=0.))then
              gradqp_x(i-il,j-jl,k-kl,4) = -gradqp_x(i-iu,j-ju,k-ku,4)
              gradqp_y(i-il,j-jl,k-kl,4) = -gradqp_y(i-iu,j-ju,k-ku,4)
              gradqp_z(i-il,j-jl,k-kl,4) = -gradqp_z(i-iu,j-ju,k-ku,4)

            end if
            !parallel component of gradient
            do l=1,n_grad
              dot = (gradqp_x(i-iu,j-ju,k-ku,l)*nx) + (gradqp_y(i-iu,j-ju,k-ku,l)*ny) + (gradqp_z(i-iu,j-ju,k-ku,l)*nz)
              gradqp_x(i-il,j-jl,k-kl,l) = gradqp_x(i-il,j-jl,k-kl,l) + (gradqp_x(i-iu,j-ju,k-ku,l) - dot*nx)
              gradqp_y(i-il,j-jl,k-kl,l) = gradqp_y(i-il,j-jl,k-kl,l) + (gradqp_y(i-iu,j-ju,k-ku,l) - dot*ny)
              gradqp_z(i-il,j-jl,k-kl,l) = gradqp_z(i-il,j-jl,k-kl,l) + (gradqp_z(i-iu,j-ju,k-ku,l) - dot*nz)
            end do
          end do
        end do
      end do

    end subroutine apply_gradient_bc_face
end module gradients
