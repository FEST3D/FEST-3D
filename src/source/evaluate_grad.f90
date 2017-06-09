module evaluate_grad
  !----------------------------------------------------------------------------
  !170698 - jatinder Pal Singh Sandhu
  ! Aim : general subroutine to computer cell center gradient
  !----------------------------------------------------------------------------

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx

  use global_vars, only : xn
  use global_vars, only : yn
  use global_vars, only : zn
  use global_vars, only : xA, yA, zA    !face area
  use global_vars, only : volume

  use global_vars, only : gm
  use global_vars, only : R_gas
  use global_vars, only : density
  use global_vars, only : pressure
  use utils      , only:   alloc
  use utils      , only: dealloc
  use utils      , only: dmsg
  use utils      , only: turbulence_read_error
  use string

  implicit none
  private

  public :: compute_gradient_G
  public :: compute_gradient_T

  contains

  subroutine compute_gradient_G(grad, var, dir)
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
        nx(1:imx  ,1:jmx-1,1:kmx-1) => xn(:,:,:,1)
        ny(1:imx-1,1:jmx  ,1:kmx-1) => yn(:,:,:,1)
        nz(1:imx-1,1:jmx-1,1:kmx  ) => zn(:,:,:,1)
      case('y')
        nx(1:imx  ,1:jmx-1,1:kmx-1) => xn(:,:,:,2)
        ny(1:imx-1,1:jmx  ,1:kmx-1) => yn(:,:,:,2)
        nz(1:imx-1,1:jmx-1,1:kmx  ) => zn(:,:,:,2)
      case('z')
        nx(1:imx  ,1:jmx-1,1:kmx-1) => xn(:,:,:,3)
        ny(1:imx-1,1:jmx  ,1:kmx-1) => yn(:,:,:,3)
        nz(1:imx-1,1:jmx-1,1:kmx  ) => zn(:,:,:,3)
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
                       )/(2*volume(i,j,k))
        end do
      end do
    end do

  end subroutine compute_gradient_G

  subroutine compute_gradient_T(grad, dir)

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
        nx(1:imx  ,1:jmx-1,1:kmx-1) => xn(:,:,:,1)
        ny(1:imx-1,1:jmx  ,1:kmx-1) => yn(:,:,:,1)
        nz(1:imx-1,1:jmx-1,1:kmx  ) => zn(:,:,:,1)
      case('y')
        nx(1:imx  ,1:jmx-1,1:kmx-1) => xn(:,:,:,2)
        ny(1:imx-1,1:jmx  ,1:kmx-1) => yn(:,:,:,2)
        nz(1:imx-1,1:jmx-1,1:kmx  ) => zn(:,:,:,2)
      case('z')
        nx(1:imx  ,1:jmx-1,1:kmx-1) => xn(:,:,:,3)
        ny(1:imx-1,1:jmx  ,1:kmx-1) => yn(:,:,:,3)
        nz(1:imx-1,1:jmx-1,1:kmx  ) => zn(:,:,:,3)
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
                       )/(2*volume(i,j,k))
        end do
      end do
    end do

  end subroutine compute_gradient_T


end module evaluate_grad
