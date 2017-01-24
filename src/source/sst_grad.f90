module sst_grad

  use global_vars  , only : imx
  use global_vars  , only : jmx
  use global_vars  , only : kmx
  use global_vars  , only : xnx
  use global_vars  , only : xny
  use global_vars  , only : xnz
  use global_vars  , only : ynx
  use global_vars  , only : yny
  use global_vars  , only : ynz
  use global_vars  , only : znx
  use global_vars  , only : zny
  use global_vars  , only : znz
  use global_vars  , only : xA
  use global_vars  , only : yA
  use global_vars  , only : zA
  use global_vars  , only : volume
  use global_vars  , only : gradtk_x
  use global_vars  , only : gradtk_y
  use global_vars  , only : gradtk_z
  use global_vars  , only : gradtw_x
  use global_vars  , only : gradtw_y
  use global_vars  , only : gradtw_z

  use face_interpolant , only : x_tk_left
  use face_interpolant , only : y_tk_left
  use face_interpolant , only : z_tk_left
  use face_interpolant , only : x_tk_right
  use face_interpolant , only : y_tk_right
  use face_interpolant , only : z_tk_right
  use face_interpolant , only : x_tw_left
  use face_interpolant , only : y_tw_left
  use face_interpolant , only : z_tw_left
  use face_interpolant , only : x_tw_right
  use face_interpolant , only : y_tw_right
  use face_interpolant , only : z_tw_right

  use utils , only : alloc
  use utils , only : dealloc
  use utils , only : dmsg
  use string
  implicit none
  private

  public :: setup_sst_grad
  public :: destroy_sst_grad
  public :: calculate_sst_grad


  contains

    subroutine setup_sst_grad()
      implicit none
      character(len=*), parameter :: err = "ERROR: Unable to allocate memory to "

      call dmsg(1, 'sst_grad', 'setup_sst_grad')

      call alloc(gradtk_x, 0,imx, 0,jmx, 0,kmx, errmsg=err//'Gradtk_x')
      call alloc(gradtk_y, 0,imx, 0,jmx, 0,kmx, errmsg=err//'Gradtk_y')
      call alloc(gradtk_z, 0,imx, 0,jmx, 0,kmx, errmsg=err//'Gradtk_z')

      call alloc(gradtw_x, 0,imx, 0,jmx, 0,kmx, errmsg=err//'Gradtw_x')
      call alloc(gradtw_y, 0,imx, 0,jmx, 0,kmx, errmsg=err//'Gradtw_y')
      call alloc(gradtw_z, 0,imx, 0,jmx, 0,kmx, errmsg=err//'Gradtw_z')

    end subroutine setup_sst_grad

    subroutine destroy_sst_grad()
      implicit none

      call dmsg(1, 'sst_grad', 'destroy_sst_grad')

      call dealloc(gradtk_x)
      call dealloc(gradtk_y)
      call dealloc(gradtk_z)

      call dealloc(gradtw_x)
      call dealloc(gradtw_y)
      call dealloc(gradtw_z)

    end subroutine destroy_sst_grad

    subroutine calculate_sst_grad()
      implicit none
      integer :: i,j,k
      real :: tk_face
      real :: tw_face

      call dmsg(1, "sst_graf", "calculate_sst_grad")

      gradtk_x=0.
      gradtk_y=0.
      gradtk_z=0.

      gradtw_x=0.
      gradtw_y=0.
      gradtw_z=0.

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
          ! turbulent kinectic energy gradient
          ! x_face(i, j, k)
          tk_face = (x_tk_left(i,j,k)+ x_tk_right(i,j,k))
          gradtk_x(i, j, k) = gradtk_x(i, j, k) - &
                             (tk_face * xnx(i, j, k) * xA(i, j, k))
          gradtk_y(i, j, k) = gradtk_y(i, j, k) - &
                             (tk_face * xny(i, j, k) * xA(i, j, k))
          gradtk_z(i, j, k) = gradtk_z(i, j, k) - &
                             (tk_face * xnz(i, j, k) * xA(i, j, k))

          ! x_face(i+1, j, k)
          tk_face = (x_tk_left(i+1,j,k)+ x_tk_right(i+1,j,k))
          gradtk_x(i, j, k) = gradtk_x(i, j, k) + &
                             (tk_face * xnx(i+1, j, k) * xA(i+1, j, k))
          gradtk_y(i, j, k) = gradtk_y(i, j, k) + &
                             (tk_face * xny(i+1, j, k) * xA(i+1, j, k))
          gradtk_z(i, j, k) = gradtk_z(i, j, k) + &
                             (tk_face * xnz(i+1, j, k) * xA(i+1, j, k))

          ! y_face(i, j, k)
          tk_face = (y_tk_left(i,j,k)+ y_tk_right(i,j,k))
          gradtk_x(i, j, k) = gradtk_x(i, j, k) - &
                             (tk_face * ynx(i, j, k) * yA(i, j, k))
          gradtk_y(i, j, k) = gradtk_y(i, j, k) - &
                             (tk_face * yny(i, j, k) * yA(i, j, k))
          gradtk_z(i, j, k) = gradtk_z(i, j, k) - &
                             (tk_face * ynz(i, j, k) * yA(i, j, k))

          ! y_face(i, j+1, k)
          tk_face = (y_tk_left(i,j+1,k)+ y_tk_right(i,j+1,k))
          gradtk_x(i, j, k) = gradtk_x(i, j, k) + &
                             (tk_face * ynx(i, j+1, k) * yA(i, j+1, k))
          gradtk_y(i, j, k) = gradtk_y(i, j, k) + &
                             (tk_face * yny(i, j+1, k) * yA(i, j+1, k))
          gradtk_z(i, j, k) = gradtk_z(i, j, k) + &
                             (tk_face * ynz(i, j+1, k) * yA(i, j+1, k))

          ! z_face(i, j, k)
          tk_face = (z_tk_left(i,j,k)+ z_tk_right(i,j,k))
          gradtk_x(i, j, k) = gradtk_x(i, j, k) - &
                             (tk_face * znx(i, j, k) * zA(i, j, k))
          gradtk_y(i, j, k) = gradtk_y(i, j, k) - &
                             (tk_face * zny(i, j, k) * zA(i, j, k))
          gradtk_z(i, j, k) = gradtk_z(i, j, k) - &
                             (tk_face * znz(i, j, k) * zA(i, j, k))

          ! z_face(i, j, k+1)
          tk_face = (z_tk_left(i,j,k+1)+ z_tk_right(i,j,k+1))
          gradtk_x(i, j, k) = gradtk_x(i, j, k) + &
                             (tk_face * znx(i, j, k+1) * zA(i, j, k+1))
          gradtk_y(i, j, k) = gradtk_y(i, j, k) + &
                             (tk_face * zny(i, j, k+1) * zA(i, j, k+1))
          gradtk_z(i, j, k) = gradtk_z(i, j, k) + &
                             (tk_face * znz(i, j, k+1) * zA(i, j, k+1))

          ! Factor of volume
          gradtk_x(i, j, k) = gradtk_x(i, j, k) * 0.5 / volume(i, j, k)
          gradtk_y(i, j, k) = gradtk_y(i, j, k) * 0.5 / volume(i, j, k)
          gradtk_z(i, j, k) = gradtk_z(i, j, k) * 0.5 / volume(i, j, k)

          !omega gradient

          ! x_face(i, j, k)
          tw_face = (x_tw_left(i,j,k)+ x_tw_right(i,j,k))
          gradtw_x(i, j, k) = gradtw_x(i, j, k) - &
                             (tw_face * xnx(i, j, k) * xA(i, j, k))
          gradtw_y(i, j, k) = gradtw_y(i, j, k) - &
                             (tw_face * xny(i, j, k) * xA(i, j, k))
          gradtw_z(i, j, k) = gradtw_z(i, j, k) - &
                             (tw_face * xnz(i, j, k) * xA(i, j, k))

          ! x_face(i+1, j, k)
          tw_face = (x_tw_left(i+1,j,k)+ x_tw_right(i+1,j,k))
          gradtw_x(i, j, k) = gradtw_x(i, j, k) + &
                             (tw_face * xnx(i+1, j, k) * xA(i+1, j, k))
          gradtw_y(i, j, k) = gradtw_y(i, j, k) + &
                             (tw_face * xny(i+1, j, k) * xA(i+1, j, k))
          gradtw_z(i, j, k) = gradtw_z(i, j, k) + &
                             (tw_face * xnz(i+1, j, k) * xA(i+1, j, k))

          ! y_face(i, j, k)
          tw_face = (y_tw_left(i,j,k)+ y_tw_right(i,j,k))
          gradtw_x(i, j, k) = gradtw_x(i, j, k) - &
                             (tw_face * ynx(i, j, k) * yA(i, j, k))
          gradtw_y(i, j, k) = gradtw_y(i, j, k) - &
                             (tw_face * yny(i, j, k) * yA(i, j, k))
          gradtw_z(i, j, k) = gradtw_z(i, j, k) - &
                             (tw_face * ynz(i, j, k) * yA(i, j, k))

          ! y_face(i, j+1, k)
          tw_face = (y_tw_left(i,j+1,k)+ y_tw_right(i,j+1,k))
          gradtw_x(i, j, k) = gradtw_x(i, j, k) + &
                             (tw_face * ynx(i, j+1, k) * yA(i, j+1, k))
          gradtw_y(i, j, k) = gradtw_y(i, j, k) + &
                             (tw_face * yny(i, j+1, k) * yA(i, j+1, k))
          gradtw_z(i, j, k) = gradtw_z(i, j, k) + &
                             (tw_face * ynz(i, j+1, k) * yA(i, j+1, k))

          ! z_face(i, j, k)
          tw_face = (z_tw_left(i,j,k)+ z_tw_right(i,j,k))
          gradtw_x(i, j, k) = gradtw_x(i, j, k) - &
                             (tw_face * znx(i, j, k) * zA(i, j, k))
          gradtw_y(i, j, k) = gradtw_y(i, j, k) - &
                             (tw_face * zny(i, j, k) * zA(i, j, k))
          gradtw_z(i, j, k) = gradtw_z(i, j, k) - &
                             (tw_face * znz(i, j, k) * zA(i, j, k))

          ! z_face(i, j, k+1)
          tw_face = (z_tw_left(i,j,k+1)+ z_tw_right(i,j,k+1))
          gradtw_x(i, j, k) = gradtw_x(i, j, k) + &
                             (tw_face * znx(i, j, k+1) * zA(i, j, k+1))
          gradtw_y(i, j, k) = gradtw_y(i, j, k) + &
                             (tw_face * zny(i, j, k+1) * zA(i, j, k+1))
          gradtw_z(i, j, k) = gradtw_z(i, j, k) + &
                             (tw_face * znz(i, j, k+1) * zA(i, j, k+1))

          ! Factor of volume
          gradtw_x(i, j, k) = gradtw_x(i, j, k) * 0.5 / volume(i, j, k)
          gradtw_y(i, j, k) = gradtw_y(i, j, k) * 0.5 / volume(i, j, k)
          gradtw_z(i, j, k) = gradtw_z(i, j, k) * 0.5 / volume(i, j, k)
        end do
      end do
    end do

    end subroutine calculate_sst_grad


end module sst_grad
