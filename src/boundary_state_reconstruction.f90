module boundary_state_reconstruction
  use utils,                only: dmsg
  use grid,                 only: imx, jmx, kmx
  use state,                only: qp
  use state,                only: n_var
  use face_interpolant,     only: x_qp_left, x_qp_right
  use face_interpolant,     only: y_qp_left, y_qp_right
  use face_interpolant,     only: z_qp_left, z_qp_right
  use boundary_conditions,  only: bc_imn, bc_imx
  use boundary_conditions,  only: bc_jmn, bc_jmx
  use boundary_conditions,  only: bc_kmn, bc_kmx
  use boundary_conditions,  only: BC_INTERFACE

  implicit none
  private

  public :: reconstruct_boundary_state

  contains

    subroutine reconstruct_boundary_state(interpolant)

      implicit none
      character(len=*), intent(in)  :: interpolant
      call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state')
      if(interpolant /='none')then
        if(bc_imn(1,1) /= BC_INTERFACE)then
          call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state', 'imin')
          call reconstruct_imin()
        end if
        if(bc_imx(1,1) /= BC_INTERFACE)then
          call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state', 'imax')
          call reconstruct_imax()
        end if
        if(bc_jmn(1,1) /= BC_INTERFACE)then
          call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state', 'jmin')
          call reconstruct_jmin()
        end if
        if(bc_jmx(1,1) /= BC_INTERFACE)then
          call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state', 'jmax')
          call reconstruct_jmax()
        end if
        if(bc_kmn(1,1) /= BC_INTERFACE)then
          call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state', 'kmin')
          call reconstruct_kmin()
        end if
        if(bc_kmx(1,1) /= BC_INTERFACE)then
        call dmsg(1,'boundary_state_reconstruction', 'reconstruct_boundary_state', 'kmax')
          call reconstruct_kmax()
        end if
      end if

    end subroutine reconstruct_boundary_state


    subroutine reconstruct_imin()

      implicit none
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

      phi = 1.0
      kappa = -1.0

      do l = 1, n_var
       do k = 1, kmx - 1
        do j = 1, jmx - 1
         do i = 1, 1 
          fd = qp(i+1, j, k, l) - qp(i, j, k, l)
          bd = qp(i, j, k, l) - qp(i-1, j, k, l)
          r = fd / bd
          psi1 = max(0., min(2*r, (2 + r)/3., 2.))
          r = bd / fd
          psi2 = max(0., min(2*r, (2 + r)/3., 2.))

          x_qp_left(i, j, k, l) = qp(i-1, j, k, l)!0.5 * (qp(i-1, j, k, l) + qp(i, j, k, l))
          x_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
              (((1.+kappa) * psi1 * bd) + ((1.-kappa) * psi2 * fd))
         end do
        end do
       end do
      end do

    end subroutine reconstruct_imin


    subroutine reconstruct_imax()

      implicit none
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

      phi = 1.0
      kappa = -1.0

      do l = 1, n_var
       do k = 1, kmx - 1
        do j = 1, jmx - 1
         do i = imx-1, imx-1 
          fd = qp(i+1, j, k, l) - qp(i, j, k, l)
          bd = qp(i, j, k, l) - qp(i-1, j, k, l)
          r = fd / bd
          psi1 = max(0., min(2*r, (2 + r)/3., 2.))
          r = bd / fd
          psi2 = max(0., min(2*r, (2 + r)/3., 2.))

          x_qp_left(i+1, j, k, l) = qp(i, j, k, l) + 0.25*phi* &
              (((1.-kappa) * psi1 * bd) + ((1.+kappa) * psi2 * fd))
          x_qp_right(i+1, j, k, l) = qp(i+1, j, k, l) !0.5 * (qp(i+1, j, k, l) + qp(i, j, k, l))
         end do
        end do
       end do
      end do

    end subroutine reconstruct_imax


    subroutine reconstruct_jmin()

      implicit none
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

      phi = 1.0
      kappa = -1.0

      do l = 1, n_var
       do k = 1, kmx - 1
        do j = 1, 1
         do i = 1, imx - 1
          fd = qp(i, j+1, k, l) - qp(i, j, k, l)
          bd = qp(i, j, k, l) - qp(i, j-1, k, l)
          r = fd / bd
          psi1 = max(0., min(2*r, (2 + r)/3., 2.))
          r = bd / fd
          psi2 = max(0., min(2*r, (2 + r)/3., 2.))

          y_qp_left(i, j, k, l) = 0.5 * (qp(i, j, k, l) + qp(i, j-1, k, l))
          y_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
              (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
         end do
        end do
       end do
      end do

    end subroutine reconstruct_jmin


    subroutine reconstruct_jmax()

      implicit none
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

      phi = 1.0
      kappa = -1.0

      do l = 1, n_var
       do k = 1, kmx - 1
        do j = jmx-1, jmx-1
         do i = 1, imx - 1
          fd = qp(i, j+1, k, l) - qp(i, j, k, l)
          bd = qp(i, j, k, l) - qp(i, j-1, k, l)
          r = fd / bd
          psi1 = max(0., min(2*r, (2 + r)/3., 2.))
          r = bd / fd
          psi2 = max(0., min(2*r, (2 + r)/3., 2.))

          y_qp_left(i, j+1, k, l) = qp(i, j, k, l) + 0.25*phi* &
              (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
          y_qp_right(i, j+1, k, l) = 0.5 * (qp(i, j, k, l) + qp(i, j+1, k, l))
         end do
        end do
       end do
      end do

    end subroutine reconstruct_jmax


    subroutine reconstruct_kmin()

      implicit none
      real :: psi1, psi2, fd, bd, r
      integer :: i, j, k, l
      real :: kappa, phi
      
      phi = 1.0
      kappa = -1.0

      do k = 1, 1
       do l = 1, n_var
        do j = 1, jmx - 1
         do i = 1, imx - 1
          fd = qp(i, j, k+1, l) - qp(i, j, k, l)
          bd = qp(i, j, k, l) - qp(i, j, k-1, l)
          r = fd / bd
          psi1 = max(0., min(2*r, (2 + r)/3., 2.))
          r = bd / fd
          psi2 = max(0., min(2*r, (2 + r)/3., 2.))

          z_qp_left(i, j, k, l) = 0.5 * (qp(i, j, k, l) + qp(i, j, k-1, l))
          z_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
              (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
         end do
        end do
       end do
      end do

    end subroutine reconstruct_kmin


    subroutine reconstruct_kmax()

      implicit none
      real :: psi1, psi2, fd, bd, r
      integer :: i, j, k, l
      real :: kappa, phi
    
      phi = 1.0
      kappa = -1.0

      do k = kmx-1, kmx-1
       do l = 1, n_var
        do j = 1, jmx - 1
         do i = 1, imx - 1
          fd = qp(i, j, k+1, l) - qp(i, j, k, l)
          bd = qp(i, j, k, l) - qp(i, j, k-1, l)
          r = fd / bd
          psi1 = max(0., min(2*r, (2 + r)/3., 2.))
          r = bd / fd
          psi2 = max(0., min(2*r, (2 + r)/3., 2.))

          z_qp_left(i, j, k+1, l) = qp(i, j, k, l) + 0.25*phi* &
              (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
          z_qp_right(i, j, k+1, l) = 0.5 * (qp(i, j, k, l) + qp(i, j, k+1, l))
         end do
        end do
       end do
      end do

    end subroutine reconstruct_kmax

end module boundary_state_reconstruction
