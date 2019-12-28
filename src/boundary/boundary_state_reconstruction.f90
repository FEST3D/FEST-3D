  !< Reconstruct the boundary face in case of 4th and 5th order methods
module boundary_state_reconstruction
  !< Reconstruct the boundary face in case of 4th and 5th order higher order
  !< face state reconstruction method. Since the limited information
  !< is available at the boundaries, the boundary face is limiter to 
  !< 3rd order accurate and is reconstructed using MUSCL Scheme even when
  !< rest of the domain is using WENO or PPM
#include "../debug.h"
#include "../error.h"

   use vartypes
  use face_interpolant,     only: x_qp_left, x_qp_right
  use face_interpolant,     only: y_qp_left, y_qp_right
  use face_interpolant,     only: z_qp_left, z_qp_right

  implicit none
  private

  integer :: ppm_flag=0
  !< Flag to check if reconstruction is required
  integer :: switch_L=1
  !< Limiter switch
  integer :: imx, jmx, kmx, n_var
  public :: reconstruct_boundary_state

  contains

    subroutine reconstruct_boundary_state(qp, control, scheme, bc, dims)
      !< Call reconstruction based on the flag and boundary condition

      implicit none
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc

      DebugCall('recons_boundary_state')

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      n_var = control%n_var
      if (scheme%interpolant == 'ppm' .or. scheme%interpolant=='weno' .or. scheme%interpolant=='weno_NM') ppm_flag=1
      if (bc%imin_id==-7 .or. bc%jmin_id==-7 .or. bc%kmin_id==-7) ppm_flag=1
      if (bc%imax_id==-7 .or. bc%jmax_id==-7 .or. bc%kmax_id==-7) ppm_flag=1
      if(scheme%interpolant /='none')then
        if(bc%imin_id<0 .and. bc%imin_id/=-10)then
          DebugCall('recons_bndry_state: imin')
          call reconstruct_imin(qp, scheme, bc)
        end if
        if(bc%imax_id<0 .and. bc%imax_id/=-10)then
          DebugCall('recons_bndry_state: imax')
          call reconstruct_imax(qp, scheme, bc)
        end if
        if(bc%jmin_id<0 .and. bc%jmin_id/=-10)then
          DebugCall('recons_bndry_state: jmin')
          call reconstruct_jmin(qp, scheme, bc)
        end if
        if(bc%jmax_id<0 .and. bc%jmax_id/=-10)then
          DebugCall('recons_bndry_state: jmax')
          call reconstruct_jmax(qp, scheme, bc)
        end if
        if(bc%kmin_id<0 .and. bc%kmin_id/=-10)then
          DebugCall('recons_bndry_state: kmin')
          call reconstruct_kmin(qp, scheme, bc)
        end if
        if(bc%kmax_id<0 .and. bc%kmax_id/=-10)then
        DebugCall('recons_bndry_state: kmax')
          call reconstruct_kmax(qp, scheme, bc)
        end if
      end if

    end subroutine reconstruct_boundary_state


    subroutine reconstruct_imin(qp, scheme, bc)
      !< Reconstruct state at the IMIN boundary face with MUSCL scheme

      implicit none
      real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var), intent(in) :: qp
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc
      integer :: i, j, k, l
      real(wp) :: psi1, psi2, fd, bd, r
      real(wp) :: kappa, phi

      phi = 1.0
      kappa = 1./3.
      switch_L=scheme%ilimiter_switch

      if (ppm_flag==1) then
        do l = 1, n_var
          if(l>=6) switch_L=scheme%itlimiter_switch
         do k = 1, kmx - 1
          do j = 1, jmx - 1
           do i = 1, 1 

            ! reconstruct first cell faces for ppm scheme
              fd = qp(i+1, j, k, l) - qp(i  , j, k, l)
              bd = qp(i  , j, k, l) - qp(i-1, j, k, l)

              r = fd / bd
              psi1 = max(0., min(2*r, (2 + r)/3., 2.))
              psi1 = (1 - (1 - psi1)*switch_L )
              r = bd / fd
              psi2 = max(0., min(2*r, (2 + r)/3., 2.))
              psi2 = (1 - (1 - psi2)*switch_L )

              ! right state of firsrt interior cell
              x_qp_left(i+1, j, k, l) = qp(i, j, k, l) + 0.25*phi* &
                  (((1.-kappa) * psi1 * bd) + ((1.+kappa) * psi2 * fd))

              ! left face of first interior cell
              x_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                  (((1.+kappa) * psi1 * bd) + ((1.-kappa) * psi2 * fd))
              
           end do
          end do
         end do
        end do
      end if
      if(bc%imin_id==-8 .or. bc%imin_id==-9)then
         x_qp_left(1,1:jmx-1,1:kmx-1,1:n_var) = qp(0,1:jmx-1,1:kmx-1,1:n_var) 
        x_qp_right(1,1:jmx-1,1:kmx-1,1:n_var) = qp(0,1:jmx-1,1:kmx-1,1:n_var) 
      else
        ! right face of first ghost cell
        x_qp_left(1,1:jmx-1,1:kmx-1,1:n_var) = 0.5*(qp(0,1:jmx-1,1:kmx-1,1:n_var)&
                                                   +qp(1,1:jmx-1,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_imin


    subroutine reconstruct_imax(qp, scheme, bc)
      !< Reconstruct state at the IMAX boundary face with MUSCL scheme

      implicit none
      real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var), intent(in) :: qp
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc
      integer :: i, j, k, l
      real(wp) :: psi1, psi2, fd, bd, r
      real(wp) :: kappa, phi

      phi = 1.0
      kappa = 1./3.
      switch_L=scheme%ilimiter_switch

      if (ppm_flag==1) then
        do l = 1, n_var
          if(l>=6) switch_L=scheme%itlimiter_switch
         do k = 1, kmx - 1
          do j = 1, jmx - 1
           do i = imx-1, imx-1 

             fd = qp(i+1, j, k, l) - qp(i  , j, k, l)
             bd = qp(i  , j, k, l) - qp(i-1, j, k, l)

             r = fd / bd
             psi1 = max(0., min(2*r, (2 + r)/3., 2.))
             psi1 = (1 - (1 - psi1)*switch_L )
             r = bd / fd
             psi2 = max(0., min(2*r, (2 + r)/3., 2.))
             psi2 = (1 - (1 - psi2)*switch_L )

             ! right face of last interior cell
             x_qp_left(i+1, j, k, l) = qp(i, j, k, l) + 0.25*phi* &
                 (((1.-kappa) * psi1 * bd) + ((1.+kappa) * psi2 * fd))

             ! left face of last interior cell
             x_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                 (((1.+kappa) * psi1 * bd) + ((1.-kappa) * psi2 * fd))
           end do
          end do
         end do
        end do
      end if
      if(bc%imax_id==-8 .or. bc%imax_id==-9)then
         x_qp_left(imx,1:jmx-1,1:kmx-1,1:n_var) = qp(imx,1:jmx-1,1:kmx-1,1:n_var) 
        x_qp_right(imx,1:jmx-1,1:kmx-1,1:n_var) = qp(imx,1:jmx-1,1:kmx-1,1:n_var) 
      else
        x_qp_right(imx,1:jmx-1,1:kmx-1,1:n_var)  = 0.5*(qp(imx-1,1:jmx-1,1:kmx-1,1:n_var)&
                                                       +qp(imx  ,1:jmx-1,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_imax


    subroutine reconstruct_jmin(qp, scheme, bc)
      !< Reconstruct state at the JMIN boundary face with MUSCL scheme

      implicit none
      real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var), intent(in) :: qp
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc
      integer :: i, j, k, l
      real(wp) :: psi1, psi2, fd, bd, r
      real(wp) :: kappa, phi

      phi = 1.0
      kappa = 1./3.
      switch_L=scheme%jlimiter_switch

      if (ppm_flag==1) then
        do l = 1, n_var
          if(l>=6) switch_L=scheme%jtlimiter_switch
         do k = 1, kmx - 1
          do j = 1, 1
           do i = 1, imx - 1

              fd = qp(i, j+1, k, l) - qp(i, j, k, l)
              bd = qp(i, j, k, l) - qp(i, j-1, k, l)

              r = fd / bd
              psi1 = max(0., min(2*r, (2 + r)/3., 2.))
              psi1 = (1 - (1 - psi1)*switch_L )
              r = bd / fd
              psi2 = max(0., min(2*r, (2 + r)/3., 2.))
              psi2 = (1 - (1 - psi2)*switch_L )

              ! right face of first j cell
              y_qp_left(i, j+1, k, l) = qp(i, j, k, l) + 0.25*phi* &
                  (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))

              ! left face of first j cell
              y_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                  (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
           end do
          end do
         end do
        end do
      end if
      if(bc%jmin_id==-8 .or. bc%jmin_id==-9)then
         y_qp_left(1:imx-1,1,1:kmx-1,1:n_var) = qp(1:imx-1,0,1:kmx-1,1:n_var) 
        y_qp_right(1:imx-1,1,1:kmx-1,1:n_var) = qp(1:imx-1,0,1:kmx-1,1:n_var) 
      else
         y_qp_left(1:imx-1,1,1:kmx-1,1:n_var) = 0.5*(qp(1:imx-1,0,1:kmx-1,1:n_var)&
                                                    +qp(1:imx-1,1,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_jmin


    subroutine reconstruct_jmax(qp, scheme, bc)
      !< Reconstruct state at the JMAX boundary face with MUSCL scheme

      implicit none
      real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var), intent(in) :: qp
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc
      integer :: i, j, k, l
      real(wp) :: psi1, psi2, fd, bd, r
      real(wp) :: kappa, phi

      phi = 1.0
      kappa = 1./3.
      switch_L=scheme%jlimiter_switch

      if (ppm_flag==1) then
        do l = 1, n_var
          if(l>=6) switch_L=scheme%jtlimiter_switch
         do k = 1, kmx - 1
          do j = jmx-1, jmx-1
           do i = 1, imx - 1

              fd = qp(i, j+1, k, l) - qp(i, j, k, l)
              bd = qp(i, j, k, l) - qp(i, j-1, k, l)
              r = fd / bd
              psi1 = max(0., min(2*r, (2 + r)/3., 2.))
              psi1 = (1 - (1 - psi1)*switch_L )
              r = bd / fd
              psi2 = max(0., min(2*r, (2 + r)/3., 2.))
              psi2 = (1 - (1 - psi2)*switch_L )

              ! right face of last j cell
              y_qp_left(i, j+1, k, l) = qp(i, j, k, l) + 0.25*phi* &
                  (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
            
              ! left face of last j cell
              y_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                  (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
           end do
          end do
         end do
        end do
      end if
      if(bc%jmax_id==-8 .or. bc%jmax_id==-9)then
         y_qp_left(1:imx-1,jmx,1:kmx-1,1:n_var) = qp(1:imx-1,jmx,1:kmx-1,1:n_var) 
        y_qp_right(1:imx-1,jmx,1:kmx-1,1:n_var) = qp(1:imx-1,jmx,1:kmx-1,1:n_var) 
      else
        y_qp_right(1:imx-1,jmx,1:kmx-1,1:n_var)  = 0.5*(qp(1:imx-1,jmx-1,1:kmx-1,1:n_var)& 
                                                       +qp(1:imx-1,jmx  ,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_jmax


    subroutine reconstruct_kmin(qp, scheme, bc)
      !< Reconstruct state at the KMIN boundary face with MUSCL scheme

      implicit none
      real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var), intent(in) :: qp
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc
      real(wp) :: psi1, psi2, fd, bd, r
      integer :: i, j, k, l
      real(wp) :: kappa, phi
      
      phi = 1.0
      kappa = 1./3.
      switch_L=scheme%klimiter_switch

      if (ppm_flag==1) then
        do k = 1, 1
         do l = 1, n_var
           if(l>=6) switch_L=scheme%ktlimiter_switch
           if(l<6) switch_L=scheme%klimiter_switch
          do j = 1, jmx - 1
           do i = 1, imx - 1

              fd = qp(i, j, k+1, l) - qp(i, j, k, l)
              bd = qp(i, j, k, l) - qp(i, j, k-1, l)

              r = fd / bd
              psi1 = max(0., min(2*r, (2 + r)/3., 2.))
              psi1 = (1 - (1 - psi1)*switch_L )
              r = bd / fd
              psi2 = max(0., min(2*r, (2 + r)/3., 2.))
              psi2 = (1 - (1 - psi2)*switch_L )
              
              ! right face of first k cell
              z_qp_left(i, j, k+1, l) = qp(i, j, k, l) + 0.25*phi* &
                  (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))

              ! left face of first k cell
              z_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                  (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
           end do
          end do
         end do
        end do
      end if
      if(bc%kmin_id==-8 .or. bc%kmin_id==-9)then
         z_qp_left(1:imx-1,1:jmx-1,1,1:n_var) = qp(1:imx-1,1:jmx-1,0,1:n_var) 
        z_qp_right(1:imx-1,1:jmx-1,1,1:n_var) = qp(1:imx-1,1:jmx-1,0,1:n_var) 
      else
         z_qp_left(1:imx-1,1:jmx-1,1,1:n_var) = 0.5*(qp(1:imx-1,1:jmx-1,0,1:n_var)&
                                                    +qp(1:imx-1,1:jmx-1,1,1:n_var))
      end if

    end subroutine reconstruct_kmin


    subroutine reconstruct_kmax(qp, scheme, bc)
      !< Reconstruct state at the KMAX boundary face with MUSCL scheme

      implicit none
      real(wp), dimension(-2:imx+2, -2:jmx+2, -2:kmx+2, 1:n_var), intent(in) :: qp
      type(schemetype), intent(in) :: scheme
      type(boundarytype), intent(in) :: bc
      real(wp) :: psi1, psi2, fd, bd, r
      integer :: i, j, k, l
      real(wp) :: kappa, phi
    
      phi = 1.0
      kappa = 1./3.
      switch_L=scheme%klimiter_switch

      do k = kmx-1, kmx-1
       do l = 1, n_var
         if(l>=6) switch_L=scheme%ktlimiter_switch
         if(l<6) switch_L=scheme%klimiter_switch
        do j = 1, jmx - 1
         do i = 1, imx - 1
          ! left face of kmx ghost cell
          z_qp_right(i, j, k+1, l) = 0.5 * (qp(i, j, k, l) + qp(i, j, k+1, l))

          if (ppm_flag==1) then

            fd = qp(i, j, k+1, l) - qp(i, j, k, l)
            bd = qp(i, j, k, l) - qp(i, j, k-1, l)

            r = fd / bd
            psi1 = max(0., min(2*r, (2 + r)/3., 2.))
            psi1 = (1 - (1 - psi1)*switch_L )
            r = bd / fd
            psi2 = max(0., min(2*r, (2 + r)/3., 2.))
            psi2 = (1 - (1 - psi2)*switch_L )

            ! right face of last k interior cell
            z_qp_left(i, j, k+1, l) = qp(i, j, k, l) + 0.25*phi* &
                (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))

            ! left face of last k cell
            z_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
         end if
         end do
        end do
       end do
      end do
      if(bc%kmax_id==-8 .or. bc%kmax_id==-9)then
         z_qp_left(1:imx-1,1:jmx-1,kmx,1:n_var) = qp(1:imx-1,1:jmx-1,kmx,1:n_var) 
        z_qp_right(1:imx-1,1:jmx-1,kmx,1:n_var) = qp(1:imx-1,1:jmx-1,kmx,1:n_var) 
      else
        z_qp_right(1:imx-1,1:jmx-1,kmx,1:n_var) = 0.5*(qp(1:imx-1,1:jmx-1,kmx-1,1:n_var)&
                                                      +qp(1:imx-1,1:jmx-1,kmx  ,1:n_var))
      end if

    end subroutine reconstruct_kmax

end module boundary_state_reconstruction
