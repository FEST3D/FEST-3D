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
  use global_vars,          only: imin_id
  use global_vars,          only: jmin_id
  use global_vars,          only: kmin_id
  use global_vars,          only: imax_id
  use global_vars,          only: jmax_id
  use global_vars,          only: kmax_id

  use global_vars,          only: qp
!  use global_vars,          only: n_var
!  use global_vars,          only: ilimiter_switch
!  use global_vars,          only: jlimiter_switch
!  use global_vars,          only: klimiter_switch
!  use global_vars,          only: itlimiter_switch
!  use global_vars,          only: jtlimiter_switch
!  use global_vars,          only: ktlimiter_switch

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

    subroutine reconstruct_boundary_state(control, scheme, dims)
      !< Call reconstruction based on the flag and boundary condition

      implicit none
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims

      DebugCall('recons_boundary_state')

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      n_var = control%n_var
      if (scheme%interpolant == 'ppm' .or. scheme%interpolant=='weno' .or. scheme%interpolant=='weno_NM') ppm_flag=1
      if (imin_id==-7 .or. jmin_id==-7 .or. kmin_id==-7) ppm_flag=1
      if (imax_id==-7 .or. jmax_id==-7 .or. kmax_id==-7) ppm_flag=1
      if(scheme%interpolant /='none')then
        if(imin_id<0 .and. imin_id/=-10)then
          DebugCall('recons_bndry_state: imin')
          call reconstruct_imin(scheme)
        end if
        if(imax_id<0 .and. imax_id/=-10)then
          DebugCall('recons_bndry_state: imax')
          call reconstruct_imax(scheme)
        end if
        if(jmin_id<0 .and. jmin_id/=-10)then
          DebugCall('recons_bndry_state: jmin')
          call reconstruct_jmin(scheme)
        end if
        if(jmax_id<0 .and. jmax_id/=-10)then
          DebugCall('recons_bndry_state: jmax')
          call reconstruct_jmax(scheme)
        end if
        if(kmin_id<0 .and. kmin_id/=-10)then
          DebugCall('recons_bndry_state: kmin')
          call reconstruct_kmin(scheme)
        end if
        if(kmax_id<0 .and. kmax_id/=-10)then
        DebugCall('recons_bndry_state: kmax')
          call reconstruct_kmax(scheme)
        end if
      end if

    end subroutine reconstruct_boundary_state


    subroutine reconstruct_imin(scheme)
      !< Reconstruct state at the IMIN boundary face with MUSCL scheme

      implicit none
      type(schemetype), intent(in) :: scheme
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

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
      if(imin_id==-8 .or. imin_id==-9)then
         x_qp_left(1,1:jmx-1,1:kmx-1,1:n_var) = qp(0,1:jmx-1,1:kmx-1,1:n_var) 
        x_qp_right(1,1:jmx-1,1:kmx-1,1:n_var) = qp(0,1:jmx-1,1:kmx-1,1:n_var) 
      else
        ! right face of first ghost cell
        x_qp_left(1,1:jmx-1,1:kmx-1,1:n_var) = 0.5*(qp(0,1:jmx-1,1:kmx-1,1:n_var)&
                                                   +qp(1,1:jmx-1,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_imin


    subroutine reconstruct_imax(scheme)
      !< Reconstruct state at the IMAX boundary face with MUSCL scheme

      implicit none
      type(schemetype), intent(in) :: scheme
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

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
      if(imax_id==-8 .or. imax_id==-9)then
         x_qp_left(imx,1:jmx-1,1:kmx-1,1:n_var) = qp(imx,1:jmx-1,1:kmx-1,1:n_var) 
        x_qp_right(imx,1:jmx-1,1:kmx-1,1:n_var) = qp(imx,1:jmx-1,1:kmx-1,1:n_var) 
      else
        x_qp_right(imx,1:jmx-1,1:kmx-1,1:n_var)  = 0.5*(qp(imx-1,1:jmx-1,1:kmx-1,1:n_var)&
                                                       +qp(imx  ,1:jmx-1,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_imax


    subroutine reconstruct_jmin(scheme)
      !< Reconstruct state at the JMIN boundary face with MUSCL scheme

      implicit none
      type(schemetype), intent(in) :: scheme
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

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
      if(jmin_id==-8 .or. jmin_id==-9)then
         y_qp_left(1:imx-1,1,1:kmx-1,1:n_var) = qp(1:imx-1,0,1:kmx-1,1:n_var) 
        y_qp_right(1:imx-1,1,1:kmx-1,1:n_var) = qp(1:imx-1,0,1:kmx-1,1:n_var) 
      else
         y_qp_left(1:imx-1,1,1:kmx-1,1:n_var) = 0.5*(qp(1:imx-1,0,1:kmx-1,1:n_var)&
                                                    +qp(1:imx-1,1,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_jmin


    subroutine reconstruct_jmax(scheme)
      !< Reconstruct state at the JMAX boundary face with MUSCL scheme

      implicit none
      type(schemetype), intent(in) :: scheme
      integer :: i, j, k, l
      real :: psi1, psi2, fd, bd, r
      real :: kappa, phi

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
      if(jmax_id==-8 .or. jmax_id==-9)then
         y_qp_left(1:imx-1,jmx,1:kmx-1,1:n_var) = qp(1:imx-1,jmx,1:kmx-1,1:n_var) 
        y_qp_right(1:imx-1,jmx,1:kmx-1,1:n_var) = qp(1:imx-1,jmx,1:kmx-1,1:n_var) 
      else
        y_qp_right(1:imx-1,jmx,1:kmx-1,1:n_var)  = 0.5*(qp(1:imx-1,jmx-1,1:kmx-1,1:n_var)& 
                                                       +qp(1:imx-1,jmx  ,1:kmx-1,1:n_var))
      end if

    end subroutine reconstruct_jmax


    subroutine reconstruct_kmin(scheme)
      !< Reconstruct state at the KMIN boundary face with MUSCL scheme

      implicit none
      type(schemetype), intent(in) :: scheme
      real :: psi1, psi2, fd, bd, r
      integer :: i, j, k, l
      real :: kappa, phi
      
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
      if(kmin_id==-8 .or. kmin_id==-9)then
         z_qp_left(1:imx-1,1:jmx-1,1,1:n_var) = qp(1:imx-1,1:jmx-1,0,1:n_var) 
        z_qp_right(1:imx-1,1:jmx-1,1,1:n_var) = qp(1:imx-1,1:jmx-1,0,1:n_var) 
      else
         z_qp_left(1:imx-1,1:jmx-1,1,1:n_var) = 0.5*(qp(1:imx-1,1:jmx-1,0,1:n_var)&
                                                    +qp(1:imx-1,1:jmx-1,1,1:n_var))
      end if

    end subroutine reconstruct_kmin


    subroutine reconstruct_kmax(scheme)
      !< Reconstruct state at the KMAX boundary face with MUSCL scheme

      implicit none
      type(schemetype), intent(in) :: scheme
      real :: psi1, psi2, fd, bd, r
      integer :: i, j, k, l
      real :: kappa, phi
    
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
      if(kmax_id==-8 .or. kmax_id==-9)then
         z_qp_left(1:imx-1,1:jmx-1,kmx,1:n_var) = qp(1:imx-1,1:jmx-1,kmx,1:n_var) 
        z_qp_right(1:imx-1,1:jmx-1,kmx,1:n_var) = qp(1:imx-1,1:jmx-1,kmx,1:n_var) 
      else
        z_qp_right(1:imx-1,1:jmx-1,kmx,1:n_var) = 0.5*(qp(1:imx-1,1:jmx-1,kmx-1,1:n_var)&
                                                      +qp(1:imx-1,1:jmx-1,kmx  ,1:n_var))
      end if

    end subroutine reconstruct_kmax

end module boundary_state_reconstruction
