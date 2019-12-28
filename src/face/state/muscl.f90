    !< Higher order face state reconstruction method: MUSCL. 
module muscl
    !<
    !<Reference:Hirsch, C., Numerical computation of internal 
    !<and external flows: The fundamentals of computational fluid 
    !<dynamics, Elsevier, 2007
    !<
    !< MUSCL (Monotone Upwing Schemes for Scalar Conservation Laws is
    !< a scheme which replaces the piecewise constant approximation by
    !< reconstructing the states at the left and right side of each face.
    !< This is a one parameter upwind scheme which results in at most 3rd
    !< order accuracy.
    !
    ! The MUSCL scheme alone creates non-physical oscillations near 
    ! discontinuities like shocks. Hence, MUSCL is combined with
    ! some TVD (Total Variation Diminishing) to reduce such oscillations.
    ! TVD schemes also ensure that no new extrema of the state variables
    ! is created at the faces.
    !-----------------------------------------------------------------

#include "../../error.h"
#include "../../debug.h"

    use vartypes
    implicit none
    private

    real(wp) :: phi=1.0, kappa=1./3.
    integer :: switch_L=1
    !< Limiter switch 

    ! Public members
    public :: compute_muscl_states

    contains
        
        subroutine pressure_based_switching(qp, f_qp_left, f_qp_right, pdif,  flags,  flow, dims)
          !< Pressure based switching. 
          !< User x,y, or z for I,J,or K face respectively
          !----------------------------------------------

            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            integer, dimension(3), intent(in) :: flags
            !< flags for direction switch
            real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2), 1-flags(3):dims%kmx-1+2*flags(3), 1:dims%n_var), intent(inout) :: f_qp_left, f_qp_right
            !< primitive variable at cell faces
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(inout) :: pdif
            !< pressure difference 
            integer :: i, j, k, i_end, j_end, k_end
            integer :: i_f, j_f, k_f  ! Flags to determine face direction
            real(wp) :: pd2

            DebugCall('pressure_based_switching')


            i_f = flags(1)
            j_f = flags(2)
            k_f = flags(3)
            i_end = dims%imx - 1 +i_f
            j_end = dims%jmx - 1 +j_f
            k_end = dims%kmx - 1 +k_f

            ! i_end and j_end denote number of faces
            ! Total number of cells including ghost_cells is
            ! (i_end+1) * j_end for xi faces and i_end*(j_end+1) for
            ! eta faces. 

            ! Loop over cells (physical)
            do k = 1, dims%kmx - 1
             do j = 1, dims%jmx - 1
              do i = 1, dims%imx - 1
                pd2 = abs(qp(i + i_f*1, j + j_f*1, k + k_f*1, 5) - &  !pressure
                          qp(i - i_f*1, j - j_f*1, k - k_f*1, 5))
                pdif(i, j, k) = 1 - (pd2/(pd2 + flow%pressure_inf))
              end do
             end do
            end do

            ! Update at ghost cells
            pdif((1-i_f):(1-i_f)*(dims%imx-1), (1-j_f):(1-j_f)*(dims%jmx-1), &
                 (1-k_f):(1-k_f)*(dims%kmx-1)) = &
                pdif(1:dims%imx-1 - i_f*(dims%imx-2), 1:dims%jmx-1 - j_f*(dims%jmx-2), &
                     1:dims%kmx-1 - k_f*(dims%kmx-2))
                    
            pdif((1-i_f*(-dims%imx+1)):(i_f)+(dims%imx-1), (1-j_f*(-dims%jmx+1)):(j_f)+(dims%jmx-1), &
                 (1-k_f*(-dims%kmx+1)):(k_f)+(dims%kmx-1)) = &
                pdif(1+i_f*(dims%imx-2):dims%imx-1 , 1+j_f*(dims%jmx-2):dims%jmx-1, &
                     1+k_f*(dims%kmx-2):dims%kmx-1)
             

            ! Loop over faces
            do k = 1, dims%kmx - (1 - k_f)            
             do j = 1, dims%jmx - (1 - j_f)
              do i = 1, dims%imx - (1 - i_f)
                f_qp_left(i, j, k, :) = qp(i - i_f*1, j - j_f*1, k - k_f*1, :) + (&
                    pdif(i - i_f*1, j - j_f*1, k - k_f*1) * ( &
                    f_qp_left(i, j, k, :) - qp(i - i_f*1, j - j_f*1, k - k_f*1, :)))

                f_qp_right(i, j, k, :) = qp(i, j, k, :) - (&
                    pdif(i, j, k) * ( &
                    qp(i, j, k, :) - f_qp_right(i, j, k, :)))
              end do
             end do
            end do

        end subroutine pressure_based_switching


        subroutine compute_face_state(qp, f_qp_left, f_qp_right, flags, lam_switch, turb_switch, dims)
          !< Subroutine to calculate state at the face, generalized for
          !< all direction : I,J, and K.
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            integer, dimension(3), intent(in) :: flags
            !< Flags for direction switch
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
            !< Store primitive variable at cell center
            real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2), 1-flags(3):dims%kmx-1+2*flags(3), 1:dims%n_var), intent(inout) :: f_qp_left, f_qp_right
            !< primitive variable at cell faces
            integer, intent(in) :: lam_switch
            !< Limiter switch for laminar variables
            integer, intent(in) :: turb_switch
            !< Limiter switch for turbulent variables
            integer :: i, j, k, l
            !< integer used for DO loop
            integer :: ii, jj, kk
            !< Variable for ALFA family limiter
            real(wp) :: alpha
            !< Flags to determine face direction
            real(wp) :: psi1, psi2
            !< limiters
            real(wp) :: fd
            !< forward difference
            real(wp) :: bd
            !< backward difference
            real(wp) :: r
            !< ratio of differences

            DebugCall('compute_face_state')


            alpha = 2./3. !Koren limiter 
            phi = 1.0
            kappa = 1./3.
            switch_L=lam_switch

            ii = flags(1)
            jj = flags(2)
            kk = flags(3)

            do l = 1, dims%n_var
              if(l>=6)then
                switch_L=turb_switch
              end if
             do k = 1-kk, dims%kmx - 1 + kk
              do j = 1-jj, dims%jmx - 1 + jj
               do i = 1-ii, dims%imx - 1 + ii
                ! Cell based
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                fd = qp(i+ii, j+jj, k+kk, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i-ii, j-jj, k-kk, l)

                r = fd / max(bd,1e-10)
                psi1 = max(0., min(2*r, alpha*(r-1.0) + 1.0, 2.))  !alpha limiter
!                psi1 = max(0., min(2*r,1.), min(r,2.))    ! superbee
!                psi1 = ((r*r) + r)/((r*r) + 1.0)          ! Van-Albda 
!                psi1 = (abs(r) + r)/(abs(r) + 1.0)          ! Van-Leer

                r = bd / max(fd, 1e-10)
                psi2 = max(0., min(2*r, alpha*(r-1.0) + 1.0, 2.))
!                psi2 = max(0., min(2*r,1.), min(r,2.))
!                psi2 = ((r*r) + r)/((r*r) + 1.0)          ! Van-Albda 
!                psi2 = (abs(r) + r)/(abs(r) + 1.0)          ! Van-Leer

                psi1 = (1 - (1 - psi1)*switch_L )
                psi2 = (1 - (1 - psi2)*switch_L )

                f_qp_left(i+ii, j+jj, k+kk, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1.-kappa) * psi1 * bd) + ((1.+kappa) * psi2 * fd))
                f_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1.+kappa) * psi1 * bd) + ((1.-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do


        end subroutine compute_face_state
        
        
        subroutine compute_muscl_states(qp, x_qp_l, x_qp_r, y_qp_l, y_qp_r, z_qp_l, z_qp_r, pdif, scheme, flow, dims)
            !< Implement MUSCL scheme to get left and right states at
            !< each face. The computation is done through all cells
            !< and first level ghost cells
            !---------------------------------------------------------
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            real(wp), dimension(0:dims%imx+1,1:dims%jmx-1,1:dims%kmx-1,1:dims%n_var), intent(inout) :: x_qp_l, x_qp_r
            !< Store primitive state at the I-face 
            real(wp), dimension(1:dims%imx-1,0:dims%jmx+1,1:dims%kmx-1,1:dims%n_var), intent(inout) :: y_qp_l, y_qp_r
            !< Store primitive state at the J-face 
            real(wp), dimension(1:dims%imx-1,1:dims%jmx-1,0:dims%kmx+1,1:dims%n_var), intent(inout) :: z_qp_l, z_qp_r
            !< Store primitive state at the K-face 
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(inout) :: pdif
            !< pressure difference
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            type(flowtype), intent(in) ::flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            integer, dimension(3) :: flags
            !< Flags for direction

            
            flags=(/1,0,0/)
            call compute_face_state(qp, x_qp_l, x_qp_r, flags, scheme%ilimiter_switch, scheme%itlimiter_switch, dims)
            if(scheme%iPB_switch==1)then
              call pressure_based_switching(qp, x_qp_l, x_qp_r, pdif, flags, flow, dims)
            end if
            flags=(/0,1,0/)
            call compute_face_state(qp, y_qp_l, y_qp_r, flags, scheme%jlimiter_switch, scheme%jtlimiter_switch, dims)
            if(scheme%jPB_switch==1)then
              call pressure_based_switching(qp, y_qp_l, y_qp_r, pdif, flags, flow, dims)
            end if
            flags=(/0,0,1/)
            call compute_face_state(qp, z_qp_l, z_qp_r, flags, scheme%klimiter_switch, scheme%ktlimiter_switch, dims)
            if(scheme%kPB_switch==1)then
              call pressure_based_switching(qp, z_qp_l, z_qp_r, pdif, flags, flow, dims)
            end if

        end subroutine compute_muscl_states

end module muscl
