    !< Higher order face state reconstruction method:PPM
module ppm
    !<
    !<Reference: Colella, P. and Woodward, P.R., The piecewise 
    !<parabolic method (PPM) for gas-dynamical simulations, Journal
    !<of computational physics, vol. 54, no. 1, pp.174-201, 1984
    !-------------------------------------------------------------------

    use vartypes
#include "../../debug.h"
#include "../../error.h"

    implicit none
    private

    ! Public members
    public :: compute_ppm_states

    contains

        subroutine compute_face_estimates(qp, f_qp_left, f_qp_right, flags, dims)
          !< Subroutine to calculate state at the face, generalized for

            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            integer, dimension(3), intent(in) :: flags
            !< Flags for direction switch
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2),&
            1-flags(3):dims%kmx-1+2*flags(3),1:dims%n_var), intent(inout) :: f_qp_left, f_qp_right
            !< primitive state at faces
            integer :: i, j, k 
            integer :: i_f, j_f, k_f ! Flags to determine face direction

            DebugCall('compute_face_estimates')
            
            i_f = flags(1)
            j_f = flags(2)
            k_f = flags(3)

            ! Interior faces
            do k = (1 - k_f), dims%kmx - 1 + 2*k_f
             do j = (1 - j_f), dims%jmx - 1 + 2*j_f
              do i = (1 - i_f), dims%imx - 1 + 2*i_f
                f_qp_left(i, j, k, :) = (7. * (qp(i, j, k, :) + &
                    qp(i - i_f, j - j_f, k - k_f, :)) - (qp(i + i_f, j + j_f, k + k_f, :) + &
                    qp(i - 2*i_f, j - 2*j_f, k - 2*k_f, :))) / 12.
              end do
             end do
            end do
            f_qp_right= f_qp_left

        end subroutine compute_face_estimates

        subroutine remove_extrema(qp, f_qp_left, f_qp_right, flags, dims)
          !< Remove extrema from the state estimated. 
          !< Limiting the value in case of PPM

            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            integer, dimension(3), intent(in) :: flags
            !< flags for direction switch
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2),&
            1-flags(3):dims%kmx-1+2*flags(3),1:dims%n_var), intent(inout) :: f_qp_left, f_qp_right
            !< primitve state variable at faces
            integer :: i, j, k, l
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            real(wp) :: dqrl, dq6

            DebugCall('remove_extrema')
            
            i_f = flags(1)
            j_f = flags(2)
            k_f = flags(3)
            
            ! Loop over cells (including ghost cells)
            do l = 1, dims%n_var            
             do k = 1 - k_f, dims%kmx - 1 + k_f
              do j = 1 - j_f, dims%jmx - 1 + j_f
               do i = 1 - i_f, dims%imx - 1 + i_f
                if ((f_qp_left(i+i_f, j+j_f, k+k_f, l) - qp(i, j, k, l)) * &
                    (qp(i, j, k, l) - f_qp_right(i, j, k, l)) <= 0) then
                    f_qp_left(i+i_f, j+j_f, k+k_f, l) = qp(i, j, k, l)
                    f_qp_right(i, j, k, l) = qp(i, j, k, l)
                else      
                    dqrl = f_qp_left(i+i_f, j+j_f, k+k_f, l) - f_qp_right(i, j, k, l)
                    dq6 = 6. * (qp(i, j, k, l) - 0.5*(f_qp_left(i+i_f, j+j_f, k+k_f, l) + &
                                                      f_qp_right(i, j, k, l)))
                    if (dqrl * dq6 > dqrl*dqrl) then
                        f_qp_right(i, j, k, l) = 3.*qp(i, j, k, l) - &
                                                 2.*f_qp_left(i+i_f, j+j_f, k+k_f, l)
                    else if (-dqrl*dqrl > dqrl * dq6) then
                        f_qp_left(i+i_f, j+j_f, k+k_f, l) = 3.*qp(i, j, k, l) - &
                                                 2.*f_qp_right(i, j, k, l)
                    end if
                end if
               end do
              end do 
             end do
            end do

        end subroutine remove_extrema

        subroutine pressure_based_switching(qp, f_qp_left, f_qp_right, pdif, flags, flow, dims)
          !< Pressure based switching. 
          !< User x,y, or z for I,J,or K face respectively
          !----------------------------------------------

            implicit none
            type(extent), intent(in) :: dims
            !< Extent of the domain:imx,jmx,kmx
            integer, dimension(3), intent(in) :: flags
            !< flags for direction switch
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            !< Store primitive variable at cell center
            real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2),&
            1-flags(3):dims%kmx-1+2*flags(3),1:dims%n_var), intent(inout) :: f_qp_left, f_qp_right
            !< primitive state at faces
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(inout) :: pdif
            !< pressure difference
            type(flowtype), intent(in) :: flow
            ! Character can be x or y or z
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
                pd2 = abs(qp(i + i_f*1, j + j_f*1, k + k_f*1, 5) - &
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
            pdif(((dims%imx-1)*i_f)+1:dims%imx-1+i_f, &
                 ((dims%jmx-1)*j_f)+1:dims%jmx-1+j_f, &
                 ((dims%kmx-1)*k_f)+1:dims%kmx-1+k_f) &
                                        =   &
                pdif(i_f*(dims%imx-2)+1:dims%imx-1, &
                     j_f*(dims%jmx-2)+1:dims%jmx-1, &
                     k_f*(dims%kmx-2)+1:dims%kmx-1)

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
        
     
        subroutine compute_ppm_states(qp, x_qp_l, x_qp_r, y_qp_l, y_qp_r, z_qp_l, z_qp_r, pdif, scheme, flow, dims)
          !< Call PPM face-state reconstruction for each face
          !< with optional call for remove extrema based on
          !< input limter switch and call pressure based switching
          !< based on input pressure based switch

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
            type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            integer, dimension(3) :: flags
            !< flags for different directions

            flags=(/1,0,0/)
            call compute_face_estimates(qp, x_qp_l, x_qp_r, flags, dims)
            if(scheme%ilimiter_switch==1)then
              call remove_extrema(qp, x_qp_l, x_qp_r, flags, dims)
            end if
            if (scheme%iPB_switch==1)then
              call pressure_based_switching(qp, x_qp_l, x_qp_r, pdif, flags, flow, dims)
            end if

            flags=(/0,1,0/)
            call compute_face_estimates(qp, y_qp_l, y_qp_r, flags, dims)
            if(scheme%jlimiter_switch==1)then
              call remove_extrema(qp, y_qp_l, y_qp_r, flags, dims)
            end if
            if (scheme%jPB_switch==1)then
              call pressure_based_switching(qp, y_qp_l, y_qp_r, pdif, flags, flow, dims)
            end if

            flags=(/0,0,1/)
            call compute_face_estimates(qp, z_qp_l, z_qp_r, flags, dims)
            if(scheme%klimiter_switch==1)then
              call remove_extrema(qp, z_qp_l, z_qp_r, flags, dims)
            end if
            if (scheme%kPB_switch==1)then
              call pressure_based_switching(qp, z_qp_l, z_qp_r, pdif, flags, flow, dims)
            end if

        end subroutine compute_ppm_states


end module ppm
