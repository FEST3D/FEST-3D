    !< Flux splitting scheme: AUSM
module ausm
    !<
    !< Reference: Liou, M.S. and Steffen Jr, C.J., A new flux splitting scheme, 
    !< Journal of Computational physics, vol. 107, no. 1, pp.23-39, 1993
    !-------------------------------------------------------------------
#include "../../../debug.h"
#include "../../../error.h"
     use vartypes
!    use global_vars, only : make_F_flux_zero
!    use global_vars, only : make_G_flux_zero
!    use global_vars, only : make_H_flux_zero
    use face_interpolant, only: x_qp_left, x_qp_right 
    use face_interpolant, only: y_qp_left, y_qp_right
    use face_interpolant, only:  z_qp_left, z_qp_right

    implicit none
    private
    public :: compute_fluxes
    
    contains

        subroutine compute_flux(Flux, f_dir, faces, flags, flow, bc, dims)
          !< A generalized subroutine to calculate
          !< flux through the input-argument direction, :x,y, or z
          !< which corresponds to the I,J, or K direction respectively
          !------------------------------------------------------------

            implicit none
            integer, dimension(3), intent(in) :: flags
            type(extent), intent(in) :: dims
            type(flowtype), intent(in) :: flow
            type(boundarytype), intent(in) :: bc
            real, dimension(:, :, :, :), intent(inout) :: Flux
            !< Store fluxes throught the any(I,J,K) faces
            type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
            character, intent(in) :: f_dir
            !< Input direction for which flux are calcuated and store
            integer :: i, j, k 
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            !real, dimension(:, :, :), pointer :: fA, nx, ny, nz
            real, dimension(:,:,:,:), pointer :: f_qp_left, f_qp_right
            real, dimension(1:dims%n_var) :: F_plus, F_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            real :: temp_c

            DebugCall('compute_flux')
            i_f = flags(1)
            j_f = flags(2)
            k_f = flags(3)
            
            select case (f_dir)
                case ('x')
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                case ('y')
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                case ('z')
                    f_qp_left => z_qp_left
                    f_qp_right => z_qp_right
                case default
                    Fatal_error
            end select
            

            do k = 1, dims%kmx - 1 + k_f
             do j = 1, dims%jmx - 1 + j_f 
              do i = 1, dims%imx - 1 + i_f
                sound_speed_avg = 0.5 * (sqrt(flow%gm * f_qp_left(i, j, k,5) / &
                                            f_qp_left(i, j, k,1) ) + &
                                          sqrt(flow%gm * f_qp_right(i, j, k,5) / &
                                            f_qp_right(i, j, k,1) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = f_qp_left(i, j, k,2) * faces(i, j, k)%nx + &
                                     f_qp_left(i, j, k,3) * faces(i, j, k)%ny + &
                                     f_qp_left(i, j, k,4) * faces(i, j, k)%nz
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! Compute '-' direction quantities
                face_normal_speeds = f_qp_right(i, j, k, 2) * faces(i, j, k)%nx + &
                                     f_qp_right(i, j, k, 3) * faces(i, j, k)%ny + &
                                     f_qp_right(i, j, k, 4) * faces(i, j, k)%nz
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = -0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)
                
                ! AUSM modification             
                temp_c = c_plus + c_minus
                c_plus = max(0., temp_c)
                c_minus = min(0., temp_c)



                ! F plus mass flux
                F_plus(1) = f_qp_left(i, j, k,1) * sound_speed_avg * c_plus
                ! F minus mass flux
                F_minus(1) = f_qp_right(i, j, k,1) * sound_speed_avg * c_minus
                F_plus(1)  = F_plus(1) *(i_f*bc%make_F_flux_zero(i) &
                                       + j_f*bc%make_G_flux_zero(j) &
                                       + k_f*bc%make_H_flux_zero(k))
                F_minus(1) = F_minus(1)*(i_f*bc%make_F_flux_zero(i) &
                                       + j_f*bc%make_G_flux_zero(j) &
                                       + k_f*bc%make_H_flux_zero(k))

                ! Construct other fluxes in terms of the F mass flux
                F_plus(2) = (F_plus(1) * f_qp_left(i, j, k,2)) + &
                            (scrD_plus * f_qp_left(i, j, k,5) * faces(i, j, k)%nx)
                F_plus(3) = (F_plus(1) * f_qp_left(i, j, k,3)) + &
                            (scrD_plus * f_qp_left(i, j, k,5) * faces(i, j, k)%ny)
                F_plus(4) = (F_plus(1) * f_qp_left(i, j, k,4)) + &
                            (scrD_plus * f_qp_left(i, j, k,5) * faces(i, j, k)%nz)
                F_plus(5) = F_plus(1) * &
                            ((0.5 * (f_qp_left(i, j, k,2) ** 2. + &
                                     f_qp_left(i, j, k,3) ** 2. + &
                                     f_qp_left(i, j, k,4) ** 2.)) + &
                            ((flow%gm / (flow%gm - 1.)) * f_qp_left(i, j, k,5) / &
                             f_qp_left(i, j, k,1)))


                
                ! Construct other fluxes in terms of the F mass flux
                F_minus(2) = (F_minus(1) * f_qp_right(i, j, k,2)) + &
                             (scrD_minus * f_qp_right(i, j, k,5) * faces(i, j, k)%nx)
                F_minus(3) = (F_minus(1) * f_qp_right(i, j, k,3)) + &
                             (scrD_minus * f_qp_right(i, j, k,5) * faces(i, j, k)%ny)
                F_minus(4) = (F_minus(1) * f_qp_right(i, j, k,4)) + &
                             (scrD_minus * f_qp_right(i, j, k,5) * faces(i, j, k)%nz)
                F_minus(5) = F_minus(1) * &
                            ((0.5 * (f_qp_right(i, j, k,2) ** 2. + &
                                     f_qp_right(i, j, k,3) ** 2. + &
                                     f_qp_right(i, j, k,4) ** 2.)) + &
                            ((flow%gm / (flow%gm - 1.)) * f_qp_right(i, j, k,5) / &
                             f_qp_right(i, j, k,1)))
         
                !turbulent fluxes
                if(dims%n_var>5) then
                  F_plus(6:)  = F_Plus(1)  * f_qp_left(i,j,k,6:)
                  F_minus(6:) = F_minus(1) * f_qp_right(i,j,k,6:)
                end if

                ! Multiply in the face areas
                F_plus(:) = F_plus(:) * faces(i, j, k)%A
                F_minus(:) = F_minus(:) * faces(i, j, k)%A

                ! Get the total flux for a face
                Flux(i, j, k, :) = F_plus(:) + F_minus(:)
              end do
             end do
            end do 


        end subroutine compute_flux

        subroutine compute_fluxes(F,G,H, Ifaces, Jfaces, Kfaces, flow, bc, dims)
        !subroutine compute_fluxes(F,G,H, flow, dims)
          !< Call to compute fluxes throught faces in each direction
            
            implicit none
            type(extent), intent(in) :: dims
            type(flowtype), intent(in) :: flow
            real, dimension(:, :, :, :), intent(inout) :: F
            !< Store fluxes throught the I faces
            real, dimension(:, :, :, :), intent(inout) :: G
            !< Store fluxes throught the J faces
            real, dimension(:, :, :, :), intent(inout) :: H
            !< Store fluxes throught the K faces
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces 
            type(boundarytype), intent(in) :: bc
            integer, dimension(3) :: flags

            
            DebugCall('compute_fluxes')

            flags=(/1,0,0/)
            call compute_flux(F, 'x', Ifaces, flags, flow, bc, dims)
            if (any(isnan(F))) then
              Fatal_error
            end if    

            flags=(/0,1,0/)
            call compute_flux(G, 'y', Jfaces, flags, flow, bc, dims)
            if (any(isnan(G))) then 
              Fatal_error
            end if    
            
            if(dims%kmx==2) then
              H = 0.
            else
              flags=(/0,0,1/)
              call compute_flux(H, 'z', Kfaces, flags, flow, bc, dims)
            end if
            if (any(isnan(H))) then
              Fatal_error
            end if

        end subroutine compute_fluxes


end module ausm
