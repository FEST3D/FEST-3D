    !< Flux splitting scheme: AUSM+
module van_leer
    !<
    !< Reference: Van Leer, B., Flux-vector splitting for the Euler equation. 
    !< In Upwind and High-Resolution Schemes, Springer, Berlin, Heidelberg, pp. 80-89, 1997
    !-------------------------------------------------------------------
#include "../../../debug.h"
#include "../../../error.h"
    use vartypes
!    use utils, only: alloc
!    use global_vars, only : xnx, xny, xnz !face unit normal x
!    use global_vars, only : ynx, yny, ynz !face unit normal y
!    use global_vars, only : znx, zny, znz !face unit normal z
!    use global_vars, only : xA, yA, zA    !face area
    use global_vars, only : process_id
    use global_vars, only : make_F_flux_zero
    use global_vars, only : make_G_flux_zero
    use global_vars, only : make_H_flux_zero
    use face_interpolant, only: x_qp_left, x_qp_right 
    use face_interpolant, only: y_qp_left, y_qp_right
    use face_interpolant, only:  z_qp_left, z_qp_right
    implicit none
    private
!
!    real, public, dimension(:, :, :, :), allocatable, target :: F, G, H, residue
!    !< Array to store the flux through I, J, and K direction face and residue in the cell
!    real, dimension(:, :, :, :), pointer :: flux_p
!    !< A general flux pointer
!
!    integer :: imx, jmx, kmx, n_var
!
!    ! Public members
!    public :: setup_scheme
!    public :: destroy_scheme
    public :: compute_fluxes
!    public :: get_residue
    
    contains

!        subroutine setup_scheme(control, dims)
!          !< Allocate memory to the flux variables
!
!            implicit none
!            type(controltype), intent(in) :: control
!            type(extent), intent(in) :: dims
!
!            imx = dims%imx
!            jmx = dims%jmx
!            kmx = dims%kmx
!
!            n_var = control%n_var
!
!            DebugCall('setup_scheme')
!
!            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'F - van_leer.')
!            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'G - van_leer.')
!            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'H - van_leer.')
!            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'residue - van_leer.')
!
!        end subroutine setup_scheme

!        subroutine destroy_scheme()
!          !< Deallocate memory
!
!            implicit none
!
!            DebugCall('destroy_scheme')
!            
!            call dealloc(F)
!            call dealloc(G)
!            call dealloc(H)
!
!        end subroutine destroy_scheme
!
        subroutine compute_flux(Flux, f_dir, faces, flags, flow, dims)
          !< A generalized subroutine to calculate
          !< flux through the input-argument direction, :x,y, or z
          !< which corresponds to the I,J, or K direction respectively
          !------------------------------------------------------------

            implicit none
            integer, dimension(3), intent(in) :: flags
            type(extent), intent(in) :: dims
            type(flowtype), intent(in) :: flow
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

            DebugCall('compute_flux')
            i_f = flags(1)
            j_f = flags(2)
            k_f = flags(3)
            
            select case (f_dir)
                case ('x')
                    !i_f = 1
                    !j_f = 0
                    !k_f = 0
                    !!flux_p => F
                    !fA => xA
                    !nx => xnx
                    !ny => xny
                    !nz => xnz
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                case ('y')
                    !i_f = 0
                    !j_f = 1
                    !k_f = 0
                    !!flux_p => G
                    !fA => yA
                    !nx => ynx
                    !ny => yny
                    !nz => ynz
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                case ('z')
                    !i_f = 0
                    !j_f = 0
                    !k_f = 1
                    !!flux_p => H
                    !fA => zA
                    !nx => znx
                    !ny => zny
                    !nz => znz
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

                ! First construct the F plus mass flux
                F_plus(1) = f_qp_left(i, j, k,1) * sound_speed_avg * c_plus
                ! First construct the F minus mass flux
                F_minus(1) = f_qp_right(i, j, k,1) * sound_speed_avg * c_minus
                F_plus(1)  = F_plus(1) *(i_f*make_F_flux_zero(i) &
                                       + j_f*make_G_flux_zero(j) &
                                       + k_f*make_H_flux_zero(k))
                F_minus(1) = F_minus(1)*(i_f*make_F_flux_zero(i) &
                                       + j_f*make_G_flux_zero(j) &
                                       + k_f*make_H_flux_zero(k))


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


        subroutine compute_fluxes(F,G,H, Ifaces, Jfaces, Kfaces, flow, dims)
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
            integer, dimension(3) :: flags

            
            DebugCall('compute_fluxes')

            flags=(/1,0,0/)
            call compute_flux(F, 'x', Ifaces, flags, flow, dims)
            if (any(isnan(F))) then
              Fatal_error
            end if    

            flags=(/0,1,0/)
            call compute_flux(G, 'y', Jfaces, flags, flow, dims)
            if (any(isnan(G))) then 
              Fatal_error
            end if    
            
            if(dims%kmx==2) then
              H = 0.
            else
              flags=(/0,0,1/)
              call compute_flux(H, 'z', Kfaces, flags, flow, dims)
            end if
            if (any(isnan(H))) then
              Fatal_error
            end if

        end subroutine compute_fluxes
!
!        subroutine get_residue()
!            !< Compute the residue using the Van-Leer scheme
!            !-----------------------------------------------------------
!
!            implicit none
!            
!            integer :: i, j, k, l
!
!            DebugCall('compute_residue')
!
!            do l = 1, n_var
!             do k = 1, kmx - 1
!              do j = 1, jmx - 1
!               do i = 1, imx - 1
!               residue(i, j, k, l) = F(i+1, j, k, l) - F(i, j, k, l) &
!                                   + G(i, j+1, k, l) - G(i, j, k, l) &
!                                   + H(i, j, k+1, l) - H(i, j, k, l)
!               end do
!              end do
!             end do
!            end do
!        
!        end subroutine get_residue

end module van_leer
