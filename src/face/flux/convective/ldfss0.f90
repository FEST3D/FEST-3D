    !< Flux-splitting scheme: LDFSS
module ldfss0
    !<
    !< Reference: Edwards, J.R., A low-diffusion flux-splitting scheme 
    !< for Navier-Stokes calculations, Computers & Fluids, vol. 26,
    !< no. 6, pp.635-659, 1997
    !-------------------------------------------------------------------

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx

    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z
    use global_vars, only : xA, yA, zA    !face area

    use global_vars, only : gm
    use global_vars, only : n_var
    use global_vars, only : turbulence
    use global_vars, only : make_F_flux_zero
    use global_vars, only : make_G_flux_zero
    use global_vars, only : make_H_flux_zero

    use utils, only: alloc, dealloc, dmsg
    use face_interpolant, only: x_qp_left, x_qp_right 
    use face_interpolant, only: y_qp_left, y_qp_right
    use face_interpolant, only:  z_qp_left, z_qp_right

    implicit none
    private

    real, public, dimension(:, :, :, :), allocatable, target :: F
    !< Store fluxes throught the I faces
    real, public, dimension(:, :, :, :), allocatable, target :: G
    !< Store fluxes throught the J faces
    real, public, dimension(:, :, :, :), allocatable, target :: H
    !< Store fluxes throught the K faces
    real, public, dimension(:, :, :, :), allocatable, target :: residue
    !< Store residue at each cell-center
    real, dimension(:, :, :, :), pointer :: flux_p
    !< Pointer/alias for the either F, G, or H

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: get_residue
    
    contains

        subroutine setup_scheme()
          !< Allocate memory to the flux variables

            implicit none

            call dmsg(1, 'ldfss0', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - ldfss0.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - ldfss0.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - ldfss0.')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - ldfss0.')

        end subroutine setup_scheme

        subroutine destroy_scheme()
          !< Deallocate memory

            implicit none

            call dmsg(1, 'ldfss0', 'destroy_scheme')
            
            call dealloc(F)
            call dealloc(G)
            call dealloc(H)

        end subroutine destroy_scheme

        subroutine compute_flux(f_dir)
          !< A generalized subroutine to calculate
          !< flux through the input direction, :x,y, or z
          !< which corresponds to the I,J, or K direction respectively
          !------------------------------------------------------------

            implicit none
            character, intent(in) :: f_dir
            !< Input direction for which flux are calcuated and store
            integer :: i, j, k 
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            real, dimension(:, :, :), pointer :: fA, nx, ny, nz
            real, dimension(:,:,:,:), pointer :: f_qp_left, f_qp_right
            real, dimension(1:n_var) :: F_plus, F_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            real :: M_ldfss, M_plus_ldfss, M_minus_ldfss

            call dmsg(1, 'ldfss0', 'compute_flux')
            
            select case (f_dir)
                case ('x')
                    i_f = 1
                    j_f = 0
                    k_f = 0
                    flux_p => F
                    fA => xA
                    nx => xnx
                    ny => xny
                    nz => xnz
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                case ('y')
                    i_f = 0
                    j_f = 1
                    k_f = 0
                    flux_p => G
                    fA => yA
                    nx => ynx
                    ny => yny
                    nz => ynz
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                case ('z')
                    i_f = 0
                    j_f = 0
                    k_f = 1
                    flux_p => H
                    fA => zA
                    nx => znx
                    ny => zny
                    nz => znz
                    f_qp_left => z_qp_left
                    f_qp_right => z_qp_right
                case default
                    call dmsg(5, 'ldfss0', 'compute_flux', &
                            'Direction not recognised')
                    stop
            end select

            do k = 1, kmx - 1 + k_f
             do j = 1, jmx - 1 + j_f 
              do i = 1, imx - 1 + i_f
                sound_speed_avg = 0.5 * (sqrt(gm * f_qp_left(i, j, k,5) / &
                                            f_qp_left(i, j, k,1) ) + &
                                          sqrt(gm * f_qp_right(i, j, k,5) / &
                                            f_qp_right(i, j, k,1) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = f_qp_left(i, j, k,2) * nx(i, j, k) + &
                                     f_qp_left(i, j, k,3) * ny(i, j, k) + &
                                     f_qp_left(i, j, k,4) * nz(i, j, k)
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
                face_normal_speeds = f_qp_right(i, j, k, 2) * nx(i, j, k) + &
                                     f_qp_right(i, j, k, 3) * ny(i, j, k) + &
                                     f_qp_right(i, j, k, 4) * nz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = -0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)
                
                ! LDFSS0 modification             
                M_ldfss = 0.25 * beta_left * beta_right * &
                    (sqrt((M_perp_left ** 2 + M_perp_right ** 2) * 0.5) &
                     - 1) ** 2
                M_plus_ldfss = M_ldfss * &
                        (1 - (f_qp_left(i, j, k,5) - f_qp_right(i, j, k,5)) / &
                            (2 * f_qp_left(i, j, k,1) * (sound_speed_avg ** 2)))
                M_minus_ldfss = M_ldfss * &
                        (1 - (f_qp_left(i, j, k,5) - f_qp_right(i, j, k,5)) / &
                            (2 * f_qp_right(i, j, k,1) * (sound_speed_avg ** 2)))
                c_plus = c_plus - M_plus_ldfss
                c_minus = c_minus + M_minus_ldfss

                ! F plus mass flux
                F_plus(1) = f_qp_left(i, j, k,1) * sound_speed_avg * c_plus
                ! F minus mass flux
                F_minus(1) = f_qp_right(i, j, k,1) * sound_speed_avg * c_minus
                F_plus(1)  = F_plus(1) *(i_f*make_F_flux_zero(i) &
                                       + j_f*make_G_flux_zero(j) &
                                       + k_f*make_H_flux_zero(k))
                F_minus(1) = F_minus(1)*(i_f*make_F_flux_zero(i) &
                                       + j_f*make_G_flux_zero(j) &
                                       + k_f*make_H_flux_zero(k))


                ! Construct other fluxes in terms of the F mass flux
                F_plus(2) = (F_plus(1) * f_qp_left(i, j, k,2)) + &
                            (scrD_plus * f_qp_left(i, j, k,5) * nx(i, j, k))
                F_plus(3) = (F_plus(1) * f_qp_left(i, j, k,3)) + &
                            (scrD_plus * f_qp_left(i, j, k,5) * ny(i, j, k))
                F_plus(4) = (F_plus(1) * f_qp_left(i, j, k,4)) + &
                            (scrD_plus * f_qp_left(i, j, k,5) * nz(i, j, k))
                F_plus(5) = F_plus(1) * &
                            ((0.5 * (f_qp_left(i, j, k,2) ** 2. + &
                                     f_qp_left(i, j, k,3) ** 2. + &
                                     f_qp_left(i, j, k,4) ** 2.)) + &
                            ((gm / (gm - 1.)) * f_qp_left(i, j, k,5) / &
                             f_qp_left(i, j, k,1)))


                
                ! Construct other fluxes in terms of the F mass flux
                F_minus(2) = (F_minus(1) * f_qp_right(i, j, k,2)) + &
                             (scrD_minus * f_qp_right(i, j, k,5) * nx(i, j, k))
                F_minus(3) = (F_minus(1) * f_qp_right(i, j, k,3)) + &
                             (scrD_minus * f_qp_right(i, j, k,5) * ny(i, j, k))
                F_minus(4) = (F_minus(1) * f_qp_right(i, j, k,4)) + &
                             (scrD_minus * f_qp_right(i, j, k,5) * nz(i, j, k))
                F_minus(5) = F_minus(1) * &
                            ((0.5 * (f_qp_right(i, j, k,2) ** 2. + &
                                     f_qp_right(i, j, k,3) ** 2. + &
                                     f_qp_right(i, j, k,4) ** 2.)) + &
                            ((gm / (gm - 1.)) * f_qp_right(i, j, k,5) / &
                             f_qp_right(i, j, k,1)))
         
                !turbulent fluxes
                if(n_var>5) then
                  F_plus(6:)  = F_Plus(1)  * f_qp_left(i,j,k,6:)
                  F_minus(6:) = F_minus(1) * f_qp_right(i,j,k,6:)
                end if

                ! Multiply in the face areas
                F_plus(:) = F_plus(:) * fA(i, j, k)
                F_minus(:) = F_minus(:) * fA(i, j, k)

                ! Get the total flux for a face
                flux_p(i, j, k, :) = F_plus(:) + F_minus(:)
              end do
             end do
            end do 

        end subroutine compute_flux

        subroutine compute_fluxes()
          !< Call to compute fluxes throught faces in each direction
            
            implicit none
            
            call dmsg(1, 'ldfss0', 'compute_fluxes')

            call compute_flux('x')
            if (any(isnan(F))) then
                call dmsg(5, 'ldfss0', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_flux('y')
            if (any(isnan(G))) then 
                call dmsg(5, 'ldfss0', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            call compute_flux('z')
            if (any(isnan(H))) then
                call dmsg(5, 'ldfss0', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        subroutine get_residue()
            !< Compute the residue for the ldfss0 scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l

            call dmsg(1, 'ldfss0', 'compute_residue')

            do l = 1, n_var
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 1, imx - 1
               residue(i, j, k, l) = F(i+1, j, k, l) - F(i, j, k, l) &
                                   + G(i, j+1, k, l) - G(i, j, k, l) &
                                   + H(i, j, k+1, l) - H(i, j, k, l)
               end do
              end do
             end do
            end do
        
        end subroutine get_residue

end module ldfss0
