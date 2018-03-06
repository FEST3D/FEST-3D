module ausm
    !-------------------------------------------------------------------
    ! The ausm scheme is a type of flux-splitting scheme
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
    use global_vars, only : density
    use global_vars, only : turbulence
    use global_vars, only : process_id
    use global_vars, only : current_iter
    use global_vars, only : max_iters
    use global_vars, only : make_F_flux_zero
    use global_vars, only : make_G_flux_zero
    use global_vars, only : make_H_flux_zero

    use utils, only: alloc, dealloc, dmsg
    use face_interpolant, only: x_qp_left, x_qp_right, y_qp_left, y_qp_right, &
                z_qp_left, z_qp_right, &
            x_density_left, x_x_speed_left, x_y_speed_left, x_z_speed_left, &
                x_pressure_left, &
            x_density_right, x_x_speed_right, x_y_speed_right, x_z_speed_right, &
                x_pressure_right, &
            y_density_left, y_x_speed_left, y_y_speed_left, y_z_speed_left, &
                y_pressure_left, &
            y_density_right, y_x_speed_right, y_y_speed_right, y_z_speed_right, &
                y_pressure_right, &
            z_density_left, z_x_speed_left, z_y_speed_left, z_z_speed_left, &
                z_pressure_left, &
            z_density_right, z_x_speed_right, z_y_speed_right, z_z_speed_right, &
                z_pressure_right
    include "turbulence_models/include/ausm/import_module.inc" !for turbulent variables

    implicit none
    private

    real, public, dimension(:, :, :, :), allocatable, target :: F, G, H, residue
    real, dimension(:, :, :, :), pointer :: flux_p

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: get_residue
    
    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ausm', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - ausm.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - ausm.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - ausm.')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - ausm.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ausm', 'destroy_scheme')
            
            call dealloc(F)
            call dealloc(G)
            call dealloc(H)

        end subroutine destroy_scheme

        subroutine compute_flux(f_dir)

            implicit none
            character, intent(in) :: f_dir
            integer :: i, j, k 
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            real, dimension(:, :, :), pointer :: fA, nx, ny, nz, &
                f_x_speed_left, f_x_speed_right, &
                f_y_speed_left, f_y_speed_right, &
                f_z_speed_left, f_z_speed_right, &
                f_density_left, f_density_right, &
                f_pressure_left, f_pressure_right
            real, dimension(1:n_var) :: F_plus, F_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            real :: temp_c

            !include compute_flux_variable and select.inc before select
            !as it contains variables deceleration
            include "turbulence_models/include/ausm/compute_flux_var.inc"
            include "turbulence_models/include/ausm/compute_flux_select.inc"

            call dmsg(1, 'ausm', 'compute_flux')
            
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
                    f_x_speed_left => x_x_speed_left
                    f_x_speed_right => x_x_speed_right
                    f_y_speed_left => x_y_speed_left
                    f_y_speed_right => x_y_speed_right
                    f_z_speed_left => x_z_speed_left
                    f_z_speed_right => x_z_speed_right
                    f_density_left => x_density_left
                    f_density_right => x_density_right
                    f_pressure_left => x_pressure_left
                    f_pressure_right => x_pressure_right
                case ('y')
                    i_f = 0
                    j_f = 1
                    k_f = 0
                    flux_p => G
                    fA => yA
                    nx => ynx
                    ny => yny
                    nz => ynz
                    f_x_speed_left => y_x_speed_left
                    f_x_speed_right => y_x_speed_right
                    f_y_speed_left => y_y_speed_left
                    f_y_speed_right => y_y_speed_right
                    f_z_speed_left => y_z_speed_left
                    f_z_speed_right => y_z_speed_right
                    f_density_left => y_density_left
                    f_density_right => y_density_right
                    f_pressure_left => y_pressure_left
                    f_pressure_right => y_pressure_right
                case ('z')
                    i_f = 0
                    j_f = 0
                    k_f = 1
                    flux_p => H
                    fA => zA
                    nx => znx
                    ny => zny
                    nz => znz
                    f_x_speed_left => z_x_speed_left
                    f_x_speed_right => z_x_speed_right
                    f_y_speed_left => z_y_speed_left
                    f_y_speed_right => z_y_speed_right
                    f_z_speed_left => z_z_speed_left
                    f_z_speed_right => z_z_speed_right
                    f_density_left => z_density_left
                    f_density_right => z_density_right
                    f_pressure_left => z_pressure_left
                    f_pressure_right => z_pressure_right
                case default
                    call dmsg(5, 'ausm', 'compute_flux', &
                            'Direction not recognised')
                    stop
            end select
            

            do k = 1, kmx - 1 + k_f
             do j = 1, jmx - 1 + j_f 
              do i = 1, imx - 1 + i_f
                sound_speed_avg = 0.5 * (sqrt(gm * f_pressure_left(i, j, k) / &
                                            f_density_left(i, j, k) ) + &
                                          sqrt(gm * f_pressure_right(i, j, k) / &
                                            f_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = f_x_speed_left(i, j, k) * nx(i, j, k) + &
                                     f_y_speed_left(i, j, k) * ny(i, j, k) + &
                                     f_z_speed_left(i, j, k) * nz(i, j, k)
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
                face_normal_speeds = f_x_speed_right(i, j, k) * nx(i, j, k) + &
                                     f_y_speed_right(i, j, k) * ny(i, j, k) + &
                                     f_z_speed_right(i, j, k) * nz(i, j, k)
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
                F_plus(1) = f_density_left(i, j, k) * sound_speed_avg * c_plus
                ! F minus mass flux
                F_minus(1) = f_density_right(i, j, k) * sound_speed_avg * c_minus
                F_plus(1)  = F_plus(1) *(i_f*make_F_flux_zero(i) &
                                       + j_f*make_G_flux_zero(j) &
                                       + k_f*make_H_flux_zero(k))
                F_minus(1) = F_minus(1)*(i_f*make_F_flux_zero(i) &
                                       + j_f*make_G_flux_zero(j) &
                                       + k_f*make_H_flux_zero(k))

                ! Construct other fluxes in terms of the F mass flux
                F_plus(2) = (F_plus(1) * f_x_speed_left(i, j, k)) + &
                            (scrD_plus * f_pressure_left(i, j, k) * nx(i, j, k))
                F_plus(3) = (F_plus(1) * f_y_speed_left(i, j, k)) + &
                            (scrD_plus * f_pressure_left(i, j, k) * ny(i, j, k))
                F_plus(4) = (F_plus(1) * f_z_speed_left(i, j, k)) + &
                            (scrD_plus * f_pressure_left(i, j, k) * nz(i, j, k))
                F_plus(5) = F_plus(1) * &
                            ((0.5 * (f_x_speed_left(i, j, k) ** 2. + &
                                     f_y_speed_left(i, j, k) ** 2. + &
                                     f_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * f_pressure_left(i, j, k) / &
                             f_density_left(i, j, k)))


                
                ! Construct other fluxes in terms of the F mass flux
                F_minus(2) = (F_minus(1) * f_x_speed_right(i, j, k)) + &
                             (scrD_minus * f_pressure_right(i, j, k) * nx(i, j, k))
                F_minus(3) = (F_minus(1) * f_y_speed_right(i, j, k)) + &
                             (scrD_minus * f_pressure_right(i, j, k) * ny(i, j, k))
                F_minus(4) = (F_minus(1) * f_z_speed_right(i, j, k)) + &
                             (scrD_minus * f_pressure_right(i, j, k) * nz(i, j, k))
                F_minus(5) = F_minus(1) * &
                            ((0.5 * (f_x_speed_right(i, j, k) ** 2. + &
                                     f_y_speed_right(i, j, k) ** 2. + &
                                     f_z_speed_right(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * f_pressure_right(i, j, k) / &
                             f_density_right(i, j, k)))
         
                !turbulent fluxes
                include "turbulence_models/include/ausm/Fcompute_flux.inc"

                ! Multiply in the face areas
                F_plus(:) = F_plus(:) * fA(i, j, k)
                F_minus(:) = F_minus(:) * fA(i, j, k)

                ! Get the total flux for a face
                flux_p(i, j, k, :) = F_plus(:) + F_minus(:)
                !write(*, '(2(i1.1,x),5(E20.10,2x))') i, j, flux_p(i,j,k,:)
              end do
             end do
            end do 

            nullify(fA)
            nullify(nx)
            nullify(ny)
            nullify(nz)
            nullify(f_x_speed_left )
            nullify(f_x_speed_right)
            nullify(f_y_speed_left )
            nullify(f_y_speed_right)
            nullify(f_z_speed_left )
            nullify(f_z_speed_right)
            nullify(f_density_left )
            nullify(f_density_right)
            nullify(f_pressure_left)
            nullify(f_pressure_right)

        end subroutine compute_flux

        subroutine compute_fluxes()
            
            implicit none
            
            call dmsg(1, 'ausm', 'compute_fluxes')

            call compute_flux('x')
            if (any(isnan(F))) then
                call dmsg(5, 'ausm', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_flux('y')
            if (any(isnan(G))) then 
                call dmsg(5, 'ausm', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            !if(kmx==2) then
            !  H = 0.
            !else
              call compute_flux('z')
            !end if
            if (any(isnan(H))) then
                call dmsg(5, 'ausm', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        subroutine get_residue()
            !-----------------------------------------------------------
            ! Compute the residue for the AUSM scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l

            call dmsg(1, 'ausm', 'compute_residue')

            do l = 1, n_var
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 1, imx - 1
               residue(i, j, k, l) = (F(i+1, j, k, l) - F(i, j, k, l)) &
                                   + (G(i, j+1, k, l) - G(i, j, k, l)) &
                                   + (H(i, j, k+1, l) - H(i, j, k, l))
               !print*, i, j, l, F(i+1, j,k,l), F(i,j,k,l), G(i,j+1,k,l)&
               !       , G(i,j,k,l), H(i,j,k+1,l), H(i,j,k,l)
               end do
              end do
             end do
            end do
        
        end subroutine get_residue

end module ausm
