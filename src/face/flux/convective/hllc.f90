module hllc
    !-------------------------------------------------------------------
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
    use global_vars, only : process_id
    use global_vars, only : current_iter
    use global_vars, only : max_iters
    use global_vars, only : imin_id
    use global_vars, only : imax_id
    use global_vars, only : jmin_id
    use global_vars, only : jmax_id
    use global_vars, only : kmin_id
    use global_vars, only : kmax_id
    use global_vars, only : merror

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
    include "turbulence_models/include/hllc/import_module.inc" !for turbulent variables

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

            call dmsg(1, 'hllc', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - hllc.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - hllc.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - hllc.')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - hllc.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'hllc', 'destroy_scheme')
            
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
            integer :: id
            real :: f_sound_left
            real :: f_sound_right
            real :: rho_sound_left
            real :: rho_sound_right
            real :: rho_speed_left
            real :: rho_speed_right

            !include compute_flux_variable and select.inc before select
            !as it contains variables deceleration
            include "turbulence_models/include/hllc/compute_flux_var.inc"

            call dmsg(1, 'hllc', 'compute_flux')
            
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
                    call dmsg(5, 'hllc', 'compute_flux', &
                            'Direction not recognised')
                    stop
            end select
            

            do k = 1, kmx - 1 + k_f
             do j = 1, jmx - 1 + j_f 
              do i = 1, imx - 1 + i_f

               f_sound_left = sqrt(gm*f_pressure_left(i,j,k)/f_density_left(i,j,k))
               f_sound_right = sqrt(gm*f_pressure_right(i,j,k)/f_density_right(i,j,k))
               f_face_normal_speed_left = f_x_speed_left(i,j,k)*nx &
                                        + f_y_speed_left(i,j,k)*ny &
                                        + f_z_speed_left(i,j,k)*nz
               f_face_normal_speed_right= f_x_speed_right(i,j,k)*nx &
                                        + f_y_speed_right(i,j,k)*ny &
                                        + f_z_speed_right(i,j,k)*nz
               f_face_normal_speed_avg = (sqrt(f_density_left(i,j,k))*face_normal_speed_left &
                                       +  sqrt(f_density_right(i,j,k))*face_normal_speed_right)&
                                       /( sqrt(f_density_left(i,j,k))+sqrt(f_density_right(i,j,k)))
               lamda1uf      = (f_x_speed(i,j,k)**2 + f_y_speed(i,j,k)**2 + f_z_speed(i,j,k)**2)



                ! F plus mass flux
                F_plus(1) = f_density_left(i, j, k) * sound_speed_avg * c_plus
                ! F minus mass flux
                F_minus(1) = f_density_right(i, j, k) * sound_speed_avg * c_minus
                include "mass_flux.inc"

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
                include "turbulence_models/include/hllc/Fcompute_flux.inc"

                ! Multiply in the face areas
                F_plus(:) = F_plus(:) * fA(i, j, k)
                F_minus(:) = F_minus(:) * fA(i, j, k)

                ! Get the total flux for a face
                flux_p(i, j, k, :) = F_plus(:) + F_minus(:)
                !if(process_id==3 .and. f_dir=='x' ) print*, i,j,k, f_pressure_left(i,j,k), f_pressure_right(i,j,k), f_density_left(i,j,k), f_density_right(i,j,k), sound_speed_avg
              end do
             end do
            end do 

        end subroutine compute_flux

        subroutine compute_fluxes()
            
            implicit none
            
            call dmsg(1, 'hllc', 'compute_fluxes')

            call compute_flux('x')
            if (any(isnan(F))) then
                call dmsg(5, 'hllc', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_flux('y')
            if (any(isnan(G))) then 
                call dmsg(5, 'hllc', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            if(kmx==2) then
              H = 0.
            else
              call compute_flux('z')
            end if
            if (any(isnan(H))) then
                call dmsg(5, 'hllc', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        subroutine get_residue()
            !-----------------------------------------------------------
            ! Compute the residue for the hllc scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l

            call dmsg(1, 'hllc', 'compute_residue')

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

end module hllc
