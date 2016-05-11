module van_leer
    !-------------------------------------------------------------------
    ! The Van-Leer scheme is a type of flux-splitting scheme
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx
!                   sphere_indices, n_sph_ind
    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz, xA, yA, zA
    use state, only: gm, n_var
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

    implicit none
    private

    real, public, dimension(:, :, :, :), allocatable, target :: F, G, H
    real, dimension(:, :, :, :), pointer :: flux_p

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: get_residue
    
    !TODO: Viscous: Change to single subroutine for all directions  

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'van_leer', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - van_leer.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - van_leer.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - van_leer.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'van_leer', 'destroy_scheme')
            
            call dealloc(F)
            call dealloc(G)
            call dealloc(H)

        end subroutine destroy_scheme

!       subroutine compute_flux(f_dir)

!           implicit none
!           character, intent(in) :: f_dir
!           integer :: i, j, k 
!           integer :: i_f, j_f, k_f ! Flags to determine face direction
!           real, dimension(:, :, :), pointer :: fA, nx, ny, nz, &
!               f_x_speed_left, f_x_speed_right, &
!               f_y_speed_left, f_y_speed_right, &
!               f_z_speed_left, f_z_speed_right, &
!               f_density_left, f_density_right, &
!               f_pressure_left, f_pressure_right

!           call dmsg(1, 'ppm', 'compute_face_estimates')
!           
!           select case (f_dir)
!               case ('x')
!                   i_f = 1
!                   j_f = 0
!                   k_f = 0
!                   flux_p => F
!                   fA => xA
!                   nx => xnx
!                   ny => xny
!                   nz => xnz
!                   f_x_speed_left => x_x_speed_left
!                   f_x_speed_right => x_x_speed_right
!                   f_y_speed_left => x_y_speed_left
!                   f_y_speed_right => x_y_speed_right
!                   f_z_speed_left => x_z_speed_left
!                   f_z_speed_right => x_z_speed_right
!                   f_density_left => x_density_left
!                   f_density_right => x_density_right
!                   f_pressure_left => x_pressure_left
!                   f_pressure_right => x_pressure_right
!               case ('y')
!                   i_f = 0
!                   j_f = 1
!                   k_f = 0
!                   flux_p => G
!                   fA => yA
!                   nx => ynx
!                   ny => yny
!                   nz => ynz
!                   f_x_speed_left => y_x_speed_left
!                   f_x_speed_right => y_x_speed_right
!                   f_y_speed_left => y_y_speed_left
!                   f_y_speed_right => y_y_speed_right
!                   f_z_speed_left => y_z_speed_left
!                   f_z_speed_right => y_z_speed_right
!                   f_density_left => y_density_left
!                   f_density_right => y_density_right
!                   f_pressure_left => y_pressure_left
!                   f_pressure_right => y_pressure_right
!               case ('z')
!                   i_f = 0
!                   j_f = 0
!                   k_f = 1
!                   flux_p => H
!                   fA => zA
!                   nx => znx
!                   ny => zny
!                   nz => znz
!                   f_x_speed_left => z_x_speed_left
!                   f_x_speed_right => z_x_speed_right
!                   f_y_speed_left => z_y_speed_left
!                   f_y_speed_right => z_y_speed_right
!                   f_z_speed_left => z_z_speed_left
!                   f_z_speed_right => z_z_speed_right
!                   f_density_left => z_density_left
!                   f_density_right => z_density_right
!                   f_pressure_left => z_pressure_left
!                   f_pressure_right => z_pressure_right
!               case default
!                   call dmsg(5, 'van_leer', 'compute_flux', &
!                           'Direction not recognised')
!                   stop
!           end select

!           real, dimension(1:n_var) :: F_plus, F_minus
!           real :: M_perp_left, M_perp_right
!           real :: alpha_plus, alpha_minus
!           real :: beta_left, beta_right
!           real :: M_plus, M_minus
!           real :: D_plus, D_minus
!           real :: c_plus, c_minus
!           real :: scrD_plus, scrD_minus
!           real :: sound_speed_avg, face_normal_speeds
!          !real :: sound_speed_left, sound_speed_right
!           integer :: i , j, k
!           
!           do k = 1, kmx - 1 + k_f
!            do j = 1, jmx - 1 + j_f 
!             do i = 1, imx - 1 + i_f
!               sound_speed_avg = 0.5 * (sqrt(gm * x_pressure_left(i, j, k) / &
!                                           x_density_left(i, j, k) ) + &
!                                         sqrt(gm * x_pressure_right(i, j, k) / &
!                                           x_density_right(i, j, k) ) )
!               
!               ! Compute '+' direction quantities

!               face_normal_speeds = x_x_speed_left(i, j, k) * xnx(i, j, k) + &
!                                    x_y_speed_left(i, j, k) * xny(i, j, k) + &
!                                    x_z_speed_left(i, j, k) * xnz(i, j, k)
!               M_perp_left = face_normal_speeds / sound_speed_avg
!               alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
!               beta_left = -max(0, 1 - floor(abs(M_perp_left)))
!               M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
!               D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
!               c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
!                         beta_left * M_plus
!               scrD_plus = (alpha_plus * (1. + beta_left)) - &
!                       (beta_left * D_plus)

!               ! First construct the F mass flux
!               F_plus(1) = x_density_left(i, j, k) * sound_speed_avg * c_plus
!               
!               ! Is convective flux zero anywhere?
!               if (check_if_sphere_wall(i, j, k)) then
!                   F_plus(1) = 0.
!               end if


!               ! Construct other fluxes in terms of the F mass flux
!               F_plus(2) = (F_plus(1) * x_x_speed_left(i, j, k)) + &
!                           (scrD_plus * x_pressure_left(i, j, k) * xnx(i, j, k))
!               F_plus(3) = (F_plus(1) * x_y_speed_left(i, j, k)) + &
!                           (scrD_plus * x_pressure_left(i, j, k) * xny(i, j, k))
!               F_plus(4) = (F_plus(1) * x_z_speed_left(i, j, k)) + &
!                           (scrD_plus * x_pressure_left(i, j, k) * xnz(i, j, k))
!               F_plus(5) = F_plus(1) * &
!                           ((0.5 * (x_x_speed_left(i, j, k) ** 2. + &
!                                    x_y_speed_left(i, j, k) ** 2. + &
!                                    x_z_speed_left(i, j, k) ** 2.)) + &
!                           ((gm / (gm - 1.)) * x_pressure_left(i, j, k) / &
!                            x_density_left(i, j, k)))

!               ! Multiply in the face areas
!               F_plus(1) = F_plus(1) * xA(i, j, k)
!               F_plus(2) = F_plus(2) * xA(i, j, k)
!               F_plus(3) = F_plus(3) * xA(i, j, k)
!               F_plus(4) = F_plus(4) * xA(i, j, k)
!               F_plus(5) = F_plus(5) * xA(i, j, k)

!               ! Compute '-' direction quantities

!               face_normal_speeds = x_x_speed_right(i, j, k) * xnx(i, j, k) + &
!                                    x_y_speed_right(i, j, k) * xny(i, j, k) + &
!                                    x_z_speed_right(i, j, k) * xnz(i, j, k)
!               M_perp_right = face_normal_speeds / sound_speed_avg
!               alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
!               beta_right = -max(0, 1 - floor(abs(M_perp_right)))
!               M_minus = -0.25 * ((1. - M_perp_right) ** 2.)
!               D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
!               c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
!                         beta_right * M_minus
!               scrD_minus = (alpha_minus * (1. + beta_right)) - &
!                            (beta_right * D_minus)

!               ! First construct the F mass flux
!               F_minus(1) = x_density_right(i, j, k) * sound_speed_avg * c_minus
!               
!               ! Is convective flux zero anywhere?
!               if (check_if_sphere_wall(i, j, k)) then
!                   F_minus(1) = 0.
!               end if

!               ! Construct other fluxes in terms of the F mass flux
!               F_minus(2) = (F_minus(1) * x_x_speed_right(i, j, k)) + &
!                            (scrD_minus * x_pressure_right(i, j, k) * xnx(i, j, k))
!               F_minus(3) = (F_minus(1) * x_y_speed_right(i, j, k)) + &
!                            (scrD_minus * x_pressure_right(i, j, k) * xny(i, j, k))
!               F_minus(4) = (F_minus(1) * x_z_speed_right(i, j, k)) + &
!                            (scrD_minus * x_pressure_right(i, j, k) * xnz(i, j, k))
!               F_minus(5) = F_minus(1) * &
!                           ((0.5 * (x_x_speed_right(i, j, k) ** 2. + &
!                                    x_y_speed_right(i, j, k) ** 2. + &
!                                    x_z_speed_right(i, j, k) ** 2.)) + &
!                           ((gm / (gm - 1.)) * x_pressure_right(i, j, k) / &
!                            x_density_right(i, j, k)))
!        
!               ! Multiply in the face areas
!               F_minus(1) = F_minus(1) * xA(i, j, k)
!               F_minus(2) = F_minus(2) * xA(i, j, k)
!               F_minus(3) = F_minus(3) * xA(i, j, k)
!               F_minus(4) = F_minus(4) * xA(i, j, k)
!               F_minus(5) = F_minus(5) * xA(i, j, k)

!               ! Get the total flux for a face
!               F(i, j, k, :) = F_plus(:) + F_minus(:)
!             end do
!            end do
!           end do 
!       end subroutine compute_flux(f_dir)

        subroutine compute_F_flux()
            !-----------------------------------------------------------
            ! Compute the flux 'F' along the xi faces
            !-----------------------------------------------------------
            
            implicit none

            real, dimension(1:n_var) :: F_plus, F_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            integer :: i , j, k
            
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx
                sound_speed_avg = 0.5 * (sqrt(gm * x_pressure_left(i, j, k) / &
                                            x_density_left(i, j, k) ) + &
                                          sqrt(gm * x_pressure_right(i, j, k) / &
                                            x_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities

                face_normal_speeds = x_x_speed_left(i, j, k) * xnx(i, j, k) + &
                                     x_y_speed_left(i, j, k) * xny(i, j, k) + &
                                     x_z_speed_left(i, j, k) * xnz(i, j, k)
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! First construct the F mass flux
                F_plus(1) = x_density_left(i, j, k) * sound_speed_avg * c_plus
                
                ! Construct other fluxes in terms of the F mass flux
                F_plus(2) = (F_plus(1) * x_x_speed_left(i, j, k)) + &
                            (scrD_plus * x_pressure_left(i, j, k) * xnx(i, j, k))
                F_plus(3) = (F_plus(1) * x_y_speed_left(i, j, k)) + &
                            (scrD_plus * x_pressure_left(i, j, k) * xny(i, j, k))
                F_plus(4) = (F_plus(1) * x_z_speed_left(i, j, k)) + &
                            (scrD_plus * x_pressure_left(i, j, k) * xnz(i, j, k))
                F_plus(5) = F_plus(1) * &
                            ((0.5 * (x_x_speed_left(i, j, k) ** 2. + &
                                     x_y_speed_left(i, j, k) ** 2. + &
                                     x_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * x_pressure_left(i, j, k) / &
                             x_density_left(i, j, k)))

                ! Multiply in the face areas
                F_plus(1) = F_plus(1) * xA(i, j, k)
                F_plus(2) = F_plus(2) * xA(i, j, k)
                F_plus(3) = F_plus(3) * xA(i, j, k)
                F_plus(4) = F_plus(4) * xA(i, j, k)
                F_plus(5) = F_plus(5) * xA(i, j, k)

                ! Compute '-' direction quantities

                face_normal_speeds = x_x_speed_right(i, j, k) * xnx(i, j, k) + &
                                     x_y_speed_right(i, j, k) * xny(i, j, k) + &
                                     x_z_speed_right(i, j, k) * xnz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = -0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)

                ! First construct the F mass flux
                F_minus(1) = x_density_right(i, j, k) * sound_speed_avg * c_minus
                
                ! Construct other fluxes in terms of the F mass flux
                F_minus(2) = (F_minus(1) * x_x_speed_right(i, j, k)) + &
                             (scrD_minus * x_pressure_right(i, j, k) * xnx(i, j, k))
                F_minus(3) = (F_minus(1) * x_y_speed_right(i, j, k)) + &
                             (scrD_minus * x_pressure_right(i, j, k) * xny(i, j, k))
                F_minus(4) = (F_minus(1) * x_z_speed_right(i, j, k)) + &
                             (scrD_minus * x_pressure_right(i, j, k) * xnz(i, j, k))
                F_minus(5) = F_minus(1) * &
                            ((0.5 * (x_x_speed_right(i, j, k) ** 2. + &
                                     x_y_speed_right(i, j, k) ** 2. + &
                                     x_z_speed_right(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * x_pressure_right(i, j, k) / &
                             x_density_right(i, j, k)))
         
                ! Multiply in the face areas
                F_minus(1) = F_minus(1) * xA(i, j, k)
                F_minus(2) = F_minus(2) * xA(i, j, k)
                F_minus(3) = F_minus(3) * xA(i, j, k)
                F_minus(4) = F_minus(4) * xA(i, j, k)
                F_minus(5) = F_minus(5) * xA(i, j, k)

                ! Get the total flux for a face
                F(i, j, k, :) = F_plus(:) + F_minus(:)
              end do
             end do
            end do 

        end subroutine compute_F_flux

        subroutine compute_G_flux()
            !-----------------------------------------------------------
            ! Compute the flux 'G' along the eta faces
            !-----------------------------------------------------------
            
            implicit none

            real, dimension(1:n_var) :: G_plus, G_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            integer :: i , j, k
            
            do k = 1, kmx - 1
             do j = 1, jmx
              do i = 1, imx - 1
                sound_speed_avg = 0.5 * (sqrt(gm * y_pressure_left(i, j, k) / &
                                            y_density_left(i, j, k) ) + &
                                          sqrt(gm * y_pressure_right(i, j, k) / &
                                            y_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = y_x_speed_left(i, j, k) * ynx(i, j, k) + &
                                     y_y_speed_left(i, j, k) * yny(i, j, k) + &
                                     y_z_speed_left(i, j, k) * ynz(i, j, k)
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! First construct the F mass flux
                G_plus(1) = y_density_left(i, j, k) * sound_speed_avg * c_plus
                
                ! Construct other fluxes in terms of the F mass flux
                G_plus(2) = (G_plus(1) * y_x_speed_left(i, j, k)) + &
                            (scrD_plus * y_pressure_left(i, j, k) * ynx(i, j, k))
                G_plus(3) = (G_plus(1) * y_y_speed_left(i, j, k)) + &
                            (scrD_plus * y_pressure_left(i, j, k) * yny(i, j, k))
                G_plus(4) = (G_plus(1) * y_z_speed_left(i, j, k)) + &
                            (scrD_plus * y_pressure_left(i, j, k) * ynz(i, j, k))
                G_plus(5) = G_plus(1) * &
                            ((0.5 * (y_x_speed_left(i, j, k) ** 2. + &
                                     y_y_speed_left(i, j, k) ** 2. + &
                                     y_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * y_pressure_left(i, j, k) / &
                             y_density_left(i, j, k)))
         
                ! Multiply in the face areas
                G_plus(1) = G_plus(1) * yA(i, j, k)
                G_plus(2) = G_plus(2) * yA(i, j, k)
                G_plus(3) = G_plus(3) * yA(i, j, k)
                G_plus(4) = G_plus(4) * yA(i, j, k)
                G_plus(5) = G_plus(5) * yA(i, j, k)

                ! Compute '-' direction quantities

                face_normal_speeds = y_x_speed_right(i, j, k) * ynx(i, j, k) + &
                                     y_y_speed_right(i, j, k) * yny(i, j, k) + &
                                     y_z_speed_right(i, j, k) * ynz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = - 0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)

                ! First construct the G mass flux
                G_minus(1) = y_density_right(i, j, k) * sound_speed_avg * c_minus
                
                ! Construct other fluxes in terms of the G mass flux
                G_minus(2) = (G_minus(1) * y_x_speed_right(i, j, k)) + &
                             (scrD_minus * y_pressure_right(i, j, k) * ynx(i, j, k))
                G_minus(3) = (G_minus(1) * y_y_speed_right(i, j, k)) + &
                             (scrD_minus * y_pressure_right(i, j, k) * yny(i, j, k))
                G_minus(4) = (G_minus(1) * y_z_speed_right(i, j, k)) + &
                             (scrD_minus * y_pressure_right(i, j, k) * ynz(i, j, k))
                G_minus(5) = G_minus(1) * &
                             ((0.5 * (y_x_speed_right(i, j, k) ** 2. + &
                                      y_y_speed_right(i, j, k) ** 2. + &
                                      y_z_speed_right(i, j, k) ** 2.)) + &
                             ((gm / (gm - 1.)) * y_pressure_right(i, j, k) / &
                               y_density_right(i, j, k)))
         
                ! Multiply in the face areas
                G_minus(1) = G_minus(1) * yA(i, j, k)
                G_minus(2) = G_minus(2) * yA(i, j, k)
                G_minus(3) = G_minus(3) * yA(i, j, k)
                G_minus(4) = G_minus(4) * yA(i, j, k)
                G_minus(5) = G_minus(5) * yA(i, j, k)

                ! Get the total flux for a face
                G(i, j, k, :) = G_plus(:) + G_minus(:)
              end do
             end do
            end do 

          ! print *, 'Checking if G_conv flux is zero'
          ! print *, G(4, jmx, 2, 1)
          ! print *, G(20, 1, 1, 1)

        end subroutine compute_G_flux

        subroutine compute_H_flux()
            !-----------------------------------------------------------
            ! Compute the flux 'G' along the zeta faces
            !-----------------------------------------------------------
            
            implicit none

            real, dimension(1:n_var) :: H_plus, H_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            integer :: i , j, k
            
            do k = 1, kmx
             do j = 1, jmx - 1
              do i = 1, imx - 1
                sound_speed_avg = 0.5 * (sqrt(gm * z_pressure_left(i, j, k) / &
                                            z_density_left(i, j, k) ) + &
                                          sqrt(gm * z_pressure_right(i, j, k) / &
                                            z_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = z_x_speed_left(i, j, k) * znx(i, j, k) + &
                                     z_y_speed_left(i, j, k) * zny(i, j, k) + &
                                     z_z_speed_left(i, j, k) * znz(i, j, k)
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! First construct the H mass flux
                H_plus(1) = z_density_left(i, j, k) * sound_speed_avg * c_plus
                
                ! Construct other fluxes in terms of the H mass flux
                H_plus(2) = (H_plus(1) * z_x_speed_left(i, j, k)) + &
                            (scrD_plus * z_pressure_left(i, j, k) * znx(i, j, k))
                H_plus(3) = (H_plus(1) * z_y_speed_left(i, j, k)) + &
                            (scrD_plus * z_pressure_left(i, j, k) * zny(i, j, k))
                H_plus(4) = (H_plus(1) * z_z_speed_left(i, j, k)) + &
                            (scrD_plus * z_pressure_left(i, j, k) * znz(i, j, k))
                H_plus(5) = H_plus(1) * &
                            ((0.5 * (z_x_speed_left(i, j, k) ** 2. + & 
                                     z_y_speed_left(i, j, k) ** 2. + &
                                     z_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * z_pressure_left(i, j, k) / &
                             z_density_left(i, j, k)))
         
                ! Multiply in the face areas
                H_plus(1) = H_plus(1) * zA(i, j, k)
                H_plus(2) = H_plus(2) * zA(i, j, k)
                H_plus(3) = H_plus(3) * zA(i, j, k)
                H_plus(4) = H_plus(4) * zA(i, j, k)
                H_plus(5) = H_plus(5) * zA(i, j, k)

                ! Compute '-' direction quantities

                face_normal_speeds = z_x_speed_right(i, j, k) * znx(i, j, k) + &
                                     z_y_speed_right(i, j, k) * zny(i, j, k) + &
                                     z_z_speed_right(i, j, k) * znz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = - 0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)

                ! First construct the F mass flux
                H_minus(1) = z_density_right(i, j, k) * sound_speed_avg * c_minus
                
                ! Construct other fluxes in terms of the F mass flux
                H_minus(2) = (H_minus(1) * z_x_speed_right(i, j, k)) + &
                             (scrD_minus * z_pressure_right(i, j, k) * znx(i, j, k))
                H_minus(3) = (H_minus(1) * z_y_speed_right(i, j, k)) + &
                             (scrD_minus * z_pressure_right(i, j, k) * zny(i, j, k))
                H_minus(4) = (H_minus(1) * z_z_speed_right(i, j, k)) + &
                             (scrD_minus * z_pressure_right(i, j, k) * znz(i, j, k))
                H_minus(5) = H_minus(1) * &
                            ((0.5 * (z_x_speed_right(i, j, k) ** 2. + &
                                     z_y_speed_right(i, j, k) ** 2. + &
                                     z_z_speed_right(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * z_pressure_right(i, j, k) / &
                             z_density_right(i, j, k)))
         
                ! Multiply in the face areas
                H_minus(1) = H_minus(1) * zA(i, j, k)
                H_minus(2) = H_minus(2) * zA(i, j, k)
                H_minus(3) = H_minus(3) * zA(i, j, k)
                H_minus(4) = H_minus(4) * zA(i, j, k)
                H_minus(5) = H_minus(5) * zA(i, j, k)

                ! Get the total flux for a face
                H(i, j, k, :) = H_plus(:) + H_minus(:)
              end do
             end do
            end do 
            
        end subroutine compute_H_flux

        subroutine compute_fluxes()
            
            implicit none
            
            call dmsg(1, 'van_leer', 'compute_fluxes')

            call compute_F_flux()
            if (any(isnan(F))) then
                call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_G_flux()
            if (any(isnan(G))) then 
                call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            call compute_H_flux()
            if (any(isnan(H))) then
                call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        function get_residue() result(residue)
            !-----------------------------------------------------------
            ! Compute the residue using the Van-Leer scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l
            real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, 1:n_var) :: residue

            call dmsg(1, 'van_leer', 'compute_residue')

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
        
        end function get_residue

end module van_leer
