    !< Flux splitting scheme: AUSM+
module ausmP
    !<
    !< Reference: Liou, M. S., “A sequel to AUSM: AUSM+,” 
    !< Journal of Computational Physics, vol. 129, pp. 364–382, 1996
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

            call dmsg(1, 'AUSM+', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - AUSM+.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - AUSM+.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - AUSM+.')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - AUSM+.')

        end subroutine setup_scheme

        subroutine destroy_scheme()
          !< Deallocate memory

            implicit none

            call dmsg(1, 'AUSM+', 'destroy_scheme')
            
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
            !< Integer for DO loop
            integer :: i_f, j_f, k_f 
            !< Flags to determine face direction
            real, dimension(:, :, :), pointer :: fA, nx, ny, nz
            !< Pointer to the face area and normal
            real, dimension(:,:,:,:), pointer :: f_qp_left, f_qp_right
            real, dimension(1:n_var) :: F_plus
            !< Right flux through the face
            real, dimension(1:n_var) ::F_minus
            !< Left flux through  the face
            real :: pbar
            real :: mass
            real :: HL, HR
            !< Enthalpy
            real :: uL, uR
            !< X-component of velocity
            real :: vL, vR
            !< Y-component of velocity
            real :: wL, wR
            !< Z-component of velocity
            real :: pL, pR
            !< Pressure
            real :: rL, rR
            !< Density
            real :: cL, cR
            !< Speed sound left/right
            real :: C
            !< Speed of sound at face
            real :: ML, MR
            !< Mach number left/right
            real :: VnL, VnR
            !< Face normal velocity left/right
            real :: betaL, betaR
            real :: alphaL, alphaR
            real :: FmL, FmR
            real :: Mface
            real :: Cs


            call dmsg(1, 'AUSM+', 'compute_flux '//trim(f_dir))
            
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
                    call dmsg(5, 'AUSM+', 'compute_flux', &
                            'Direction not recognised')
                    stop
            end select
            

            do k = 1, kmx - 1 + k_f
             do j = 1, jmx - 1 + j_f 
              do i = 1, imx - 1 + i_f

                ! -- primitve face state assignment --
                ! ---- left face quantities ----
                rL = f_qp_left(i,j,k,1)
                uL = f_qp_left(i,j,k,2)
                vL = f_qp_left(i,j,k,3)
                wL = f_qp_left(i,j,k,4)
                pL = f_qp_left(i,j,k,5)

                ! ---- right face quantities ----
                rR = f_qp_right(i,j,k,1)
                uR = f_qp_right(i,j,k,2)
                vR = f_qp_right(i,j,k,3)
                wR = f_qp_right(i,j,k,4)
                pR = f_qp_right(i,j,k,5)

                !-- calculated quntaties --
                ! ---- total enthalpy ----
                HL = (0.5*(uL*uL + vL*vL + wL*wL)) + ((gm/(gm - 1.))*pL/rL)
                HR = (0.5*(uR*uR + vR*vR + wR*wR)) + ((gm/(gm - 1.))*pR/rR)

                ! ---- face normal velocity ----
                VnL = uL*nx(i, j, k) + vL*ny(i, j, k) + wL*nz(i, j, k)
                VnR = uR*nx(i, j, k) + vR*ny(i, j, k) + wR*nz(i, j, k)

                ! ---- speed of sound ----
                cs = sqrt(2.0*(gm-1.0)*(0.5*(HL + HR))/(gm+1.0))
                cL = cs*cs/(max(cs, abs(VnL)))
                cR = cs*cs/(max(cs, abs(VnR)))
                C  = min(cL, CR)

                ! ---- Mach at face ----
                ML = VnL/C
                MR = VnR/C

                ! ---- switch for supersonic flow ----
                alphaL= max(0, 1-floor(abs(ML)))
                alphaR= max(0, 1-floor(abs(MR)))

                
                ! Compute '+' direction quantities
                FmL   = (0.5*(1.0+sign(1.0,ML))*(1.0-alphaL)*ML) + alphaL*0.25*((1.0+ML)**2)
                betaL = (0.5*(1.0+sign(1.0,ML))*(1.0-alphaL))    + alphaL*0.25*((1.0+ML)**2) * (2.0 - ML)

                ! Compute '-' direction quantities
                FmR   = (0.5*(1.0-sign(1.0,MR))*(1.0-alphaR)*MR) - alphaR*0.25*((1.0-MR)**2)
                betaR = (0.5*(1.0-sign(1.0,MR))*(1.0-alphaR))    + alphaR*0.25*((1.0-MR)**2)*(2.0 + MR)

                !AUSM+modification
                ! Compute '+' direction quantities
                FmL   = FmL   + alphaL*0.1250*((ML**2-1.0)**2)
                betaL = betaL + alphaL*0.1875*((ML**2-1.0)**2)*ML

                ! Compute '-' direction quantities
                FmR   = FmR   - alphaR*0.1250*((MR**2-1.0)**2)
                betaR = betaR - alphaR*0.1875*((MR**2-1.0)**2)*MR

                
                ! mass coefficient
                Mface = FmL + FmR
                ! -- Pressure coeffient--
                pbar = betaL*pL + betaR*pR


                ! -- mass --
                if(Mface>0.0)then
                    mass = Mface*C*rL
                else
                    mass = Mface*c*rR
                end if

                mass = mass *(i_f*make_F_flux_zero(i) &
                            + j_f*make_G_flux_zero(j) &
                            + k_f*make_H_flux_zero(k))


                ! F plus mass flux
                ! Construct other fluxes in terms of the F mass flux
                F_plus(1) = 0.5*(mass + abs(mass))
                F_plus(2) = (F_plus(1) * uL)
                F_plus(3) = (F_plus(1) * vL)
                F_plus(4) = (F_plus(1) * wL)
                F_plus(5) = (F_plus(1) * HL)
                
                ! F minus mass flux
                ! Construct other fluxes in terms of the F mass flux
                F_minus(1) = 0.5*(mass - abs(mass))
                F_minus(2) = (F_minus(1) * uR)
                F_minus(3) = (F_minus(1) * vR)
                F_minus(4) = (F_minus(1) * wR)
                F_minus(5) = (F_minus(1) * HR)

                !! -- Turbulence variables mass flux --
                if(n_var>5) then
                  F_plus(6:)  = F_Plus(1)  * f_qp_left(i,j,k,6:)
                  F_minus(6:) = F_minus(1) * f_qp_right(i,j,k,6:)
                end if

                ! total flux
                flux_p(i, j, k, :) = F_plus(:) + F_minus(:)

                ! Get the total flux for a face
                ! -- Pressure flux addition --
                flux_p(i, j, K, 2) = flux_p(i, j, k, 2) + (pbar * nx(i, j, k))
                flux_p(i, j, K, 3) = flux_p(i, j, k, 3) + (pbar * ny(i, j, k))
                flux_p(i, j, K, 4) = flux_p(i, j, k, 4) + (pbar * nz(i, j, k))

                flux_p(i, j, k, :) = flux_p(i, j, k, :)*fA(i,j,k)
              end do
             end do
            end do 

        end subroutine compute_flux

        subroutine compute_fluxes()
          !< Call to compute fluxes throught faces in each direction
            
            implicit none
            
            call dmsg(1, 'AUSM+', 'compute_fluxes')

            call compute_flux('x')
            if (any(isnan(F))) then
                call dmsg(5, 'AUSM+', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_flux('y')
            if (any(isnan(G))) then 
                call dmsg(5, 'AUSM+', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            if(kmx==2) then
              H = 0.
            else
              call compute_flux('z')
            end if
            if (any(isnan(H))) then
                call dmsg(5, 'AUSM+', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        subroutine get_residue()
            !< Compute the residue for the AUSM+ scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l

            call dmsg(1, 'AUSM+', 'compute_residue')

            do l = 1, n_var
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 1, imx - 1
               residue(i, j, k, l) = (F(i+1, j, k, l) - F(i, j, k, l)) &
                                   + (G(i, j+1, k, l) - G(i, j, k, l)) &
                                   + (H(i, j, k+1, l) - H(i, j, k, l))
               end do
              end do
             end do
            end do
        
        end subroutine get_residue

end module ausmP
