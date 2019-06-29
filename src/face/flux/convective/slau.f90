    !< Flux splitting scheme: SLAU
module slau
    !< Shima, E., and Kitamura, K., “Parameter-Free Simple
    !< Low-Dissipation AUSM-Family Scheme for All Speeds,” 
    !< AIAA Journal, vol. 49, pp. 1693–1709, 2011
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

            call dmsg(1, 'slau', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - slau.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - slau.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - slau.')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - slau.')

        end subroutine setup_scheme

        subroutine destroy_scheme()
          !< Deallocate memory

            implicit none

            call dmsg(1, 'slau', 'destroy_scheme')
            
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
            real :: xi
            real :: vnabs
            real :: delp
            real :: delrho
            real :: fnG
            real :: pbar
            real :: Mcap
            real :: vtface
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
            real :: VnabsL, VnabsR

            call dmsg(1, 'slau', 'compute_flux '//trim(f_dir))
            
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
                    call dmsg(5, 'slau', 'compute_flux', &
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

                ! ---- speed of sound ----
                cL = sqrt(gm*pL/rL)
                cR = sqrt(gm*pR/rR)
                C  = 0.5*(cL + cR)

                ! ---- delta quantities ----
                delp   = pR-pL!pL - pR
                delrho = rR-rL!rL - rR

                ! ---- face normal velocity ----
                VnL = uL*nx(i, j, k) + vL*ny(i, j, k) + wL*nz(i, j, k)
                VnR = uR*nx(i, j, k) + vR*ny(i, j, k) + wR*nz(i, j, k)

                ! ---- Mach at face ----
                ML = VnL/C
                MR = VnR/C

                ! ---- switch for supersonic flow ----
                alphaL= max(0.0, 1.0-floor(abs(ML)))
                alphaR= max(0.0, 1.0-floor(abs(MR)))
                !Above two line of code is eqvivalent to following code
                    !if(abs(ML)>=1.0) then
                    !  alphaL = 0.0
                    !else
                    !  alphaL=1.0
                    !end if
                    !if(abs(MR)>=1.0) then
                    !  alphaR=0.0
                    !else
                    !  alphaR=1.0
                    !end if

                ! -- pressure factor --
                betaL = (1.0-alphaL)*0.5*(1.0+sign(1.0,ML)) + (alphaL)*0.25*(2.0-ML)*((ML+1.0)**2)
                betaR = (1.0-alphaR)*0.5*(1.0-sign(1.0,MR)) + (alphaR)*0.25*(2.0+MR)*((MR-1.0)**2)
                
                ! -- xi calculation --
                vtface = sqrt(0.5*((uL*uL) + (vL*vL) + (wL*wL) + (uR*uR) + (vR*vR) + (wR*wR)))
                Mcap   = min(1.0, vtface/C)
                Xi     = (1.0 - Mcap)**2

                ! -- |Vn| --
                Vnabs = (rL *abs(VnL) + rR*abs(VnR))/(rL + rR)
                
                ! -- function G --
                fnG = -1.0*max(min(ML,0.0),-1.0)*min(max(MR,0.0),1.0)

                ! -- Pressure --
                pbar = 0.5*((pL+pR) + (betaL-betaR)*(pL-pR) + (1.0-xi)*(betaL+betaR-1.0)*(pL+pR))

                ! -- mass --
                !mass = 0.5*((rL*VnL + rR*VnR - Vnabs*delrho)*(1.0-fnG) - (Xi*delp/C))
                VnabsL = (1.0 - fnG)*Vnabs + fnG*abs(VnL)
                VnabsR = (1.0 - fnG)*Vnabs + fnG*abs(VnR)
                mass = 0.5*((rL*(VnL+VnabsL) + rR*(VnR-VnabsR)) - (Xi*delp/C))
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

                ! -- Turbulence variables mass flux --
                if(n_var>5) then
                  F_plus(6:)  = F_Plus(1)  * f_qp_left(i,j,k,6:)
                  F_minus(6:) = F_minus(1) * f_qp_right(i,j,k,6:)
                end if

                ! Get the total flux for a face
                flux_p(i, j, k, :) = F_plus(:) + F_minus(:)

                ! -- Pressure flux addition --
                flux_p(i, j, K, 2) = flux_p(i, j, k, 2) + (pbar * nx(i, j, k))
                flux_p(i, j, K, 3) = flux_p(i, j, k, 3) + (pbar * ny(i, j, k))
                flux_p(i, j, K, 4) = flux_p(i, j, k, 4) + (pbar * nz(i, j, k))

                flux_P(i, j, k, :) = flux_p(i, j, k, :) * fA(i, j, k)

              end do
             end do
            end do 

        end subroutine compute_flux

        subroutine compute_fluxes()
          !< Call to compute fluxes throught faces in each direction

            
            implicit none
            
            call dmsg(1, 'slau', 'compute_fluxes')

            call compute_flux('x')
            if (any(isnan(F))) then
                call dmsg(5, 'slau', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_flux('y')
            if (any(isnan(G))) then 
                call dmsg(5, 'slau', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            if(kmx==2) then
              H = 0.
            else
              call compute_flux('z')
            end if
            if (any(isnan(H))) then
                call dmsg(5, 'slau', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        subroutine get_residue()
            !< Compute the residue for the slau scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l

            call dmsg(1, 'slau', 'compute_residue')

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

end module slau
