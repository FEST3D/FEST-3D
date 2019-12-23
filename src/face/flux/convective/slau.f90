    !< Flux splitting scheme: SLAU
module slau
    !< Shima, E., and Kitamura, K., “Parameter-Free Simple
    !< Low-Dissipation AUSM-Family Scheme for All Speeds,” 
    !< AIAA Journal, vol. 49, pp. 1693–1709, 2011
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
          !< flux through the input direction, :x,y, or z
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
            !< Integer for DO loop
            integer :: i_f, j_f, k_f 
            !< Flags to determine face direction
            real, dimension(:,:,:,:), pointer :: f_qp_left, f_qp_right
            real, dimension(1:dims%n_var) :: F_plus
            !< Right flux through the face
            real, dimension(1:dims%n_var) ::F_minus
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

            DebugCall('compute_flux '//trim(f_dir))
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
                HL = (0.5*(uL*uL + vL*vL + wL*wL)) + ((flow%gm/(flow%gm - 1.))*pL/rL)
                HR = (0.5*(uR*uR + vR*vR + wR*wR)) + ((flow%gm/(flow%gm - 1.))*pR/rR)

                ! ---- speed of sound ----
                cL = sqrt(flow%gm*pL/rL)
                cR = sqrt(flow%gm*pR/rR)
                C  = 0.5*(cL + cR)

                ! ---- delta quantities ----
                delp   = pR-pL!pL - pR
                delrho = rR-rL!rL - rR

                ! ---- face normal velocity ----
                VnL = uL*faces(i, j, k)%nx + vL*faces(i, j, k)%ny + wL*faces(i, j, k)%nz
                VnR = uR*faces(i, j, k)%nx + vR*faces(i, j, k)%ny + wR*faces(i, j, k)%nz

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
                mass = mass *(i_f*bc%make_F_flux_zero(i) &
                            + j_f*bc%make_G_flux_zero(j) &
                            + k_f*bc%make_H_flux_zero(k))


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
                if(dims%n_var>5) then
                  F_plus(6:)  = F_Plus(1)  * f_qp_left(i,j,k,6:)
                  F_minus(6:) = F_minus(1) * f_qp_right(i,j,k,6:)
                end if

                ! Get the total flux for a face
                Flux(i, j, k, :) = F_plus(:) + F_minus(:)

                ! -- Pressure flux addition --
                Flux(i, j, K, 2) = Flux(i, j, k, 2) + (pbar * faces(i, j, k)%nx)
                Flux(i, j, K, 3) = Flux(i, j, k, 3) + (pbar * faces(i, j, k)%ny)
                Flux(i, j, K, 4) = Flux(i, j, k, 4) + (pbar * faces(i, j, k)%nz)

                Flux(i, j, k, :) = Flux(i, j, k, :)*faces(i,j,k)%A

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
            call compute_flux(F, 'x', Ifaces, flags, flow, bc,dims)
            if (any(isnan(F))) then
              Fatal_error
            end if    

            flags=(/0,1,0/)
            call compute_flux(G, 'y', Jfaces, flags, flow, bc,dims)
            if (any(isnan(G))) then 
              Fatal_error
            end if    
            
            if(dims%kmx==2) then
              H = 0.
            else
              flags=(/0,0,1/)
              call compute_flux(H, 'z', Kfaces, flags, flow, bc,dims)
            end if
            if (any(isnan(H))) then
              Fatal_error
            end if

        end subroutine compute_fluxes



end module slau
