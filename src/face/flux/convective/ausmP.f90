    !< Flux splitting scheme: AUSM+
module ausmP
    !<
    !< Reference: Liou, M. S., “A sequel to AUSM: AUSM+,” 
    !< Journal of Computational Physics, vol. 129, pp. 364–382, 1996
    !-------------------------------------------------------------------
#include "../../../debug.h"
#include "../../../error.h"    
    use vartypes
!    use global_vars, only : xnx, xny, xnz !face unit normal x
!    use global_vars, only : ynx, yny, ynz !face unit normal y
!    use global_vars, only : znx, zny, znz !face unit normal z
!    use global_vars, only : xA, yA, zA    !face area
    use global_vars, only : process_id
    use global_vars, only : make_F_flux_zero
    use global_vars, only : make_G_flux_zero
    use global_vars, only : make_H_flux_zero

!    use utils, only: alloc
    use face_interpolant, only: x_qp_left, x_qp_right 
    use face_interpolant, only: y_qp_left, y_qp_right
    use face_interpolant, only:  z_qp_left, z_qp_right


    implicit none
    private

!    real, public, dimension(:, :, :, :), allocatable, target :: F
!    !< Store fluxes throught the I faces
!    real, public, dimension(:, :, :, :), allocatable, target :: G
!    !< Store fluxes throught the J faces
!    real, public, dimension(:, :, :, :), allocatable, target :: H
!    !< Store fluxes throught the K faces
!    real, public, dimension(:, :, :, :), allocatable, target :: residue
!    !< Store residue at each cell-center
!    real, dimension(:, :, :, :), pointer :: flux_p
!    !< Pointer/alias for the either F, G, or H

!    integer :: imx, jmx, kmx, n_var
    ! Public members
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
!
!            DebugCall('setup_scheme')
!
!            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'F - AUSM+.')
!            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'G - AUSM+.')
!            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'H - AUSM+.')
!            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'residue - AUSM+.')
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

        subroutine compute_flux(Flux, f_dir, faces, flags, flow, dims)
          !< A generalized subroutine to calculate
          !< flux through the input direction, :x,y, or z
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
            !< Integer for DO loop
            integer :: i_f, j_f, k_f 
            !< Flags to determine face direction
            !real, dimension(:, :, :), pointer :: fA, nx, ny, nz
            !< Pointer to the face area and normal
            real, dimension(:,:,:,:), pointer :: f_qp_left, f_qp_right
            real, dimension(1:dims%n_var) :: F_plus
            !< Right flux through the face
            real, dimension(1:dims%n_var) ::F_minus
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


            DebugCall('compute_flux '//trim(f_dir))
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

                ! ---- face normal velocity ----
                VnL = uL*faces(i, j, k)%nx + vL*faces(i, j, k)%ny + wL*faces(i, j, k)%nz
                VnR = uR*faces(i, j, k)%nx + vR*faces(i, j, k)%ny + wR*faces(i, j, k)%nz

                ! ---- speed of sound ----
                cs = sqrt(2.0*(flow%gm-1.0)*(0.5*(HL + HR))/(flow%gm+1.0))
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
                if(dims%n_var>5) then
                  F_plus(6:)  = F_Plus(1)  * f_qp_left(i,j,k,6:)
                  F_minus(6:) = F_minus(1) * f_qp_right(i,j,k,6:)
                end if

                ! total flux
                Flux(i, j, k, :) = F_plus(:) + F_minus(:)

                ! Get the total flux for a face
                ! -- Pressure flux addition --
                Flux(i, j, K, 2) = Flux(i, j, k, 2) + (pbar * faces(i, j, k)%nx)
                Flux(i, j, K, 3) = Flux(i, j, k, 3) + (pbar * faces(i, j, k)%ny)
                Flux(i, j, K, 4) = Flux(i, j, k, 4) + (pbar * faces(i, j, k)%nz)

                Flux(i, j, k, :) = Flux(i, j, k, :)*faces(i,j,k)%A
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

!        subroutine get_residue()
!            !< Compute the residue for the AUSM+ scheme
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
!               residue(i, j, k, l) = (F(i+1, j, k, l) - F(i, j, k, l)) &
!                                   + (G(i, j+1, k, l) - G(i, j, k, l)) &
!                                   + (H(i, j, k+1, l) - H(i, j, k, l))
!               end do
!              end do
!             end do
!            end do
!        
!        end subroutine get_residue

end module ausmP
