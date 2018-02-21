module ausmUP
    !-------------------------------------------------------------------
    ! The ausmUP scheme is a type of flux-splitting scheme
    !-------------------------------------------------------------------
    
    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx

    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z
    use global_vars, only : xA, yA, zA    !face area

    use global_vars, only : MInf
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

            call dmsg(1, 'ausmUP', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - ausmUP.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - ausmUP.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - ausmUP.')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - ausmUP.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ausmUP', 'destroy_scheme')
            
            call dealloc(F)
            call dealloc(G)
            call dealloc(H)

        end subroutine destroy_scheme

        subroutine compute_flux(f_dir)

            implicit none
            character, intent(in) :: f_dir
            integer :: i, j, k 
            integer :: i_f, j_f, k_f ! Flags to determine face direction
            real, dimension(:, :, :), pointer :: fA, nx, ny, nz
            real, dimension(:,:,:,:), pointer :: f_qp_left, f_qp_right
            real, dimension(1:n_var) :: F_plus, F_minus
            real :: xi
            real :: vnabs
            real :: delp
            real :: delrho
            real :: fnG
            real :: pbar
            real :: Mcap
            real :: vtface
            real :: mass
            real :: HL, HR !enthalpy
            real :: uL, uR
            real :: vL, vR
            real :: wL, wR
            real :: pL, pR
            real :: rL, rR
            real :: cL, cR
            real :: C, Cstar2
            real :: ML, MR
            real :: VnL, VnR
            real :: betaL, betaR
            real :: alphaL, alphaR
            real :: VnabsL, VnabsR
            real :: fna
            real :: Mo
            real :: MF
            real :: MB2
            real :: Pu
            real :: fnML, fnMR
            real :: alfa
            real, parameter :: Ku = 0.75
            real, parameter :: Kp = 0.25
            real, parameter :: sigma = 1.0
            real :: Mp

            call dmsg(1, 'ausmUP', 'compute_flux '//trim(f_dir))
            
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
                    call dmsg(5, 'ausmUP', 'compute_flux', &
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
                Cstar2 = 1.0*(gm-1.0)*(HL + HR)/(gm+1.0)
                cL = Cstar2/(max(sqrt(Cstar2),VnL))
                cR = Cstar2/(max(sqrt(Cstar2),-VnR))
                C  = min(cL,cR)

                ! ---- Mach at face ----
                ML = VnL/C
                MR = VnR/C
                MB2= 0.5*((VnL*VnL) + (VnR*VnR))/(C*C)


                ! -- function at face --
                Mo = sqrt(min(1.0, max(MB2, MInf*MInf)))
                fna = Mo*(2.0 - Mo)
                alfa = 3.0*(-4.0 + 5.0*fna*fna)/16.0

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
                betaL = (1.0-alphaL)*0.5*(1.0+sign(1.0,ML)) + (alphaL)*0.25*(2.0-ML)*((ML+1.0)**2) &
                                                            + (alphaL)*alfa*ML*((ML*ML)-1.0)**2
                betaR = (1.0-alphaR)*0.5*(1.0-sign(1.0,MR)) + (alphaR)*0.25*(2.0+MR)*((MR-1.0)**2) &
                                                            - (alphaR)*alfa*MR*((MR*MR)-1.0)**2
                Pu    = -Ku*betaL*betaR*(rL+rR)*fna*C*(VnR-VnL)
                
                ! -- mass flux factor --
                fnML =((1.0-alphaL)*0.5*(1.0+sign(1.0,ML))*ML) + (alphaL)*((0.25*(ML+1.0)**2)+(0.125*((ML*ML)-1.0)**2))
                fnMR =((1.0-alphaR)*0.5*(1.0-sign(1.0,MR))*MR) - (alphaR)*((0.25*(MR-1.0)**2)-(0.125*((MR*MR)-1.0)**2))
                Mp   = -2.0*Kp*max(1.0-sigma*MB2,0.0)*(pR-pL)/(fna*(rL+rR)*C*C)

                ! -- Pressure --
                pbar = betaL*pL + betaR*pR + Pu

                ! -- mass --
                MF = fnML + fnMR + Mp
                mass = 0.5*C*(rL*(MF+abs(MF)) + rR*(MF-abs(MF))) 
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

                if(i_f==1 .and. (i==1 .or. i==2) .and. j==1 .and. k==1)then
                  print*, process_id, flux_p(i,j,k,:)
                end if

              end do
             end do
            end do 

        end subroutine compute_flux

        subroutine compute_fluxes()
            
            implicit none
            
            call dmsg(1, 'ausmUP', 'compute_fluxes')

            call compute_flux('x')
            if (any(isnan(F))) then
                call dmsg(5, 'ausmUP', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_flux('y')
            if (any(isnan(G))) then 
                call dmsg(5, 'ausmUP', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            if(kmx==2) then
              H = 0.
            else
              call compute_flux('z')
            end if
            if (any(isnan(H))) then
                call dmsg(5, 'ausmUP', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        subroutine get_residue()
            !-----------------------------------------------------------
            ! Compute the residue for the ausmUP scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l

            call dmsg(1, 'ausmUP', 'compute_residue')

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

end module ausmUP
