  !< Time integration module
module update
  !< This module march the solution is time.
  use vartypes
  use global_kkl , only : cphi1
  use global_kkl , only : cphi2
  use global_kkl , only : fphi
  use global_kkl , only : eta
  use global_kkl , only : cd1
  use global_kkl , only : cmu
  use global_sst , only : beta1
  use global_sst , only : beta2
  use global_sst , only : bstar
  use global_sst , only : sst_F1
  use global_sa  , only : cb1
  use global_sa  , only : cw1
  use global_sa  , only : cw2
  use global_sa  , only : cw3
  use global_sa  , only : cv1
  use global_sa  , only : kappa_sa
  use gradients  ,only :   gradu_x
  use gradients  ,only :   gradu_y
  use gradients  ,only :   gradu_z
  use gradients  ,only :   gradv_x
  use gradients  ,only :   gradv_y
  use gradients  ,only :   gradv_z
  use gradients  ,only :   gradw_x
  use gradients  ,only :   gradw_y
  use gradients  ,only :   gradw_z
  use wall_dist, only : dist
  use viscosity, only : mu
  use viscosity, only : mu_t

  use utils, only: alloc

  !subroutine for residual calculation
  use interface1,                      only: apply_interface
  use bc_primitive,                   only: populate_ghost_primitive
  use face_interpolant,               only: compute_face_interpolant
  use boundary_state_reconstruction,  only: reconstruct_boundary_state
  use scheme,                         only: compute_fluxes
  use gradients,                      only: evaluate_all_gradients
  use viscosity                      ,only: calculate_viscosity
  use viscous,                        only: compute_viscous_fluxes
  use scheme,                         only: compute_residue
  use source,                         only: add_source_term_residue

  use time,                           only : compute_time_step
  !--- sst implicit update ---!
  use global_sst, only : sst_F1
  use global_sst, only : sigma_k1
  use global_sst, only : sigma_k2
  use global_sst, only : sigma_w1
  use global_sst, only : sigma_w2

  use plusgs     , only : update_with_plusgs
  use plusgs     , only : setup_plusgs
  use lusgs     , only : update_with_lusgs
  use lusgs     , only : setup_lusgs
#include "debug.h"
#include "error.h"
    private

    real(wp), dimension(:,:,:,:), allocatable :: U_store
    !< Array to store the intermediate solution
    real(wp), dimension(:,:,:,:), allocatable :: R_store
    !< Array to store the intermediate Residue
    real(wp), dimension(:,:,:,:), allocatable, target :: aux
    !< Array to store some auxilary intermediate variables
    real(wp), dimension(:)      , allocatable :: u1
    !< Variable array old for each cell center
    real(wp), dimension(:)      , allocatable :: u2
    !< Variable array new for each cell center
    real(wp), dimension(:)      , allocatable :: R
    !< Residue array for each cell center
    integer :: imx, jmx, kmx, n_var

    ! Public methods
    public :: setup_update
    public :: get_next_solution

    contains


      subroutine setup_update(control, scheme,flow, dims)
        !< Allocate memory to variables required based 
        !< on the time-integration scheme.
        implicit none
        type(controltype), intent(in) :: control
        !< Control parameters
        type(schemetype), intent(in) :: scheme
        !< finite-volume Schemes
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx

        imx = dims%imx
        jmx = dims%jmx
        kmx = dims%kmx

        n_var = control%n_var

        call alloc(u1,1,n_var)
        call alloc(u2,1,n_var)
        call alloc(R ,1,n_var)
        call alloc(aux,-2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)

        select case (scheme%time_step_accuracy)
          case ("none")
            ! Do nothing
            continue
          case ("RK2", "RK4")
            call alloc(U_store,-2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)
            call alloc(R_store, 1,imx-1, 1,jmx-1, 1,kmx-1,1,n_var)
          case ("TVDRK2", "TVDRK3")
            call alloc(U_store,-2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)
          case ("implicit")
            call setup_lusgs(control, scheme, flow, dims)
          case ("plusgs")
            call setup_plusgs(control, scheme, flow, dims)
          case default
            Fatal_error
        end select

      end subroutine setup_update


    subroutine get_next_solution(qp, Temp, residue, delta_t, cells, F,G,H, Ifaces, Jfaces, Kfaces, control, scheme, flow, bc, dims)
        !< Get solution at next time-step using scheme
        !< given in the input file.
        implicit none
        type(controltype), intent(in) :: control
        !< Control parameters
        type(schemetype), intent(in) :: scheme
        !< finite-volume Schemes
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(boundarytype), intent(in) :: bc
        !< boundary conditions and fixed values
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx
        real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var), intent(inout):: qp
        !< Store primitive variable at cell center
        real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(inout):: Temp
        !< Store Temperature variable at cell center
        real(wp), dimension(:, :, :, :), intent(inout)  :: residue
        !< Store residue at each cell-center
        real(wp), dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(inout)  :: delta_t
        !< Local time increment value at each cell center
        !< Store residue at each cell-center
        real(wp), dimension(:, :, :, :), intent(inout) :: F
        !< Store fluxes throught the I faces
        real(wp), dimension(:, :, :, :), intent(inout) :: G
        !< Store fluxes throught the J faces
        real(wp), dimension(:, :, :, :), intent(inout) :: H
        !< Store fluxes throught the K faces
        type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
        !< Input cell quantities: volume
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Store face quantites for I faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Store face quantites for J faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Store face quantites for K faces 
        real(wp) :: CFL 
        CFL = control%CFL
        !finding the updated Temperature using ideal gas law
        !T=P/(R_gas*Rho)
        Temp = qp(:,:,:,5)/(flow%R_gas*qp(:,:,:,1) )
        select case (trim(scheme%time_step_accuracy))
            case ("none")
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1. ,1., .FALSE.) 
            case ("RK4")
              R_store=0.
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 0.5  , 1., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 0.5  , 2., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1.0  , 2., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1./6., 1., .TRUE. , R_store, U_store) 
            case("RK2")
              R_store=0.
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 0.5  , 1., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 0.5  , 1., .TRUE., R_store, U_store) 
            case ("TVDRK3")
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1.0  , 1.) 
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1.0  , 1.) 
              qp = 0.75*U_store + 0.25*qp
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1.0  , 1.) 
              qp = (1./3.)*U_store + (2./3.)*qp
            case ("TVDRK2")
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1.0  , 1.) 
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call update_with(qp, residue, delta_t, cells, scheme, flow, "conservative", 1.0  , 1.) 
              qp = 0.5*U_store + 0.5*qp
            case ("implicit")
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with_lusgs(qp,residue, delta_t, cells,Ifaces,Jfaces,Kfaces, scheme, dims)
            case ("plusgs")
              call get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
              call compute_time_step(qp, delta_t, CFL, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with_plusgs(qp,delta_t, cells,Ifaces,Jfaces,Kfaces, residue, scheme, dims)
            case default
              Fatal_error
        end select
      end subroutine get_next_solution

      subroutine update_with(qp, residue, delta_t, cells, scheme, flow, type, time_factor, store_factor, use, Rn, un)
        !< A generalized scheme to updat the solution explicitly using
        !< any RK method and even first order euler explicit.
        implicit none
        real(wp), dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(inout), target:: qp
        !< Store primitive variable at cell center
        type(celltype), dimension(-2:imx+2,-2:jmx+2,-2:kmx+2), intent(in) :: cells
        !< Cell center quantities: volume
        real(wp), dimension(1:imx-1,1:jmx-1,1:kmx-1), intent(in) :: delta_t
        !< Local time increment value at each cell center
        real(wp), dimension(:, :, :, :), intent(in)  :: residue
        !< Store residue at each cell-center
        type(schemetype), intent(in) :: scheme
        !< finite-volume Schemes
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        character(len=*), intent(in) :: type
        real(wp), intent(in), optional :: time_factor ! time factor
        real(wp), intent(in), optional :: store_factor
        logical, intent(in), optional :: use
        real(wp), dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in), optional, target :: un
        real(wp), dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(inout), optional :: Rn
        real(wp)               :: TF = 1.0 !time factor
        real(wp)               :: SF = 1.0!store factor
        Logical               :: TU = .FALSE. !to use or not
        real(wp), dimension(:,:,:,:), pointer :: Quse
        integer :: i,j,k
        real(wp) :: KE=0.
        real(wp) :: beta

        !sa variables
        real(wp) :: vort
        real(wp) :: fv1
        real(wp) :: fv2
        real(wp) :: fw
        real(wp) :: g
        real(wp) :: scap
        real(wp) :: rsa
        real(wp) :: kd2
        real(wp) :: xi 
        real(wp) :: mass_residue
        real(wp) :: x_mom_residue, y_mom_residue, z_mom_residue
        real(wp) :: energy_residue
        real(wp) :: TKE_residue, Omega_residue, kl_residue

        if(present(time_factor)) TF=time_factor
        if(present(store_factor)) SF=store_factor
        if(present(use)) TU=use
        !check if user want to update from particular solution
        if(present(un))then
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
        else
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
        end if

        select case(type)
          case('primitive')

            !update primitive variable
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1

                  mass_residue = residue(i,j,k,1)
                  x_mom_residue = residue(i,j,k,2)
                  y_mom_residue = residue(i,j,k,3)
                  z_mom_residue = residue(i,j,k,4)
                  energy_residue = residue(i,j,k,5)
            
                  u1(1:n_var) = Quse(i,j,k,1:n_var)
            
                  ! finding primitive residue
                  R(1) = mass_residue
                  R(2) = -1*(u1(2)/u1(1))*mass_residue + x_mom_residue/u1(1)
                  R(3) = -1*(u1(3)/u1(1))*mass_residue + y_mom_residue/u1(1)
                  R(4) = -1*(u1(4)/u1(1))*mass_residue + z_mom_residue/u1(1)
                  R(5) = 0.5*(flow%gm-1.)*(sum(u1(2:4)**2)*mass_residue) &
                         -(flow%gm-1.)*u1(2)*x_mom_residue               &
                         -(flow%gm-1.)*u1(3)*y_mom_residue               &
                         -(flow%gm-1.)*u1(4)*z_mom_residue               &
                         +(flow%gm-1.)*energy_residue
            
                  select case(scheme%turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sst', 'sst2003')
                      TKE_residue = residue(i,j,k,6)
                      omega_residue = residue(i,j,k,7)
                      beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
                      R(5) = R(5) - (flow%gm-1.)*TKE_residue
                      R(6) = -(u1(6)/u1(1))*mass_residue&
                             +(1./(1.+bstar*u1(6)*delta_t(i,j,k)))*TKE_residue/u1(1)
                      R(7) = -(u1(7)/u1(1))*mass_residue&
                             +(1./(1.+2.*beta*u1(6)*delta_t(i,j,k)))*omega_residue/u1(1)
                    case('kkl')
                      TKE_residue = residue(i,j,k,6)
                      kl_residue  = residue(i,j,k,7)
                      eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
                      fphi = (1+cd1*eta)/(1+eta**4)
                      R(5) = R(5) - (flow%gm-1.)*TKE_residue
                      R(6) = -(u1(6)/u1(1))*mass_residue&
                             + (1./(1.+((2.5*((cmu**0.75)*u1(1)*(u1(6)**1.5)/max(u1(7),1.e-20))&
                             +(2*mu(i,j,k)/dist(i,j,k)**2))*delta_t(i,j,k))))*TKE_residue/u1(1)
                      R(7) = -(u1(7)/u1(1))*mass_residue&
                             +(1./(1.+(6*mu(i,j,k)*fphi/dist(i,j,k)**2)*delta_t(i,j,k)))*kl_residue/u1(1)
                    case DEFAULT
                      Fatal_error
                  end select
            
                        
                 !check if user want to store residue
                  if(present(Rn)) then
                    Rn(i,j,k,1:n_var) = Rn(i,j,k,1:n_var) + SF*R(1:n_var)
                    if(TU) R(:) = Rn(i,j,k,:)
                  end if
                 
            
                 !update
                 u2(:) = u1(:) - R(:)*(TF*delta_t(i,j,k)/cells(i,j,k)%volume)
            
                  !check solution for non pyhysical results
                  if((u2(1) < 0.) .or. (u2(5)) < 0.)then
                    Fatal_error
                  else !update
                    qp(i,j,k,1:5) = u2(1:5)
                    select case(trim(scheme%turbulence))
                     case('sst', 'sst2003', 'kkl')
                       if(u2(6)>0.) qp(i,j,k,6) = u2(6)
                       if(u2(7)>0.) qp(i,j,k,7) = u2(7)
                     case DEFAULT
                       ! do nothing
                       continue
                    end select
                  end if
                end do
              end do
            end do
            
          case('conservative')
            !include "update_conservative.inc"

            !update conservative variable
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1

                  ! getting conservative variable
                  u1(1)  = Quse(i,j,k,1)
                  u1(2:) = Quse(i,j,k,2:)*u1(1)
                  select case(scheme%turbulence)
                    case('sst', 'sst2003', 'kkl')
                      KE = 0.0!u1(6)
                    case('sa','saBC')
                      KE=0.0
                    case DEFAULT
                      KE = 0.
                  end select
                  u1(5) = (u1(5)/(flow%gm-1.) + 0.5*sum(u1(2:4)**2))/u1(1) + KE

                 ! get R
                  R(1:n_var) = residue(i,j,k,1:n_var) 
                  ! point implicit destruction term
                  select case(trim(scheme%turbulence))
                    case('none')
                      !do nothing
                      continue
                    case('sst', 'sst2003')
                      beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
                      R(6) = R(6)/(1+(beta*qp(i,j,k,7)*delta_t(i,j,k)))
                      R(7) = R(7)/(1+(2*beta*qp(i,j,k,7)*delta_t(i,j,k)))
                    case('kkl')
                      eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
                      fphi = (1+cd1*eta)/(1+eta**4)
                      R(6) = R(6)/(1.+((2.5*((cmu**0.75)*sqrt(u1(1))*(u1(6)**1.5)/max(u1(7),1.e-20))&
                             +(2*mu(i,j,k)/(dist(i,j,k)**2)))*delta_t(i,j,k)))
                      R(7) = R(7)/(1.+(6*mu(i,j,k)*fphi/(dist(i,j,k)**2))*delta_t(i,j,k))
                    case('sa', 'saBC')
                      vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                                      + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                                      + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                                       )&
                                 )
                      kd2  = (kappa_sa*dist(i,j,k))**2
                      xi   = U1(6)*qp(i,j,k,1)/mu(i,j,k)
                      fv1  = xi**3/(xi**3 + cv1**3)
                      fv2  = 1.0 - xi/(1 + xi*fv1)
                      scap = vort + U1(6)*fv2/(kd2)
                      rsa    = min(U1(6)/(Scap*kd2), 10.0)
                      g    = rsa + cw2*(rsa**6 - rsa)
                      fw   = g*( (1.0+cw3**6)/(g**6+cw3**6) )**(1.0/6.0)
                      R(6) = R(6)/(1.+((-1.0*u1(1)*cb1*scap)+(2.0*u1(1)*cw1*fw*u1(6)/(dist(i,j,k)**2)))*delta_t(i,j,k))
                    case DEFAULT
                      Fatal_error
                  end select

                 !check if user want to store residue
                 if(present(Rn)) then
                   Rn(i,j,k,1:n_var) = Rn(i,j,k,1:n_var) + SF*R(1:n_var)
                   if(TU) R(:) = Rn(i,j,k,:)
                 end if

                 !update
                 u2(1:n_var) = u1(1:n_var) - R(1:n_var)*(TF*delta_t(i,j,k)/cells(i,j,k)%volume)

                ! getting primitve variable back variable
                  u2(1)  = u2(1)
                  u2(2:) = u2(2:)/u2(1)
                  select case(scheme%turbulence)
                    case('sst', 'sst2003', 'kkl')
                      KE = 0.0!u2(6)
                    case('sa', 'saBC')
                      !u2(6) = u2(6)*u2(1)
                      KE=0.0
                    case DEFAULT
                      KE = 0.
                  end select
                  u2(5) = (flow%gm-1.)*u2(1)*(u2(5) - (0.5*sum(u2(2:4)**2)) - KE)

                  !check solution for non pyhysical results
                  if((u2(1) < 0.) .or. (u2(5)) < 0. .or. any(isnan(u2)))then
                    print*, u2(:)
                    print*, "R: ", R
                    print*, "old ", U1
                    Fatal_error
                  else !update
                    qp(i,j,k,1:5) = u2(1:5)
                    select case(trim(scheme%turbulence))
                     case('sst', 'sst2003', 'kkl')
                       if(u2(6)>=0.) then
                         qp(i,j,k,6) = u2(6)
                       else
                       !  qp(i,j,k,6) = tk_inf
                       !  qp(i,j,k,6) = (max(qp(i-1,j,k,6),0.) + max(qp(i+1,j,k,6),0.) &
                       !                +max(qp(i,j-1,k,6),0.) + max(qp(i,j+1,k,6),0.) &
                       !                )/4
                       !  qp(i,j,k,6) = 1.e-3*maxval(qp(i-1:i+1,j-1:j+1,k-1:k+1,6))
                       end if
                       if(u2(7)>=0.) then
                        qp(i,j,k,7) = u2(7)
                       else
                       !  qp(i,j,k,7) = tkl_inf
                       !  qp(i,j,k,7) = (max(qp(i-1,j,k,7),0.) + max(qp(i+1,j,k,7),0.) &
                       !                +max(qp(i,j-1,k,7),0.) + max(qp(i,j+1,k,7),0.) &
                       !                )/4
                       end if
                     case('sa', 'saBC')
                       qp(i,j,k,6) = max(u2(6), 1.e-12)
                     case DEFAULT
                       ! do nothing
                       continue
                    end select
                  end if
                  !print*, i,j, R(1:n_var)

                end do
              end do
            end do

          case DEFAULT
            Fatal_error
        end select

      end subroutine update_with


      !subroutine get_total_conservative_Residue(qp, Temp, residue, F,G,H, control, scheme, flow, dims)
      subroutine get_total_conservative_Residue(qp, Temp, cells, residue, F,G,H, Ifaces,Jfaces,Kfaces, control, scheme, flow, bc, dims)
        !< For each iteration it apply boundary conditions,
        !< use higher order method to reconstruct state at
        !< face, evalute fluxes at each face, calculate 
        !< inviscid residual, and introuduce additional 
        !< residual due to  viscosity, turbulence and source
        !< terms.
        implicit none
        type(controltype), intent(in) :: control
        !< Control parameters
        type(schemetype), intent(in) :: scheme
        !< finite-volume Schemes
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx
        type(boundarytype), intent(in) :: bc
        !< boundary conditions and fixed values
        real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var), intent(inout):: qp
        !< Store primitive variable at cell center
        real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in):: Temp
        !< Store Temperature variable at cell center
        real(wp), dimension(:, :, :, :), intent(inout)  :: residue
        !< Store residue at each cell-center
        real(wp), dimension(:, :, :, :), intent(inout) :: F
        !< Store fluxes throught the I faces
        real(wp), dimension(:, :, :, :), intent(inout) :: G
        !< Store fluxes throught the J faces
        real(wp), dimension(:, :, :, :), intent(inout) :: H
        !< Store fluxes throught the K faces
        type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
        !< Input cell quantities: volume
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Input varaible which stores I faces' area and unit normal
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Input varaible which stores J faces' area and unit normal
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Input varaible which stores K faces' area and unit normal

        call apply_interface(qp, control, bc, dims)
        call populate_ghost_primitive(qp, Ifaces, Jfaces, Kfaces, control, scheme, flow, bc, dims)
        call compute_face_interpolant(qp, cells, scheme, flow, dims)
        call reconstruct_boundary_state(qp, control, scheme, bc, dims)
        call compute_fluxes(F,G,H,Ifaces,Jfaces,Kfaces,scheme, flow, bc, dims)
        if (flow%mu_ref /= 0.0) then
          call evaluate_all_gradients(qp,Temp,cells,Ifaces,Jfaces,Kfaces,scheme,bc,dims)
          call calculate_viscosity(qp, scheme, flow, bc, dims)
          call compute_viscous_fluxes(F, G, H, qp, cells, Ifaces,Jfaces,Kfaces,scheme, flow, dims)
        end if
        call compute_residue(residue, F,G,H,dims)
        call add_source_term_residue(qp, residue, cells, Ifaces,Jfaces,Kfaces,scheme, flow, dims)

      end subroutine get_total_conservative_Residue

end module update
