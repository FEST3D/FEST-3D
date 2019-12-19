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
  use global_vars, only : volume
!  use global_vars, only : qp
!  use global_vars, only : density
!  use global_vars, only : x_speed
!  use global_vars, only : y_speed
!  use global_vars, only : z_speed
!  use global_vars, only : pressure
  use global_vars, only : dist
  use global_vars, only : mu
  use global_vars, only : mu_t
  use global_vars, only : delta_t
  use global_vars, only : process_id

  use global_vars, only: F_p
  use global_vars, only: G_p
  use global_vars, only: H_p
  use global_vars, only: mass_residue
  use global_vars, only: x_mom_residue
  use global_vars, only: y_mom_residue
  use global_vars, only: z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only: TKE_residue
  use global_vars, only: omega_residue
  use global_vars, only: kl_residue
  use global_vars, only: residue

  use geometry   , only: CellCenter

  use utils, only: alloc
  use utils, only:  dealloc 

  !subroutine for residual calculation
  use interface1,                      only: apply_interface
  use bc_primitive,                   only: populate_ghost_primitive
  use face_interpolant,               only: compute_face_interpolant
  use boundary_state_reconstruction,  only: reconstruct_boundary_state
  use scheme,                         only: compute_fluxes
  use gradients,                      only: evaluate_all_gradients
  use viscosity                      ,only: calculate_viscosity
  use viscous,                        only: compute_viscous_fluxes
!  use turbulent_fluxes,               only: compute_turbulent_fluxes
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
!  use plusgs     , only : destroy_plusgs
  use lusgs     , only : update_with_lusgs
  use lusgs     , only : setup_lusgs
!  use lusgs     , only : destroy_lusgs
#include "debug.h"
#include "error.h"
#include "mpi.inc"
    private

    real, dimension(:,:,:,:), allocatable :: U_store
    !< Array to store the intermediate solution
    real, dimension(:,:,:,:), allocatable :: R_store
    !< Array to store the intermediate Residue
    real, dimension(:,:,:,:), allocatable, target :: aux
    !< Array to store some auxilary intermediate variables
    real, dimension(:)      , allocatable :: u1
    !< Variable array old for each cell center
    real, dimension(:)      , allocatable :: u2
    !< Variable array new for each cell center
    real, dimension(:)      , allocatable :: R
    !< Residue array for each cell center
    integer :: imx, jmx, kmx, n_var

    ! Public methods
    public :: setup_update
!    public :: destroy_update
    public :: get_next_solution

    contains


      subroutine setup_update(control, scheme,flow, dims)
        !< Allocate memory to variables required based 
        !< on the time-integration scheme.
        implicit none
        type(controltype), intent(in) :: control
        type(schemetype), intent(in) :: scheme
        type(flowtype), intent(in) :: flow
        type(extent), intent(in) :: dims

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

!
!      subroutine destroy_update()
!        !< Dellocate memory from all variables
!        implicit none
!
!        select case (time_step_accuracy)
!          case ("none")
!            ! Do nothing
!            continue
!          case ("RK2", "RK4")
!            call dealloc(U_store)
!            call dealloc(R_store)
!          case ("TVDRK2","TVDRK3")
!            call dealloc(U_store)
!          case ("implicit")
!            call destroy_lusgs()
!          case ("plusgs")
!            call destroy_plusgs()
!          case default
!            Fatal_error
!        end select
!        call dealloc(u1)
!        call dealloc(u2)
!        call dealloc(R)
!        call dealloc(aux)
!
!      end subroutine destroy_update


      subroutine get_next_solution(qp, Temp, control, scheme, flow, dims)
        !< Get solution at next time-step using scheme
        !< given in the input file.
        implicit none
        type(controltype), intent(in) :: control
        type(schemetype), intent(in) :: scheme
        type(flowtype), intent(in) :: flow
        type(extent), intent(in) :: dims
        real, dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var), intent(inout):: qp
        real, dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(inout):: Temp
        real :: CFL 
        CFL = control%CFL
        !finding the updated Temperature using ideal gas law
        !T=P/(R_gas*Rho)
        Temp = qp(:,:,:,5)/(flow%R_gas*qp(:,:,:,1) )
        select case (trim(scheme%time_step_accuracy))
            case ("none")
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, scheme, flow, "conservative", 1. ,1., .FALSE.) 
            case ("RK4")
              R_store=0.
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, scheme, flow, "conservative", 0.5  , 1., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 0.5  , 2., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 1.0  , 2., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 1./6., 1., .TRUE. , R_store, U_store) 
            case("RK2")
              R_store=0.
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, scheme, flow, "conservative", 0.5  , 1., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 0.5  , 1., .TRUE., R_store, U_store) 
            case ("TVDRK3")
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, scheme, flow, "conservative", 1.0  , 1.) 
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 1.0  , 1.) 
              qp = 0.75*U_store + 0.25*qp
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 1.0  , 1.) 
              qp = (1./3.)*U_store + (2./3.)*qp
            case ("TVDRK2")
              U_store = qp
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with(qp, scheme, flow, "conservative", 1.0  , 1.) 
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call update_with(qp, scheme, flow, "conservative", 1.0  , 1.) 
              qp = 0.5*U_store + 0.5*qp
            case ("implicit")
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with_lusgs(qp, scheme, dims)
            case ("plusgs")
              call get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
              call compute_time_step(qp, CFL, scheme, flow, dims) ! has to be after get_..._Residue()
              call update_with_plusgs(qp, scheme, dims)
            case default
              Fatal_error
        end select
      end subroutine get_next_solution

      subroutine update_with(qp, scheme, flow, type, time_factor, store_factor, use, Rn, un)
        !< A generalized scheme to updat the solution explicitly using
        !< any RK method and even first order euler explicit.
        implicit none
        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(inout), target:: qp
        type(schemetype), intent(in) :: scheme
        type(flowtype), intent(in) :: flow
        character(len=*), intent(in) :: type
        real, intent(in), optional :: time_factor ! time factor
        real, intent(in), optional :: store_factor
        logical, intent(in), optional :: use
        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in), optional, target :: un
        real, dimension(1:imx-1,1:jmx-1,1:kmx-1,1:n_var), intent(inout), optional :: Rn
        real               :: TF = 1.0 !time factor
        real               :: SF = 1.0!store factor
        Logical               :: TU = .FALSE. !to use or not
        real, dimension(:,:,:,:), pointer :: Quse
        integer :: i,j,k
        real :: KE=0.
        real :: beta

        !sa variables
        real :: vort
        real :: fv1
        real :: fv2
        real :: fw
        real :: g
        real :: scap
        real :: rsa
        real :: kd2
        real :: xi 

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
            
                  u1(1:n_var) = Quse(i,j,k,1:n_var)
            
                  ! finding primitive residue
                  R(1) = mass_residue(i,j,k)
                  R(2) = -1*(u1(2)/u1(1))*mass_residue(i,j,k) + x_mom_residue(i,j,k)/u1(1)
                  R(3) = -1*(u1(3)/u1(1))*mass_residue(i,j,k) + y_mom_residue(i,j,k)/u1(1)
                  R(4) = -1*(u1(4)/u1(1))*mass_residue(i,j,k) + z_mom_residue(i,j,k)/u1(1)
                  R(5) = 0.5*(flow%gm-1.)*(sum(u1(2:4)**2)*mass_residue(i,j,k)) &
                         -(flow%gm-1.)*u1(2)*x_mom_residue(i,j,k)               &
                         -(flow%gm-1.)*u1(3)*y_mom_residue(i,j,k)               &
                         -(flow%gm-1.)*u1(4)*z_mom_residue(i,j,k)               &
                         +(flow%gm-1.)*energy_residue(i,j,k)
            
                  select case(scheme%turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sst', 'sst2003')
                      beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
                      R(5) = R(5) - (flow%gm-1.)*TKE_residue(i,j,k)
                      R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                             +(1./(1.+bstar*u1(6)*delta_t(i,j,k)))*TKE_residue(i,j,k)/u1(1)
                      R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                             +(1./(1.+2.*beta*u1(6)*delta_t(i,j,k)))*omega_residue(i,j,k)/u1(1)
                    case('kkl')
                      eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
                      fphi = (1+cd1*eta)/(1+eta**4)
                      R(5) = R(5) - (flow%gm-1.)*TKE_residue(i,j,k)
                      R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                             + (1./(1.+((2.5*((cmu**0.75)*u1(1)*(u1(6)**1.5)/max(u1(7),1.e-20))&
                             +(2*mu(i,j,k)/dist(i,j,k)**2))*delta_t(i,j,k))))*TKE_residue(i,j,k)/u1(1)
                      R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                             +(1./(1.+(6*mu(i,j,k)*fphi/dist(i,j,k)**2)*delta_t(i,j,k)))*kl_residue(i,j,k)/u1(1)
                    case DEFAULT
                      Fatal_error
                  end select
            
                        
                 !check if user want to store residue
                  if(present(Rn)) then
                    Rn(i,j,k,1:n_var) = Rn(i,j,k,1:n_var) + SF*R(1:n_var)
                    if(TU) R(:) = Rn(i,j,k,:)
                  end if
                 
            
                 !update
                 u2(:) = u1(:) - R(:)*(TF*delta_t(i,j,k)/volume(i,j,k))
            
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
                 u2(1:n_var) = u1(1:n_var) - R(1:n_var)*(TF*delta_t(i,j,k)/volume(i,j,k))

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


      subroutine get_total_conservative_Residue(qp, Temp, control, scheme, flow, dims)
        !< For each iteration it apply boundary conditions,
        !< use higher order method to reconstruct state at
        !< face, evalute fluxes at each face, calculate 
        !< inviscid residual, and introuduce additional 
        !< residual due to  viscosity, turbulence and source
        !< terms.
        implicit none
        type(controltype), intent(in) :: control
        type(schemetype), intent(in) :: scheme
        type(flowtype), intent(in) :: flow
        type(extent), intent(in) :: dims
        real, dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var), intent(inout):: qp
        real, dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in):: Temp

        call apply_interface(qp, control, dims)
        call populate_ghost_primitive(qp, control, scheme, flow, dims)
        call compute_face_interpolant(qp, scheme, flow)
        call reconstruct_boundary_state(qp, control, scheme, dims)
        call compute_fluxes(scheme, flow)
        if (flow%mu_ref /= 0.0) then
          call evaluate_all_gradients(qp,Temp,scheme)
          call calculate_viscosity(qp, scheme, flow, dims)
          call compute_viscous_fluxes(F_p, G_p, H_p, qp, control, scheme, flow, dims)
        end if
        call compute_residue(scheme)
        call add_source_term_residue(qp, control, scheme, flow, dims)

      end subroutine get_total_conservative_Residue
!
!
!      subroutine update_laminar_variables_primitive(time_factor, store_factor, use, tostore, Rn, un)
!        !< Update first five primitive variables in time
!        implicit none
!        real,    intent(in), optional :: time_factor
!        real,    intent(in), optional :: store_factor
!        logical, intent(in), optional :: use
!        logical, intent(in), optional :: tostore
!        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
!        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn
!
!        !local variables
!        real                              :: TF       = 1.0     !time factor
!        real                              :: SF       = 1.0     !store factor
!        real                              :: SFU      = 1.0     !store factor to use inside ijk loop
!        integer                           :: R_switch = 0       !R_store use switch based on TU
!        Logical                           :: TU       = .FALSE. !to use R_store or not
!        Logical                           :: TS       = .FALSE. !to store R_store or not
!        real, dimension(:,:,:,:), pointer :: Quse
!        real, dimension(:,:,:,:), pointer :: Rstore
!        real, dimension(5)                :: Res
!        real                              :: t1 !temp variable 1/density
!        integer                           :: i,j,k
!
!        if(present(time_factor)) TF=time_factor
!        if(present(store_factor)) SF=store_factor
!        if(present(use)) TU=use
!        if(present(tostore)) TS=tostore
!        if(TU)then
!          R_switch=1
!          SFU=SF
!        else
!          R_switch=0
!          SFU=1.0
!        end if
!
!        !check if user want to update from particular solution
!        if(present(un))then
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
!        else
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
!        end if
!        if(present(Rn))then
!          Rstore=>Rn
!        else
!          Rstore=>aux!making it point to junk values 
!          R_switch=0 !but never used as switch is zero
!        end if
!
!
!        !--- Start update---!
!        do k = 1,kmx-1
!          do j = 1,jmx-1
!            do i = 1,imx-1
!        
!              u1(1:5) = Quse(i,j,k,1:5)
!              t1      = 1.0/u1(1)
!              Res(1:5)= Residue(i,j,k,1:5)
!        
!              ! finding primitive residue
!              R(1) = Res(1)
!              R(2) = (-1*(u1(2)*Res(1)) + Res(2))*t1
!              R(3) = (-1*(u1(3)*Res(1)) + Res(3))*t1
!              R(4) = (-1*(u1(4)*Res(1)) + Res(4))*t1
!              R(5) = 0.5*(gm-1.)*(u1(2)*u1(2)+u1(3)*u1(3)+u1(4)*u1(4))*Res(1) &
!                     -(gm-1.)*u1(2)*Res(2)               &
!                     -(gm-1.)*u1(3)*Res(3)               &
!                     -(gm-1.)*u1(4)*Res(4)               &
!                     +(gm-1.)*Res(5)
!        
!                    
!             !store residue
!              aux(i,j,k,1:5)=R(1:5)
!             
!        
!             !update
!             u2(1:5) = u1(1:5) - (SFU*R(1:5)+R_switch*Rstore(i,j,k,1:5))&
!                                *(TF*delta_t(i,j,k)/volume(i,j,k))
!        
!             qp(i,j,k,1:5) = u2(1:5)
!            end do
!          end do
!        end do
!        if(present(Rn) .and. TS) then
!          Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) = Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) &
!                                      + SF*aux(1:imx-1,1:jmx-1,1:kmx-1,1:5)
!        end if
!        
!        if(any(qp(:,:,:,1)<0.) .and. any(qp(:,:,:,5)<0))then
!          Fatal_error
!        end if
!        nullify(Rstore)
!        nullify(Quse)
!
!      end subroutine update_laminar_variables_primitive
!
!      subroutine update_laminar_variables_conservative(time_factor, store_factor, use, tostore, Rn, un)
!        !< Update first five conservative variables in time
!        implicit none
!        real,    intent(in), optional :: time_factor
!        real,    intent(in), optional :: store_factor
!        logical, intent(in), optional :: use
!        logical, intent(in), optional :: tostore
!        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
!        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn
!
!        !local variables
!        real                              :: TF       = 1.0     !time factor
!        real                              :: SF       = 1.0     !store factor
!        real                              :: SFU      = 1.0     !store factor to use inside ijk loop
!        integer                           :: R_switch = 0       !R_store use switch based on TU
!        Logical                           :: TU       = .FALSE. !to use R_store or not
!        Logical                           :: TS       = .FALSE. !to store R_store or not
!        real, dimension(:,:,:,:), pointer :: Quse
!        real, dimension(:,:,:,:), pointer :: Rstore
!        integer                           :: i,j,k
!
!        if(present(time_factor)) TF=time_factor
!        if(present(store_factor)) SF=store_factor
!        if(present(use)) TU=use
!        if(present(tostore)) TS=tostore
!        if(TU)then
!          R_switch=1
!          SFU=SF
!        else
!          R_switch=0
!          SFU=1.0
!        end if
!
!        !check if user want to update from particular solution
!        if(present(un))then
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
!        else
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
!        end if
!        if(present(Rn))then
!          Rstore=>Rn
!        else
!          Rstore=>aux!making it point to junk values 
!          R_switch=0 !but never used as switch is zero
!        end if
!
!
!        !--- Start update---!
!
!        do k = 1,kmx-1
!          do j = 1,jmx-1
!            do i = 1,imx-1
!
!              ! getting conservative variable
!              u1(1)  = Quse(i,j,k,1)
!              u1(2:) = Quse(i,j,k,2:)*u1(1)
!              u1(5) = (u1(5)/(gm-1.) + 0.5*sum(u1(2:4)**2))/u1(1)
!
!             ! get R
!              R(1:5) = residue(i,j,k,1:5)                                                              
!
!             !store conservative variables
!             aux(i,j,k,1:n_var)=u1(1:n_var)
!
!             !update
!             u2(1:5) = u1(1:5) - (SFU*R(1:5)+R_switch*Rstore(i,j,k,1:5))&
!                                *(TF*delta_t(i,j,k)/volume(i,j,k))
!        
!
!              ! getting primitve variable back variable
!              u2(1)  = u2(1)
!              u2(2:) = u2(2:)/u2(1)
!              u2(5) = (gm-1.)*u2(1)*(u2(5) - (0.5*sum(u2(2:4)**2)) )
!
!              qp(i,j,k,1:5) = u2(1:5)
!
!            end do
!          end do
!        end do
!
!        if(present(Rn) .and. TS) then
!          Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) = Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) &
!                                      + SF*aux(1:imx-1,1:jmx-1,1:kmx-1,1:5)
!        end if
!        
!        if(any(qp(:,:,:,1)<0.) .and. any(qp(:,:,:,5)<0))then
!          Fatal_error
!        end if
!        nullify(Rstore)
!        nullify(Quse)
!
!      end subroutine update_laminar_variables_conservative
!
!
!      subroutine update_turbulent_variables_primitive(time_factor, store_factor, use, tostore, Rn, un)
!        !< Update primitive turbulence variables in time
!        implicit none
!        !arguments
!        real,    intent(in), optional :: time_factor ! time factor
!        real,    intent(in), optional :: store_factor
!        logical, intent(in), optional :: use
!        logical, intent(in), optional :: tostore
!        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
!        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn
!
!        !local variables
!        real                              :: TF       = 1.0     !time factor
!        real                              :: SF       = 1.0     !store factor
!        real                              :: SFU      = 1.0     !store factor to use in ijk loop
!        integer                           :: R_switch = 0       !R_store use switch based on TU
!        Logical                           :: TU       = .FALSE. !to use R_store or not
!        Logical                           :: TS       = .FALSE. !to store R_store or not
!        real, dimension(:,:,:,:), pointer :: Quse
!        real, dimension(:,:,:,:), pointer :: Rstore
!        integer                           :: i,j,k
!        real                              :: beta
!
!        if(present(time_factor)) TF=time_factor
!        if(present(store_factor)) SF=store_factor
!        if(present(use)) TU=use
!        if(present(tostore)) TS=tostore
!        if(TU)then
!          R_switch=1
!          SFU=SF
!        else
!          R_switch=0
!          SFU=1.0
!        end if
!
!        !check if user want to update from particular solution
!        if(present(un))then
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
!        else
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
!        end if
!        if(present(Rn))then
!          Rstore=>Rn
!        else
!          Rstore=>aux!making it point to junk values 
!          R_switch=0 !but never used as switch is zero
!        end if
!
!
!        !--- Start update---!
!
!        select case(turbulence)
!          case('none')
!            !do nothing
!            continue
!          
!          case('sst', 'sst2003')
!            do k = 1,kmx-1
!              do j = 1,jmx-1
!                do i = 1,imx-1
!                  !update primitive variable
!                  u1(6:n_var) = Quse(i,j,k,6:n_var)
!                  beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
!                  R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
!                         +(1./(1.+bstar*u1(6)*delta_t(i,j,k)))*TKE_residue(i,j,k)/u1(1)
!                  R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
!                         +(1./(1.+2.*beta*u1(6)*delta_t(i,j,k)))*omega_residue(i,j,k)/u1(1)
!                  !store residue
!                  aux(i,j,k,6:n_var)=R(6:n_var)
!                  !update
!                  qp(i,j,k,6:n_var) = u1(6:n_var) - (SFU*R(6:n_var)+R_switch*Rstore(i,j,k,6:n_var))&
!                                              *(TF*delta_t(i,j,k)/volume(i,j,k))
!                end do
!              end do
!            end do
!          case('kkl')
!            do k = 1,kmx-1
!              do j = 1,jmx-1
!                do i = 1,imx-1
!                  u1(6:n_var) = Quse(i,j,k,6:n_var)
!                  eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
!                  fphi = (1+cd1*eta)/(1+eta**4)
!                  R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
!                         + (1./(1.+((2.5*((cmu**0.75)*u1(1)*(u1(6)**1.5)/max(u1(7),1.e-20))&
!                         +(2*mu(i,j,k)/dist(i,j,k)**2))*delta_t(i,j,k))))*TKE_residue(i,j,k)/u1(1)
!                  R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
!                         +(1./(1.+(6*mu(i,j,k)*fphi/dist(i,j,k)**2)*delta_t(i,j,k)))*kl_residue(i,j,k)/u1(1)
!                  !store residue
!                  aux(i,j,k,6:n_var)=R(6:n_var)
!                  !update
!                  qp(i,j,k,6:n_var) = u1(6:n_var) - (SFU*R(6:n_var)+R_switch*Rstore(i,j,k,6:n_var))&
!                                              *(TF*delta_t(i,j,k)/volume(i,j,k))
!                  qp(i,j,k,6:n_var) = max(1e-10,qp(i,j,k,6:n_var))
!                end do
!              end do
!            end do
!          case DEFAULT
!            Fatal_error
!        end select
!
!        !--- end update ---!
!        
!        if(present(Rn) .and. TS) then
!          Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) = Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) &
!                                          + SF*aux(1:imx-1,1:jmx-1,1:kmx-1,6:n_var)
!        end if
!        
!        !if(any(qp(:,:,:,6)<0.) .and. any(qp(:,:,:,n_var)<0))then
!        !  Fatal_error
!        !end if
!        nullify(Rstore)
!
!      end subroutine update_turbulent_variables_primitive
!
!      subroutine update_turbulent_variables_conservative(time_factor, store_factor, use, tostore, Rn, un)
!        !< Update conservative turbulence variables in time
!        implicit none
!        !arguments
!        real,    intent(in), optional :: time_factor
!        real,    intent(in), optional :: store_factor
!        logical, intent(in), optional :: use
!        logical, intent(in), optional :: tostore
!        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
!        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn
!
!        !local variables
!        real                              :: TF       = 1.0     !time factor
!        real                              :: SF       = 1.0     !store factor
!        real                              :: SFU      = 1.0     !store factor to use inside ijk loop
!        integer                           :: R_switch = 0       !R_store use switch based on TU
!        Logical                           :: TU       = .FALSE. !to use R_store or not
!        Logical                           :: TS       = .FALSE. !to store R_store or not
!        real, dimension(:,:,:,:), pointer :: Quse
!        real, dimension(:,:,:,:), pointer :: Rstore
!        real                              :: beta
!        integer                           :: i,j,k
!
!        if(present(time_factor)) TF=time_factor
!        if(present(store_factor)) SF=store_factor
!        if(present(use)) TU=use
!        if(present(tostore)) TS=tostore
!        if(TU)then
!          R_switch=1
!          SFU=SF
!        else
!          R_switch=0
!          SFU=1.0
!        end if
!
!        !check if user want to update from particular solution
!        if(present(un))then
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
!        else
!          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
!        end if
!        if(present(Rn))then
!          Rstore=>Rn
!        else
!          Rstore=>aux!making it point to junk values 
!          R_switch=0 !but never used as switch is zero
!        end if
!
!
!        !--- Start update---!
!
!        select case(turbulence)
!          case('none')
!            !do nothing
!            continue
!          
!          case('sst', 'sst2003')
!            do k = 1,kmx-1
!              do j = 1,jmx-1
!                do i = 1,imx-1
!                  u1(1:n_var) = aux(i,j,k,1:n_var)
!                  !get R
!                  R(6:n_var) = residue(i,j,k,6:n_var) 
!                  !point implicit
!                  beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
!                  R(6) = R(6)/(1+(beta*qp(i,j,k,7)*delta_t(i,j,k)))
!                  R(7) = R(7)/(1+(2*beta*qp(i,j,k,7)*delta_t(i,j,k)))
!                  !update
!                  u2(6:n_var) = u1(6:n_var) - (SFU*R(6:n_var)&
!                                              +R_switch*Rstore(i,j,k,6:n_var))&
!                                             *(TF*delta_t(i,j,k)/volume(i,j,k))
!                  if(u2(6)>0.) qp(i,j,k,6) = u2(6)/aux(i,j,k,1)
!                  if(u2(7)>0.) qp(i,j,k,7) = u2(7)/aux(i,j,k,1)
!                end do
!              end do
!            end do
!          case('kkl')
!            do k = 1,kmx-1
!              do j = 1,jmx-1
!                do i = 1,imx-1
!                  u1(6:n_var) = aux(i,j,k,6:n_var)
!                  !get R
!                  R(6:n_var) = residue(i,j,k,6:n_var) 
!                  !point implicit
!                  eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
!                  fphi = (1+cd1*eta)/(1+eta**4)
!                  R(6) = R(6)/(1.+((2.5*((cmu**0.75)*sqrt(u1(1))*(u1(6)**1.5)/max(u1(7),1.e-20))&
!                         +(2*mu(i,j,k)/(dist(i,j,k)**2)))*delta_t(i,j,k)))
!                  R(7) = R(7)/(1.+(6*mu(i,j,k)*fphi/(dist(i,j,k)**2))*delta_t(i,j,k))
!                  !update
!                  u2(6:n_var) = u1(6:n_var) - (R(6:n_var)+(R_switch*SF*Rstore(i,j,k,6:n_var)))&
!                                             *(TF*delta_t(i,j,k)/volume(i,j,k))
!                  qp(i,j,k,6:n_var) = max(1.e-10, u2(6:n_var)/aux(i,j,k,1))
!                end do
!              end do
!            end do
!          case DEFAULT
!            Fatal_error
!        end select
!
!        !--- end update ---!
!        
!        if(present(Rn)) then
!          Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) = Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) &
!                                          + SF*residue(1:imx-1,1:jmx-1,1:kmx-1,6:n_var)
!        end if
!        
!        if(any(qp(1:imx-1,1:jmx-1,1:kmx-1,6:n_var)<0.))then
!          Fatal_error
!        end if
!        nullify(Rstore)
!
!      end subroutine update_turbulent_variables_conservative
!
end module update
