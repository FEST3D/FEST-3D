module update
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
  use global_vars  ,only :   gradu_x
  use global_vars  ,only :   gradu_y
  use global_vars  ,only :   gradu_z
  use global_vars  ,only :   gradv_x
  use global_vars  ,only :   gradv_y
  use global_vars  ,only :   gradv_z
  use global_vars  ,only :   gradw_x
  use global_vars  ,only :   gradw_y
  use global_vars  ,only :   gradw_z
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : R_gas
  use global_vars, only : Pr
  use global_vars, only : tPr

  use global_vars, only : volume
  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
    
  use global_vars, only : n_var
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
  use global_vars, only : gm
  use global_vars, only : sst_n_var
  use global_vars, only : qp
  use global_vars, only : qp_inf
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : pressure
  use global_vars, only : dist
  use global_vars, only : mu
  use global_vars, only : mu_t
  use global_vars, only : tk_inf
  use global_vars, only : tkl_inf

  use global_vars, only : time_stepping_method
  use global_vars, only : time_step_accuracy
  use global_vars, only : global_time_step
  use global_vars, only : delta_t
  use global_vars, only : turbulence
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
  use global_vars, only: mu_ref

  use geometry   , only: CellCenter

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL

  use string

  !subroutine for residual calculation
  use face_interpolant,               only: interpolant
  use global_vars,                    only : mu_ref
  use interface,                      only: apply_interface
  use bc_primitive,                   only: populate_ghost_primitive
  use face_interpolant,               only: compute_face_interpolant
  use boundary_state_reconstruction,  only: reconstruct_boundary_state
  use scheme,                         only: compute_fluxes
  use summon_grad_evaluation,         only: evaluate_all_gradients
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

  use lusgs     , only : update_with_lusgs
  use lusgs     , only : setup_lusgs
  use lusgs     , only : destroy_lusgs
#include "error.inc"
#include "mpi.inc"
    private

    real, dimension(:,:,:,:), allocatable :: U_store
    real, dimension(:,:,:,:), allocatable :: R_store
    real, dimension(:,:,:,:), allocatable, target :: aux
    real, dimension(:)      , allocatable :: u1
    real, dimension(:)      , allocatable :: u2
    real, dimension(:)      , allocatable :: R
    real :: eps=0.05
    real, dimension(:,:,:,:), allocatable :: delQ
    real, dimension(:,:,:,:), allocatable :: delQstar

    ! Public methods
    public :: setup_update
    public :: destroy_update
    public :: get_next_solution

    contains


      subroutine setup_update()
        implicit none

        call alloc(u1,1,n_var)
        call alloc(u2,1,n_var)
        call alloc(R ,1,n_var)
        call alloc(aux,-2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)

        select case (time_step_accuracy)
          case ("none")
            ! Do nothing
            continue
          case ("RK2", "RK4")
            call alloc(U_store,-2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)
            call alloc(R_store, 1,imx-1, 1,jmx-1, 1,kmx-1,1,n_var)
          case ("TVDRK2", "TVDRK3", "TVDRK4")
            call alloc(U_store,-2,imx+2,-2,jmx+2,-2,kmx+2,1,n_var)
          case ("implicit")
            call setup_lusgs()
            call alloc(delQ, 0, imx, 0, jmx, 0, kmx, 1, n_var)
            call alloc(delQstar, 0, imx, 0, jmx, 0, kmx, 1, n_var)
          case default
            Fatal_error
        end select

      end subroutine setup_update


      subroutine destroy_update()
        implicit none

        select case (time_step_accuracy)
          case ("none")
            ! Do nothing
            continue
          case ("RK2", "RK4")
            call dealloc(U_store)
            call dealloc(R_store)
          case ("TVDRK2","TVDRK3", "TVDRK4")
            call dealloc(U_store)
          case ("implicit")
            call destroy_lusgs()
            call dealloc(delQ)
            call dealloc(delQstar)
          case default
            Fatal_error
        end select
        call dealloc(u1)
        call dealloc(u2)
        call dealloc(R)
        call dealloc(aux)

      end subroutine destroy_update


      subroutine get_next_solution()
        implicit none
        select case (time_step_accuracy)
            case ("none")
              call get_total_conservative_Residue()
              call compute_time_step() ! has to be after get_..._Residue()
              call update_with("conservative", 1. ,1., .FALSE.) 
              !call update_with("primitive", 1. ,1., .FALSE.) 
              !call update_laminar_variables_primitive(1. ,1., .FALSE.) 
              !if(turbulence/='none')then
              !  call update_turbulent_variables_primitive(1. ,1., .FALSE.) 
              !end if
            case ("RK4")
              R_store=0.
              U_store = qp
              call get_total_conservative_Residue()
              call compute_time_step()
              call update_with("conservative", 0.5  , 1., .FALSE., R_store, U_store) 
              !call update_laminar_variables_primitive(0.5, 1., .FALSE., .True., R_store, U_store) 
              !call update_turbulent_variables_primitive(0.5, 1., .FALSE., .True., R_store, U_store) 
              call get_total_conservative_Residue()
              call update_with("conservative", 0.5  , 2., .FALSE., R_store, U_store) 
              !call update_laminar_variables_primitive(0.5, 2., .FALSE., .True., R_store, U_store) 
              !call update_turbulent_variables_primitive(0.5, 2., .FALSE., .True., R_store, U_store) 
              call get_total_conservative_Residue()
              call update_with("conservative", 1.0  , 2., .FALSE., R_store, U_store) 
              !call update_laminar_variables_primitive(1.0, 2., .FALSE., .True., R_store, U_store) 
              !call update_turbulent_variables_primitive(1.0, 2., .FALSE., .True., R_store, U_store) 
              call get_total_conservative_Residue()
              call update_with("conservative", 1./6., 1., .TRUE. , R_store, U_store) 
              !call update_laminar_variables_primitive(1./6., 1., .TRUE., .FALSE., R_store, U_store) 
              !call update_turbulent_variables_primitive(1./6., 1., .TRUE., .FALSE., R_store, U_store) 
            case("RK2")
              R_store=0.
              U_store = qp
              call get_total_conservative_Residue()
              call compute_time_step()
              call update_with("conservative", 0.5  , 1., .FALSE., R_store, U_store) 
              call get_total_conservative_Residue()
              call update_with("conservative", 0.5  , 1., .TRUE., R_store, U_store) 
            case ("TVDRK3")
              U_store = qp
              call get_total_conservative_Residue()
              call compute_time_step()
              call update_with("conservative", 1.0  , 1.) 
              call get_total_conservative_Residue()
              call update_with("conservative", 1.0  , 1.) 
              qp = 0.75*U_store + 0.25*qp
              call get_total_conservative_Residue()
              call update_with("conservative", 1.0  , 1.) 
              qp = (1./3.)*U_store + (2./3.)*qp
            case ("TVDRK2")
              U_store = qp
              call get_total_conservative_Residue()
              call compute_time_step()
              call update_with("conservative", 1.0  , 1.) 
              call get_total_conservative_Residue()
              call update_with("conservative", 1.0  , 1.) 
              qp = 0.5*U_store + 0.5*qp
            case ("TVDRK4")
              U_store = qp
              call get_total_conservative_Residue()
              call compute_time_step()
              call update_laminar_variables_conservative(0.25, un=U_store) 
              if(turbulence/='none')then
                call update_turbulent_variables_conservative(1.0, un=U_store)
              end if
              call get_total_conservative_Residue()
              call update_laminar_variables_conservative(0.333333, un=U_store) 
              !if(turbulence/='none')then
              !  call update_turbulent_variables_conservative(0.333333, un=U_store)
              !end if
              call get_total_conservative_Residue()
              call update_laminar_variables_conservative(0.5, un=U_store) 
              !if(turbulence/='none')then
              !  call update_turbulent_variables_conservative(0.5, un=U_store)
              !end if
              call get_total_conservative_Residue()
              call update_laminar_variables_conservative(1.00, un=U_store) 
              !if(turbulence/='none')then
              !  call update_turbulent_variables_conservative(1.00, un=U_store)
              !end if
            case ("implicit")
              call get_total_conservative_Residue()
              call compute_time_step() ! has to be after get_..._Residue()
              call update_with_lusgs()
              !if(mu_ref/=0.0 .and. trim(turbulence)=='none') then
              !  call matrix_free_implicit_update_viscous()
              !elseif(mu_ref/=0.0 .and. trim(turbulence)=='sst') then
              !  call matrix_free_implicit_update_SST()
              !else
              ! call matrix_free_implicit_update()
              !end if
            case default
              Fatal_error
        end select
      end subroutine get_next_solution

      subroutine update_with(type, time_factor, store_factor, use, Rn, un)
        implicit none
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
            !include "update_primitive.inc"

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
                  R(5) = 0.5*(gm-1.)*(sum(u1(2:4)**2)*mass_residue(i,j,k)) &
                         -(gm-1.)*u1(2)*x_mom_residue(i,j,k)               &
                         -(gm-1.)*u1(3)*y_mom_residue(i,j,k)               &
                         -(gm-1.)*u1(4)*z_mom_residue(i,j,k)               &
                         +(gm-1.)*energy_residue(i,j,k)
            
                  select case(turbulence)
                    case('none')
                      !do nothing
                      continue
                    case('sst')
                      beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
                      R(5) = R(5) - (gm-1.)*TKE_residue(i,j,k)
                      R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                             +(1./(1.+bstar*u1(6)*delta_t(i,j,k)))*TKE_residue(i,j,k)/u1(1)
                      R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                             +(1./(1.+2.*beta*u1(6)*delta_t(i,j,k)))*omega_residue(i,j,k)/u1(1)
                    case('kkl')
                      eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
                      fphi = (1+cd1*eta)/(1+eta**4)
                      R(5) = R(5) - (gm-1.)*TKE_residue(i,j,k)
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
                    select case(trim(turbulence))
                     case('sst', 'kkl')
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
                  select case(turbulence)
                    case('sst', 'kkl')
                      KE = 0.0!u1(6)
                    case('sa','saBC')
                      KE=0.0
                    case DEFAULT
                      KE = 0.
                  end select
                  u1(5) = (u1(5)/(gm-1.) + 0.5*sum(u1(2:4)**2))/u1(1) + KE

                 ! get R
                  R(1:n_var) = residue(i,j,k,1:n_var) 
                  ! point implicit destruction term
                  select case(trim(turbulence))
                    case('none')
                      !do nothing
                      continue
                    case('sst')
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
                      xi   = U1(6)*density(i,j,k)/mu(i,j,k)
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
                  select case(turbulence)
                    case('sst', 'kkl')
                      KE = 0.0!u2(6)
                    case('sa', 'saBC')
                      !u2(6) = u2(6)*u2(1)
                      KE=0.0
                    case DEFAULT
                      KE = 0.
                  end select
                  u2(5) = (gm-1.)*u2(1)*(u2(5) - (0.5*sum(u2(2:4)**2)) - KE)

                  !check solution for non pyhysical results
                  if((u2(1) < 0.) .or. (u2(5)) < 0. .or. any(isnan(u2)))then
                    print*, u2(:)
                    print*, "R: ", R
                    print*, "old ", U1
                    Fatal_error
                  else !update
                    qp(i,j,k,1:5) = u2(1:5)
                    select case(trim(turbulence))
                     case('sst', 'kkl')
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


      subroutine get_total_conservative_Residue()
        implicit none

!        call send_recv(3) ! parallel call-argument:no of layers 
        call apply_interface()
        call populate_ghost_primitive()
        call compute_face_interpolant()
        call reconstruct_boundary_state(interpolant)
        call compute_fluxes()
        if (mu_ref /= 0.0) then
          call evaluate_all_gradients()
          call calculate_viscosity()
          call compute_viscous_fluxes(F_p, G_p, H_p)
!          call compute_turbulent_fluxes(F_p, G_p, H_p)
        end if
        call compute_residue()
        call add_source_term_residue()

      end subroutine get_total_conservative_Residue


      subroutine update_laminar_variables_primitive(time_factor, store_factor, use, tostore, Rn, un)
        implicit none
        real,    intent(in), optional :: time_factor
        real,    intent(in), optional :: store_factor
        logical, intent(in), optional :: use
        logical, intent(in), optional :: tostore
        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn

        !local variables
        real                              :: TF       = 1.0     !time factor
        real                              :: SF       = 1.0     !store factor
        real                              :: SFU      = 1.0     !store factor to use inside ijk loop
        integer                           :: R_switch = 0       !R_store use switch based on TU
        Logical                           :: TU       = .FALSE. !to use R_store or not
        Logical                           :: TS       = .FALSE. !to store R_store or not
        real, dimension(:,:,:,:), pointer :: Quse
        real, dimension(:,:,:,:), pointer :: Rstore
        real, dimension(5)                :: Res
        real                              :: t1 !temp variable 1/density
        integer                           :: i,j,k

        if(present(time_factor)) TF=time_factor
        if(present(store_factor)) SF=store_factor
        if(present(use)) TU=use
        if(present(tostore)) TS=tostore
        if(TU)then
          R_switch=1
          SFU=SF
        else
          R_switch=0
          SFU=1.0
        end if

        !check if user want to update from particular solution
        if(present(un))then
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
        else
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
        end if
        if(present(Rn))then
          Rstore=>Rn
        else
          Rstore=>aux!making it point to junk values 
          R_switch=0 !but never used as switch is zero
        end if


        !--- Start update---!
        do k = 1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
        
              u1(1:5) = Quse(i,j,k,1:5)
              t1      = 1.0/u1(1)
              Res(1:5)= Residue(i,j,k,1:5)
        
              ! finding primitive residue
              R(1) = Res(1)
              R(2) = (-1*(u1(2)*Res(1)) + Res(2))*t1
              R(3) = (-1*(u1(3)*Res(1)) + Res(3))*t1
              R(4) = (-1*(u1(4)*Res(1)) + Res(4))*t1
              R(5) = 0.5*(gm-1.)*(u1(2)*u1(2)+u1(3)*u1(3)+u1(4)*u1(4))*Res(1) &
                     -(gm-1.)*u1(2)*Res(2)               &
                     -(gm-1.)*u1(3)*Res(3)               &
                     -(gm-1.)*u1(4)*Res(4)               &
                     +(gm-1.)*Res(5)
        
                    
             !store residue
              aux(i,j,k,1:5)=R(1:5)
             
        
             !update
             u2(1:5) = u1(1:5) - (SFU*R(1:5)+R_switch*Rstore(i,j,k,1:5))&
                                *(TF*delta_t(i,j,k)/volume(i,j,k))
        
             qp(i,j,k,1:5) = u2(1:5)
            end do
          end do
        end do
        if(present(Rn) .and. TS) then
          Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) = Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) &
                                      + SF*aux(1:imx-1,1:jmx-1,1:kmx-1,1:5)
        end if
        
        if(any(qp(:,:,:,1)<0.) .and. any(qp(:,:,:,5)<0))then
          Fatal_error
        end if
        nullify(Rstore)
        nullify(Quse)

      end subroutine update_laminar_variables_primitive

      subroutine update_laminar_variables_conservative(time_factor, store_factor, use, tostore, Rn, un)
        implicit none
        real,    intent(in), optional :: time_factor
        real,    intent(in), optional :: store_factor
        logical, intent(in), optional :: use
        logical, intent(in), optional :: tostore
        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn

        !local variables
        real                              :: TF       = 1.0     !time factor
        real                              :: SF       = 1.0     !store factor
        real                              :: SFU      = 1.0     !store factor to use inside ijk loop
        integer                           :: R_switch = 0       !R_store use switch based on TU
        Logical                           :: TU       = .FALSE. !to use R_store or not
        Logical                           :: TS       = .FALSE. !to store R_store or not
        real, dimension(:,:,:,:), pointer :: Quse
        real, dimension(:,:,:,:), pointer :: Rstore
        integer                           :: i,j,k

        if(present(time_factor)) TF=time_factor
        if(present(store_factor)) SF=store_factor
        if(present(use)) TU=use
        if(present(tostore)) TS=tostore
        if(TU)then
          R_switch=1
          SFU=SF
        else
          R_switch=0
          SFU=1.0
        end if

        !check if user want to update from particular solution
        if(present(un))then
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
        else
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
        end if
        if(present(Rn))then
          Rstore=>Rn
        else
          Rstore=>aux!making it point to junk values 
          R_switch=0 !but never used as switch is zero
        end if


        !--- Start update---!

        do k = 1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1

              ! getting conservative variable
              u1(1)  = Quse(i,j,k,1)
              u1(2:) = Quse(i,j,k,2:)*u1(1)
              u1(5) = (u1(5)/(gm-1.) + 0.5*sum(u1(2:4)**2))/u1(1)

             ! get R
              R(1:5) = residue(i,j,k,1:5)                                                              

             !store conservative variables
             aux(i,j,k,1:n_var)=u1(1:n_var)

             !update
             u2(1:5) = u1(1:5) - (SFU*R(1:5)+R_switch*Rstore(i,j,k,1:5))&
                                *(TF*delta_t(i,j,k)/volume(i,j,k))
        

              ! getting primitve variable back variable
              u2(1)  = u2(1)
              u2(2:) = u2(2:)/u2(1)
              u2(5) = (gm-1.)*u2(1)*(u2(5) - (0.5*sum(u2(2:4)**2)) )

              qp(i,j,k,1:5) = u2(1:5)

            end do
          end do
        end do

        if(present(Rn) .and. TS) then
          Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) = Rn(1:imx-1,1:jmx-1,1:kmx-1,1:5) &
                                      + SF*aux(1:imx-1,1:jmx-1,1:kmx-1,1:5)
        end if
        
        if(any(qp(:,:,:,1)<0.) .and. any(qp(:,:,:,5)<0))then
          Fatal_error
        end if
        nullify(Rstore)
        nullify(Quse)

      end subroutine update_laminar_variables_conservative


      subroutine update_turbulent_variables_primitive(time_factor, store_factor, use, tostore, Rn, un)
        implicit none
        !arguments
        real,    intent(in), optional :: time_factor ! time factor
        real,    intent(in), optional :: store_factor
        logical, intent(in), optional :: use
        logical, intent(in), optional :: tostore
        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn

        !local variables
        real                              :: TF       = 1.0     !time factor
        real                              :: SF       = 1.0     !store factor
        real                              :: SFU      = 1.0     !store factor to use in ijk loop
        integer                           :: R_switch = 0       !R_store use switch based on TU
        Logical                           :: TU       = .FALSE. !to use R_store or not
        Logical                           :: TS       = .FALSE. !to store R_store or not
        real, dimension(:,:,:,:), pointer :: Quse
        real, dimension(:,:,:,:), pointer :: Rstore
        integer                           :: i,j,k
        real                              :: beta

        if(present(time_factor)) TF=time_factor
        if(present(store_factor)) SF=store_factor
        if(present(use)) TU=use
        if(present(tostore)) TS=tostore
        if(TU)then
          R_switch=1
          SFU=SF
        else
          R_switch=0
          SFU=1.0
        end if

        !check if user want to update from particular solution
        if(present(un))then
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
        else
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
        end if
        if(present(Rn))then
          Rstore=>Rn
        else
          Rstore=>aux!making it point to junk values 
          R_switch=0 !but never used as switch is zero
        end if


        !--- Start update---!

        select case(turbulence)
          case('none')
            !do nothing
            continue
          
          case('sst')
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1
                  !update primitive variable
                  u1(6:n_var) = Quse(i,j,k,6:n_var)
                  beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
                  R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                         +(1./(1.+bstar*u1(6)*delta_t(i,j,k)))*TKE_residue(i,j,k)/u1(1)
                  R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                         +(1./(1.+2.*beta*u1(6)*delta_t(i,j,k)))*omega_residue(i,j,k)/u1(1)
                  !store residue
                  aux(i,j,k,6:n_var)=R(6:n_var)
                  !update
                  qp(i,j,k,6:n_var) = u1(6:n_var) - (SFU*R(6:n_var)+R_switch*Rstore(i,j,k,6:n_var))&
                                              *(TF*delta_t(i,j,k)/volume(i,j,k))
                end do
              end do
            end do
          case('kkl')
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1
                  u1(6:n_var) = Quse(i,j,k,6:n_var)
                  eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
                  fphi = (1+cd1*eta)/(1+eta**4)
                  R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                         + (1./(1.+((2.5*((cmu**0.75)*u1(1)*(u1(6)**1.5)/max(u1(7),1.e-20))&
                         +(2*mu(i,j,k)/dist(i,j,k)**2))*delta_t(i,j,k))))*TKE_residue(i,j,k)/u1(1)
                  R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                         +(1./(1.+(6*mu(i,j,k)*fphi/dist(i,j,k)**2)*delta_t(i,j,k)))*kl_residue(i,j,k)/u1(1)
                  !store residue
                  aux(i,j,k,6:n_var)=R(6:n_var)
                  !update
                  qp(i,j,k,6:n_var) = u1(6:n_var) - (SFU*R(6:n_var)+R_switch*Rstore(i,j,k,6:n_var))&
                                              *(TF*delta_t(i,j,k)/volume(i,j,k))
                  qp(i,j,k,6:n_var) = max(1e-10,qp(i,j,k,6:n_var))
                end do
              end do
            end do
          case DEFAULT
            Fatal_error
        end select

        !--- end update ---!
        
        if(present(Rn) .and. TS) then
          Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) = Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) &
                                          + SF*aux(1:imx-1,1:jmx-1,1:kmx-1,6:n_var)
        end if
        
        !if(any(qp(:,:,:,6)<0.) .and. any(qp(:,:,:,n_var)<0))then
        !  Fatal_error
        !end if
        nullify(Rstore)

      end subroutine update_turbulent_variables_primitive

      subroutine update_turbulent_variables_conservative(time_factor, store_factor, use, tostore, Rn, un)
        implicit none
        !arguments
        real,    intent(in), optional :: time_factor
        real,    intent(in), optional :: store_factor
        logical, intent(in), optional :: use
        logical, intent(in), optional :: tostore
        real, dimension(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var), intent(in)   , optional, target :: un
        real, dimension( 1:imx-1, 1:jmx-1, 1:kmx-1,1:n_var), intent(inout), optional, target :: Rn

        !local variables
        real                              :: TF       = 1.0     !time factor
        real                              :: SF       = 1.0     !store factor
        real                              :: SFU      = 1.0     !store factor to use inside ijk loop
        integer                           :: R_switch = 0       !R_store use switch based on TU
        Logical                           :: TU       = .FALSE. !to use R_store or not
        Logical                           :: TS       = .FALSE. !to store R_store or not
        real, dimension(:,:,:,:), pointer :: Quse
        real, dimension(:,:,:,:), pointer :: Rstore
        real                              :: beta
        integer                           :: i,j,k

        if(present(time_factor)) TF=time_factor
        if(present(store_factor)) SF=store_factor
        if(present(use)) TU=use
        if(present(tostore)) TS=tostore
        if(TU)then
          R_switch=1
          SFU=SF
        else
          R_switch=0
          SFU=1.0
        end if

        !check if user want to update from particular solution
        if(present(un))then
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>un(:,:,:,:)
        else
          Quse(-2:imx+2,-2:jmx+2,-2:kmx+2,1:n_var)=>qp(:,:,:,:)
        end if
        if(present(Rn))then
          Rstore=>Rn
        else
          Rstore=>aux!making it point to junk values 
          R_switch=0 !but never used as switch is zero
        end if


        !--- Start update---!

        select case(turbulence)
          case('none')
            !do nothing
            continue
          
          case('sst')
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1
                  u1(1:n_var) = aux(i,j,k,1:n_var)
                  !get R
                  R(6:n_var) = residue(i,j,k,6:n_var) 
                  !point implicit
                  beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
                  R(6) = R(6)/(1+(beta*qp(i,j,k,7)*delta_t(i,j,k)))
                  R(7) = R(7)/(1+(2*beta*qp(i,j,k,7)*delta_t(i,j,k)))
                  !update
                  u2(6:n_var) = u1(6:n_var) - (SFU*R(6:n_var)&
                                              +R_switch*Rstore(i,j,k,6:n_var))&
                                             *(TF*delta_t(i,j,k)/volume(i,j,k))
                  if(u2(6)>0.) qp(i,j,k,6) = u2(6)/aux(i,j,k,1)
                  if(u2(7)>0.) qp(i,j,k,7) = u2(7)/aux(i,j,k,1)
                end do
              end do
            end do
          case('kkl')
            do k = 1,kmx-1
              do j = 1,jmx-1
                do i = 1,imx-1
                  u1(6:n_var) = aux(i,j,k,6:n_var)
                  !get R
                  R(6:n_var) = residue(i,j,k,6:n_var) 
                  !point implicit
                  eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
                  fphi = (1+cd1*eta)/(1+eta**4)
                  R(6) = R(6)/(1.+((2.5*((cmu**0.75)*sqrt(u1(1))*(u1(6)**1.5)/max(u1(7),1.e-20))&
                         +(2*mu(i,j,k)/(dist(i,j,k)**2)))*delta_t(i,j,k)))
                  R(7) = R(7)/(1.+(6*mu(i,j,k)*fphi/(dist(i,j,k)**2))*delta_t(i,j,k))
                  !update
                  u2(6:n_var) = u1(6:n_var) - (R(6:n_var)+(R_switch*SF*Rstore(i,j,k,6:n_var)))&
                                             *(TF*delta_t(i,j,k)/volume(i,j,k))
                  qp(i,j,k,6:n_var) = max(1.e-10, u2(6:n_var)/aux(i,j,k,1))
                end do
              end do
            end do
          case DEFAULT
            Fatal_error
        end select

        !--- end update ---!
        
        if(present(Rn)) then
          Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) = Rn(1:imx-1,1:jmx-1,1:kmx-1,6:n_var) &
                                          + SF*residue(1:imx-1,1:jmx-1,1:kmx-1,6:n_var)
        end if
        
        if(any(qp(1:imx-1,1:jmx-1,1:kmx-1,6:n_var)<0.))then
          Fatal_error
        end if
        nullify(Rstore)

      end subroutine update_turbulent_variables_conservative


      subroutine matrix_free_implicit_update()
        !----------------------
        !Reference:
        !AIAA 2000-0927
        !----------------------
        implicit none
        integer :: i,j,k
        real, dimension(1:5)     :: deltaU
        real                     :: D
        real, dimension(1:5)     :: conservativeQ
        real, dimension(1:5)     :: OldIminusFlux
        real, dimension(1:5)     :: OldJminusFlux
        real, dimension(1:5)     :: OldKminusFlux
        real, dimension(1:5)     :: NewIminusFlux
        real, dimension(1:5)     :: NewJminusFlux
        real, dimension(1:5)     :: NewKminusFlux
        real, dimension(1:5)     :: DelIminusFlux
        real, dimension(1:5)     :: DelJminusFlux
        real, dimension(1:5)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:n_var) :: State0
        real, dimension(1:n_var) :: state
        real, dimension(1:n_var) :: stateChange


        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              State0        =       qp(i  ,j,k,1:n_var)
              state         =       qp(i-1,j,k,1:n_var)
              stateChange   = delQstar(i-1,j,k,1:n_var)
              NewIminusFlux = flux(state, stateChange, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k))
              stateChange   = 0.0
              OldIminusFlux = flux(state, stateChange, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k))
       LambdaTimesArea(1)=SpectralRadius(State0,state, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k))

              state         =       qp(i,j-1,k,1:n_var)
              stateChange   = delQstar(i,j-1,k,1:n_var)
              NewJminusFlux = flux(state, stateChange, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k))
              stateChange   = 0.0
              OldJminusFlux = flux(state, stateChange, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k))
       LambdaTimesArea(2)=SpectralRadius(State0,state, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k))

              state         =       qp(i,j,k-1,1:n_var)
              stateChange   = delQstar(i,j,k-1,1:n_var)
              NewKminusFlux = flux(state, stateChange, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k))
              stateChange   = 0.0
              OldKminusFlux = flux(state, stateChange, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k))
       LambdaTimesArea(3)=SpectralRadius(State0,state, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k))

              state         = qp(i+1,j,k,1:n_var)
       LambdaTimesArea(4)=SpectralRadius(State0,state, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k))
              state         = qp(i,j+1,k,1:n_var)
       LambdaTimesArea(5)=SpectralRadius(State0,state, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k))
              state         = qp(i,j,k+1,1:n_var)
       LambdaTimesArea(6)=SpectralRadius(State0,state, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1))


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              deltaU(1:5) = -residue(i,j,k,1:5) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQstar(i-1,j,k,1:5)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQstar(i,j-1,k,1:5)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQstar(i,j,k-1,1:5)) )

              delQstar(i,j,k,1:5) = deltaU(1:5)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              State0        =   qp(i  ,j,k,1:n_var)
              state         =   qp(i+1,j,k,1:n_var)
              stateChange   = delQ(i+1,j,k,1:n_var)
              NewIminusFlux = flux(state, stateChange, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k))
              stateChange   = 0.0
              OldIminusFlux = flux(state, stateChange, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k))
       LambdaTimesArea(1)=SpectralRadius(State0,state, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k))

              state         =   qp(i,j+1,k,1:n_var)
              stateChange   = delQ(i,j+1,k,1:n_var)
              NewJminusFlux = flux(state, stateChange, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k))
              stateChange   = 0.0
              OldJminusFlux = flux(state, stateChange, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k))
       LambdaTimesArea(2)=SpectralRadius(State0,state, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k))

              state         =   qp(i,j,k+1,1:n_var)
              stateChange   = delQ(i,j,k+1,1:n_var)
              NewKminusFlux = flux(state, stateChange, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1))
              stateChange   = 0.0
              OldKminusFlux = flux(state, stateChange, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1))
       LambdaTimesArea(3)=SpectralRadius(State0,state, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1))

              state         = qp(i-1,j,k,1:n_var)
       LambdaTimesArea(4)=SpectralRadius(State0,state, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k))
              state         = qp(i,j-1,k,1:n_var)
       LambdaTimesArea(5)=SpectralRadius(State0,state, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k))
              state         = qp(i,j,k-1,1:n_var)
       LambdaTimesArea(6)=SpectralRadius(State0,state, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k))


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              delQ(i,j,k,1:5) = delQstar(i,j,k,1:5) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQ(i+1,j,k,1:5)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j+1,k,1:5)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k+1,1:5)) )/D

            end do
          end do
        end do
        
        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              
              ! add new change into conservative solution
              conservativeQ(1:5) = conservativeQ(1:5) + delQ(i,j,k,1:5)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
            end do
          end do
        end do
      end subroutine matrix_free_implicit_update


      function flux(q,du,nx, ny, nz, Area)
        implicit none
        real                    , intent(in) :: Area
        real                    , intent(in) :: nx
        real                    , intent(in) :: ny
        real                    , intent(in) :: nz
        real, dimension(1:n_var), intent(in) :: q
        real, dimension(1:n_var), intent(in) :: du
        real, dimension(1:n_var)             :: flux
        real, dimension(1:n_var)             :: U ! conservative variables
        real, dimension(1:n_var)             :: w ! new primitive variables

        real    :: HalfRhoUsquare
        real    :: RhoHt
        real    :: FaceNormalVelocity


        ! find conservative variable
        U(1) = q(1)
        U(2) = q(1) * q(2)
        U(3) = q(1) * q(3)
        U(4) = q(1) * q(4)
        U(5) = ( q(5) / (gm-1.0) ) + ( 0.5 * q(1) * sum(q(2:4)**2)  )

        U(1:5) = U(1:5) + du(1:5)

        W(1) = U(1)
        W(2) = U(2) / U(1)
        W(3) = U(3) / U(1)
        W(4) = U(4) / U(1)
        W(5) = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )

        FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)

        flux(1) =   W(1) * FaceNormalVelocity

        flux(2) = ( W(2) * flux(1) ) + ( W(5) * nx )

        flux(3) = ( W(3) * flux(1) ) + ( W(5) * ny )

        flux(4) = ( W(4) * flux(1) ) + ( W(5) * nz )

        HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
        RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
        flux(5)        = RhoHt * FaceNormalVelocity
        !flux(5)        = ( U(5) + W(5) ) * FaceNormalVelocity
        flux           = flux * Area

      end function flux

      function SpectralRadius(q1, q2, nx, ny, nz, Area)
        real, dimension(1:n_var), intent(in) :: q1
        real, dimension(1:n_var), intent(in) :: q2
        real                    , intent(in) :: nx
        real                    , intent(in) :: ny
        real                    , intent(in) :: nz
        real                    , intent(in) :: Area
        real                                 :: SpectralRadius
        real                                 :: NormalSpeed
        real                                 :: SpeedOfSound

        ! in state vector q (2-4) are the cell center velocity
        NormalSpeed = 0.5 * ( ( ( q1(2) + q2(2) ) * nx ) &
                            + ( ( q1(3) + q2(3) ) * ny ) &
                            + ( ( q1(4) + q2(4) ) * nz ) &
                            )
        NormalSpeed = abs(NormalSpeed)

        !SpeedOfSound = sqrt(gm*(q1(5) + q2(5))/(q1(1) + q2(1)))
        SpeedOfSound = 0.5*( sqrt(gm*q1(5)/q1(1)) + sqrt(gm*q2(5)/q2(1)) )

        SpectralRadius = ( NormalSpeed + SpeedOfSound ) * Area

      end function SpectralRadius


      !----------------------------------------------------
      ! viscous flow implicit update
      !----------------------------------------------------

      subroutine matrix_free_implicit_update_viscous()
        !----------------------
        !Reference:
        !AIAA 2000-0927
        !----------------------
        implicit none
        integer :: i,j,k
        real, dimension(1:5)     :: deltaU
        real                     :: D
        real, dimension(1:5)     :: conservativeQ
        real, dimension(1:5)     :: OldIminusFlux
        real, dimension(1:5)     :: OldJminusFlux
        real, dimension(1:5)     :: OldKminusFlux
        real, dimension(1:5)     :: NewIminusFlux
        real, dimension(1:5)     :: NewJminusFlux
        real, dimension(1:5)     :: NewKminusFlux
        real, dimension(1:5)     :: DelIminusFlux
        real, dimension(1:5)     :: DelJminusFlux
        real, dimension(1:5)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:n_var) :: State0
        real, dimension(1:n_var) :: state
        real, dimension(1:n_var) :: stateChange


        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              State0        =       qp(i  ,j,k,1:n_var)
              state         =       qp(i-1,j,k,1:n_var)
              stateChange   = delQstar(i-1,j,k,1:n_var)
              NewIminusFlux = vflux(state, stateChange, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k),&
                                    0.5*(mu(i-1,j,k) + mu(i,j,k)), 0.5*(volume(i-1,j,k)+volume(i,j,k)), state0)
              stateChange   = 0.0
              OldIminusFlux = vflux(state, stateChange, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k),&
                                   0.5*(mu(i-1,j,k) + mu(i,j,k)), 0.5*(volume(i-1,j,k)+volume(i,j,k)), state0)
       LambdaTimesArea(1)=vSpectralRadius(State0,state, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k)&
                                     , (mu(i-1,j,k) + mu(i,j,k)), CellCenter(i-1,j,k,:), CellCenter(i,j,k,:))

              state         =       qp(i,j-1,k,1:n_var)
              stateChange   = delQstar(i,j-1,k,1:n_var)
              NewJminusFlux = vflux(state, stateChange, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k),&
                                   0.5*(mu(i,j-1,k) + mu(i,j,k)), 0.5*(volume(i,j-1,k)+volume(i,j,k)), state0)
              stateChange   = 0.0
              OldJminusFlux = vflux(state, stateChange, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k)&
                                   ,0.5*(mu(i,j-1,k) + mu(i,j,k)) ,0.5*(volume(i,j-1,k)+volume(i,j,k)), state0)
       LambdaTimesArea(2)=vSpectralRadius(State0,state, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k)&
                                     , (mu(i,j-1,k) + mu(i,j,k)), CellCenter(i,j-1,k,:), CellCenter(i,j,k,:))

              state         =       qp(i,j,k-1,1:n_var)
              stateChange   = delQstar(i,j,k-1,1:n_var)
              NewKminusFlux = vflux(state, stateChange, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k)&
                                    ,0.5*(mu(i,j,k-1) + mu(i,j,k)),0.5*(volume(i,j,k-1)+volume(i,j,k)), state0)
              stateChange   = 0.0
              OldKminusFlux = vflux(state, stateChange, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k)&
                                    ,0.5*(mu(i,j,k-1) + mu(i,j,k)),0.5*(volume(i,j,k-1)+volume(i,j,k)), state0)
       LambdaTimesArea(3)=vSpectralRadius(State0,state, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k)&
                                     , (mu(i,j,k-1) + mu(i,j,k)), CellCenter(i,j,k-1,:), CellCenter(i,j,k,:))

              state         = qp(i+1,j,k,1:n_var)
       LambdaTimesArea(4)=vSpectralRadius(State0,state, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k)&
                                     , (mu(i+1,j,k) + mu(i,j,k)), CellCenter(i+1,j,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j+1,k,1:n_var)
       LambdaTimesArea(5)=vSpectralRadius(State0,state, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k)&
                                     , (mu(i,j+1,k) + mu(i,j,k)), CellCenter(i,j+1,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j,k+1,1:n_var)
       LambdaTimesArea(6)=vSpectralRadius(State0,state, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1)&
                                     , (mu(i,j,k+1) + mu(i,j,k)), CellCenter(i,j,k+1,:), CellCenter(i,j,k,:))


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              deltaU(1:5) = -residue(i,j,k,1:5) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQstar(i-1,j,k,1:5)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQstar(i,j-1,k,1:5)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQstar(i,j,k-1,1:5)) )

              delQstar(i,j,k,1:5) = deltaU(1:5)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              State0        =   qp(i  ,j,k,1:n_var)
              state         =   qp(i+1,j,k,1:n_var)
              stateChange   = delQ(i+1,j,k,1:n_var)
              NewIminusFlux = vflux(state, stateChange, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k)&
                                    ,0.5*(mu(i+1,j,k) + mu(i,j,k)),0.5*(volume(i+1,j,k)+volume(i,j,k)), state0)
              stateChange   = 0.0
              OldIminusFlux = vflux(state, stateChange, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k)&
                                    ,0.5*(mu(i+1,j,k) + mu(i,j,k)),0.5*(volume(i+1,j,k)+volume(i,j,k)), state0)
       LambdaTimesArea(1)=vSpectralRadius(State0,state, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k)&
                                     , (mu(i+1,j,k) + mu(i,j,k)), CellCenter(i+1,j,k,:), CellCenter(i,j,k,:))

              state         =   qp(i,j+1,k,1:n_var)
              stateChange   = delQ(i,j+1,k,1:n_var)
              NewJminusFlux = vflux(state, stateChange, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k)&
                                    ,0.5*(mu(i,j+1,k) + mu(i,j,k)),0.5*(volume(i,j+1,k)+volume(i,j,k)), state0)
              stateChange   = 0.0
              OldJminusFlux = vflux(state, stateChange, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k)&
                                    ,0.5*(mu(i,j+1,k) + mu(i,j,k)),0.5*(volume(i,j+1,k)+volume(i,j,k)), state0)
       LambdaTimesArea(2)=vSpectralRadius(State0,state, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k)&
                                     , (mu(i,j+1,k) + mu(i,j,k)), CellCenter(i,j+1,k,:), CellCenter(i,j,k,:))

              state         =   qp(i,j,k+1,1:n_var)
              stateChange   = delQ(i,j,k+1,1:n_var)
              NewKminusFlux = vflux(state, stateChange, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1)&
                                    ,0.5*(mu(i,j,k+1) + mu(i,j,k)),0.5*(volume(i,j,k+1)+volume(i,j,k)), state0)
              stateChange   = 0.0
              OldKminusFlux = vflux(state, stateChange, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1)&
                                    ,0.5*(mu(i,j,k+1) + mu(i,j,k)),0.5*(volume(i,j,k+1)+volume(i,j,k)), state0)
       LambdaTimesArea(3)=vSpectralRadius(State0,state, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1)&
                                     , (mu(i,j,k+1) + mu(i,j,k)), CellCenter(i,j,k+1,:), CellCenter(i,j,k,:))

              state         = qp(i-1,j,k,1:n_var)
       LambdaTimesArea(4)=vSpectralRadius(State0,state, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k)&
                                     , (mu(i-1,j,k) + mu(i,j,k)), CellCenter(i-1,j,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j-1,k,1:n_var)
       LambdaTimesArea(5)=vSpectralRadius(State0,state, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k)&
                                     , (mu(i,j-1,k) + mu(i,j,k)), CellCenter(i,j-1,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j,k-1,1:n_var)
       LambdaTimesArea(6)=vSpectralRadius(State0,state, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k)&
                                     , (mu(i,j,k-1) + mu(i,j,k)), CellCenter(i,j,k-1,:), CellCenter(i,j,k,:))


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)

              delQ(i,j,k,1:5) = delQstar(i,j,k,1:5) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQ(i+1,j,k,1:5)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j+1,k,1:5)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k+1,1:5)) )/D

            end do
          end do
        end do
        
        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              
              ! add new change into conservative solution
              conservativeQ(1:5) = conservativeQ(1:5) + delQ(i,j,k,1:5)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
            end do
          end do
        end do
      end subroutine matrix_free_implicit_update_viscous


      function vflux(q,du,nx, ny, nz, Area, mu, volume, q2)
        implicit none
        real                    , intent(in) :: Area
        real                    , intent(in) :: nx
        real                    , intent(in) :: ny
        real                    , intent(in) :: nz
        real                    , intent(in) :: volume
        real                    , intent(in) :: mu
        real, dimension(1:n_var), intent(in) :: q
        real, dimension(1:n_var), intent(in) :: q2
        real, dimension(1:n_var), intent(in) :: du
        real, dimension(1:n_var)             :: vflux
        real, dimension(1:n_var)             :: U ! conservative variables
        real, dimension(1:n_var)             :: w ! new primitive variables

        real    :: dudx
        real    :: dudy
        real    :: dudz
        real    :: dvdx
        real    :: dvdy
        real    :: dvdz
        real    :: dwdx
        real    :: dwdy
        real    :: dwdz
        real    :: dTdx
        real    :: dTdy
        real    :: dTdz
        real    :: T1, T2
        real    :: uface
        real    :: vface
        real    :: wface
        real    :: trace
        real    :: Tauxx
        real    :: Tauyy
        real    :: Tauzz
        real    :: Tauxy
        real    :: Tauxz
        real    :: Tauyz
        real    :: Qx
        real    :: Qy
        real    :: Qz
        real    :: HalfRhoUsquare
        real    :: RhoHt
        real    :: K_heat
        real    :: FaceNormalVelocity


        ! find conservative variable
        U(1) = q(1)
        U(2) = q(1) * q(2)
        U(3) = q(1) * q(3)
        U(4) = q(1) * q(4)
        U(5) = ( q(5) / (gm-1.0) ) + ( 0.5 * q(1) * sum(q(2:4)**2)  )

        U(1:5) = U(1:5) + du(1:5)

        W(1) = U(1)
        W(2) = U(2) / U(1)
        W(3) = U(3) / U(1)
        W(4) = U(4) / U(1)
        W(5) = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )

        FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)

        vflux(1) =   W(1) * FaceNormalVelocity

        vflux(2) = ( W(2) * vflux(1) ) + ( W(5) * nx )

        vflux(3) = ( W(3) * vflux(1) ) + ( W(5) * ny )

        vflux(4) = ( W(4) * vflux(1) ) + ( W(5) * nz )

        HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
        RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
        vflux(5)        = RhoHt * FaceNormalVelocity
        !vflux(5)        = ( U(5) + W(5) ) * FaceNormalVelocity

        ! viscous terms
        T1    =  W(5) / ( W(1) * R_gas )
        T2    =  q2(5) / ( q2(1) * R_gas )
        dTdx  =  ( T2    - T1   ) * nx * Area / Volume
        dTdy  =  ( T2    - T1   ) * ny * Area / Volume
        dTdz  =  ( T2    - T1   ) * nz * Area / Volume
        dudx  =  ( q2(2) - W(2) ) * nx * Area / Volume
        dudy  =  ( q2(2) - W(2) ) * ny * Area / Volume
        dudz  =  ( q2(2) - W(2) ) * nz * Area / Volume
        dvdx  =  ( q2(3) - W(3) ) * nx * Area / Volume
        dvdy  =  ( q2(3) - W(3) ) * ny * Area / Volume
        dvdz  =  ( q2(3) - W(3) ) * nz * Area / Volume
        dwdx  =  ( q2(4) - W(4) ) * nx * Area / Volume
        dwdy  =  ( q2(4) - W(4) ) * ny * Area / Volume
        dwdz  =  ( q2(4) - W(4) ) * nz * Area / Volume

        trace = dudx + dvdy + dwdz
        Tauxx =  2. * mu * (dudx - trace/3.0)
        Tauyy =  2. * mu * (dvdy - trace/3.0)
        Tauzz =  2. * mu * (dwdz - trace/3.0)
        Tauxy = mu * (dvdx + dudy)
        Tauxz = mu * (dwdx + dudz)
        Tauyz = mu * (dwdy + dvdz)

        K_heat = ( mu / Pr ) * gm * R_gas / ( gm - 1.0 )
        Qx = K_heat*dTdx
        Qy = K_heat*dTdy
        Qz = K_heat*dTdz

        vflux(2) = vflux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
        vflux(3) = vflux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
        vflux(4) = vflux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
        uface = 0.5 * ( W(2) + q2(2) )
        vface = 0.5 * ( W(3) + q2(3) )
        wface = 0.5 * ( W(4) + q2(4) )

        vflux(5) = vflux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
        vflux(5) = vflux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
        vflux(5) = vflux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz
        vflux    = vflux * Area

      end function vflux

      function vSpectralRadius(q1, q2, nx, ny, nz, Area, mu, c1, c2)
        real, dimension(1:n_var), intent(in) :: q1
        real, dimension(1:n_var), intent(in) :: q2
        real                    , intent(in) :: nx
        real                    , intent(in) :: ny
        real                    , intent(in) :: nz
        real                    , intent(in) :: Area
        real                    , intent(in) :: mu
        real, dimension(1:3)    , intent(in) :: c1
        real, dimension(1:3)    , intent(in) :: c2
        real                                 :: vSpectralRadius
        real                                 :: NormalSpeed
        real                                 :: SpeedOfSound
        real                                 :: vis

        ! in state vector q (2-4) are the cell center velocity
        NormalSpeed = 0.5 * ( ( ( q1(2) + q2(2) ) * nx ) &
                            + ( ( q1(3) + q2(3) ) * ny ) &
                            + ( ( q1(4) + q2(4) ) * nz ) &
                            )
        NormalSpeed = abs(NormalSpeed)

        !SpeedOfSound = sqrt(gm*(q1(5) + q2(5))/(q1(1) + q2(1)))
        SpeedOfSound = 0.5*( sqrt(gm*q1(5)/q1(1)) + sqrt(gm*q2(5)/q2(1)) )

        ! visocus part
       ! vis =  max(4./3.) * gm * mu / ( Pr * 0.5*(q1(1) + q2(1)) * abs((c1(1)-c2(1))*nx + (c1(2)-c2(2))*ny +(c1(3)-c2(3))*nz) )
        vis =  max(4./3.,gm) * mu / ( Pr * 0.5*(q1(1) + q2(1)) &
          * sqrt((c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 +(c1(3)-c2(3))**2) )
        vSpectralRadius = ( NormalSpeed + SpeedOfSound + vis) * Area

      end function vSpectralRadius

      !---------------------------------------------------
      ! SST implicit update
      !---------------------------------------------------

      subroutine matrix_free_implicit_update_SST()
        !----------------------
        !Reference:
        !AIAA 2000-0927
        !----------------------
        implicit none
        integer :: i,j,k
        real, dimension(1:n_var)     :: deltaU
        real, dimension(1:n_var)     :: D
        real, dimension(1:n_var)     :: conservativeQ
        real, dimension(1:n_var)     :: OldIminusFlux
        real, dimension(1:n_var)     :: OldJminusFlux
        real, dimension(1:n_var)     :: OldKminusFlux
        real, dimension(1:n_var)     :: NewIminusFlux
        real, dimension(1:n_var)     :: NewJminusFlux
        real, dimension(1:n_var)     :: NewKminusFlux
        real, dimension(1:n_var)     :: DelIminusFlux
        real, dimension(1:n_var)     :: DelJminusFlux
        real, dimension(1:n_var)     :: DelKminusFlux
        real, dimension(1:6)     :: LambdaTimesArea
        real, dimension(1:n_var) :: State0
        real, dimension(1:n_var) :: state
        real, dimension(1:n_var) :: stateChange
        real, dimension(1:6,1:3) :: fn
        real, dimension(1:6)     :: Vl
        real, dimension(1:6)     :: mm
        real, dimension(1:6)     :: tm
        real, dimension(1:6)     :: fA
        real, dimension(1:6)     :: F1
        real                     :: beta


        !intialize delQ
        delQstar = 0.0

        !forward sweep
        do k=1,kmx-1
          do j=1,jmx-1
            do i=1,imx-1
              fn(1,1)       = -xnx(i,j,k)  ! face normal x direction
              fn(1,2)       = -xny(i,j,k)  ! face normal y direction
              fn(1,3)       = -xnz(i,j,k)  ! face normal z direction
              fn(3,1)       = -ynx(i,j,k)
              fn(3,2)       = -yny(i,j,k)
              fn(3,3)       = -ynz(i,j,k)
              fn(5,1)       = -znx(i,j,k)
              fn(5,2)       = -zny(i,j,k)
              fn(5,3)       = -znz(i,j,k)
              fA(1)         =   xA(i,j,k)  ! face area 
              fA(3)         =   yA(i,j,k)  ! face area 
              fA(5)         =   zA(i,j,k)  ! face area 
              mm(1)         =  (mu(i,j,k) + mu(i-1,j,k))*0.5    ! molecular viscosity
              mm(3)         =  (mu(i,j,k) + mu(i,j-1,k))*0.5    ! molecular viscosity
              mm(5)         =  (mu(i,j,k) + mu(i,j,k-1))*0.5    ! molecular viscosity
              tm(1)         = (mu_t(i,j,k) + mu_t(i-1,j,k))*0.5 ! turbulenct viscoisty
              tm(3)         = (mu_t(i,j,k) + mu_t(i,j-1,k))*0.5 ! turbulenct viscoisty
              tm(5)         = (mu_t(i,j,k) + mu_t(i,j,k-1))*0.5 ! turbulenct viscoisty
              Vl(1)          =  (volume(i,j,k) + volume(i-1,j,k))*0.5 !volume at face
              Vl(3)          =  (volume(i,j,k) + volume(i,j-1,k))*0.5 !volume at face
              Vl(5)          =  (volume(i,j,k) + volume(i,j,k-1))*0.5 !volume at face
              F1(1)         = (sst_F1(i,j,k) + sst_F1(i-1,j,k))*0.5
              F1(3)         = (sst_F1(i,j,k) + sst_F1(i,j-1,k))*0.5
              F1(5)         = (sst_F1(i,j,k) + sst_F1(i,j,k-1))*0.5
              State0        =       qp(i  ,j,k,1:n_var)
              state         =       qp(i-1,j,k,1:n_var)
              stateChange   = delQstar(i-1,j,k,1:n_var)

              NewIminusFlux = SSTflux(state, stateChange, fn(1,1), fn(1,2), fn(1,3), fA(1), mm(1), tm(1), F1(1), Vl(1), state0)
              stateChange   = 0.0
              OldIminusFlux = SSTflux(state, stateChange, fn(1,1), fn(1,2), fn(1,3), fA(1), mm(1), tm(1), F1(1), Vl(1), state0)
       LambdaTimesArea(1)=SSTSpectralRadius(State0,state, fn(1,1), fn(1,2), fn(1,3), fA(1), mm(1), tm(1), CellCenter(i-1,j,k,:), CellCenter(i,j,k,:))

              state         =       qp(i,j-1,k,1:n_var)
              stateChange   = delQstar(i,j-1,k,1:n_var)
              NewJminusFlux = SSTflux(state, stateChange, fn(3,1), fn(3,2), fn(3,3), fA(3), mm(3), tm(3), F1(3), Vl(3), state0)
              stateChange   = 0.0
              OldJminusFlux = SSTflux(state, stateChange, fn(3,1), fn(3,2), fn(3,3), fA(3), mm(3), tm(3), F1(3), Vl(3), state0)
       LambdaTimesArea(2)=SSTSpectralRadius(State0,state, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k)&
                        , (mu(i,j-1,k) + mu(i,j,k)), (mu_t(i,j-1,k) + mu_t(i,j,k)), CellCenter(i,j-1,k,:), CellCenter(i,j,k,:))

              state         =       qp(i,j,k-1,1:n_var)
              stateChange   = delQstar(i,j,k-1,1:n_var)
              NewKminusFlux = SSTflux(state, stateChange, fn(5,1), fn(5,2), fn(5,3), fA(5), mm(5), tm(5), F1(5), Vl(5), state0)
              stateChange   = 0.0
              OldKminusFlux = SSTflux(state, stateChange, fn(5,1), fn(5,2), fn(5,3), fA(5), mm(5), tm(5), F1(5), Vl(5), state0)
       LambdaTimesArea(3)=SSTSpectralRadius(State0,state, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k)&
                       , (mu(i,j,k-1) + mu(i,j,k)) , (mu_t(i,j,k-1) + mu_t(i,j,k)), CellCenter(i,j,k-1,:), CellCenter(i,j,k,:))

              state         = qp(i+1,j,k,1:n_var)
       LambdaTimesArea(4)=SSTSpectralRadius(State0,state, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k)&
                      , (mu(i+1,j,k) + mu(i,j,k)), (mu_t(i+1,j,k) + mu_t(i,j,k)), CellCenter(i+1,j,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j+1,k,1:n_var)
       LambdaTimesArea(5)=SSTSpectralRadius(State0,state, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k)&
                      , (mu(i,j+1,k) + mu(i,j,k)), (mu_t(i,j+1,k) + mu_t(i,j,k)), CellCenter(i,j+1,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j,k+1,1:n_var)
       LambdaTimesArea(6)=SSTSpectralRadius(State0,state, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1)&
                      , (mu(i,j,k+1) + mu(i,j,k)), (mu_t(i,j,k+1) + mu_t(i,j,k)), CellCenter(i,j,k+1,:), CellCenter(i,j,k,:))


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux


              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
              D(6) = D(6) + bstar*qp(i,j,k,1)*qp(i,j,k,6)*volume(i,j,k)
              D(7) = D(7) + beta*qp(i,j,k,1)*qp(i,j,k,7)*volume(i,j,k)

              deltaU(1:n_var) = -residue(i,j,k,1:n_var) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQstar(i-1,j,k,1:n_var)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQstar(i,j-1,k,1:n_var)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQstar(i,j,k-1,1:n_var)) )

              delQstar(i,j,k,1:n_var) = deltaU(1:n_var)/D
            end do
          end do
        end do

        delQ=0.0
        !backward sweep
            do i=imx-1,1,-1
          do j=jmx-1,1,-1
        do k=kmx-1,1,-1
              fn(2,1)       =  xnx(i+1,j,k)
              fn(2,2)       =  xny(i+1,j,k)
              fn(2,3)       =  xnz(i+1,j,k)
              fn(4,1)       =  ynx(i,j+1,k)
              fn(4,2)       =  yny(i,j+1,k)
              fn(4,3)       =  ynz(i,j+1,k)
              fn(6,1)       =  znx(i,j,k+1)
              fn(6,2)       =  zny(i,j,k+1)
              fn(6,3)       =  znz(i,j,k+1)
              fA(2)         =   xA(i+1,j,k)  ! face area 
              fA(4)         =   yA(i,j+1,k)  ! face area 
              fA(6)         =   zA(i,j,k+1)  ! face area 
              mm(2)         =  (mu(i,j,k) + mu(i+1,j,k))*0.5    ! molecular viscosity
              mm(4)         =  (mu(i,j,k) + mu(i,j+1,k))*0.5    ! molecular viscosity
              mm(6)         =  (mu(i,j,k) + mu(i,j,k+1))*0.5    ! molecular viscosity
              tm(2)         = (mu_t(i,j,k) + mu_t(i+1,j,k))*0.5 ! turbulenct viscoisty
              tm(4)         = (mu_t(i,j,k) + mu_t(i,j+1,k))*0.5 ! turbulenct viscoisty
              tm(6)         = (mu_t(i,j,k) + mu_t(i,j,k+1))*0.5 ! turbulenct viscoisty
              Vl(2)          =  (volume(i,j,k) + volume(i+1,j,k))*0.5 !volume at face
              Vl(4)          =  (volume(i,j,k) + volume(i,j+1,k))*0.5 !volume at face
              Vl(6)          =  (volume(i,j,k) + volume(i,j,k+1))*0.5 !volume at face
              F1(2)         = (sst_F1(i,j,k) + sst_F1(i+1,j,k))*0.5
              F1(4)         = (sst_F1(i,j,k) + sst_F1(i,j+1,k))*0.5
              F1(6)         = (sst_F1(i,j,k) + sst_F1(i,j,k+1))*0.5
              State0        =   qp(i  ,j,k,1:n_var)
              state         =   qp(i+1,j,k,1:n_var)
              stateChange   = delQ(i+1,j,k,1:n_var)
              NewIminusFlux = SSTflux(state, stateChange, fn(2,1), fn(2,2), fn(2,3), fA(2), mm(2), tm(2), F1(2), Vl(2), state0)
              stateChange   = 0.0
              OldIminusFlux = SSTflux(state, stateChange, fn(2,1), fn(2,2), fn(2,3), fA(2), mm(2), tm(2), F1(2), Vl(2), state0)
       LambdaTimesArea(1)=SSTSpectralRadius(State0,state, xnx(i+1,j,k), xny(i+1,j,k), xnz(i+1,j,k), xA(i+1,j,k)&
                          , (mu(i+1,j,k) + mu(i,j,k)), (mu_t(i+1,j,k) + mu_t(i,j,k)), CellCenter(i+1,j,k,:), CellCenter(i,j,k,:))

              state         =   qp(i,j+1,k,1:n_var)
              stateChange   = delQ(i,j+1,k,1:n_var)
              NewJminusFlux = SSTflux(state, stateChange, fn(4,1), fn(4,2), fn(4,3), fA(4), mm(4), tm(4), F1(4), Vl(4), state0)
              stateChange   = 0.0
              OldJminusFlux = SSTflux(state, stateChange, fn(4,1), fn(4,2), fn(4,3), fA(4), mm(4), tm(4), F1(4), Vl(4), state0)
       LambdaTimesArea(2)=SSTSpectralRadius(State0,state, ynx(i,j+1,k), yny(i,j+1,k), ynz(i,j+1,k), yA(i,j+1,k)&
                          , (mu(i,j+1,k) + mu(i,j,k)), (mu_t(i,j+1,k) + mu_t(i,j,k)), CellCenter(i,j+1,k,:), CellCenter(i,j,k,:))

              state         =   qp(i,j,k+1,1:n_var)
              stateChange   = delQ(i,j,k+1,1:n_var)
              NewKminusFlux = SSTflux(state, stateChange, fn(6,1), fn(6,2), fn(6,3), fA(6), mm(6), tm(6), F1(6), Vl(6), state0)
              stateChange   = 0.0
              OldKminusFlux = SSTflux(state, stateChange,  fn(6,1), fn(6,2), fn(6,3), fA(6), mm(6), tm(6), F1(6), Vl(6), state0)
       LambdaTimesArea(3)=SSTSpectralRadius(State0,state, znx(i,j,k+1), zny(i,j,k+1), znz(i,j,k+1), zA(i,j,k+1)&
                       , (mu(i,j,k+1) + mu(i,j,k)) , (mu_t(i,j,k+1) + mu_t(i,j,k)), CellCenter(i,j,k+1,:), CellCenter(i,j,k,:))

              state         = qp(i-1,j,k,1:n_var)
       LambdaTimesArea(4)=SSTSpectralRadius(State0,state, -xnx(i,j,k), -xny(i,j,k), -xnz(i,j,k), xA(i,j,k)&
                       , (mu(i-1,j,k) + mu(i,j,k)),(mu_t(i-1,j,k) + mu_t(i,j,k)), CellCenter(i-1,j,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j-1,k,1:n_var)
       LambdaTimesArea(5)=SSTSpectralRadius(State0,state, -ynx(i,j,k), -yny(i,j,k), -ynz(i,j,k), yA(i,j,k)&
                       , (mu(i,j-1,k) + mu(i,j,k)), (mu_t(i,j-1,k) + mu_t(i,j,k)), CellCenter(i,j-1,k,:), CellCenter(i,j,k,:))
              state         = qp(i,j,k-1,1:n_var)
       LambdaTimesArea(6)=SSTSpectralRadius(State0,state, -znx(i,j,k), -zny(i,j,k), -znz(i,j,k), zA(i,j,k)&
                       , (mu(i,j,k-1) + mu(i,j,k)), (mu_t(i,j,k-1) + mu_t(i,j,k)), CellCenter(i,j,k-1,:), CellCenter(i,j,k,:))


              ! multiply above flux with area to get correct values
              DelIminusFlux =  NewIminusFlux - OldIminusFlux
              DelJminusFlux =  NewJminusFlux - OldJminusFlux
              DelKminusFlux =  NewKminusFlux - OldKminusFlux

              D = (volume(i,j,k)/delta_t(i,j,k)) + 0.5*SUM(LambdaTimesArea)
              beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
              D(6) = D(6) + bstar*qp(i,j,k,1)*qp(i,j,k,6)*volume(i,j,k)
              D(7) = D(7) + beta*qp(i,j,k,1)*qp(i,j,k,7)*volume(i,j,k)

              delQ(i,j,k,1:n_var) = delQstar(i,j,k,1:n_var) &
                - 0.5*((DelIminusFlux - LambdaTimesArea(1)*delQ(i+1,j,k,1:n_var)) &
                     + (DelJminusFlux - LambdaTimesArea(2)*delQ(i,j+1,k,1:n_var)) &
                     + (DelKminusFlux - LambdaTimesArea(3)*delQ(i,j,k+1,1:n_var)) )/D

            end do
          end do
        end do
        
        do k=1,kmx-1
          do j = 1,jmx-1
            do i = 1,imx-1
              conservativeQ(1) = qp(i,j,k,1)
              conservativeQ(2) = qp(i,j,k,1) * qp(i,j,k,2)
              conservativeQ(3) = qp(i,j,k,1) * qp(i,j,k,3)
              conservativeQ(4) = qp(i,j,k,1) * qp(i,j,k,4)
              conservativeQ(5) = (qp(i,j,k,5) / (gm-1.0)) + ( 0.5 * qp(i,j,k,1) * sum( qp(i,j,k,2:4)**2) )
              conservativeQ(6) = qp(i,j,k,1) * qp(i,j,k,6)
              conservativeQ(7) = qp(i,j,k,1) * qp(i,j,k,7)
              
              ! add new change into conservative solution
              conservativeQ(1:n_var) = conservativeQ(1:n_var) + delQ(i,j,k,1:n_var)

              ! convert back conservative to primitive
              qp(i,j,k,1) = conservativeQ(1)
              qp(i,j,k,2) = conservativeQ(2) / conservativeQ(1)
              qp(i,j,k,3) = conservativeQ(3) / conservativeQ(1)
              qp(i,j,k,4) = conservativeQ(4) / conservativeQ(1)
              qp(i,j,k,5) = (gm-1.0) * ( conservativeQ(5) - (0.5 * sum(conservativeQ(2:4)**2) / conservativeQ(1)) )
              qp(i,j,k,6) = conservativeQ(6) / conservativeQ(1)
              qp(i,j,k,7) = conservativeQ(7) / conservativeQ(1)
            end do
          end do
        end do
      end subroutine matrix_free_implicit_update_SST

      function SSTflux(q,du,nx, ny, nz, Area, mmu, tmu, F1, volume, q2)
        implicit none
        real                    , intent(in) :: Area
        real                    , intent(in) :: nx
        real                    , intent(in) :: ny
        real                    , intent(in) :: nz
        real                    , intent(in) :: volume
        real                    , intent(in) :: mmu
        real                    , intent(in) :: tmu
        real                    , intent(in) :: F1
        real, dimension(1:n_var), intent(in) :: q
        real, dimension(1:n_var), intent(in) :: q2
        real, dimension(1:n_var), intent(in) :: du
        real, dimension(1:n_var)             :: SSTflux
        real, dimension(1:n_var)             :: U ! conservative variables
        real, dimension(1:n_var)             :: w ! new primitive variables

        real    :: dudx
        real    :: dudy
        real    :: dudz
        real    :: dvdx
        real    :: dvdy
        real    :: dvdz
        real    :: dwdx
        real    :: dwdy
        real    :: dwdz
        real    :: dTdx
        real    :: dTdy
        real    :: dTdz
        real    :: dtkdx
        real    :: dtkdy
        real    :: dtkdz
        real    :: dtwdx
        real    :: dtwdy
        real    :: dtwdz
        real    :: T1, T2
        real    :: uface
        real    :: vface
        real    :: wface
        real    :: trace
        real    :: Tauxx
        real    :: Tauyy
        real    :: Tauzz
        real    :: Tauxy
        real    :: Tauxz
        real    :: Tauyz
        real    :: Qx
        real    :: Qy
        real    :: Qz
        real    :: HalfRhoUsquare
        real    :: RhoHt
        real    :: K_heat
        real    :: FaceNormalVelocity
        real    :: mu
        real    :: sigma_k
        real    :: sigma_w


        ! find conservative variable
        U(1) = q(1)
        U(2) = q(1) * q(2)
        U(3) = q(1) * q(3)
        U(4) = q(1) * q(4)
        U(5) = ( q(5) / (gm-1.0) ) + ( 0.5 * q(1) * sum(q(2:4)**2)  )
        U(6) = q(1) * q(6)
        U(7) = q(1) * q(7)

        U(1:n_var) = U(1:n_var) + du(1:n_var)
        

        W(1) = U(1)
        W(2) = U(2) / U(1)
        W(3) = U(3) / U(1)
        W(4) = U(4) / U(1)
        W(5) = (gm-1.0) * ( U(5) - ( 0.5 * SUM(U(2:4)**2) / U(1) ) )
        W(6) = U(6) / U(1)
        W(7) = U(7) / U(1)

        FaceNormalVelocity = (W(2) * nx) + (W(3) * ny) + (W(4) * nz)

        SSTflux(1) =   W(1) * FaceNormalVelocity

        SSTflux(2) = ( W(2) * SSTflux(1) ) + ( W(5) * nx )

        SSTflux(3) = ( W(3) * SSTflux(1) ) + ( W(5) * ny )

        SSTflux(4) = ( W(4) * SSTflux(1) ) + ( W(5) * nz )

        HalfRhoUsquare = 0.5 * W(1) * ( W(2)*W(2) + W(3)*W(3) + W(4)*W(4) )
        RhoHt          = ( (gm/(gm-1.0)) * W(5) ) + HalfRhoUsquare
        SSTflux(5)        = RhoHt * FaceNormalVelocity
        !SSTflux(5)        = ( U(5) + W(5) ) * FaceNormalVelocity

        SSTflux(6) = ( W(6) * SSTflux(1) )

        SSTflux(7) = ( W(7) * SSTflux(1) )


        ! viscous terms
        mu = mmu + tmu
        T1    =  W(5) / ( W(1) * R_gas )
        T2    =  q2(5) / ( q2(1) * R_gas )
        dTdx  =  ( T2    - T1   ) * nx * Area / Volume
        dTdy  =  ( T2    - T1   ) * ny * Area / Volume
        dTdz  =  ( T2    - T1   ) * nz * Area / Volume
        dudx  =  ( q2(2) - W(2) ) * nx * Area / Volume
        dudy  =  ( q2(2) - W(2) ) * ny * Area / Volume
        dudz  =  ( q2(2) - W(2) ) * nz * Area / Volume
        dvdx  =  ( q2(3) - W(3) ) * nx * Area / Volume
        dvdy  =  ( q2(3) - W(3) ) * ny * Area / Volume
        dvdz  =  ( q2(3) - W(3) ) * nz * Area / Volume
        dwdx  =  ( q2(4) - W(4) ) * nx * Area / Volume
        dwdy  =  ( q2(4) - W(4) ) * ny * Area / Volume
        dwdz  =  ( q2(4) - W(4) ) * nz * Area / Volume
        dtkdx  =  ( q2(6) - W(6) ) * nx * Area / Volume
        dtkdy  =  ( q2(6) - W(6) ) * ny * Area / Volume
        dtkdz  =  ( q2(6) - W(6) ) * nz * Area / Volume
        dtwdx  =  ( q2(7) - W(7) ) * nx * Area / Volume
        dtwdy  =  ( q2(7) - W(7) ) * ny * Area / Volume
        dtwdz  =  ( q2(7) - W(7) ) * nz * Area / Volume

        trace = dudx + dvdy + dwdz
        Tauxx =  2. * mu * (dudx - trace/3.0)
        Tauyy =  2. * mu * (dvdy - trace/3.0)
        Tauzz =  2. * mu * (dwdz - trace/3.0)
        Tauxy = mu * (dvdx + dudy)
        Tauxz = mu * (dwdx + dudz)
        Tauyz = mu * (dwdy + dvdz)

        K_heat = ( mmu / Pr  + tmu/tpr) * gm * R_gas / ( gm - 1.0 )
        Qx = K_heat*dTdx
        Qy = K_heat*dTdy
        Qz = K_heat*dTdz

        SSTflux(2) = SSTflux(2) - ( Tauxx * nx + Tauxy * ny + Tauxz * nz )
        SSTflux(3) = SSTflux(3) - ( Tauxy * nx + Tauyy * ny + Tauyz * nz )
        SSTflux(4) = SSTflux(4) - ( Tauxz * nx + Tauyz * ny + Tauzz * nz )
        uface = 0.5 * ( W(2) + q2(2) )
        vface = 0.5 * ( W(3) + q2(3) )
        wface = 0.5 * ( W(4) + q2(4) )

        SSTflux(5) = SSTflux(5) - ( Tauxx * uface + Tauxy * vface + Tauxz * wface + Qx ) * nx
        SSTflux(5) = SSTflux(5) - ( Tauxy * uface + Tauyy * vface + Tauyz * wface + Qy ) * ny
        SSTflux(5) = SSTflux(5) - ( Tauxz * uface + Tauyz * vface + Tauzz * wface + Qz ) * nz

        sigma_k = sigma_k1*F1 + sigma_k2*(1.0 - F1)
        sigma_w = sigma_w1*F1 + sigma_w2*(1.0 - F1)

        SSTflux(6) = SSTflux(6) + (mmu + sigma_k*tmu)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)
        SSTflux(7) = SSTflux(7) + (mmu + sigma_w*tmu)*(dtwdx*nx + dtwdy*ny + dtwdz*nz)

        SSTflux    = SSTflux * Area

      end function SSTflux

      function SSTSpectralRadius(q1, q2, nx, ny, nz, Area, mm, tm, c1, c2)
        real, dimension(1:n_var), intent(in) :: q1
        real, dimension(1:n_var), intent(in) :: q2
        real                    , intent(in) :: nx
        real                    , intent(in) :: ny
        real                    , intent(in) :: nz
        real                    , intent(in) :: Area
        real                    , intent(in) :: mm
        real                    , intent(in) :: tm
        real, dimension(1:3)    , intent(in) :: c1
        real, dimension(1:3)    , intent(in) :: c2
        real                                 :: SSTSpectralRadius
        real                                 :: NormalSpeed
        real                                 :: SpeedOfSound
        real                                 :: vis
        real                                 :: mu

        ! in state vector q (2-4) are the cell center velocity
        NormalSpeed = 0.5 * ( ( ( q1(2) + q2(2) ) * nx ) &
                            + ( ( q1(3) + q2(3) ) * ny ) &
                            + ( ( q1(4) + q2(4) ) * nz ) &
                            )
        NormalSpeed = abs(NormalSpeed)

        !SpeedOfSound = sqrt(gm*(q1(5) + q2(5))/(q1(1) + q2(1)))
        SpeedOfSound = 0.5*( sqrt(gm*q1(5)/q1(1)) + sqrt(gm*q2(5)/q2(1)) )

        ! visocus part
       ! vis =  max(4./3.) * gm * mu / ( Pr * 0.5*(q1(1) + q2(1)) * abs((c1(1)-c2(1))*nx + (c1(2)-c2(2))*ny +(c1(3)-c2(3))*nz) )
       mu = mm/Pr + tm/tPr
        vis =  max(4./3.,gm) * mu / ( Pr * 0.5*(q1(1) + q2(1)) &
          * sqrt((c1(1)-c2(1))**2 + (c1(2)-c2(2))**2 +(c1(3)-c2(3))**2) )
        SSTSpectralRadius = ( NormalSpeed + SpeedOfSound + vis) * Area

      end function SSTSpectralRadius


end module update
