module scheme

#include "../../../debug.h"
#include "../../../error.h"
    use vartypes
    use global, only: SCHEME_NAME_LENGTH

!    use global_vars, only : imx
!    use global_vars, only : jmx
!    use global_vars, only : kmx

    use global_vars, only : mu_ref
!    use global_vars, only : n_var
    use global_vars, only : F_p
    use global_vars, only : G_p
    use global_vars, only : H_p
    use global_vars, only : mass_residue
    use global_vars, only : x_mom_residue
    use global_vars, only : y_mom_residue
    use global_vars, only : z_mom_residue
    use global_vars, only : energy_residue
    use global_vars, only : TKE_residue
    use global_vars, only : omega_residue
    use global_vars, only : KL_residue
    use global_vars, only : dissipation_residue
    use global_vars, only : tv_residue
    use global_vars, only : residue
    use global_vars, only : turbulence

    use global_vars, only : scheme_name

    use utils, only: alloc, dealloc
    use face_interpolant, only: setup_interpolant_scheme, &
            destroy_interpolant_scheme
    use van_leer, only: &
            setup_scheme_van_leer => setup_scheme, &
            destroy_scheme_van_leer => destroy_scheme, &
            compute_fluxes_van_leer => compute_fluxes, &
            get_residue_van_leer => get_residue, &
            F_van_leer => F, &
            G_van_leer => G, &
            H_van_leer => H, &
            residue_van_leer => residue
    use ausm, only: &
            setup_scheme_ausm => setup_scheme, &
            destroy_scheme_ausm => destroy_scheme, &
            compute_fluxes_ausm => compute_fluxes, &
            get_residue_ausm => get_residue, &
            F_ausm => F, &
            G_ausm => G, &
            H_ausm => H, &
            residue_ausm => residue
    use ausmP, only: &
            setup_scheme_ausmP => setup_scheme, &
            destroy_scheme_ausmP => destroy_scheme, &
            compute_fluxes_ausmP => compute_fluxes, &
            get_residue_ausmP => get_residue, &
            F_ausmP => F, &
            G_ausmP => G, &
            H_ausmP => H, &
            residue_ausmP => residue
    use ausmUP, only: &
            setup_scheme_ausmUP => setup_scheme, &
            destroy_scheme_ausmUP => destroy_scheme, &
            compute_fluxes_ausmUP => compute_fluxes, &
            get_residue_ausmUP => get_residue, &
            F_ausmUP => F, &
            G_ausmUP => G, &
            H_ausmUP => H, &
            residue_ausmUP => residue
    use slau, only: &
            setup_scheme_slau => setup_scheme, &
            destroy_scheme_slau => destroy_scheme, &
            compute_fluxes_slau => compute_fluxes, &
            get_residue_slau => get_residue, &
            F_slau => F, &
            G_slau => G, &
            H_slau => H, &
            residue_slau => residue
    use ldfss0, only: &
            setup_scheme_ldfss0=> setup_scheme, &
            destroy_scheme_ldfss0=> destroy_scheme, &
            compute_fluxes_ldfss0 => compute_fluxes, &
            get_residue_ldfss0 => get_residue, &
            F_ldfss0 => F, &
            G_ldfss0 => G, &
            H_ldfss0 => H, &
            residue_ldfss0 => residue
!   use hlle, only: &
!           setup_scheme_hlle => setup_scheme, &
!           destroy_scheme_hlle => destroy_scheme, &
!           get_residue_hlle => get_residue
    include "turbulence_models/include/scheme/import_module.inc"


    implicit none
    integer :: imx, jmx, kmx, n_var
    private

    include "turbulence_models/include/scheme/variable_deceleration.inc" 

    type :: residual
      real  ::  residue
    end type
    type :: fluxLamComp
      real ::  mass
      real ::  xmom
      real ::  ymom
      real ::  zmom
      real ::  energy
    end type fluxLamComp

    type, extends (fluxLamComp) :: fluxkwComp
      real :: tke
      real :: omega
    end type fluxkwComp

    type, extends (fluxLamComp) :: fluxSAComp
      real :: tv
    end type fluxSAComp

    type, extends (fluxkwComp) :: fluxkwGammaComp
      real :: tgm
    end type fluxkwGammaComp
      


    ! Public members
    public :: scheme_name
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: compute_residue
    public :: residue

    contains

        subroutine setup_scheme(control, dims)
            implicit none
            type(controltype), intent(in) :: control
            type(extent), intent(in) :: dims

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            n_var = control%n_var

            call setup_interpolant_scheme(control, dims)

            select case (scheme_name)
                case ("van_leer")
                    call setup_scheme_van_leer(control, dims)
                    F_p => F_van_leer
                    G_p => G_van_leer
                    H_p => H_van_leer
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 5)
                    select case (turbulence)

                        case ("none")
                            !include nothing
                            continue
                        case ("sst", "sst2003", "kw" , "des-sst")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 7)
                        case ("kkl")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 7)
                        case ("sa", "saBC")
                          tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
                        case ("ke")
                                  TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
                          dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 7)
                        case DEFAULT
                            Fatal_error

                    end select
                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_van_leer(:,:,:,:)
                case ("ausm")
                    call setup_scheme_ausm(control, dims)
                    F_p => F_ausm
                    G_p => G_ausm
                    H_p => H_ausm
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 5)
                    select case (turbulence)
                        case ("none")
                            !include nothing
                            continue
                        case ("sst", "sst2003", "kw" , "des-sst")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 7)
                        case ("kkl")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 7)
                        case ("sa", "saBC")
                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
                        case ("ke")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 7)
                        case DEFAULT
                            Fatal_error
                    end select
                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ausm(:,:,:,:)
                case ("ausmP")
                    call setup_scheme_ausmP(control, dims)
                    F_p => F_ausmP
                    G_p => G_ausmP
                    H_p => H_ausmP
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 5)
                    select case (turbulence)
                        case ("none")
                            !include nothing
                            continue
                        case ("sst", "sst2003", "kw" , "des-sst")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 7)
                        case ("kkl")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 7)
                        case ("sa", "saBC")
                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
                        case ("ke")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 7)
                        case DEFAULT
                            Fatal_error
                    end select
                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ausmP(:,:,:,:)
                case ("ausmUP")
                    call setup_scheme_ausmUP(control, dims)
                    F_p => F_ausmUP
                    G_p => G_ausmUP
                    H_p => H_ausmUP
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 5)
                    select case (turbulence)
                        case ("none")
                            !include nothing
                            continue
                        case ("sst", "sst2003", "kw" , "des-sst")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 7)
                        case ("kkl")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 7)
                        case ("sa", "saBC")
                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
                        case ("ke")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 7)
                        case DEFAULT
                            Fatal_error
                    end select
                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ausmUP(:,:,:,:)
                case ("slau")
                    call setup_scheme_slau(control, dims)
                    F_p => F_slau
                    G_p => G_slau
                    H_p => H_slau
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 5)
                    select case (turbulence)
                        case ("none")
                            !include nothing
                            continue
                        case ("sst", "sst2003", "kw" , "des-sst")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 7)
                        case ("kkl")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 7)
                        case ("sa", "saBC")
                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
                        case ("ke")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 7)
                        case DEFAULT
                            Fatal_error
                    end select
                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_slau(:,:,:,:)
                case ("ldfss0")
                    call setup_scheme_ldfss0(control, dims)
                    F_p => F_ldfss0
                    G_p => G_ldfss0
                    H_p => H_ldfss0
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 5)
                    select case (turbulence)
                        case ("none")
                            !include nothing
                            continue
                        case ("sst", "sst2003", "kw" , "des-sst")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 7)
                        case ("kkl")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 7)
                        case ("sa", "saBC")
                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
                        case ("ke")
                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 7)
                        case DEFAULT
                            Fatal_error
                    end select
                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ldfss0(:,:,:,:)
!               case ("hlle")
!                   call setup_scheme_hlle()
                case default
                    Fatal_error
            end select

        end subroutine setup_scheme

        subroutine deallocate_memory()

            implicit none

            nullify(F_p)
            nullify(G_p)
            nullify(H_p)
            nullify(mass_residue)
            nullify(x_mom_residue)
            nullify(y_mom_residue)
            nullify(z_mom_residue)
            nullify(energy_residue)
            nullify(residue)
            nullify(TKE_residue)
            nullify(omega_residue)
            nullify(KL_residue)
            nullify(tv_residue)
            nullify(dissipation_residue)

        end subroutine deallocate_memory

        subroutine destroy_scheme()

            implicit none

            select case (scheme_name)
                case ("van_leer")
                    call destroy_scheme_van_leer
                case ("ausm")
                    call destroy_scheme_ausm()
                case ("ausmP")
                    call destroy_scheme_ausmP()
                case ("ausmUP")
                    call destroy_scheme_ausmUP()
                case ("slau")
                    call destroy_scheme_slau()
                case ("ldfss0")
                    call destroy_scheme_ldfss0()
!               case ("hlle")
!                   call destroy_scheme_hlle()
                case default
                    Fatal_error
            end select
            
            call destroy_interpolant_scheme()
            call deallocate_memory()

        end subroutine destroy_scheme

        subroutine compute_fluxes
        
            implicit none

            select case (scheme_name)
                case ("van_leer")
                    call compute_fluxes_van_leer()
                case ("ausm")
                    call compute_fluxes_ausm()
                case ("ausmP")
                    call compute_fluxes_ausmP()
                case ("ausmUP")
                    call compute_fluxes_ausmUP()
                case ("slau")
                    call compute_fluxes_slau()
                case ("ldfss0")
                    call compute_fluxes_ldfss0()
!               case ("hlle")
!                   call compute_fluxes_hlle()
                case default
                    Fatal_error
            end select
            
        end subroutine compute_fluxes

        subroutine compute_residue()
            
            implicit none
            
            select case (scheme_name)
                case ("van_leer")
                    call get_residue_van_leer()
                case ("ausm")
                    call get_residue_ausm()
                case ("ausmP")
                    call get_residue_ausmP()
                case ("ausmUP")
                    call get_residue_ausmUP()
                case ("slau")
                    call get_residue_slau()
                case ("ldfss0")
                    call get_residue_ldfss0()
!               case ("hlle")
!                   call get_residue_hlle()
                case default
                    Fatal_error
            end select
        
        end subroutine compute_residue

end module scheme
