module scheme

#include "../../../debug.h"
#include "../../../error.h"
    use vartypes
!    use global_vars, only : F_p
!    use global_vars, only : G_p
!    use global_vars, only : H_p
!    use global_vars, only : mass_residue
!    use global_vars, only : x_mom_residue
!    use global_vars, only : y_mom_residue
!    use global_vars, only : z_mom_residue
!    use global_vars, only : energy_residue
!    use global_vars, only : TKE_residue
!    use global_vars, only : omega_residue
!    use global_vars, only : KL_residue
!    use global_vars, only : dissipation_residue
!    use global_vars, only : tv_residue
!    use global_vars, only : residue
    use utils, only: alloc
    use face_interpolant, only: setup_interpolant_scheme
!    use van_leer, only: &
!            setup_scheme_van_leer => setup_scheme, &
!            compute_fluxes_van_leer => compute_fluxes, &
!            get_residue_van_leer => get_residue, &
!            F_van_leer => F, &
!            G_van_leer => G, &
!            H_van_leer => H, &
!            residue_van_leer => residue
!            !destroy_scheme_van_leer => destroy_scheme, &
!    use ausm, only: &
!            setup_scheme_ausm => setup_scheme, &
!            compute_fluxes_ausm => compute_fluxes, &
!            get_residue_ausm => get_residue, &
!            F_ausm => F, &
!            G_ausm => G, &
!            H_ausm => H, &
!            residue_ausm => residue
!            !destroy_scheme_ausm => destroy_scheme, &
!    use ausmP, only: &
!            setup_scheme_ausmP => setup_scheme, &
!            compute_fluxes_ausmP => compute_fluxes, &
!            get_residue_ausmP => get_residue, &
!            F_ausmP => F, &
!            G_ausmP => G, &
!            H_ausmP => H, &
!            residue_ausmP => residue
!            !destroy_scheme_ausmP => destroy_scheme, &
!    use ausmUP, only: &
!            setup_scheme_ausmUP => setup_scheme, &
!            compute_fluxes_ausmUP => compute_fluxes, &
!            get_residue_ausmUP => get_residue, &
!            F_ausmUP => F, &
!            G_ausmUP => G, &
!            H_ausmUP => H, &
!            residue_ausmUP => residue
!            !destroy_scheme_ausmUP => destroy_scheme, &
!    use slau, only: &
!            setup_scheme_slau => setup_scheme, &
!            compute_fluxes_slau => compute_fluxes, &
!            get_residue_slau => get_residue, &
!            F_slau => F, &
!            G_slau => G, &
!            H_slau => H, &
!            residue_slau => residue
!            !destroy_scheme_slau => destroy_scheme, &
!    use ldfss0, only: &
!            setup_scheme_ldfss0=> setup_scheme, &
!            compute_fluxes_ldfss0 => compute_fluxes, &
!            get_residue_ldfss0 => get_residue, &
!            F_ldfss0 => F, &
!            G_ldfss0 => G, &
!            H_ldfss0 => H, &
!            residue_ldfss0 => residue
!            !destroy_scheme_ldfss0=> destroy_scheme, &
!
    use ausmP,    only: compute_fluxes_ausmP  => compute_fluxes
    use ausmUP,   only: compute_fluxes_ausmUP => compute_fluxes
    use slau,     only: compute_fluxes_slau   => compute_fluxes
    use ausm,     only: compute_fluxes_ausm   => compute_fluxes
    use ldfss0,   only: compute_fluxes_ldfss0   => compute_fluxes
    use van_leer, only: compute_fluxes_van_leer => compute_fluxes

    implicit none
    integer :: imx, jmx, kmx, n_var
    private

    ! Public members
    public :: setup_scheme
!    public :: destroy_scheme
    public :: compute_fluxes
    public :: compute_residue
!    public :: residue

    contains

        subroutine setup_scheme(residue, F,G,H, control, scheme, dims)
            implicit none
            type(controltype), intent(in) :: control
            type(schemetype), intent(in) :: scheme
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real, dimension(:, :, :, :), allocatable, intent(out), target :: residue
            !< Store residue at each cell-center
            real, dimension(:, :, :, :), allocatable, intent(out) :: F
            !< Store fluxes throught the I faces
            real, dimension(:, :, :, :), allocatable, intent(out) :: G
            !< Store fluxes throught the J faces
            real, dimension(:, :, :, :), allocatable, intent(out) :: H
            !< Store fluxes throught the K faces

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            n_var = control%n_var

            call setup_interpolant_scheme(control, scheme, dims)

!            select case (scheme%scheme_name)
!                case ("van_leer")
!                    call setup_scheme_van_leer(control, dims)
!                    F_p => F_van_leer
!                    G_p => G_van_leer
!                    H_p => H_van_leer
!                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 1)
!                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 2)
!                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 3)
!                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 4)
!                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 5)
!                    select case (scheme%turbulence)
!
!                        case ("none")
!                            !include nothing
!                            continue
!                        case ("sst", "sst2003", "kw" , "des-sst")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
!                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 7)
!                        case ("kkl")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
!                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 7)
!                        case ("sa", "saBC")
!                          tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
!                        case ("ke")
!                                  TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 6)
!                          dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 7)
!                        case DEFAULT
!                            Fatal_error
!
!                    end select
!                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_van_leer(:,:,:,:)
!                case ("ausm")
!                    call setup_scheme_ausm(control, dims)
!                    F_p => F_ausm
!                    G_p => G_ausm
!                    H_p => H_ausm
!                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 1)
!                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 2)
!                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 3)
!                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 4)
!                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 5)
!                    select case (scheme%turbulence)
!                        case ("none")
!                            !include nothing
!                            continue
!                        case ("sst", "sst2003", "kw" , "des-sst")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
!                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 7)
!                        case ("kkl")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
!                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 7)
!                        case ("sa", "saBC")
!                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
!                        case ("ke")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 6)
!                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 7)
!                        case DEFAULT
!                            Fatal_error
!                    end select
!                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ausm(:,:,:,:)
!                case ("ausmP")
!                    call setup_scheme_ausmP(control, dims)
!                    F_p => F_ausmP
!                    G_p => G_ausmP
!                    H_p => H_ausmP
!                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 1)
!                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 2)
!                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 3)
!                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 4)
!                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 5)
!                    select case (scheme%turbulence)
!                        case ("none")
!                            !include nothing
!                            continue
!                        case ("sst", "sst2003", "kw" , "des-sst")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
!                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 7)
!                        case ("kkl")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
!                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 7)
!                        case ("sa", "saBC")
!                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
!                        case ("ke")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 6)
!                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmP(:, :, :, 7)
!                        case DEFAULT
!                            Fatal_error
!                    end select
!                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ausmP(:,:,:,:)
!                case ("ausmUP")
!                    call setup_scheme_ausmUP(control, dims)
!                    F_p => F_ausmUP
!                    G_p => G_ausmUP
!                    H_p => H_ausmUP
!                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 1)
!                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 2)
!                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 3)
!                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 4)
!                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 5)
!                    select case (scheme%turbulence)
!                        case ("none")
!                            !include nothing
!                            continue
!                        case ("sst", "sst2003", "kw" , "des-sst")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
!                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 7)
!                        case ("kkl")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
!                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 7)
!                        case ("sa", "saBC")
!                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
!                        case ("ke")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 6)
!                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausmUP(:, :, :, 7)
!                        case DEFAULT
!                            Fatal_error
!                    end select
!                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ausmUP(:,:,:,:)
!                case ("slau")
!                    call setup_scheme_slau(control, dims)
!                    F_p => F_slau
!                    G_p => G_slau
!                    H_p => H_slau
!                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 1)
!                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 2)
!                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 3)
!                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 4)
!                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 5)
!                    select case (scheme%turbulence)
!                        case ("none")
!                            !include nothing
!                            continue
!                        case ("sst", "sst2003", "kw" , "des-sst")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
!                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 7)
!                        case ("kkl")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
!                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 7)
!                        case ("sa", "saBC")
!                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
!                        case ("ke")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 6)
!                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_slau(:, :, :, 7)
!                        case DEFAULT
!                            Fatal_error
!                    end select
!                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_slau(:,:,:,:)
!                case ("ldfss0")
!                    call setup_scheme_ldfss0(control, dims)
!                    F_p => F_ldfss0
!                    G_p => G_ldfss0
!                    H_p => H_ldfss0
!                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 1)
!                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 2)
!                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 3)
!                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 4)
!                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 5)
!                    select case (scheme%turbulence)
!                        case ("none")
!                            !include nothing
!                            continue
!                        case ("sst", "sst2003", "kw" , "des-sst")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
!                        omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 7)
!                        case ("kkl")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
!                           KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 7)
!                        case ("sa", "saBC")
!                           tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
!                        case ("ke")
!                          TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 6)
!                  dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 7)
!                        case DEFAULT
!                            Fatal_error
!                    end select
!                    residue(1:imx-1,1:jmx-1,1:kmx-1,1:n_var)=>residue_ldfss0(:,:,:,:)
!!               case ("hlle")
!!                   call setup_scheme_hlle()
!                case default
!                    Fatal_error
!            end select
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
       !     mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 1)
       !     x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 2)
       !     y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 3)
       !     z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 4)
       !     energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 5)
       !     select case (scheme%turbulence)
       !         case ("none")
       !             !include nothing
       !             continue
       !         case ("sst", "sst2003", "kw" , "des-sst")
       !           TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 6)
       !         omega_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 7)
       !         case ("kkl")
       !           TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 6)
       !            KL_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 7)
       !         case ("sa", "saBC")
       !            tv_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 6)
       !         case ("ke")
       !           TKE_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 6)
       !   dissipation_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 7)
       !         case DEFAULT
       !             Fatal_error
       !     end select

        end subroutine setup_scheme

!        subroutine deallocate_memory()
!
!            implicit none
!
!            nullify(F_p)
!            nullify(G_p)
!            nullify(H_p)
!            nullify(mass_residue)
!            nullify(x_mom_residue)
!            nullify(y_mom_residue)
!            nullify(z_mom_residue)
!            nullify(energy_residue)
!            nullify(residue)
!            nullify(TKE_residue)
!            nullify(omega_residue)
!            nullify(KL_residue)
!            nullify(tv_residue)
!            nullify(dissipation_residue)
!
!        end subroutine deallocate_memory
!
!        subroutine destroy_scheme()
!
!            implicit none
!
!            select case (scheme%scheme_name)
!                case ("van_leer")
!                    call destroy_scheme_van_leer
!                case ("ausm")
!                    call destroy_scheme_ausm()
!                case ("ausmP")
!                    call destroy_scheme_ausmP()
!                case ("ausmUP")
!                    call destroy_scheme_ausmUP()
!                case ("slau")
!                    call destroy_scheme_slau()
!                case ("ldfss0")
!                    call destroy_scheme_ldfss0()
!!               case ("hlle")
!!                   call destroy_scheme_hlle()
!                case default
!                    Fatal_error
!            end select
!            
!            call destroy_interpolant_scheme()
!            call deallocate_memory()
!
!        end subroutine destroy_scheme

subroutine compute_fluxes(F,G,H, Ifaces, Jfaces, Kfaces, scheme, flow, dims)
        
            implicit none
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(in) :: flow
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
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

            select case (scheme%scheme_name)
                case ("van_leer")
                  call compute_fluxes_van_leer(F,G,H,Ifaces,Jfaces,Kfaces,flow,dims)
                case ("ausm")
                  call compute_fluxes_ausm(F,G,H,Ifaces,Jfaces,Kfaces,flow,dims)
                case ("ausmP")
                  call compute_fluxes_ausmP(F,G,H,Ifaces,Jfaces,Kfaces,flow,dims)
                case ("ausmUP")
                  call compute_fluxes_ausmUP(F,G,H,Ifaces,Jfaces,Kfaces,flow,dims)
                case ("slau")
                  call compute_fluxes_slau(F,G,H,Ifaces,Jfaces,Kfaces,flow,dims)
                case ("ldfss0")
                  call compute_fluxes_ldfss0(F,G,H,Ifaces,Jfaces,Kfaces,flow,dims)
                case default
                    Fatal_error
            end select
            
        end subroutine compute_fluxes

        subroutine compute_residue(residue,F,G,H,dims)
            
            implicit none
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real, dimension(:, :, :, :), intent(out)  :: residue
            !< Store residue at each cell-center
            real, dimension(:, :, :, :), intent(in) :: F
            !< Store fluxes throught the I faces
            real, dimension(:, :, :, :), intent(in) :: G
            !< Store fluxes throught the J faces
            real, dimension(:, :, :, :), intent(in) :: H
            !< Store fluxes throught the K faces
            
            integer :: i, j, k, l

            DebugCall('compute_residue')

            do l = 1, dims%n_var
             do k = 1, dims%kmx - 1
              do j = 1, dims%jmx - 1
               do i = 1, dims%imx - 1
               residue(i, j, k, l) = (F(i+1, j, k, l) - F(i, j, k, l)) &
                                   + (G(i, j+1, k, l) - G(i, j, k, l)) &
                                   + (H(i, j, k+1, l) - H(i, j, k, l))
               end do
              end do
             end do
            end do
!            select case (scheme%scheme_name)
!!                case ("van_leer")
!!                    call get_residue_van_leer()
!!                case ("ausm")
!!                    call get_residue_ausm()
!!                case ("ausmP")
!!                    call get_residue_ausmP()
!!                case ("ausmUP")
!!                    call get_residue_ausmUP()
!!                case ("slau")
!!                    call get_residue_slau()
!!                case ("ldfss0")
!!                    call get_residue_ldfss0()
!!                    Fatal_error
!            end select
        
        end subroutine compute_residue

end module scheme
