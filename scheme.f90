module scheme

    use global, only: SCHEME_NAME_LENGTH
    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx
    use state, only: n_var, mu_ref
    use face_interpolant, only: setup_interpolant_scheme, &
            destroy_interpolant_scheme
    use van_leer, only: &
            setup_scheme_van_leer => setup_scheme, &
            destroy_scheme_van_leer => destroy_scheme, &
            compute_fluxes_van_leer => compute_fluxes, &
            get_residue_van_leer => get_residue, &
            F_van_leer => F, &
            G_van_leer => G, &
            H_van_leer => H
    use viscous, only: setup_viscous, destroy_viscous
!   use ausm, only: &
!           setup_scheme_ausm => setup_scheme, &
!           destroy_scheme_ausm => destroy_scheme, &
!           get_residue_ausm => get_residue
!   use ldfss0, only: &
!           setup_scheme_ldfss0 => setup_scheme, &
!           destroy_scheme_ldfss0 => destroy_scheme, &
!           get_residue_ldfss0 => get_residue
!   use hlle, only: &
!           setup_scheme_hlle => setup_scheme, &
!           destroy_scheme_hlle => destroy_scheme, &
!           get_residue_hlle => get_residue

    implicit none
    private

    character(len=SCHEME_NAME_LENGTH) :: scheme_name
    real, public, dimension(:, :, :, :), allocatable, target :: residue
    real, public, dimension(:, :, :, :), pointer :: F_p, G_p, H_p

    ! Public members
    public :: scheme_name
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: compute_residue

    contains

        subroutine allocate_memory()
            
            implicit none
            
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for residue.')

        end subroutine allocate_memory

        subroutine setup_scheme()
            implicit none

            call allocate_memory()

            call setup_interpolant_scheme()

            select case (scheme_name)
                case ("van_leer")
                    call setup_scheme_van_leer()
                    F_p => F_van_leer
                    G_p => G_van_leer
                    H_p => H_van_leer
!               case ("ausm")
!                   call setup_scheme_ausm()
!               case ("ldfss0")
!                   call setup_scheme_ldfss0()
!               case ("hlle")
!                   call setup_scheme_hlle()
                case default
                    call dmsg(5, 'scheme', 'setup_scheme', &
                            'Scheme not recognized.')
                    stop
            end select

            if (mu_ref /= 0.0) then
                call setup_viscous()
            end if

        end subroutine setup_scheme

        subroutine deallocate_memory()
            implicit none
            call dealloc(residue)
            nullify(F_p)
            nullify(G_p)
            nullify(H_p)
        end subroutine deallocate_memory

        subroutine destroy_scheme()
            implicit none

            select case (scheme_name)
                case ("van_leer")
                    call destroy_scheme_van_leer
!               case ("ausm")
!                   call destroy_scheme_ausm()
!               case ("ldfss0")
!                   call destroy_scheme_ldfss0()
!               case ("hlle")
!                   call destroy_scheme_hlle()
                case default
                    call dmsg(5, 'scheme', 'destroy_scheme', &
                            'Scheme not recognized.')
                    stop
            end select
            
            if (mu_ref /= 0.0) then
                call destroy_viscous()
            end if

            call destroy_interpolant_scheme()
            call deallocate_memory()

        end subroutine destroy_scheme

        subroutine compute_fluxes
        
            implicit none

            select case (scheme_name)
                case ("van_leer")
                    call compute_fluxes_van_leer()
                case default
                    call dmsg(5, 'scheme', 'compute_fluxes', &
                            'Scheme not recognized.')
                    stop
            end select
            
        end subroutine compute_fluxes

        subroutine compute_residue()
            
            implicit none
            
            select case (scheme_name)
                case ("van_leer")
                    residue = get_residue_van_leer()
!               case ("ausm")
!                   residue = get_residue_ausm()
!               case ("ldfss0")
!                   residue = get_residue_ldfss0()
!               case ("hlle")
!                   residue = get_residue_hlle()
                case default
                    call dmsg(5, 'scheme', 'compute_residue', &
                            'Scheme not recognized.')
                    stop
            end select
        
        end subroutine compute_residue

end module scheme
