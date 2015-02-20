subroutine setup_scheme()
    implicit none
    select case (scheme)
        case ("van_leer")
            call setup_scheme_van_leer()
        case ("ldfss0")
            call setup_scheme_ldfss0
    end select
end subroutine setup_scheme

subroutine destroy_scheme()
    implicit none
    select case (scheme)
        case ("van_leer")
            call destroy_scheme_van_leer()
        case ("ldfss0")
            call destroy_scheme_ldfss0()
    end select
end subroutine destroy_scheme

subroutine compute_residue()
    implicit none
    select case (scheme)
        case ("van_leer")
            residue = get_residue_van_leer()
        case ("ldfss0")
            residue = get_residue_ldfss0()
    end select
end subroutine compute_residue

