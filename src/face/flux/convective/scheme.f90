 !< Inviscid flux calculation through faces
module scheme
 !< Inviscid flux calculation through faces

#include "../../../debug.h"
#include "../../../error.h"
    use vartypes
    use utils, only: alloc
    use face_interpolant, only: setup_interpolant_scheme
    use face_interpolant
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
    public :: compute_fluxes
    public :: compute_residue

    contains

        subroutine setup_scheme(residue, F,G,H, control, dims)
            implicit none
            type(controltype), intent(in) :: control
            !< Control parameters
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), allocatable, intent(out), target :: residue
            !< Store residue at each cell-center
            real(wp), dimension(:, :, :, :), allocatable, intent(out) :: F
            !< Store fluxes throught the I faces
            real(wp), dimension(:, :, :, :), allocatable, intent(out) :: G
            !< Store fluxes throught the J faces
            real(wp), dimension(:, :, :, :), allocatable, intent(out) :: H
            !< Store fluxes throught the K faces

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            n_var = control%n_var

            call setup_interpolant_scheme(dims)

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - Scheme')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - Scheme')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - Scheme')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - Scheme')

        end subroutine setup_scheme


        subroutine compute_fluxes(F,G,H, Ifaces, Jfaces, Kfaces, scheme, flow, bc, dims)
        
            implicit none
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), intent(inout) :: F
            !< Store fluxes throught the I faces
            real(wp), dimension(:, :, :, :), intent(inout) :: G
            !< Store fluxes throught the J faces
            real(wp), dimension(:, :, :, :), intent(inout) :: H
            !< Store fluxes throught the K faces
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces 
            type(boundarytype), intent(in) :: bc
            !< boundary conditions and fixed values

            select case (scheme%scheme_name)
                case ("van_leer")
                  call compute_fluxes_van_leer(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case ("ausm")
                  call compute_fluxes_ausm(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case ("ausmP")
                  call compute_fluxes_ausmP(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case ("ausmUP")
                  call compute_fluxes_ausmUP(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case ("slau")
                  call compute_fluxes_slau(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case ("ldfss0")
                  call compute_fluxes_ldfss0(F,G,H,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case default
                    Fatal_error
            end select
            
        end subroutine compute_fluxes

        subroutine compute_residue(residue,F,G,H,dims)
            
            implicit none
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), intent(out)  :: residue
            !< Store residue at each cell-center
            real(wp), dimension(:, :, :, :), intent(in) :: F
            !< Store fluxes throught the I faces
            real(wp), dimension(:, :, :, :), intent(in) :: G
            !< Store fluxes throught the J faces
            real(wp), dimension(:, :, :, :), intent(in) :: H
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
        
        end subroutine compute_residue

end module scheme
