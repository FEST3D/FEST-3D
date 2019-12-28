  !< The viscous module contains the viscous fluxes calculations 
module viscous
  !< The viscous module contains the viscous fluxes calculations
  !-----------------------------------------------------------------
#include "error.inc"
  use vartypes
  use gradients  , only : gradu_x
  use gradients  , only : gradu_y
  use gradients  , only : gradu_z
  use gradients  , only : gradv_x
  use gradients  , only : gradv_y
  use gradients  , only : gradv_z
  use gradients  , only : gradw_x
  use gradients  , only : gradw_y
  use gradients  , only : gradw_z
  use gradients  , only : gradT_x
  use gradients  , only : gradT_y
  use gradients  , only : gradT_z
  use gradients  , only : gradtk_x
  use gradients  , only : gradtk_y
  use gradients  , only : gradtk_z
  use gradients  , only : gradtw_x
  use gradients  , only : gradtw_y
  use gradients  , only : gradtw_z
  use gradients  , only : gradtkl_x
  use gradients  , only : gradtkl_y
  use gradients  , only : gradtkl_z
  use gradients  , only : gradtv_x
  use gradients  , only : gradtv_y
  use gradients  , only : gradtv_z
  use gradients  , only : gradtgm_x
  use gradients  , only : gradtgm_y
  use gradients  , only : gradtgm_z
  use viscosity, only : mu
  use viscosity, only : mu_t
  use global_sst , only : sst_F1
  use global_sst , only : sigma_k1
  use global_sst , only : sigma_k2
  use global_sst , only : sigma_w1
  use global_sst , only : sigma_w2
  use global_kkl , only : sigma_k
  use global_kkl , only : sigma_phi
  use global_sa  , only : sigma_sa
  use global_sa  , only : cb2
  use utils      , only : alloc
  implicit none
  private

  integer :: imx, jmx, kmx

  public :: compute_viscous_fluxes

  contains

    subroutine compute_viscous_fluxes(F, G, H, qp, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims)
      !< Call to all viscous flux subroutine based on 
      !< the drection and turbulence/transition model being
      !< used

        implicit none
        type(schemetype), intent(in) :: scheme
        !< finite-volume Schemes
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx
        real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
        real(wp), dimension(:, :, :, :), intent(inout) :: F, G, H
        type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
        !< Input cell quantities: volume
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Store face quantites for I faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Store face quantites for J faces 
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Store face quantites for K faces 
        integer, dimension(3) :: flagsi=(/1,0,0/)
        integer, dimension(3) :: flagsj=(/0,1,0/)
        integer, dimension(3) :: flagsk=(/0,0,1/)

        imx = dims%imx
        jmx = dims%jmx
        kmx = dims%kmx

        call compute_viscous_fluxes_laminar(F, qp, cells, Ifaces, flagsi, scheme, flow, dims)
        call compute_viscous_fluxes_laminar(G, qp, cells, Jfaces, flagsj, scheme, flow, dims)
        !if(kmx==2)then
        !  continue
        !else
          call compute_viscous_fluxes_laminar(H, qp, cells, Kfaces, flagsk,scheme, flow, dims)
        !end if
        
        select case(trim(scheme%turbulence))
          case('none')
            !do nothing
            continue

          case('sa', 'saBC')
            call compute_viscous_fluxes_sa(F, qp, cells, Ifaces, flagsi,dims)
            call compute_viscous_fluxes_sa(G, qp, cells, Jfaces, flagsj,dims)
            call compute_viscous_fluxes_sa(H, qp, cells, Kfaces, flagsk,dims)
          case('sst', 'sst2003')
            call compute_viscous_fluxes_sst(F, qp,cells,  Ifaces, flagsi,dims)
            call compute_viscous_fluxes_sst(G, qp,cells,  Jfaces, flagsj,dims)
            if(kmx==2)then
              continue
            else
              call compute_viscous_fluxes_sst(H, qp,cells, Kfaces, flagsk,dims)
            end if
          case('kkl')
            call compute_viscous_fluxes_kkl(F, qp,cells, Ifaces, flagsi,dims)
            call compute_viscous_fluxes_kkl(G, qp,cells, Jfaces, flagsj,dims)
            !if(kmx==2)then
            !  continue
            !else
              call compute_viscous_fluxes_kkl(H, qp,cells, Kfaces, flagsk,dims)
            !end if
          case DEFAULT
            Fatal_error
        end select


        select case(trim(scheme%transition))
          case('lctm2015')
            call compute_viscous_fluxes_lctm2015(F, qp,cells, Ifaces, flagsi,dims)
            call compute_viscous_fluxes_lctm2015(G, qp,cells, Jfaces, flagsj,dims)
            if(kmx==2)then
              continue
            else
              call compute_viscous_fluxes_lctm2015(H, qp,cells, Kfaces, flagsk,dims)
            end if
          case('none', 'bc')
            !do nothing
            continue
          case DEFAULT
            Fatal_error
        end select


            if (any(isnan(G))) then
              Fatal_error
            end if
            if (any(isnan(F))) then
              Fatal_error
            end if
            if (any(isnan(H))) then
              Fatal_error
            end if

    end subroutine compute_viscous_fluxes

    subroutine compute_viscous_fluxes_laminar(F, qp, cells, faces, flags, scheme, flow, dims)
      !< Compute viscous fluxes for first five Navier-Stokes equation
      implicit none
      real(wp), dimension(:, :, :, :), intent(inout) :: F 
      !< Flux array
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      integer, dimension(3), intent(in) :: flags
      !< flags for direction switch
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
     !< Store primitive variable at cell center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
      !< Face quantities: area and unit normal
      ! local variables
      real(wp) :: dudx, dudy, dudz
      real(wp) :: dvdx, dvdy, dvdz
      real(wp) :: dwdx, dwdy, dwdz
      real(wp) :: dTdx, dTdy, dTdz
      real(wp) :: normal_comp
      real(wp) :: d_LR
      real(wp) :: T_RE
      real(wp) :: T_LE
      real(wp) :: K_heat, Qx, Qy, Qz
      real(wp) :: mu_f
      real(wp) :: mut_f
      real(wp) :: total_mu
      real(wp) :: delx
      real(wp) :: dely
      real(wp) :: delz
      real(wp) :: delu
      real(wp) :: delv
      real(wp) :: delw
      real(wp) :: delT
      real(wp) :: Tau_xx
      real(wp) :: Tau_xy
      real(wp) :: Tau_xz
      real(wp) :: Tau_yx
      real(wp) :: Tau_yy
      real(wp) :: Tau_yz
      real(wp) :: Tau_zx
      real(wp) :: Tau_zy
      real(wp) :: Tau_zz
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
      real(wp) :: area
      real(wp) :: uface
      real(wp) :: vface
      real(wp) :: wface
      integer :: i, j, k
      integer :: ii, jj, kk

      ii = flags(1)
      jj = flags(2)
      kk = flags(3)

      !---------------------------------------------------------------------
      ! Calculating the fluxes at the faces
      !--------------------------------------------------------------------
      do k = 1, dims%kmx - 1 + kk
       do j = 1, dims%jmx - 1 + jj
        do i = 1, dims%imx - 1 + ii

          !--- FACE Gradients ---!
          ! Gradients at face as average of gradients at cell centres
          dudx = 0.5 * (gradu_x(i-ii, j-jj, k-kk) + gradu_x(i, j, k))
          dudy = 0.5 * (gradu_y(i-ii, j-jj, k-kk) + gradu_y(i, j, k))
          dudz = 0.5 * (gradu_z(i-ii, j-jj, k-kk) + gradu_z(i, j, k))
          dvdx = 0.5 * (gradv_x(i-ii, j-jj, k-kk) + gradv_x(i, j, k))
          dvdy = 0.5 * (gradv_y(i-ii, j-jj, k-kk) + gradv_y(i, j, k))
          dvdz = 0.5 * (gradv_z(i-ii, j-jj, k-kk) + gradv_z(i, j, k))
          dwdx = 0.5 * (gradw_x(i-ii, j-jj, k-kk) + gradw_x(i, j, k))
          dwdy = 0.5 * (gradw_y(i-ii, j-jj, k-kk) + gradw_y(i, j, k))
          dwdz = 0.5 * (gradw_z(i-ii, j-jj, k-kk) + gradw_z(i, j, k))
          dTdx = 0.5 * (gradT_x(i-ii, j-jj, k-kk) + gradT_x(i, j, k))
          dTdy = 0.5 * (gradT_y(i-ii, j-jj, k-kk) + gradT_y(i, j, k))
          dTdz = 0.5 * (gradT_z(i-ii, j-jj, k-kk) + gradT_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
          ! distance between cell center of adjacent cell for the i,j,k face
          delx = cells(i, j, k)%centerx - cells(i-ii, j-jj, k-kk)%centerx
          dely = cells(i, j, k)%centery - cells(i-ii, j-jj, k-kk)%centery
          delz = cells(i, j, k)%centerz - cells(i-ii, j-jj, k-kk)%centerz

          d_LR = sqrt(delx*delx + dely*dely + delz*delz)

          ! Finding the temperature of left and right element to the face i,j,k
          T_LE = qp(i-ii, j-jj, k-kk, 5) / (qp(i-ii, j-jj, k-kk, 1) * flow%R_gas)
          T_RE = qp(i, j, k, 5) / (qp(i, j, k, 1) * flow%R_gas)

          ! difference in state across face
          delu = qp(i, j, k, 2) - qp(i-ii, j-jj, k-kk, 2) !x_speed
          delv = qp(i, j, k, 3) - qp(i-ii, j-jj, k-kk, 3) !y_speed
          delw = qp(i, j, k, 4) - qp(i-ii, j-jj, k-kk, 4) !z_speed
          delT = T_RE - T_LE

          !normal_comp   = ( delta(phi) - (grad(phi).dot.delR) )/magnitudeR
          !new grad(phi) =  grad(phi) + correction(normal_comp.dot.delR/magnitudeR)

          normal_comp = (delu - (dudx*delx + dudy*dely + dudz*delz))/d_LR
          dudx        =  dudx + (normal_comp * delx / d_LR)
          dudy        =  dudy + (normal_comp * dely / d_LR)
          dudz        =  dudz + (normal_comp * delz / d_LR)

          normal_comp = (delv - (dvdx*delx + dvdy*dely + dvdz*delz))/d_LR
          dvdx        =  dvdx + (normal_comp * delx / d_LR)
          dvdy        =  dvdy + (normal_comp * dely / d_LR)
          dvdz        =  dvdz + (normal_comp * delz / d_LR)

          normal_comp = (delw - (dwdx*delx + dwdy*dely + dwdz*delz))/d_LR
          dwdx        = dwdx + (normal_comp * delx / d_LR)
          dwdy        = dwdy + (normal_comp * dely / d_LR)
          dwdz        = dwdz + (normal_comp * delz / d_LR)

          normal_comp = (delT - (dTdx*delx + dTdy*dely + dTdz*delz))/d_LR
          dTdx        = dTdx + (normal_comp * delx / d_LR)
          dTdy        = dTdy + (normal_comp * dely / d_LR)
          dTdz        = dTdz + (normal_comp * delz / d_LR)
          !--- end of ODD-EVEN coupling correction ---!

          mu_f  = 0.5*(mu(i-ii, j-jj, k-kk) + mu(i,j,k))
          if(trim(scheme%turbulence)/='none') then
            mut_f = 0.5*(mu_t(i-ii, j-jj, k-kk) + mu_t(i, j, k))
          else
            mut_f = 0.0 
          end if

          ! effective viscosity
          total_mu = mu_f + mut_f

          ! Using lambda = -2 * mu / 3
          ! diagonal terms of stress tensor
          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          ! off diagonal symmetrical part of stress tensor
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)
          Tau_yx = Tau_xy
          Tau_zx = Tau_xz
          Tau_zy = Tau_yz

          ! Pr: Prandtl Number and tPr: Turbulent Prandtl number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/flow%Pr + mut_f/flow%tPr)* flow%gm * flow%R_gas / (flow%gm - 1)
          Qx = K_heat*dTdx
          Qy = K_heat*dTdy
          Qz = K_heat*dTdz

          ! calling some element from memory and keep them handy for calculation
          nx    = faces(i,j,k)%nx
          ny    = faces(i,j,k)%ny
          nz    = faces(i,j,k)%nz
          area  = faces(i,j,k)%A
          uface = 0.5 * (qp(i-ii, j-jj, k-kk, 2) + qp(i, j, k, 2))
          vface = 0.5 * (qp(i-ii, j-jj, k-kk, 3) + qp(i, j, k, 3))
          wface = 0.5 * (qp(i-ii, j-jj, k-kk, 4) + qp(i, j, k, 4))

          ! adding viscous fluxes to stored convective flux
          F(i, j, k, 2) = F(i, j, k, 2) - ((Tau_xx*nx + Tau_xy*ny + Tau_xz*nz)*area)
          F(i, j, k, 3) = F(i, j, k, 3) - ((Tau_yx*nx + Tau_yy*ny + Tau_yz*nz)*area)
          F(i, j, k, 4) = F(i, j, k, 4) - ((Tau_zx*nx + Tau_zy*ny + Tau_zz*nz)*area)
         
          ! Energy flux
          ! (TijVj + Qi)ni
          F(i, j, k, 5) = F(i, j, k, 5) - (area * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + Qx)*nx) + &
                          ((Tau_yx*uface + Tau_yy*vface + Tau_yz*wface + Qy)*ny) + &
                          ((Tau_zx*uface + Tau_zy*vface + Tau_zz*wface + Qz)*nz) ) )
          
         
        end do
       end do
      end do

    end subroutine compute_viscous_fluxes_laminar


    subroutine compute_viscous_fluxes_sst(F, qp, cells, faces, flags, dims)
      !< Compute viscous fluxes for additianal equations due to SST turbulence model
      implicit none
      real(wp), dimension(:, :, :, :), intent(inout) :: F 
      !< flux array
      type(extent), intent(in) :: dims
      !< Extent of the domain: imx,jmx,kmx
      integer, dimension(3), intent(in) :: flags
      !< flags for direction swithc
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
     !< Store primitive variable at cell center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
      !< quantities on the face
      ! local variables
      real(wp) :: tkface
      real(wp) :: rhoface
      real(wp) :: normal_comp
      real(wp) :: d_LR
      real(wp) :: mu_f
      real(wp) :: mut_f
      real(wp) :: total_mu
      real(wp) :: delx
      real(wp) :: dely
      real(wp) :: delz
      real(wp) :: deltk
      real(wp) :: deltw
      real(wp) :: Tau_xx
      real(wp) :: Tau_yy
      real(wp) :: Tau_zz
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
      real(wp) :: area
      integer :: i, j, k
      integer :: ii, jj, kk
      !--- sst variable requirement ---!
      real(wp) :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
      real(wp) ::  F1
      real(wp) ::  sigma_kf
      real(wp) ::  sigma_wf

      ii = flags(1)
      jj = flags(2)
      kk = flags(3)

      !---------------------------------------------------------------------
      ! Calculating the turbulent viscous fluxes at the faces
      !--------------------------------------------------------------------
      do k = 1, dims%kmx - 1 + kk
       do j = 1, dims%jmx - 1 + jj
        do i = 1, dims%imx - 1 + ii

          !--- FACE Gradients ---!
          ! Gradients at face as average of gradients at cell centres
          dtkdx = 0.5 * (gradtk_x(i-ii, j-jj, k-kk) + gradtk_x(i, j, k))
          dtkdy = 0.5 * (gradtk_y(i-ii, j-jj, k-kk) + gradtk_y(i, j, k))
          dtkdz = 0.5 * (gradtk_z(i-ii, j-jj, k-kk) + gradtk_z(i, j, k))
          dtwdx = 0.5 * (gradtw_x(i-ii, j-jj, k-kk) + gradtw_x(i, j, k))
          dtwdy = 0.5 * (gradtw_y(i-ii, j-jj, k-kk) + gradtw_y(i, j, k))
          dtwdz = 0.5 * (gradtw_z(i-ii, j-jj, k-kk) + gradtw_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
          ! distance between cell center of adjacent cell for the i,j,k face
          delx = cells(i, j, k)%centerx - cells(i-ii, j-jj, k-kk)%centerx
          dely = cells(i, j, k)%centery - cells(i-ii, j-jj, k-kk)%centery
          delz = cells(i, j, k)%centerz - cells(i-ii, j-jj, k-kk)%centerz

          d_LR = sqrt(delx*delx + dely*dely + delz*delz)

          ! difference in state across face
          deltk = qp(i, j, k, 6) - qp(i-ii, j-jj, k-kk, 6) !TKE
          deltw = qp(i, j, k, 7) - qp(i-ii, j-jj, k-kk, 7) !Omega

          !normal_comp   = ( delta(phi) - (grad(phi).dot.delR) )/magnitudeR
          !new grad(phi) =  grad(phi) + correction(normal_comp.dot.delR/magnitudeR)

          normal_comp = (deltk - (dtkdx*delx + dtkdy*dely + dtkdz*delz))/d_LR
          dtkdx       =  dtkdx + (normal_comp * delx / d_LR)
          dtkdy       =  dtkdy + (normal_comp * dely / d_LR)
          dtkdz       =  dtkdz + (normal_comp * delz / d_LR)

          normal_comp = (deltw - (dtwdx*delx + dtwdy*dely + dtwdz*delz))/d_LR
          dtwdx       =  dtwdx + (normal_comp * delx / d_LR)
          dtwdy       =  dtwdy + (normal_comp * dely / d_LR)
          dtwdz       =  dtwdz + (normal_comp * delz / d_LR)
          !--- end of ODD-EVEN coupling correction ---!

          mu_f     = 0.5*(  mu(i-ii, j-jj, k-kk) +   mu(i, j, k))
          mut_f    = 0.5*(mu_t(i-ii, j-jj, k-kk) + mu_t(i, j, k))
          F1       = 0.5*(sst_F1(i-ii, j-jj, k-kk) + sst_F1(i, j, k))
          sigma_kf = sigma_k1*F1 + sigma_k2*(1.0 - F1)
          sigma_wf = sigma_w1*F1 + sigma_w2*(1.0 - F1)

          total_mu = mu_f + mut_f
          rhoface  = 0.5 * (qp(i-ii, j-jj, k-kk, 1) + qp(i, j, k, 1)) !Density
          tkface   = 0.5 * (qp(i-ii, j-jj, k-kk, 6) + qp(i, j, k, 6)) !TKE
          ! k in reynolds stress
          Tau_xx = -2.0*rhoface*tkface/3.0
          Tau_yy = Tau_xx
          Tau_zz = Tau_xx

          ! calling some element from memory and keep them handy for calculation
          nx    = faces(i,j,k)%nx
          ny    = faces(i,j,k)%ny
          nz    = faces(i,j,k)%nz
          area  = faces(i,j,k)%A

          ! adding viscous fluxes to stored convective flux
          F(i, j, k, 2) = F(i, j, k, 2) - (Tau_xx * nx * area)
          F(i, j, k, 3) = F(i, j, k, 3) - (Tau_yy * ny * area)
          F(i, j, k, 4) = F(i, j, k, 4) - (Tau_zz * nz * area)
          F(i, j, k, 5) = F(i, j, k, 5) - (area*((mu_f + sigma_kf*mut_f)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)))
          F(i, j, k, 6) = F(i, j, k, 6) - (area*((mu_f + sigma_kf*mut_f)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)))
          F(i, j, k, 7) = F(i, j, k, 7) - (area*((mu_f + sigma_wf*mut_f)*(dtwdx*nx + dtwdy*ny + dtwdz*nz)))
         
        end do
       end do
      end do
    end subroutine compute_viscous_fluxes_sst


    subroutine compute_viscous_fluxes_kkl(F, qp, cells, faces, flags, dims)
      !< Compute viscous fluxes for additianal equations due to k-kL turbulence model
      implicit none
      real(wp), dimension(:, :, :, :), intent(inout) :: F 
      !< Flux array
      type(extent), intent(in) :: dims
      !< Extent of the domain: imx,jmx,kmx
      integer, dimension(3), intent(in) :: flags
      !< flags for directions switch
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
     !< Store primitive variable at cell center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
      !< Face quantities: area and unit normal
      ! local variables
      real(wp) :: tkface
      real(wp) :: rhoface
      real(wp) :: normal_comp
      real(wp) :: d_LR
      real(wp) :: mu_f
      real(wp) :: mut_f
      real(wp) :: total_mu
      real(wp) :: delx
      real(wp) :: dely
      real(wp) :: delz
      real(wp) :: deltk
      real(wp) :: deltkl
      real(wp) :: Tau_xx
      real(wp) :: Tau_yy
      real(wp) :: Tau_zz
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
      real(wp) :: area
      integer :: i, j, k
      integer :: ii, jj, kk
      !--- kkl variable requirement  ---!
      real(wp) :: dtkdx, dtkdy, dtkdz
      real(wp) :: dtkldx, dtkldy, dtkldz


      ii = flags(1)
      jj = flags(2)
      kk = flags(3)


      !---------------------------------------------------------------------
      ! Calculating the turbulent viscous fluxes at the faces
      !--------------------------------------------------------------------
      do k = 1, dims%kmx - 1 + kk
       do j = 1, dims%jmx - 1 + jj
        do i = 1, dims%imx - 1 + ii

          !--- FACE Gradients ---!
          ! Gradients at face as average of gradients at cell centres
          dtkdx  = 0.5 * (gradtk_x(i-ii, j-jj, k-kk) + gradtk_x(i, j, k))
          dtkdy  = 0.5 * (gradtk_y(i-ii, j-jj, k-kk) + gradtk_y(i, j, k))
          dtkdz  = 0.5 * (gradtk_z(i-ii, j-jj, k-kk) + gradtk_z(i, j, k))
          dtkldx = 0.5 * (gradtkl_x(i-ii, j-jj, k-kk) + gradtkl_x(i, j, k))
          dtkldy = 0.5 * (gradtkl_y(i-ii, j-jj, k-kk) + gradtkl_y(i, j, k))
          dtkldz = 0.5 * (gradtkl_z(i-ii, j-jj, k-kk) + gradtkl_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
          ! distance between cell center of adjacent cell for the i,j,k face
          delx = cells(i, j, k)%centerx - cells(i-ii, j-jj, k-kk)%centerx
          dely = cells(i, j, k)%centery - cells(i-ii, j-jj, k-kk)%centery
          delz = cells(i, j, k)%centerz - cells(i-ii, j-jj, k-kk)%centerz

          d_LR = sqrt(delx*delx + dely*dely + delz*delz)

          ! difference in state across face
          deltk  = qp(i, j, k, 6) - qp(i-ii, j-jj, k-kk, 6)   !TKE
          deltkl = qp(i, j, k, 7) - qp(i-ii, j-jj, k-kk, 7)   !Kl

          !normal_comp   = ( delta(phi) - (grad(phi).dot.delR) )/magnitudeR
          !new grad(phi) =  grad(phi) + correction(normal_comp.dot.delR/magnitudeR)

          normal_comp = (deltk - (dtkdx*delx + dtkdy*dely + dtkdz*delz))/d_LR
          dtkdx       =  dtkdx + (normal_comp * delx / d_LR)
          dtkdy       =  dtkdy + (normal_comp * dely / d_LR)
          dtkdz       =  dtkdz + (normal_comp * delz / d_LR)

          normal_comp  = (deltkl - (dtkldx*delx + dtkldy*dely + dtkldz*delz))/d_LR
          dtkldx       =  dtkldx + (normal_comp * delx / d_LR)
          dtkldy       =  dtkldy + (normal_comp * dely / d_LR)
          dtkldz       =  dtkldz + (normal_comp * delz / d_LR)
          !--- end of ODD-EVEN coupling correction ---!

          mu_f     = 0.5*(  mu(i-ii, j-jj, k-kk) +  mu(i, j, k))
          mut_f    = 0.5*(mu_t(i-ii, j-jj, k-kk) + mu_t(i, j, k))

          total_mu = mu_f + mut_f
          rhoface  = 0.5 * (qp(i-ii, j-jj, k-kk, 1) + qp(i, j, k, 1))
          tkface   = 0.5 * (qp(i-ii, j-jj, k-kk, 6) + qp(i, j, k, 6))
          ! k in reynolds stress
          Tau_xx = -2.0*rhoface*tkface/3.0
          Tau_yy = Tau_xx
          Tau_zz = Tau_xx

          ! calling some element from memory and keep them handy for calculation
          nx    = faces(i,j,k)%nx
          ny    = faces(i,j,k)%ny
          nz    = faces(i,j,k)%nz
          area  = faces(i,j,k)%A

          ! adding viscous fluxes to stored convective flux
          F(i, j, k, 2) = F(i, j, k, 2) - (Tau_xx * nx * area)
          F(i, j, k, 3) = F(i, j, k, 3) - (Tau_yy * ny * area)
          F(i, j, k, 4) = F(i, j, k, 4) - (Tau_zz * nz * area)
          F(i, j, k, 5) = F(i, j, k, 5) - (area*((mu_f + sigma_k*mut_f)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)))
          F(i, j, k, 6) = F(i, j, k, 6) - (area*((mu_f + sigma_k*mut_f)*(dtkdx*nx + dtkdy*ny + dtkdz*nz)))
          F(i, j, k, 7) = F(i, j, k, 7) - (area*((mu_f + sigma_phi*mut_f)*(dtkldx*nx + dtkldy*ny + dtkldz*nz)))
         
        end do
       end do
      end do

    end subroutine compute_viscous_fluxes_kkl


    subroutine compute_viscous_fluxes_sa(F, qp, cells, faces, flags, dims)
      !< Compute viscous fluxes for additianal equations due to SA turbulence model
      implicit none
      real(wp), dimension(:, :, :, :), intent(inout) :: F 
      !< Flux array
      type(extent), intent(in) :: dims
      !< Extent of the domain: imx,jmx,kmx
      integer, dimension(3), intent(in) :: flags
      !< flags for direction switch
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
      !< Face quantities: area and unit normal
      ! local variables
      real(wp) :: rhoface
      real(wp) :: normal_comp
      real(wp) :: d_LR
      real(wp) :: mu_f
      real(wp) :: mut_f
      real(wp) :: delx
      real(wp) :: dely
      real(wp) :: delz
      real(wp) :: deltv
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
      real(wp) :: area
      integer :: i, j, k
      integer :: ii, jj, kk
      !--- sa variable requirement ---!
      real(wp) :: dtvdx, dtvdy, dtvdz

      ii = flags(1)
      jj = flags(2)
      kk = flags(3)

      !---------------------------------------------------------------------
      ! Calculating the turbulent viscous fluxes at the faces
      !--------------------------------------------------------------------
      do k = 1, dims%kmx - 1 + kk
       do j = 1, dims%jmx - 1 + jj
        do i = 1, dims%imx - 1 + ii

          !--- FACE Gradients ---!
          ! Gradients at face as average of gradients at cell centres
          dtvdx = 0.5 * (gradtv_x(i-ii, j-jj, k-kk) + gradtv_x(i, j, k))
          dtvdy = 0.5 * (gradtv_y(i-ii, j-jj, k-kk) + gradtv_y(i, j, k))
          dtvdz = 0.5 * (gradtv_z(i-ii, j-jj, k-kk) + gradtv_z(i, j, k))

          !--- For ODD-EVEN coupling error ---!
          ! distance between cell center of adjacent cell for the i,j,k face
          delx = cells(i, j, k)%centerx - cells(i-ii, j-jj, k-kk)%centerx
          dely = cells(i, j, k)%centery - cells(i-ii, j-jj, k-kk)%centery
          delz = cells(i, j, k)%centerz - cells(i-ii, j-jj, k-kk)%centerz

          d_LR = sqrt(delx*delx + dely*dely + delz*delz)

          ! difference in state across face
          deltv = qp(i, j, k, 6) - qp(i-ii, j-jj, k-kk, 6)

          !normal_comp   = ( delta(phi) - (grad(phi).dot.delR) )/magnitudeR
          !new grad(phi) =  grad(phi) + correction(normal_comp.dot.delR/magnitudeR)
          normal_comp = (deltv - (dtvdx*delx + dtvdy*dely + dtvdz*delz))/d_LR
          dtvdx       =  dtvdx + (normal_comp * delx / d_LR)
          dtvdy       =  dtvdy + (normal_comp * dely / d_LR)
          dtvdz       =  dtvdz + (normal_comp * delz / d_LR)
          !--- end of ODD-EVEN coupling correction ---!

          rhoface  = 0.5 * (qp(i-ii, j-jj, k-kk, 1) + qp(i, j, k, 1))
          mu_f     = 0.5*(mu(i-ii, j-jj, k-kk) + mu(i, j, k))
          mut_f    = 0.5*(qp(i-ii, j-jj, k-kk, 6) + qp(i, j, k, 6))*rhoface

          ! calling some element from memory and keep them handy for calculation
          nx    = faces(i,j,k)%nx
          ny    = faces(i,j,k)%ny
          nz    = faces(i,j,k)%nz
          area  = faces(i,j,k)%A

          ! adding viscous fluxes to stored convective flux
          F(i, j, k, 6) = F(i, j, k, 6) - (area*((mu_f + mut_f)*(dtvdx*nx + dtvdy*ny + dtvdz*nz)))/sigma_sa
         
        end do
       end do
      end do
    end subroutine compute_viscous_fluxes_sa


    subroutine compute_viscous_fluxes_lctm2015(F, qp, cells, faces, flags,dims)
      !< Compute viscous fluxes for additianal equations due to LCTM2015 transition model
      implicit none
      real(wp), dimension(:, :, :, :), intent(inout) :: F 
      !< Flux array
      type(extent), intent(in) :: dims
      !< Extent of the doamin:imx,jmx,kmx
      integer, dimension(3), intent(in) :: flags
      !<flags for direction switch
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
      !< Face quantities: area and unit normal
      ! local variables
      real(wp) :: rhoface
      real(wp) :: normal_comp
      real(wp) :: d_LR
      real(wp) :: mu_f
      real(wp) :: mut_f
      real(wp) :: delx
      real(wp) :: dely
      real(wp) :: delz
      real(wp) :: deltgm
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
      real(wp) :: area
      integer :: i, j, k
      integer :: ii, jj, kk
      !--- sa variable requirement ---!
      real(wp) :: dtgmdx, dtgmdy, dtgmdz

      ii = flags(1)
      jj = flags(2)
      kk = flags(3)


      !---------------------------------------------------------------------
      ! Calculating the turbulent viscous fluxes at the faces
      !--------------------------------------------------------------------
      do k = 1, dims%kmx - 1 + kk
       do j = 1, dims%jmx - 1 + jj
        do i = 1, dims%imx - 1 + ii

          !--- FACE Gradients ---!
          ! Gradients at face as average of gradients at cell centres
          dtgmdx = 0.5 * (gradtgm_x(i-ii, j-jj, k-kk) + gradtgm_x(i, j, k))
          dtgmdy = 0.5 * (gradtgm_y(i-ii, j-jj, k-kk) + gradtgm_y(i, j, k))
          dtgmdz = 0.5 * (gradtgm_z(i-ii, j-jj, k-kk) + gradtgm_z(i, j, k))

          !--- For ODD-EVEN coupling error ---!
          ! distance between cell center of adjacent cell for the i,j,k face
          delx = cells(i, j, k)%centerx - cells(i-ii, j-jj, k-kk)%centerx
          dely = cells(i, j, k)%centery - cells(i-ii, j-jj, k-kk)%centery
          delz = cells(i, j, k)%centerz - cells(i-ii, j-jj, k-kk)%centerz

          d_LR = sqrt(delx*delx + dely*dely + delz*delz)

          ! difference in state across face
          deltgm = qp(i, j, k, 8) - qp(i-ii, j-jj, k-kk, 8)

          !normal_comp   = ( delta(phi) - (grad(phi).dot.delR) )/magnitudeR
          !new grad(phi) =  grad(phi) + correction(normal_comp.dot.delR/magnitudeR)
          normal_comp = (deltgm - (dtgmdx*delx + dtgmdy*dely + dtgmdz*delz))/d_LR
          dtgmdx       =  dtgmdx + (normal_comp * delx / d_LR)
          dtgmdy       =  dtgmdy + (normal_comp * dely / d_LR)
          dtgmdz       =  dtgmdz + (normal_comp * delz / d_LR)
          !--- end of ODD-EVEN coupling correction ---!

          rhoface  = 0.5 * (qp(i-ii, j-jj, k-kk, 1) + qp(i, j, k, 1))
          mu_f     = 0.5*(mu(i-ii, j-jj, k-kk) + mu(i, j, k))
          mut_f    = 0.5*(mu_t(i-ii, j-jj, k-kk) + mu_t(i, j, k))

          ! calling some element from memory and keep them handy for calculation
          nx    = faces(i,j,k)%nx
          ny    = faces(i,j,k)%ny
          nz    = faces(i,j,k)%nz
          area  = faces(i,j,k)%A

          ! adding viscous fluxes to stored convective flux
          F(i, j, k, dims%n_var) = F(i, j, k, dims%n_var) - (area*((mu_f + mut_f)*(dtgmdx*nx + dtgmdy*ny + dtgmdz*nz)))
         
        end do
       end do
      end do
    end subroutine compute_viscous_fluxes_lctm2015


end module viscous
