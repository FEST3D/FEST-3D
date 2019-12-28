  !< Apply flow tangency boundary condition
module FT_bc
  !< Apply flow tangency boundary condition
  !--------------------------------------------
  use vartypes
  use copy_bc   , only : copy3

  implicit none
  private

  public :: flow_tangency

  contains

    subroutine flow_tangency(qp, face, Ifaces, Jfaces, Kfaces, dims)
      !< Apply flow tangency boundary condition
      implicit none
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout):: qp
      character(len=*), intent(in) :: face
      !< Face over which flow tangency condition has to be applied
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      real(wp) :: dot
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: imx, jmx, kmx

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      select case(face)
        case("imin")
          do k=1,kmx-1
            do j=1,jmx-1
              do l = 1,3
                dot = qp(l,j,k,2)*Ifaces(1,j,k)%nx + qp(l,j,k,3)*Ifaces(1,j,k)%ny + qp(l,j,k,4)*Ifaces(1,j,k)%nz
                qp(-l+1,j,k,2) = qp(l,j,k,2) - (2.0*dot*Ifaces(1,j,k)%nx)
                qp(-l+1,j,k,3) = qp(l,j,k,3) - (2.0*dot*Ifaces(1,j,k)%ny)
                qp(-l+1,j,k,4) = qp(l,j,k,4) - (2.0*dot*Ifaces(1,j,k)%nz)
              end do
            end do
          end do
        case("imax")
          do k=1,kmx-1
            do j=1,jmx-1
              do l = 1,3
                dot = qp(imx-l,j,k,2)*Ifaces(imx,j,k)%nx + qp(imx-l,j,k,3)*Ifaces(imx,j,k)%ny + qp(imx-l,j,k,4)*Ifaces(imx,j,k)%nz
                qp(imx+l-1,j,k,2) = qp(imx-l,j,k,2) - (2.0*dot*Ifaces(imx,j,k)%nx)
                qp(imx+l-1,j,k,3) = qp(imx-l,j,k,3) - (2.0*dot*Ifaces(imx,j,k)%ny)
                qp(imx+l-1,j,k,4) = qp(imx-l,j,k,4) - (2.0*dot*Ifaces(imx,j,k)%nz)
              end do
            end do
          end do
        case ("jmin")
          do k=1,kmx-1
            do i=1,imx-1
              do l =1,3
                dot = qp(i,l,k,2)*Jfaces(i,1,k)%nx + qp(i,l,k,3)*Jfaces(i,1,k)%ny + qp(i,l,k,4)*Jfaces(i,1,k)%nz
                qp(i,-l+1,k,2) = qp(i,l,k,2) - (2.0*dot*Ifaces(i,1,k)%nx)
                qp(i,-l+1,k,3) = qp(i,l,k,3) - (2.0*dot*Ifaces(i,1,k)%ny)
                qp(i,-l+1,k,4) = qp(i,l,k,4) - (2.0*dot*Ifaces(i,1,k)%nz)
              end do
            end do
          end do
        case ("jmax")
          do k=1,kmx-1
            do i=1,imx-1
              do l = 1,3
                dot = qp(i,jmx-l,k,2)*Jfaces(i,jmx,k)%nx + qp(i,jmx-l,k,3)*Jfaces(i,jmx,k)%ny + qp(i,jmx-l,k,4)*Jfaces(i,jmx,k)%nz
                qp(i,jmx+l-1,k,2) = qp(i,jmx-l,k,2) - (2.0*dot*Ifaces(i,jmx,k)%nx)
                qp(i,jmx+l-1,k,3) = qp(i,jmx-l,k,3) - (2.0*dot*Ifaces(i,jmx,k)%ny)
                qp(i,jmx+l-1,k,4) = qp(i,jmx-l,k,4) - (2.0*dot*Ifaces(i,jmx,k)%nz)
              end do
            end do
          end do
        case("kmin")
          do j=1,jmx-1
            do i=1,imx-1
              do l = 1,3
                dot = qp(i,j,l,2)*Kfaces(i,j,1)%nx + qp(i,j,l,3)*Kfaces(i,j,1)%ny + qp(i,j,l,4)*Kfaces(i,j,1)%nz
                qp(i,j,-l+1,2) = qp(i,j,l,2) - (2.0*dot*Ifaces(i,j,1)%nx)
                qp(i,j,-l+1,3) = qp(i,j,l,3) - (2.0*dot*Ifaces(i,j,1)%ny)
                qp(i,j,-l+1,4) = qp(i,j,l,4) - (2.0*dot*Ifaces(i,j,1)%nz)
              end do
            end do
          end do
        case("kmax")
          do j=1,jmx-1
            do i=1,imx-1
              do l=1,3
                dot = qp(i,j,kmx-l,2)*Kfaces(i,j,kmx)%nx + qp(i,j,kmx-l,3)*Kfaces(i,j,kmx)%ny + qp(i,j,kmx-l,4)*Kfaces(i,j,kmx)%nz
                qp(i,j,kmx+l-1,2) = qp(i,j,kmx-l,2) - (2.0*dot*Ifaces(i,j,kmx)%nx)
                qp(i,j,kmx+l-1,3) = qp(i,j,kmx-l,3) - (2.0*dot*Ifaces(i,j,kmx)%ny)
                qp(i,j,kmx+l-1,4) = qp(i,j,kmx-l,4) - (2.0*dot*Ifaces(i,j,kmx)%nz)
              end do
            end do
          end do
      end select
    end subroutine flow_tangency

end module FT_bc
