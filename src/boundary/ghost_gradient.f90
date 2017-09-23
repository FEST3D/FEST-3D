module ghost_gradients
  ! ---------------------------------------------------------
  !   # 170508 - Jatinder Pal Singh Sandhu
  !---------------------------------------------------------
  
  ! geometry
  use global_vars, only : imx, jmx, kmx
  use global_vars, only : xnx, xny, xnz
  use global_vars, only : ynx, yny, ynz
  use global_vars, only : znx, zny, znz
  use global_vars, only : xA, yA, zA
  use global_vars, only : volume

  ! state variable
  use global_vars, only : density
  use global_vars, only : pressure
  use global_vars, only : R_gas
  use global_vars, only : n_grad
  use global_vars, only : n_var

  ! state and gradients
  use global_vars, only : qp
  use global_vars, only : gradqp_x
  use global_vars, only : gradqp_y
  use global_vars, only : gradqp_z

  ! layout boundary condition id for face
  use global_vars, only : imin_id
  use global_vars, only : imax_id
  use global_vars, only : jmin_id
  use global_vars, only : jmax_id
  use global_vars, only : kmin_id
  use global_vars, only : kmax_id

  use utils,       only : dmsg
  
  implicit none

  public :: apply_gradient_bc

  contains

    subroutine apply_gradient_bc
      implicit none

      call dmsg(1,'ghost_gradients', 'apply_gradient_bc')
  !    if(imin_id<0)then
        call apply('imin')
  !    end if
  !    if(imax_id<0)then
        call apply('imax')
  !    end if
  !    if(jmin_id<0)then
        call apply('jmin')
  !    end if
  !    if(jmax_id<0)then
        call apply('jmax')
  !    end if
  !    if(kmin_id<0)then
        call apply('kmin')
  !    end if
  !    if(kmax_id<0)then
        call apply('kmax')
  !    end if
      
    end subroutine apply_gradient_bc
    



    subroutine apply(face)
      implicit none
      character(len=*) :: face
      real, dimension(n_grad) :: qp_I
      real, dimension(n_grad) :: qp_G
      real    :: T_I
      real    :: T_G 
      real    :: c_x
      real    :: c_y
      real    :: c_z
      integer :: n
      integer :: i,j,k

      !-----------------------------------------------------------
      ! gradqp_G = (qp_I - qp_G)*Area_W*unit_normal_G/(volume_G)
      ! volume_G = volume_I
      !-----------------------------------------------------------

      n = n_grad
      select case (face)
        
        case('imin')

          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = 1,1
                c_x  = xA(i,j,k)*xnx(i,j,k)/volume(i,j,k)
                c_y  = xA(i,j,k)*xny(i,j,k)/volume(i,j,k)
                c_z  = xA(i,j,k)*xnz(i,j,k)/volume(i,j,k)
                T_I  = pressure(i  ,j,k)/(R_gas*density(i  ,j,k))
                T_G  = pressure(i-1,j,k)/(R_gas*density(i-1,j,k))
                qp_I = qp(i  ,j,k,2:n_var)
                qp_G = qp(i-1,j,k,2:n_var)

                gradqp_x(i-1,j,k,:) = (qp_I - qp_G)*c_x 
                gradqp_y(i-1,j,k,:) = (qp_I - qp_G)*c_y
                gradqp_z(i-1,j,k,:) = (qp_I - qp_G)*c_z
                gradqp_x(i-1,j,k,4) = ( T_I -  T_G)*c_x
                gradqp_y(i-1,j,k,4) = ( T_I -  T_G)*c_y
                gradqp_z(i-1,j,k,4) = ( T_I -  T_G)*c_z
              end do
            end do
          end do

        case('imax')
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = imx,imx
                c_x  = xA(i,j,k)*xnx(i,j,k)/volume(i-1,j,k)
                c_y  = xA(i,j,k)*xny(i,j,k)/volume(i-1,j,k)
                c_z  = xA(i,j,k)*xnz(i,j,k)/volume(i-1,j,k)
                T_I  = pressure(i-1,j,k)/(R_gas*density(i-1,j,k))
                T_G  = pressure(i  ,j,k)/(R_gas*density(i  ,j,k))
                qp_I = qp(i-1,j,k,2:n_var)
                qp_G = qp(i  ,j,k,2:n_var)

                gradqp_x(i,j,k,:) = -(qp_I - qp_G)*c_x 
                gradqp_y(i,j,k,:) = -(qp_I - qp_G)*c_y
                gradqp_z(i,j,k,:) = -(qp_I - qp_G)*c_z
                gradqp_x(i,j,k,4) = -( T_I -  T_G)*c_x
                gradqp_y(i,j,k,4) = -( T_I -  T_G)*c_y
                gradqp_z(i,j,k,4) = -( T_I -  T_G)*c_z
              end do
            end do
          end do


        case('jmin')
          do k = 1,kmx-1
            do j = 1,1
              do i = 1,imx-1
                c_x  = yA(i,j,k)*ynx(i,j,k)/volume(i,j,k)
                c_y  = yA(i,j,k)*yny(i,j,k)/volume(i,j,k)
                c_z  = yA(i,j,k)*ynz(i,j,k)/volume(i,j,k)
                T_I  = pressure(i,j  ,k)/(R_gas*density(i,j  ,k))
                T_G  = pressure(i,j-1,k)/(R_gas*density(i,j-1,k))
                qp_I = qp(i,j  ,k,2:n_var)
                qp_G = qp(i,j-1,k,2:n_var)

                gradqp_x(i,j-1,k,:) = (qp_I - qp_G)*c_x 
                gradqp_y(i,j-1,k,:) = (qp_I - qp_G)*c_y
                gradqp_z(i,j-1,k,:) = (qp_I - qp_G)*c_z
                gradqp_x(i,j-1,k,4) = ( T_I -  T_G)*c_x
                gradqp_y(i,j-1,k,4) = ( T_I -  T_G)*c_y
                gradqp_z(i,j-1,k,4) = ( T_I -  T_G)*c_z
              end do
            end do
          end do

        case('jmax')
          do k = 1,kmx-1
            do j = jmx,jmx
              do i = 1,imx-1
                c_x  = yA(i,j,k)*ynx(i,j,k)/volume(i,j-1,k)
                c_y  = yA(i,j,k)*yny(i,j,k)/volume(i,j-1,k)
                c_z  = yA(i,j,k)*ynz(i,j,k)/volume(i,j-1,k)
                T_I  = pressure(i,j-1,k)/(R_gas*density(i,j-1,k))
                T_G  = pressure(i,j  ,k)/(R_gas*density(i,j  ,k))
                qp_I = qp(i,j-1,k,2:n_var)
                qp_G = qp(i,j  ,k,2:n_var)

                gradqp_x(i,j,k,:) = -(qp_I - qp_G)*c_x 
                gradqp_y(i,j,k,:) = -(qp_I - qp_G)*c_y
                gradqp_z(i,j,k,:) = -(qp_I - qp_G)*c_z
                gradqp_x(i,j,k,4) = -( T_I -  T_G)*c_x
                gradqp_y(i,j,k,4) = -( T_I -  T_G)*c_y
                gradqp_z(i,j,k,4) = -( T_I -  T_G)*c_z
              end do
            end do
          end do


        case('kmin')
          do k = 1,1
            do j = 1,jmx-1
              do i = 1,imx-1
                c_x  = zA(i,j,k)*znx(i,j,k)/volume(i,j,k)
                c_y  = zA(i,j,k)*zny(i,j,k)/volume(i,j,k)
                c_z  = zA(i,j,k)*znz(i,j,k)/volume(i,j,k)
                T_I  = pressure(i,j,k  )/(R_gas*density(i,j,k  ))
                T_G  = pressure(i,j,k-1)/(R_gas*density(i,j,k-1))
                qp_I = qp(i,j,k  ,2:n_var)
                qp_G = qp(i,j,k-1,2:n_var)

                gradqp_x(i,j,k-1,:) = (qp_I - qp_G)*c_x 
                gradqp_y(i,j,k-1,:) = (qp_I - qp_G)*c_y
                gradqp_z(i,j,k-1,:) = (qp_I - qp_G)*c_z
                gradqp_x(i,j,k-1,4) = ( T_I -  T_G)*c_x
                gradqp_y(i,j,k-1,4) = ( T_I -  T_G)*c_y
                gradqp_z(i,j,k-1,4) = ( T_I -  T_G)*c_z
              end do
            end do
          end do

        case('kmax')
          do k = kmx,kmx
            do j = 1,jmx-1
              do i = 1,imx-1
                c_x  = zA(i,j,k)*znx(i,j,k)/volume(i,j,k-1)
                c_y  = zA(i,j,k)*zny(i,j,k)/volume(i,j,k-1)
                c_z  = zA(i,j,k)*znz(i,j,k)/volume(i,j,k-1)
                T_I  = pressure(i,j,k-1)/(R_gas*density(i,j,k-1))
                T_G  = pressure(i,j,k  )/(R_gas*density(i,j,k  ))
                qp_I = qp(i,j,k-1,2:n_var)
                qp_G = qp(i,j,k  ,2:n_var)

                gradqp_x(i,j,k,:) = -(qp_I - qp_G)*c_x 
                gradqp_y(i,j,k,:) = -(qp_I - qp_G)*c_y
                gradqp_z(i,j,k,:) = -(qp_I - qp_G)*c_z
                gradqp_x(i,j,k,4) = -( T_I -  T_G)*c_x
                gradqp_y(i,j,k,4) = -( T_I -  T_G)*c_y
                gradqp_z(i,j,k,4) = -( T_I -  T_G)*c_z
              end do
            end do
          end do

        case DEFAULT
          print*, "Ghost gradients : Wrong face name string"

      end select
          
    end subroutine apply


 
end module ghost_gradients
