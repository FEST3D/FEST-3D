  !< Set value gradients in the ghost cells
module ghost_gradients
  !< Set value gradients in the ghost cells
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
  use global_vars, only : fixed_wall_temperature

  use utils,       only : dmsg
  
  implicit none

  public :: apply_gradient_bc

  contains

    subroutine apply_gradient_bc()
      !< call same subroutine for all the face
      implicit none

      call dmsg(1,'ghost_gradients', 'apply_gradient_bc')
      if(imin_id<0)then
        call apply('imin')
      end if
      if(imax_id<0)then
        call apply('imax')
      end if
      if(jmin_id<0)then
        call apply('jmin')
      end if
      if(jmax_id<0)then
        call apply('jmax')
      end if
      if(kmin_id<0)then
        call apply('kmin')
      end if
      if(kmax_id<0)then
        call apply('kmax')
      end if
      
    end subroutine apply_gradient_bc
    



    subroutine apply(face)
      !< Apply/set value of all gradient in the ghost cells
      !< gradqp_G = (qp_I - qp_G)*Area_W*unit_normal_G/(volume_G)
      !< volume_G = volume_I
      !-----------------------------------------------------------
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
      integer :: i,j,k,l
      real    :: nx
      real    :: ny
      real    :: nz
      real    :: dot

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
                nx   = xnx(i,j,k)
                ny   = xny(i,j,k)
                nz   = xnz(i,j,k)
                c_x  = xA(i,j,k)*nx/volume(i,j,k)
                c_y  = xA(i,j,k)*ny/volume(i,j,k)
                c_z  = xA(i,j,k)*nz/volume(i,j,k)
                T_I  = pressure(i  ,j,k)/(R_gas*density(i  ,j,k))
                T_G  = pressure(i-1,j,k)/(R_gas*density(i-1,j,k))
                qp_I = qp(i  ,j,k,2:n_var)
                qp_G = qp(i-1,j,k,2:n_var)

                ! normal component of gradient
                gradqp_x(i-1,j,k,:) = (qp_I - qp_G)*c_x 
                gradqp_y(i-1,j,k,:) = (qp_I - qp_G)*c_y
                gradqp_z(i-1,j,k,:) = (qp_I - qp_G)*c_z
                gradqp_x(i-1,j,k,4) = ( T_I -  T_G)*c_x
                gradqp_y(i-1,j,k,4) = ( T_I -  T_G)*c_y
                gradqp_z(i-1,j,k,4) = ( T_I -  T_G)*c_z
                if(imin_id==-5 .and. (fixed_wall_temperature(1)<1. .and. fixed_wall_temperature(1)>=0.))then
                  gradqp_x(i-1,j,k,4) = -gradqp_x(i,j,k,4)
                  gradqp_y(i-1,j,k,4) = -gradqp_y(i,j,k,4)
                  gradqp_z(i-1,j,k,4) = -gradqp_z(i,j,k,4)

                end if
                !parallel component of gradient
                do l=1,n
                  dot = (gradqp_x(i,j,k,l)*nx) + (gradqp_y(i,j,k,l)*ny) + (gradqp_z(i,j,k,l)*nz)
                  gradqp_x(i-1,j,k,l) = gradqp_x(i-1,j,k,l) + (gradqp_x(i,j,k,l) - dot*nx)
                  gradqp_y(i-1,j,k,l) = gradqp_y(i-1,j,k,l) + (gradqp_y(i,j,k,l) - dot*ny)
                  gradqp_z(i-1,j,k,l) = gradqp_z(i-1,j,k,l) + (gradqp_z(i,j,k,l) - dot*nz)
                end do
              end do
            end do
          end do

        case('imax')
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = imx,imx
                nx   = xnx(i,j,k)
                ny   = xny(i,j,k)
                nz   = xnz(i,j,k)
                c_x  = xA(i,j,k)*nx/volume(i-1,j,k)
                c_y  = xA(i,j,k)*ny/volume(i-1,j,k)
                c_z  = xA(i,j,k)*nz/volume(i-1,j,k)
                T_I  = pressure(i-1,j,k)/(R_gas*density(i-1,j,k))
                T_G  = pressure(i  ,j,k)/(R_gas*density(i  ,j,k))
                qp_I = qp(i-1,j,k,2:n_var)
                qp_G = qp(i  ,j,k,2:n_var)

                ! normal component of gradient
                gradqp_x(i,j,k,:) = -(qp_I - qp_G)*c_x 
                gradqp_y(i,j,k,:) = -(qp_I - qp_G)*c_y
                gradqp_z(i,j,k,:) = -(qp_I - qp_G)*c_z
                gradqp_x(i,j,k,4) = -( T_I -  T_G)*c_x
                gradqp_y(i,j,k,4) = -( T_I -  T_G)*c_y
                gradqp_z(i,j,k,4) = -( T_I -  T_G)*c_z
                if(imax_id==-5 .and. (fixed_wall_temperature(2)<1. .and. fixed_wall_temperature(2)>=0.))then
                gradqp_x(i,j,k,4) = -gradqp_x(i-1,j,k,4)
                gradqp_y(i,j,k,4) = -gradqp_y(i-1,j,k,4)
                gradqp_z(i,j,k,4) = -gradqp_z(i-1,j,k,4)
                end if
                !parallel component of gradient
                do l=1,n
                  dot = (gradqp_x(i-1,j,k,l)*nx) + (gradqp_y(i-1,j,k,l)*ny) + (gradqp_z(i-1,j,k,l)*nz)
                  gradqp_x(i,j,k,l) = gradqp_x(i,j,k,l) + (gradqp_x(i-1,j,k,l) - dot*nx)
                  gradqp_y(i,j,k,l) = gradqp_y(i,j,k,l) + (gradqp_y(i-1,j,k,l) - dot*ny)
                  gradqp_z(i,j,k,l) = gradqp_z(i,j,k,l) + (gradqp_z(i-1,j,k,l) - dot*nz)
                end do
              end do
            end do
          end do


        case('jmin')
          do k = 1,kmx-1
            do j = 1,1
              do i = 1,imx-1
                nx   = ynx(i,j,k)
                ny   = yny(i,j,k)
                nz   = ynz(i,j,k)
                c_x  = yA(i,j,k)*nx/volume(i,j,k)
                c_y  = yA(i,j,k)*ny/volume(i,j,k)
                c_z  = yA(i,j,k)*nz/volume(i,j,k)
                T_I  = pressure(i,j  ,k)/(R_gas*density(i,j  ,k))
                T_G  = pressure(i,j-1,k)/(R_gas*density(i,j-1,k))
                qp_I = qp(i,j  ,k,2:n_var)
                qp_G = qp(i,j-1,k,2:n_var)

                ! normal component of gradient
                gradqp_x(i,j-1,k,:) = (qp_I - qp_G)*c_x 
                gradqp_y(i,j-1,k,:) = (qp_I - qp_G)*c_y
                gradqp_z(i,j-1,k,:) = (qp_I - qp_G)*c_z
                gradqp_x(i,j-1,k,4) = ( T_I -  T_G)*c_x
                gradqp_y(i,j-1,k,4) = ( T_I -  T_G)*c_y
                gradqp_z(i,j-1,k,4) = ( T_I -  T_G)*c_z
                if(jmin_id==-5 .and. (fixed_wall_temperature(3)<1. .and. fixed_wall_temperature(3)>=0.))then
                gradqp_x(i,j-1,k,4) = -gradqp_x(i,j,k,4)
                gradqp_y(i,j-1,k,4) = -gradqp_y(i,j,k,4)
                gradqp_z(i,j-1,k,4) = -gradqp_z(i,j,k,4)
                end if
                !parallel component of gradient
                do l=1,n
                  dot = (gradqp_x(i,j,k,l)*nx) + (gradqp_y(i,j,k,l)*ny) + (gradqp_z(i,j,k,l)*nz)
                  gradqp_x(i,j-1,k,l) = gradqp_x(i,j-1,k,l) + (gradqp_x(i,j,k,l) - dot*nx)
                  gradqp_y(i,j-1,k,l) = gradqp_y(i,j-1,k,l) + (gradqp_y(i,j,k,l) - dot*ny)
                  gradqp_z(i,j-1,k,l) = gradqp_z(i,j-1,k,l) + (gradqp_z(i,j,k,l) - dot*nz)
                end do
              end do
            end do
          end do

        case('jmax')
          do k = 1,kmx-1
            do j = jmx,jmx
              do i = 1,imx-1
                nx   = ynx(i,j,k)
                ny   = yny(i,j,k)
                nz   = ynz(i,j,k)
                c_x  = yA(i,j,k)*nx/volume(i,j,k)
                c_y  = yA(i,j,k)*ny/volume(i,j,k)
                c_z  = yA(i,j,k)*nz/volume(i,j,k)
                T_I  = pressure(i,j-1,k)/(R_gas*density(i,j-1,k))
                T_G  = pressure(i,j  ,k)/(R_gas*density(i,j  ,k))
                qp_I = qp(i,j-1,k,2:n_var)
                qp_G = qp(i,j  ,k,2:n_var)

                ! normal component of gradient
                gradqp_x(i,j,k,:) = -(qp_I - qp_G)*c_x 
                gradqp_y(i,j,k,:) = -(qp_I - qp_G)*c_y
                gradqp_z(i,j,k,:) = -(qp_I - qp_G)*c_z
                gradqp_x(i,j,k,4) = -( T_I -  T_G)*c_x
                gradqp_y(i,j,k,4) = -( T_I -  T_G)*c_y
                gradqp_z(i,j,k,4) = -( T_I -  T_G)*c_z
                if(jmax_id==-5 .and. (fixed_wall_temperature(4)<1. .and. fixed_wall_temperature(4)>=0.))then
                gradqp_x(i,j,k,4) = -gradqp_x(i,j-1,k,4)
                gradqp_y(i,j,k,4) = -gradqp_y(i,j-1,k,4)
                gradqp_z(i,j,k,4) = -gradqp_z(i,j-1,k,4)
                end if
                !parallel component of gradient
                do l=1,n
                  dot = (gradqp_x(i,j-1,k,l)*nx) + (gradqp_y(i,j-1,k,l)*ny) + (gradqp_z(i,j-1,k,l)*nz)
                  gradqp_x(i,j,k,l) = gradqp_x(i,j,k,l) + (gradqp_x(i,j-1,k,l) - dot*nx)
                  gradqp_y(i,j,k,l) = gradqp_y(i,j,k,l) + (gradqp_y(i,j-1,k,l) - dot*ny)
                  gradqp_z(i,j,k,l) = gradqp_z(i,j,k,l) + (gradqp_z(i,j-1,k,l) - dot*nz)
                end do
              end do
            end do
          end do


        case('kmin')
          do k = 1,1
            do j = 1,jmx-1
              do i = 1,imx-1
                nx   = znx(i,j,k)
                ny   = zny(i,j,k)
                nz   = znz(i,j,k)
                c_x  = zA(i,j,k)*nx/volume(i,j,k)
                c_y  = zA(i,j,k)*ny/volume(i,j,k)
                c_z  = zA(i,j,k)*nz/volume(i,j,k)
                T_I  = pressure(i,j,k  )/(R_gas*density(i,j,k  ))
                T_G  = pressure(i,j,k-1)/(R_gas*density(i,j,k-1))
                qp_I = qp(i,j,k  ,2:n_var)
                qp_G = qp(i,j,k-1,2:n_var)

                ! normal component of gradient
                gradqp_x(i,j,k-1,:) = (qp_I - qp_G)*c_x 
                gradqp_y(i,j,k-1,:) = (qp_I - qp_G)*c_y
                gradqp_z(i,j,k-1,:) = (qp_I - qp_G)*c_z
                gradqp_x(i,j,k-1,4) = ( T_I -  T_G)*c_x
                gradqp_y(i,j,k-1,4) = ( T_I -  T_G)*c_y
                gradqp_z(i,j,k-1,4) = ( T_I -  T_G)*c_z
                if(kmin_id==-5 .and. (fixed_wall_temperature(5)<1. .and. fixed_wall_temperature(5)>=0.))then
                gradqp_x(i,j,k-1,4) = -gradqp_x(i,j,k,4)
                gradqp_y(i,j,k-1,4) = -gradqp_y(i,j,k,4)
                gradqp_z(i,j,k-1,4) = -gradqp_z(i,j,k,4)
                end if
                !parallel component of gradient
                do l=1,n
                  dot = (gradqp_x(i,j,k,l)*nx) + (gradqp_y(i,j,k,l)*ny) + (gradqp_z(i,j,k,l)*nz)
                  gradqp_x(i,j,k-1,l) = gradqp_x(i,j,k-1,l) + (gradqp_x(i,j,k,l) - dot*nx)
                  gradqp_y(i,j,k-1,l) = gradqp_y(i,j,k-1,l) + (gradqp_y(i,j,k,l) - dot*ny)
                  gradqp_z(i,j,k-1,l) = gradqp_z(i,j,k-1,l) + (gradqp_z(i,j,k,l) - dot*nz)
                end do
              end do
            end do
          end do

        case('kmax')
          do k = kmx,kmx
            do j = 1,jmx-1
              do i = 1,imx-1
                nx   = znx(i,j,k)
                ny   = zny(i,j,k)
                nz   = znz(i,j,k)
                c_x  = zA(i,j,k)*nx/volume(i,j,k)
                c_y  = zA(i,j,k)*ny/volume(i,j,k)
                c_z  = zA(i,j,k)*nz/volume(i,j,k)
                T_I  = pressure(i,j,k-1)/(R_gas*density(i,j,k-1))
                T_G  = pressure(i,j,k  )/(R_gas*density(i,j,k  ))
                qp_I = qp(i,j,k-1,2:n_var)
                qp_G = qp(i,j,k  ,2:n_var)

                ! normal component of gradient
                gradqp_x(i,j,k,:) = -(qp_I - qp_G)*c_x 
                gradqp_y(i,j,k,:) = -(qp_I - qp_G)*c_y
                gradqp_z(i,j,k,:) = -(qp_I - qp_G)*c_z
                gradqp_x(i,j,k,4) = -( T_I -  T_G)*c_x
                gradqp_y(i,j,k,4) = -( T_I -  T_G)*c_y
                gradqp_z(i,j,k,4) = -( T_I -  T_G)*c_z
                if(kmax_id==-5 .and. (fixed_wall_temperature(6)<1. .and. fixed_wall_temperature(6)>=0.))then
                gradqp_x(i,j,k,4) = -gradqp_x(i,j,k-1,4)
                gradqp_y(i,j,k,4) = -gradqp_y(i,j,k-1,4)
                gradqp_z(i,j,k,4) = -gradqp_z(i,j,k-1,4)
                end if
                !parallel component of gradient
                do l=1,n
                  dot = (gradqp_x(i,j,k-1,l)*nx) + (gradqp_y(i,j,k-1,l)*ny) + (gradqp_z(i,j,k-1,l)*nz)
                  gradqp_x(i,j,k,l) = gradqp_x(i,j,k,l) + (gradqp_x(i,j,k-1,l) - dot*nx)
                  gradqp_y(i,j,k,l) = gradqp_y(i,j,k,l) + (gradqp_y(i,j,k-1,l) - dot*ny)
                  gradqp_z(i,j,k,l) = gradqp_z(i,j,k,l) + (gradqp_z(i,j,k-1,l) - dot*nz)
                end do
              end do
            end do
          end do

        case DEFAULT
          print*, "Ghost gradients : Wrong face name string"

      end select
          
    end subroutine apply


 
end module ghost_gradients
