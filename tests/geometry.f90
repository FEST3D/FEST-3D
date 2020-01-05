  subroutine compute_i_face_area(imx, jmx, kmx, grid_x, grid_y, grid_z, xn)
      !-----------------------------------------------------------
      ! Compute face area vectors
      !
      ! The face area vectors denote the face area both in 
      ! magnitude and direction. They are placed in the face 
      ! normal variables for further calculation.
      !
      ! The face normal is given by d1 x d2, where d1 and d2 are
      ! the diagonals of a face
      !-----------------------------------------------------------
      
      implicit none
      integer , intent(in) :: imx, jmx, kmx
      real*8 , dimension(1:imx, 1:jmx, 1:kmx), intent(in) :: grid_x, grid_y, grid_z
      real*8 , dimension(1:imx, 1:jmx-1, 1:kmx-1,1:3), intent(out) :: xn

      real*8 :: d1x, d2x, d1y, d2y, d1z, d2z
      integer :: i, j, k

      !f2py intent(in) :: imx, jmx, kmx
      !f2py intent(in) :: grid_x, grid_y, grid_z
      !f2py intent(out) :: xn

      do k = 1, kmx-1
       do j = 1, jmx-1
        do i = 1, imx
          d1x = grid_x(i, j+1, k+1) - grid_x(i, j, k)
          d1y = grid_y(i, j+1, k+1) - grid_y(i, j, k)
          d1z = grid_z(i, j+1, k+1) - grid_z(i, j, k)
          d2x = grid_x(i, j, k+1) - grid_x(i, j+1, k)
          d2y = grid_y(i, j, k+1) - grid_y(i, j+1, k)
          d2z = grid_z(i, j, k+1) - grid_z(i, j+1, k)
          xn(i, j, k,1) = 0.5 * (d1y*d2z - d1z*d2y)
          xn(i, j, k,2) = 0.5 * (d1z*d2x - d1x*d2z)
          xn(i, j, k,3) = 0.5 * (d1x*d2y - d1y*d2x)
         end do
        end do
       end do

  end subroutine compute_i_face_area
    

  subroutine compute_j_face_area(imx, jmx, kmx, grid_x, grid_y, grid_z, yn)
      implicit none
      integer , intent(in) :: imx, jmx, kmx
      real*8 , dimension(1:imx, 1:jmx, 1:kmx), intent(in) :: grid_x, grid_y, grid_z
      real*8 , dimension(1:imx-1, 1:jmx, 1:kmx-1,1:3), intent(out) :: yn

      real*8 :: d1x, d2x, d1y, d2y, d1z, d2z
      integer :: i, j, k

      !f2py intent(in) :: imx, jmx, kmx
      !f2py intent(in) :: grid_x, grid_y, grid_z
      !f2py intent(out) :: yn

      do k = 1, kmx-1
       do j = 1, jmx
        do i = 1, imx-1
         d1x = grid_x(i+1, j, k+1) - grid_x(i, j, k)
         d1y = grid_y(i+1, j, k+1) - grid_y(i, j, k)
         d1z = grid_z(i+1, j, k+1) - grid_z(i, j, k)
         d2x = grid_x(i+1, j, k) - grid_x(i, j, k+1)
         d2y = grid_y(i+1, j, k) - grid_y(i, j, k+1)
         d2z = grid_z(i+1, j, k) - grid_z(i, j, k+1)
      
         yn(i, j, k,1) = 0.5 * (d1y*d2z - d1z*d2y)
         yn(i, j, k,2) = 0.5 * (d1z*d2x - d1x*d2z)
         yn(i, j, k,3) = 0.5 * (d1x*d2y - d1y*d2x)
        end do
       end do
      end do
  end subroutine compute_j_face_area
    

  subroutine compute_k_face_area(imx, jmx, kmx, grid_x, grid_y, grid_z, zn)
      implicit none
      integer , intent(in) :: imx, jmx, kmx
      real*8 , dimension(1:imx, 1:jmx, 1:kmx), intent(in) :: grid_x, grid_y, grid_z
      real*8 , dimension(1:imx-1, 1:jmx-1, 1:kmx,1:3), intent(out) :: zn

      real*8 :: d1x, d2x, d1y, d2y, d1z, d2z
      integer :: i, j, k

      !f2py intent(in) :: imx, jmx, kmx
      !f2py intent(in) :: grid_x, grid_y, grid_z
      !f2py intent(out) :: zn
      
      do k = 1, kmx
       do j = 1, jmx-1
        do i = 1, imx-1
         d1x = grid_x(i+1, j+1, k) - grid_x(i, j, k)
         d1y = grid_y(i+1, j+1, k) - grid_y(i, j, k)
         d1z = grid_z(i+1, j+1, k) - grid_z(i, j, k)
         d2x = grid_x(i, j+1, k) - grid_x(i+1, j, k)
         d2y = grid_y(i, j+1, k) - grid_y(i+1, j, k)
         d2z = grid_z(i, j+1, k) - grid_z(i+1, j, k)
      
         zn(i, j, k,1) = 0.5 * (d1y*d2z - d1z*d2y)
         zn(i, j, k,2) = 0.5 * (d1z*d2x - d1x*d2z)
         zn(i, j, k,3) = 0.5 * (d1x*d2y - d1y*d2x)
        end do
       end do
      end do
  end subroutine compute_k_face_area
    

  subroutine vol_tetrahedron(p1, p2, p3, p4, result)
      !-----------------------------------------------------------
      ! Compute the volume of a tetrahedron, given 4 points which
      ! are 1-D arrays
      ! Since we know that the determinant is to be evaluated of 
      ! a 3x3 matrix, we write the expression itself
      !-----------------------------------------------------------

      implicit none
      real*8, dimension(1:3), intent(in):: p1, p2, p3, p4
      real*8, dimension(1:3,1:3) :: A
      real*8 :: result

      !f2py intent(in) :: p1, p2, p3, p4
      !f2py intent(out) :: result

      A(1:3, 1) = p1 - p4
      A(1:3, 2) = p2 - p4
      A(1:3, 3) = p3 - p4

      result = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) + &
                        A(1,2) * (A(2,3)*A(3,1) - A(2,1)*A(3,3)) + &
                        A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      result = result / 6.                  
  
  end subroutine vol_tetrahedron
  
  subroutine vol_hexahedron(p_list, volume)
      implicit none
      real*8, intent(out) :: volume
      real*8, dimension(1:3, 1:8), intent(in) :: p_list
      real*8 :: vol_hex
      real*8 :: temp

      !f2py intent(in) :: p_list
      !f2py intent(out) :: volume
      
      vol_hex = 0.
      call vol_tetrahedron(p_list(1:3,1), p_list(1:3,5), &
                       p_list(1:3,8), p_list(1:3,6),temp)
      vol_hex = vol_hex + temp
                  call  vol_tetrahedron(p_list(1:3,7), p_list(1:3,8), &
                                       p_list(1:3,6), p_list(1:3,3), temp)
      vol_hex = vol_hex + temp
                  call  vol_tetrahedron(p_list(1:3,8), p_list(1:3,4), &
                                       p_list(1:3,1), p_list(1:3,3), temp)
      vol_hex = vol_hex + temp
                  call  vol_tetrahedron(p_list(1:3,6), p_list(1:3,1), &
                                       p_list(1:3,3), p_list(1:3,8), temp)
      vol_hex = vol_hex + temp
                  call  vol_tetrahedron(p_list(1:3,1), p_list(1:3,2), &
                                       p_list(1:3,6), p_list(1:3,3), temp)
      vol_hex = vol_hex + temp
      volume = vol_hex
      
  end subroutine vol_hexahedron
  
  subroutine compute_volume(imx, jmx, kmx, grid_x, grid_y, grid_z, volume)
      !-----------------------------------------------------------
      ! Compute the grid cell volumes
      ! Each grid is a hexahedron, whose volume is calculated by
      ! splitting it into 5 tetrahedrons, whose volume is known
      !-----------------------------------------------------------

      implicit none
      integer, intent(in) :: imx
      integer, intent(in) :: jmx
      integer, intent(in) :: kmx
      real*8 , dimension(imx, jmx, kmx), intent(in) :: grid_x
      real*8 , dimension(imx, jmx, kmx), intent(in) :: grid_y
      real*8 , dimension(imx, jmx, kmx), intent(in) :: grid_z
      real*8 , dimension(imx-1, jmx-1, kmx-1), intent(out) :: volume

      integer :: i,j,k
      real*8, dimension(1:3, 1:8) :: p_list

      !f2py intent(in) :: imx, jmx, kmx, grid_x, grid_y, grid_z
      !f2py intent(out) :: volume

      volume=1.
      do k = 1, kmx-1
          do j = 1, jmx-1
              do i = 1, imx-1
                  p_list(:, :) = 0.
                  p_list(:, 1) = (/ grid_x(i,j,k), grid_y(i,j,k), grid_z(i,j,k) /)
                  p_list(:, 2) = (/ grid_x(i+1,j,k), grid_y(i+1,j,k), grid_z(i+1,j,k) /)
                  p_list(:, 3) = (/ grid_x(i+1,j+1,k), grid_y(i+1,j+1,k), grid_z(i+1,j+1,k) /)
                  p_list(:, 4) = (/ grid_x(i,j+1,k), grid_y(i,j+1,k), grid_z(i,j+1,k) /)
                  p_list(:, 5) = (/ grid_x(i,j,k+1), grid_y(i,j,k+1), grid_z(i,j,k+1) /)
                  p_list(:, 6) = (/ grid_x(i+1,j,k+1), grid_y(i+1,j,k+1), grid_z(i+1,j,k+1) /)
                  p_list(:, 7) = (/ grid_x(i+1,j+1,k+1), grid_y(i+1,j+1,k+1), grid_z(i+1,j+1,k+1) /)
                  p_list(:, 8) = (/ grid_x(i,j+1,k+1), grid_y(i,j+1,k+1), grid_z(i,j+1,k+1) /)
                  !volume(i, j, k) = (vol_hexahedron(p_list))
                  call vol_hexahedron(p_list, volume(i,j,k))
              end do
          end do
      end do
      
  end subroutine compute_volume

  subroutine compute_centroid(imx, jmx, kmx, grid_x, grid_y, grid_z, CellCenter)
    implicit none
      integer, intent(in) :: imx
      integer, intent(in) :: jmx
      integer, intent(in) :: kmx
      real*8 , dimension(1:imx, 1:jmx, 1:kmx), intent(in) :: grid_x
      real*8 , dimension(1:imx, 1:jmx, 1:kmx), intent(in) :: grid_y
      real*8 , dimension(1:imx, 1:jmx, 1:kmx), intent(in) :: grid_z
      real*8 , dimension(1:imx-1, 1:jmx-1, 1:kmx-1,1:3), intent(out) :: CellCenter
    integer :: i,j,k
      !f2py intent(in) :: imx, jmx, kmx, grid_x, grid_y, grid_z
      !f2py intent(out) :: CellCenter
    !-----------------------------------
    ! compute cell center of all cell including ghost cells
    !------------------------------------------------------

    do k = 1, kmx-1
      do j = 1, jmx-1
        do i = 1, imx-1
          CellCenter(i,j,k,1) = 0.125 * ( grid_x(i  ,j  ,k  ) &
                                        + grid_x(i+1,j  ,k  ) &
                                        + grid_x(i+1,j+1,k  ) &
                                        + grid_x(i+1,j+1,k+1) &
                                        + grid_x(i+1,j  ,k+1) &
                                        + grid_x(i  ,j+1,k  ) &
                                        + grid_x(i  ,j+1,k+1) &
                                        + grid_x(i  ,j  ,k+1) &
                                        )

          CellCenter(i,j,k,2) = 0.125 * ( grid_y(i  ,j  ,k  ) &
                                        + grid_y(i+1,j  ,k  ) &
                                        + grid_y(i+1,j+1,k  ) &
                                        + grid_y(i+1,j+1,k+1) &
                                        + grid_y(i+1,j  ,k+1) &
                                        + grid_y(i  ,j+1,k  ) &
                                        + grid_y(i  ,j+1,k+1) &
                                        + grid_y(i  ,j  ,k+1) &
                                        )

          CellCenter(i,j,k,3) = 0.125 * ( grid_z(i  ,j  ,k  ) &
                                        + grid_z(i+1,j  ,k  ) &
                                        + grid_z(i+1,j+1,k  ) &
                                        + grid_z(i+1,j+1,k+1) &
                                        + grid_z(i+1,j  ,k+1) &
                                        + grid_z(i  ,j+1,k  ) &
                                        + grid_z(i  ,j+1,k+1) &
                                        + grid_z(i  ,j  ,k+1) &
                                        )
        end do
      end do
    end do

  end subroutine compute_centroid

