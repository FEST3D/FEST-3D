integer function test_geometry() result(r)
   use vartypes
   use grid,  only : setup_grid
   use geometry, only : setup_geometry
    
   type(nodetype), dimension(:,:,:), allocatable :: points
   type(extent) :: dims
   type(celltype), dimension(:,:,:), allocatable :: cell
   type(facetype), dimension(:,:,:), allocatable :: Iface, Jface, Kface
   type(boundarytype) :: bc
   character(len=*), parameter :: prefix="../../tests/"

   dims%imx = -4
   dims%jmx = -4
   dims%kmx = -4
   allocate(points(-2:dims%imx+3, -2:dims%jmx+3, -2:dims%kmx+3))
   points(-2,-2,-2)%x = 0.0
   points(-1,-2,-2)%x = 1.0
   points(-2,-1,-2)%x = 0.0
   points(-1,-1,-2)%x = 1.0
   points(-2,-2,-1)%x = 0.0
   points(-1,-2,-1)%x = 1.0
   points(-2,-1,-1)%x = 0.0
   points(-1,-1,-1)%x = 1.0
   points(-2,-2,-2)%y = 0.0
   points(-1,-2,-2)%y = 0.0
   points(-2,-1,-2)%y = 1.0
   points(-1,-1,-2)%y = 1.0
   points(-2,-2,-1)%y = 0.0
   points(-1,-2,-1)%y = 0.0
   points(-2,-1,-1)%y = 1.0
   points(-1,-1,-1)%y = 1.0
   points(-2,-2,-2)%z = 0.0
   points(-1,-2,-2)%z = 0.0
   points(-2,-1,-2)%z = 0.0
   points(-1,-1,-2)%z = 1.0
   points(-2,-2,-1)%z = 1.0
   points(-1,-2,-1)%z = 1.0
   points(-2,-1,-1)%z = 1.0
   points(-1,-1,-1)%z = 1.0

   r = 1
   call setup_geometry(cell, Iface, Jface, Kface, points, bc, dims)
   if( cell(-2,-2,-2)%volume==1. .and. Iface(-2,-2,-2)%nx==1 .and. Iface(-2,-2,-2)%ny==0)then
       r = 0
   end if

end function
