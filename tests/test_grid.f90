!integer function tests_test_grid() result(r)
integer function test_grid() result(r)
   use vartypes
   use grid,  only : setup_grid
   !use geometry, only: setup_geometry
    
   type(nodetype), dimension(:,:,:), allocatable :: points
   type(extent) :: dims
   character(len=*), parameter :: prefix="../../tests/"

   !type(celltype), dimension(:,:,:), allocatable :: cell
   !type(facetype), dimension(:,:,:), allocatable :: Iface, Jface, Kface

   r = 1
   call setup_grid(prefix//"grid.txt", prefix//"mapping.txt", prefix//"periodic.txt", points, dims)
   !call setup_geometry(cell, Iface, Jface, Kface, points, dims)
   if( dims%imx==5 .and.dims%jmx==2 .and. dims%kmx==2)then
     if (points(1,1,1)%x==-0.002 .and. points(5,1,1)%x==0.002)then
       r = 0
     end if
   end if

end function
