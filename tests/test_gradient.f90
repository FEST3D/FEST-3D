integer function test_gradient() result(r)
   use vartypes
   use utils,
    
   type(filetype) :: files
   type(controltype) :: control
   type(schemetype) :: scheme
   type(flowtype) :: flow
   type(nodetype), dimension(:,:,:), allocatable :: points
   type(extent) :: dims
   character(len=*), parameter :: prefix="../../tests/files/"

   r = 1
   files%gridfile = prefix//"grid.txt"
   files%mapfile = prefix//"mapping.txt"
   files%periodicfile = prefix//"periodic.txt"
   call setup_grid(files, points, control, dims)
   if( dims%imx==5 .and.dims%jmx==2 .and. dims%kmx==2)then
     if (points(1,1,1)%x==-0.002 .and. points(5,1,1)%x==0.002)then
       r = 0
     end if
   end if

end function
