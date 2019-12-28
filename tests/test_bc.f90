integer function test_bc() result(r)
   use vartypes
   use bc
    
   type(filetype)  :: files
   type(controltype) :: control
   type(schemetype) :: scheme
   type(boundarytype) :: boundary
   type(flowtype) :: flow
   type(extent) :: dims
   character(len=*), parameter :: prefix="../../tests/files/"

   dims%imx = 2
   dims%jmx = 2
   dims%kmx = 2

   dims%n_var = 5

   boundary%imin_id=1
   boundary%imax_id=1
   boundary%jmin_id=-5
   boundary%jmax_id=1
   boundary%kmin_id=1
   boundary%kmax_id=1

   files%bcfile=prefix//"bc_00.md"

   r = 1
   call setup_bc(files, scheme, flow, boundary, dims)

   print*, boundary%fixed_pressure(:)
   if(boundary%make_G_flux_zero(1)==0 .and. boundary%id(2)==1 .and. boundary%fixed_pressure(1)<101325)then
    r = 0
   end if

end function
