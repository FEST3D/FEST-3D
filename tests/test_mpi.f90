integer function test_mpi() result(r)
   use vartypes
   use layout
   use mpi
    
   integer :: ierr
   type(controltype) :: control

   r = 1
   control%process_id=-1
   control%total_process=-1
   call MPI_INIT(ierr)
   call get_process_data(control)

   if(control%process_id>=0 .and. control%total_process>=0)then
     r=0
   end if

end function
