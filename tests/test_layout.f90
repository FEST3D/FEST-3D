integer function test_layout() result(r)
   use vartypes
   use layout
   use mpi
    
   type(filetype) :: files
   type(boundarytype) :: bc
   type(controltype) :: control
   character(len=*), parameter :: prefix="../../tests/files/"

   r = 1
   control%process_id=0
   control%total_process=2
   files%layout_file=prefix//"layout.md"
   call read_layout_file(files, control,bc)

   if(files%bcfile=="system/mesh/bc/bc_00.md")then
     r=0
   end if

end function
