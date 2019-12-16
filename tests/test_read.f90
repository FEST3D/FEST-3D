integer function test_read() result(r)
   use vartypes
   use read
    
   type(extent) :: dims
   type(controltype):: control
   type(schemetype) :: scheme
   type(flowtype)   :: flow
   character(len=*), parameter :: prefix="../../tests/files"

   r = 1
   call read_input_and_controls(control, scheme, flow)
   if( control%CFL==100)then
       r = 0
   end if

end function
