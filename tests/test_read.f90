integer function test_read() result(r)
   use vartypes
   use read
    
   type(filetype) :: files
   type(controltype):: control
   type(schemetype) :: scheme
   type(flowtype)   :: flow
   character(len=*), parameter :: prefix="../../tests/files/"

   files%control_file = prefix//"control.md"
   files%scheme_file = prefix//"fvscheme.md"
   files%flow_file = prefix//"flow.md"
   files%outin_file = prefix//"output_control.md"
   files%res_control_file = prefix//"res_control.md"

   r = 1
   call read_input_and_controls(files, control, scheme, flow)
   if( control%CFL==100 .and. scheme%scheme_name=='ausm' .and. scheme%turbulence=='sst' .and. flow%T_ref==300)then
       r = 0
   end if

end function
