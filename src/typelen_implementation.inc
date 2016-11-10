function intlen(i) result(l)
    !-------------------------------------------------------------------
    ! Return the "string length of an integer"
    !
    ! The string length of an integer is defined as the length of the 
    ! string containing only that integer (computed using len_trim).
    !-------------------------------------------------------------------
    
    implicit none
    integer, intent(in) :: i
    integer :: l

    l = len_trim(tostr(i))

end function intlen

function reallen(r) result(l)
    !-------------------------------------------------------------------
    ! Return the "string length of a real"
    !
    ! The string length of a real is defined as the length of the string
    ! containing only that real (computed using len_trim).
    !
    ! Use the chfmt subroutine to change the format specifiers used to 
    ! convert the real into a string.
    !-------------------------------------------------------------------

    implicit none
    real, intent(in) :: r
    integer :: l

    l = len_trim(real_to_str(r))

end function reallen

