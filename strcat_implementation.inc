!-----------------------------------------------------------------------
! String concatenation functions
!-----------------------------------------------------------------------

function str_cat_str(s1, s2) result(a)

    implicit none
    character(len=*), intent(in) :: s1
    character(len=*), intent(in) :: s2
    character(len=MAX_STRING_LEN) :: a

    write(a, '(2A)') trim(s1), trim(s2)

end function str_cat_str

function str_cat_int(s, i) result(a)

    implicit none
    character(len=*), intent(in) :: s
    integer, intent(in) :: i
    character(len=MAX_STRING_LEN) :: a

    write(a, '(A,I0)') trim(s), i

end function str_cat_int

function int_cat_str(i, s) result(a)

    implicit none
    character(len=*), intent(in) :: s
    integer, intent(in) :: i
    character(len=MAX_STRING_LEN) :: a

    write(a, '(I0,A)') i, trim(s)

end function int_cat_str

function str_cat_real(s, r) result(a)

    implicit none
    character(len=*), intent(in) :: s
    real, intent(in) :: r
    character(len=MAX_STRING_LEN) :: a

    write(a, '(2A)') trim(s), trim(tostr(r))

end function str_cat_real

function real_cat_str(r, s) result(a)

    implicit none
    character(len=*), intent(in) :: s
    real, intent(in) :: r
    character(len=MAX_STRING_LEN) :: a

    write(a, *) trim(tostr(r)), trim(s)

end function real_cat_str

function str_cat_bool(s, b) result(a)

    implicit none
    character(len=*), intent(in) :: s
    logical, intent(in) :: b
    character(len=MAX_STRING_LEN) :: a

    write(a, '(A,L1)') trim(s), b

end function str_cat_bool

function bool_cat_str(b, s) result(a)

    implicit none
    character(len=*), intent(in) :: s
    logical, intent(in) :: b
    character(len=MAX_STRING_LEN) :: a

    write(a, '(A,L1)') b, trim(s)

end function bool_cat_str

