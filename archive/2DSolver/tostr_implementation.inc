function int_to_str(i) result(a)

    implicit none
    integer, intent(in) :: i
    character(len=MAX_STRING_LEN) :: a

    write(a, *) i

end function int_to_str

function real_to_str(r) result(a)

    implicit none
    real, intent(in) :: r
    character(len=MAX_STRING_LEN) :: a

    select case (form_)
        case ('F')
            write(a, '(F0.' + dec_ + ')') r
        case ('E')
            if (exp_ .ne. -1) then
                write(a, '(E0.' + dec_ + 'E' + exp_ + ')') r
            else
                write(a, '(E0.' + dec_ + ')') r
            end if
        case ('N')
            if (exp_ .ne. -1) then
                write(a, '(EN0.' + dec_ + 'E' + exp_ + ')') r
            else
                write(a, '(EN0.' + dec_ + ')') r
            end if
        case ('S')
            if (exp_ .ne. -1) then
                write(a, '(ES0.' + dec_ + 'E' + exp_ + ')') r
            else
                write(a, '(ES0.' + dec_ + ')') r
            end if
        case default
            stop
    end select

end function real_to_str

function bool_to_str(b) result(a)

    implicit none
    logical, intent(in) :: b
    character(len=MAX_STRING_LEN) :: a

    write(a, *) b

end function bool_to_str

