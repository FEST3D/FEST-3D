module string
    !-------------------------------------------------------------------
    ! String manipulation
    !
    ! This module implements the to_string family of methods which allow
    ! other type variables to be converted to strings.
    !
    ! This module also overloads the string concatenation operator (//) 
    ! so that strings can be concatenated with other strings, 
    ! integers, reals and logicals.
    !-------------------------------------------------------------------

    implicit none
    integer, parameter :: MAX_STRING_LEN = 256
    integer, private :: dec_ = 6
    integer, private :: exp_ = -1
    character, private :: form_ = 'F'
    logical, private :: reset_flag = .FALSE.
    !TODO: Implement a reset flag. If this flag is set to true, reset 
    !the parameters dec, exp, form after one operation (can be changed 
    !to some other integer?). This should be set by the chfmt method.

    interface tostr
        module procedure int_to_str, &
                real_to_str, &
                bool_to_str
    end interface tostr

    interface operator( + )
        module procedure str_cat_str, &
                str_cat_int, &
                int_cat_str, &
                str_cat_real, &
                real_cat_str, &
                str_cat_bool, &
                bool_cat_str
    end interface operator( + )

    interface len
        module procedure intlen, reallen
    end interface len

    contains

        include "tostr_implementation.inc"
        include "strcat_implementation.inc"
        include "typelen_implementation.inc"

        subroutine chfmt(d, e, f)
            !-----------------------------------------------------------
            ! Change format specifier for reals
            !
            ! Inputs:
            !   d -> integer, optional
            !       digits after decimal
            !   e -> integer, optional
            !       digits in exponent
            !       When this is set to -1, the exponent part is dropped
            !   f -> character, optional
            !       form
            !       Valid options: 'F' (decimal), 'E' (exponential),
            !           'S' (scientific), 'N' (engineering)
            !       If 'F' is specified, then the value for e is 
            !           overridden and set to -1.
            !
            ! If no arguments are passed, this function resets all the 
            ! parameters to their default values (as provided in the 
            ! following table. Else, only the passed arguments will be
            ! updated (the others will be left unchanged). 
            !
            ! This function is sticky; a format once set will continue 
            ! to apply till either it is changed or the program ends.
            !
            ! Default values: 
            !   d (digits after decimal) --> 6
            !   e (digits in exponent) --> -1
            !   f (form) --> 'F' (decimal)
            !
            !TODO: Add support for width also?
            !http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
            !-----------------------------------------------------------

            implicit none
            integer, intent(in), optional :: d
            integer, intent(in), optional :: e
            character, intent(in), optional :: f

            if (.not. (present(d) .or. present(e) .or. present(f))) then
                dec_ = 6
                exp_ = 0
                form_ = 'F'
            else
                if (present(d)) dec_ = d
                if (present(e)) exp_ = e
                if (present(f)) then
                    if (.not. (f .eq. 'F' .or. f .eq. 'E' .or. &
                            f .eq. 'S' .or. f .eq. 'N')) then
                        print *, 'Error: Unknown kind specified.'
                        stop
                    end if
                    form_ = f
                    if (form_ .eq. 'F') exp_ = -1
                end if
            end if

        end subroutine chfmt

        subroutine disp(s)
            !-----------------------------------------------------------
            ! Display the contents of the string
            !
            ! This function trims the string before printing it.
            !-----------------------------------------------------------
        
            implicit none
            character(len=MAX_STRING_LEN), intent(in) :: s

            print *, trim(s)

        end subroutine disp

end module string
