  !< All the bitwise operation
module bitwise
  !< To apply bitwise (and) and (or) to integer which represents the binary or oct number

    implicit none
    
    interface operator(.and.)
        module procedure int4_and_int4, &
                int8_and_int8
    end interface operator(.and.)

    interface operator(.or.)
        module procedure int4_or_int4, int4_1D_or_int4_1D, &
                int4_1D_or_int4, int4_2D_or_int4_2D, &
                int4_2D_or_int4, int8_or_int8
    end interface operator(.or.)

    interface assignment(=)
        module procedure int4_from_string, &
                int8_from_string
    end interface assignment(=)

    interface bin2int
        module procedure bin_to_int4, bin_to_int8
    end interface bin2int

    interface oct2int
        module procedure oct_to_int4, oct_to_int8
    end interface oct2int

    contains

        function int4_and_int4(a, b) result(r)
          !< bitwise 'AND' over two integer of kind 4
            implicit none
            integer(kind=4), intent(in) :: a
            integer(kind=4), intent(in) :: b
            integer(kind=4) :: r
            r = iand(a, b)
        end function int4_and_int4

        function int8_and_int8(a, b) result(r)
          !< bitwise 'AND' over two integer of kind 8
            implicit none
            integer(kind=8), intent(in) :: a
            integer(kind=8), intent(in) :: b
            integer(kind=8) :: r
            r = iand(a, b)
        end function int8_and_int8

        function int4_or_int4(a, b) result(r)
          !< bitwise 'OR' over two integer of kind 4
            implicit none
            integer(kind=4), intent(in) :: a
            integer(kind=4), intent(in) :: b
            integer(kind=4) :: r
            r = ior(a, b)
        end function int4_or_int4

        function int4_1D_or_int4_1D(a, b) result(r)
          !< bitwise 'OR' over two 1D integer array of kind 4
            implicit none
            integer(kind=4), dimension(:), intent(in) :: a
            integer(kind=4), dimension(:), intent(in) :: b
            integer(kind=4), dimension(size(a)) :: r
            if (size(a) /= size(b)) then
                print *, "Error: Sizes of arrays being 'or'ed should be the same."
                stop
            end if
            r = ior(a, b)
        end function int4_1D_or_int4_1D

        function int4_1D_or_int4(a, b) result(r)
          !< bitwise 'OR' over one 1D integer array and integer of kind 4
            implicit none
            integer(kind=4), dimension(:), intent(in) :: a
            integer(kind=4), intent(in) :: b
            integer(kind=4), dimension(size(a)) :: r
            integer :: i
            i = 1
            do while (i <= size(a))
                r(i) = ior(a(i), b)
                i = i + 1
            end do
        end function int4_1D_or_int4

        function int4_2D_or_int4_2D(a, b) result(r)
          !< bitwise 'OR' over two 2D integer array of kind 4
            implicit none
            integer(kind=4), dimension(:, :), intent(in) :: a
            integer(kind=4), dimension(:, :), intent(in) :: b
            integer(kind=4), dimension(1:2) :: na, nb
            integer(kind=4), dimension(:, :), allocatable :: r
            na = shape(a)
            nb = shape(b)
            if ((na(1) /= nb(1)) .or. (na(2) /= nb(2))) then
                print *, "Error: Sizes of arrays being 'or'ed should be the same."
                stop
            end if
            allocate(r(1:na(1), 1:na(2)))
            r = ior(a, b)
        end function int4_2D_or_int4_2D

        function int4_2D_or_int4(a, b) result(r)
          !< bitwise 'OR' over one 2D integer array and integer of kind 4
            implicit none
            integer(kind=4), dimension(:, :), intent(in) :: a
            integer(kind=4), intent(in) :: b
            integer(kind=4), dimension(1:2) :: n
            integer(kind=4), dimension(:, :), allocatable :: r
            integer :: i, j
            n = shape(a)
            allocate(r(1:n(1), 1:n(2)))
            do j = 1, n(2)
             do i = 1, n(1)
                r(i,j) = ior(a(i, j), b)
             end do
            end do
        end function int4_2D_or_int4

        function int8_or_int8(a, b) result(r)
          !< bitwise 'OR' over two integer of kind 8
            implicit none
            integer(kind=8), intent(in) :: a
            integer(kind=8), intent(in) :: b
            integer(kind=8) :: r
            r = ior(a, b)
        end function int8_or_int8

        subroutine bin_to_int4(r, binstr)
          !< String of binary number converted to integer of kind 4
            implicit none
            character(len=*) :: binstr
            integer(kind=4), intent(out) :: r
            integer :: current_digit
            integer :: i
            i = len(binstr)
            r = 0
            do while (i > 0)
                read (binstr(i:i), *) current_digit
                r = r + ((2 ** (len(binstr) - i)) * current_digit)
                i = i - 1
            end do
        end subroutine bin_to_int4
                
        subroutine bin_to_int8(r, binstr)
          !< String of binary number converted to integer of kind 8
            implicit none
            character(len=*) :: binstr
            integer(kind=8), intent(out) :: r
            integer :: current_digit
            integer :: i
            i = len(binstr)
            r = 0
            do while (i > 0)
                read (binstr(i:i), *) current_digit
                r = r + ((2 ** (len(binstr) - i)) * current_digit)
                i = i - 1
            end do
        end subroutine bin_to_int8
                
        subroutine oct_to_int4(r, octstr)
          !< String of octal number converted to integer of kind 4
            implicit none
            character(len=*) :: octstr
            integer(kind=4), intent(out) :: r
            integer :: current_digit
            integer :: i
            i = len(octstr)
            r = 0
            do while (i > 0)
                read (octstr(i:i), *) current_digit
                r = r + ((8 ** (len(octstr) - i)) * current_digit)
                i = i - 1
            end do
        end subroutine oct_to_int4

        subroutine oct_to_int8(r, octstr)
          !< String of octal number converted to integer of kind 8
            implicit none
            character(len=*) :: octstr
            integer(kind=8), intent(out) :: r
            integer :: current_digit
            integer :: i
            i = len(octstr)
            r = 0
            do while (i > 0)
                read (octstr(i:i), *) current_digit
                r = r + ((8 ** (len(octstr) - i)) * current_digit)
                i = i - 1
            end do
        end subroutine oct_to_int8

        subroutine int4_from_string(lhs, rhs)
          !< Get integer of kind 4 from the string which might contain either octal or binary number
            implicit none
            integer(kind=4), intent(out) :: lhs
            character(len=*), intent(in) :: rhs
            character(len=len(rhs)-1) :: temp
            temp = rhs(2:)
            if (rhs(1:1) == 'b') then
                call bin2int(lhs, temp)
            else if (rhs(1:1) == 'o') then
                call oct2int(lhs, temp)
            end if
        end subroutine int4_from_string

        subroutine int8_from_string(lhs, rhs)
          !< Get integer of kind 8 from the string which might contain either octal or binary number
            implicit none
            integer(kind=8), intent(out) :: lhs
            character(len=*), intent(in) :: rhs
            character(len=len(rhs)-1) :: temp
            temp = rhs(2:)
            if (rhs(1:1) == 'b') then
                call bin2int(lhs, temp)
            else if (rhs(1:1) == 'o') then
                call oct2int(lhs, temp)
            end if
        end subroutine int8_from_string

end module bitwise
