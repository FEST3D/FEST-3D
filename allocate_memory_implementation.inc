subroutine alloc_rank1_real(var, start1, stop1, errmsg)
    implicit none
    real, dimension(:), intent(inout), allocatable :: var
    integer, intent(in) :: start1
    integer, intent(in) :: stop1
    integer :: mem_stat
    character(len=*), intent(in), optional :: errmsg
    allocate(var(start1:stop1), stat=mem_stat)
    if (mem_stat /= 0) then
        if (present(errmsg)) then
            print *, errmsg
        else
            print *, 'Error: Could not allocate memory.'
        end if
        print *, 'Required extent: ', stop1 - start1 + 1
        stop
    end if
end subroutine alloc_rank1_real

subroutine alloc_rank2_real(var, start1, stop1, start2, stop2, errmsg)
    implicit none
    real, dimension(:, :), intent(inout), allocatable :: var
    integer, intent(in) :: start1, start2
    integer, intent(in) :: stop1, stop2
    integer :: mem_stat
    character(len=*), intent(in), optional :: errmsg
    allocate(var(start1:stop1, start2:stop2), stat=mem_stat)
    if (mem_stat /= 0) then
        if (present(errmsg)) then
            print *, errmsg
        else
            print *, 'Error: Could not allocate memory.'
        end if
        print *, 'Required extent: ', stop1 - start1 + 1, &
                stop2 - start2 + 1
        stop
    end if
end subroutine alloc_rank2_real

subroutine alloc_rank3_real(var, start1, stop1, start2, stop2, &
        start3, stop3, errmsg)
    implicit none
    real, dimension(:, :, :), intent(inout), allocatable :: var
    integer, intent(in) :: start1, start2, start3
    integer, intent(in) :: stop1, stop2, stop3
    integer :: mem_stat
    character(len=*), intent(in), optional :: errmsg
    allocate(var(start1:stop1, start2:stop2, start3:stop3), &
            stat=mem_stat)
    if (mem_stat /= 0) then
        if (present(errmsg)) then
            print *, errmsg
        else
            print *, 'Error: Could not allocate memory.'
        end if
        print *, 'Required extent: ', stop1 - start1 + 1, &
                stop2 - start2 + 1, stop3 - stop3 + 1
        stop
    end if
end subroutine alloc_rank3_real

subroutine alloc_rank4_real(var, start1, stop1, start2, stop2, &
        start3, stop3, start4, stop4, errmsg)
    implicit none
    real, dimension(:, :, :, :), intent(inout), allocatable :: var
    integer, intent(in) :: start1, start2, start3, start4
    integer, intent(in) :: stop1, stop2, stop3, stop4
    integer :: mem_stat
    character(len=*), intent(in), optional :: errmsg
    allocate(var(start1:stop1, start2:stop2, start3:stop3, &
            start4:stop4), stat=mem_stat)
    if (mem_stat /= 0) then
        if (present(errmsg)) then
            print *, errmsg
        else
            print *, 'Error: Could not allocate memory.'
        end if
        print *, 'Required extent: ', stop1 - start1 + 1, &
                stop2 - start2 + 1, stop3 - stop3 + 1, &
                stop4 - start4 + 1
        stop
    end if
end subroutine alloc_rank4_real

subroutine alloc_rank5_real(var, start1, stop1, start2, stop2, &
        start3, stop3, start4, stop4, start5, stop5, errmsg)
    implicit none
    real, dimension(:, :, :, :, :), intent(inout), allocatable :: var
    integer, intent(in) :: start1, start2, start3, start4, start5
    integer, intent(in) :: stop1, stop2, stop3, stop4, stop5
    integer :: mem_stat
    character(len=*), intent(in), optional :: errmsg
    allocate(var(start1:stop1, start2:stop2, start3:stop3, &
            start4:stop4, start5:stop5), stat=mem_stat)
    if (mem_stat /= 0) then
        if (present(errmsg)) then
            print *, errmsg
        else
            print *, 'Error: Could not allocate memory.'
        end if
        print *, 'Required extent: ', stop1 - start1 + 1, &
                stop2 - start2 + 1, stop3 - stop3 + 1, &
                stop4 - start4 + 1, stop5 - start5 + 1
        stop
    end if
end subroutine alloc_rank5_real

subroutine alloc_rank2_integer(var, start1, stop1, start2, stop2, errmsg)
    implicit none
    integer, dimension(:, :), intent(inout), allocatable :: var
    integer, intent(in) :: start1, start2
    integer, intent(in) :: stop1, stop2
    integer :: mem_stat
    character(len=*), intent(in), optional :: errmsg
    allocate(var(start1:stop1, start2:stop2), stat=mem_stat)
    if (mem_stat /= 0) then
        if (present(errmsg)) then
            print *, errmsg
        else
            print *, 'Error: Could not allocate memory.'
        end if
        print *, 'Required extent: ', stop1 - start1 + 1, &
                stop2 - start2 + 1
        stop
    end if
end subroutine alloc_rank2_integer
