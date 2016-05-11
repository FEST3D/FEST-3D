subroutine dealloc_rank1_real(var)
    implicit none
    real, dimension(:), intent(inout), allocatable :: var
    integer :: mem_stat
    deallocate(var, stat=mem_stat)
    if (mem_stat /= 0) then
        print *, 'Error: Could not deallocate memory.'
    end if
end subroutine dealloc_rank1_real

subroutine dealloc_rank2_real(var)
    implicit none
    real, dimension(:, :), intent(inout), allocatable :: var
    integer :: mem_stat
    deallocate(var, stat=mem_stat)
    if (mem_stat /= 0) then
        print *, 'Error: Could not deallocate memory.'
    end if
end subroutine dealloc_rank2_real

subroutine dealloc_rank3_real(var)
    implicit none
    real, dimension(:, :, :), intent(inout), allocatable :: var
    integer :: mem_stat
    deallocate(var, stat=mem_stat)
    if (mem_stat /= 0) then
        print *, 'Error: Could not deallocate memory.'
    end if
end subroutine dealloc_rank3_real

subroutine dealloc_rank4_real(var)
    implicit none
    real, dimension(:, :, :, :), intent(inout), allocatable :: var
    integer :: mem_stat
    deallocate(var, stat=mem_stat)
    if (mem_stat /= 0) then
        print *, 'Error: Could not deallocate memory.'
    end if
end subroutine dealloc_rank4_real

subroutine dealloc_rank5_real(var)
    implicit none
    real, dimension(:, :, :, :, :), intent(inout), allocatable :: var
    integer :: mem_stat
    deallocate(var, stat=mem_stat)
    if (mem_stat /= 0) then
        print *, 'Error: Could not deallocate memory.'
    end if
end subroutine dealloc_rank5_real

subroutine dealloc_rank2_integer(var)
    implicit none
    integer, dimension(:, :), intent(inout), allocatable :: var
    integer :: mem_stat
    deallocate(var, stat=mem_stat)
    if (mem_stat /= 0) then
        print *, 'Error: Could not deallocate memory.'
    end if
end subroutine dealloc_rank2_integer
