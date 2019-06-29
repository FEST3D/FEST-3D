subroutine dealloc_rank1_real(var)
  !< Deallcoate 1-Dimensional array of type: real
    implicit none
    real, dimension(:), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank1_real

subroutine dealloc_rank2_real(var)
  !< Deallcoate 2-Dimensional array of type: real
    implicit none
    real, dimension(:, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank2_real

subroutine dealloc_rank3_real(var)
  !< Deallcoate 3-Dimensional array of type: real
    implicit none
    real, dimension(:, :, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank3_real

subroutine dealloc_rank4_real(var)
  !< Deallcoate 4-Dimensional array of type: real
    implicit none
    real, dimension(:, :, :, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank4_real

subroutine dealloc_rank5_real(var)
  !< Deallcoate 5-Dimensional array of type: real
    implicit none
    real, dimension(:, :, :, :, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank5_real

subroutine dealloc_rank6_real(var)
  !< Deallcoate 6-Dimensional array of type: real
    implicit none
    real, dimension(:, :, :, :, :, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank6_real

subroutine dealloc_rank1_integer(var)
  !< Deallcoate 1-Dimensional array of type: integer
    implicit none
    integer, dimension(:), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank1_integer

subroutine dealloc_rank2_integer(var)
  !< Deallcoate 2-Dimensional array of type: integer
    implicit none
    integer, dimension(:, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank2_integer

subroutine dealloc_rank3_integer(var)
  !< Deallcoate 3-Dimensional array of type: integer
    implicit none
    integer, dimension(:, :, :), intent(inout), allocatable :: var
    !< Variable to which memeory is deallocated
    integer :: mem_stat
    !< Status of the memory deallocation process
    if(allocated(var)) then
      deallocate(var, stat=mem_stat)
      if (mem_stat /= 0) then
          print *, 'Error: Could not deallocate memory.'
      end if
    end if
end subroutine dealloc_rank3_integer

