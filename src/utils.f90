  !< Utility module to allocate, deallocate and debug message
module utils
  !< Utility module to allocate, deallocate and debug message

    implicit none
    private
    public :: alloc
    interface alloc
        module procedure alloc_rank1_real, &
                         alloc_rank2_real, &
                         alloc_rank3_real, &
                         alloc_rank4_real, &
                         alloc_rank5_real, &
                         alloc_rank6_real, &
                         alloc_rank1_integer,&
                         alloc_rank2_integer,&
                         alloc_rank3_integer
    end interface alloc


    contains

        subroutine alloc_rank1_real(var, start1, stop1, errmsg)
          !< Allcoate 1-Dimensional array of type: real
            implicit none
            real, dimension(:), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
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
          !< Allcoate 2-Dimensional array of type: real
            implicit none
            real, dimension(:, :), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
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
          !< Allcoate 3-Dimensional array of type: real
            implicit none
            real, dimension(:, :, :), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2, start3
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2, stop3
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
            allocate(var(start1:stop1, start2:stop2, start3:stop3), &
                    stat=mem_stat)
            if (mem_stat /= 0) then
                if (present(errmsg)) then
                    print *, errmsg
                else
                    print *, 'Error: Could not allocate memory.'
                end if
                print *, 'Required extent: ', stop1 - start1 + 1, &
                        stop2 - start2 + 1, stop3 - start3 + 1
                stop
            end if
        end subroutine alloc_rank3_real

        subroutine alloc_rank4_real(var, start1, stop1, start2, stop2, &
                start3, stop3, start4, stop4, errmsg)
          !< Allcoate 4-Dimensional array of type: real
            implicit none
            real, dimension(:, :, :, :), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2, start3, start4
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2, stop3, stop4
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
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
          !< Allcoate 5-Dimensional array of type: real
            implicit none
            real, dimension(:, :, :, :, :), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2, start3, start4, start5
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2, stop3, stop4, stop5
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
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

        subroutine alloc_rank6_real(var, start1, stop1, start2, stop2, &
          !< Allcoate 6-Dimensional array of type: real
                start3, stop3, start4, stop4, start5, stop5, start6, stop6, errmsg)
            implicit none
            real, dimension(:, :, :, :, :,:), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2, start3, start4, start5, start6
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2, stop3, stop4, stop5, stop6
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
            allocate(var(start1:stop1, start2:stop2, start3:stop3, &
                    start4:stop4, start5:stop5, start6:stop6), stat=mem_stat)
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
        end subroutine alloc_rank6_real

        subroutine alloc_rank1_integer(var, start1, stop1, errmsg)
          !< Allcoate 1-Dimensional array of type: integer
            implicit none
            integer, dimension(:), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
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
        end subroutine alloc_rank1_integer

        subroutine alloc_rank2_integer(var, start1, stop1, start2, stop2, errmsg)
          !< Allcoate 2-Dimensional array of type: integer
            implicit none
            integer, dimension(:, :), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
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


        subroutine alloc_rank3_integer(var, start1, stop1, start2, stop2, &
                start3, stop3, errmsg)
          !< Allcoate 3-Dimensional array of type: integer
            implicit none
            integer, dimension(:, :, :), intent(inout), allocatable :: var
            !< Variable to which memory is allocated
            integer, intent(in) :: start1, start2, start3
            !< Starting index of Var array's dimension
            integer, intent(in) :: stop1, stop2, stop3
            !< Last index of Var array's dimension
            integer :: mem_stat
            !< Status of the memory allocation process
            character(len=*), intent(in), optional :: errmsg
            !< Error message to print if mem_stat is not 0(successful)
            allocate(var(start1:stop1, start2:stop2, start3:stop3), &
                    stat=mem_stat)
            if (mem_stat /= 0) then
                if (present(errmsg)) then
                    print *, errmsg
                else
                    print *, 'Error: Could not allocate memory.'
                end if
                print *, 'Required extent: ', stop1 - start1 + 1, &
                        stop2 - start2 + 1, stop3 - start3 + 1
                stop
            end if
        end subroutine alloc_rank3_integer

end module utils
