module utils

    implicit none
    private
    integer, public :: DEBUG_LEVEL = 1

    public :: alloc
    interface alloc
        module procedure alloc_rank1_real, &
                         alloc_rank2_real, &
                         alloc_rank3_real, &
                         alloc_rank4_real, &
                         alloc_rank5_real, &
                         alloc_rank2_integer
    end interface alloc

    public :: dealloc
    interface dealloc
        module procedure dealloc_rank1_real, &
                         dealloc_rank2_real, &
                         dealloc_rank3_real, &
                         dealloc_rank4_real, &
                         dealloc_rank5_real, &
                         dealloc_rank2_integer
    end interface dealloc

    public :: dmsg

    contains

        include "allocate_memory_implementation.inc"
        include "deallocate_memory_implementation.inc"

        subroutine dmsg(level, prog, method, msg)
            !---------------------------------------------------------------
            ! Print a DEBUG message
            !
            ! Input arguments:
            !   level -> integer
            !       the message's debug level
            !   prog -> character
            !       module / program name
            !   method -> character
            !       subroutine / function name
            !   msg -> character
            !       message
            !---------------------------------------------------------------

            implicit none
            character(len=*), optional :: prog
            character(len=*), optional :: method
            character(len=*), optional :: msg
            integer :: level
            
            if (level < DEBUG_LEVEL) then
                ! Don't print the message
                return
            end if

            if (.not. present(prog) .and. .not. present(method) .and. &
                    .not. present(msg)) then
                print *, 'Please provide atleast one of the following:'
                print *, '- Module / program name'
                print *, '- Subroutine / function name'
                print *, '- A custom message'
                stop
            end if

            print '(A7,I1.1,A2,A,A2,A,A1)', 'Debug: ', level, ' (', &
                    prog, ', ', method, ')'

            if (present(msg)) then
                print *, msg
            end if

        end subroutine dmsg

end module utils
