module utils

    use global_vars, only : process_id
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
    public :: turbulence_read_error

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
            character(len=256)         :: ifmsg   
            integer :: level
           
!            if (process_id == 0) then
            if (level < DEBUG_LEVEL) then
                ! Don't print the message
                return
            end if


            ifmsg = ""
            if (present(msg)) then
              ifmsg = " >--> "//trim(msg)
            end if

            if (.not. present(prog) .and. .not. present(method) .and. &
                    .not. present(msg)) then
                print *, 'Please provide atleast one of the following:'
                print *, '- Module / program name'
                print *, '- Subroutine / function name'
                print *, '- A custom message'
                stop
            end if

            print '(A7,I1.1,A,I2,A2,A,A2,A,A,A1)', 'Debug: ', level," id - ", process_id, ' (', &
              trim(prog), ', ', trim(method), trim(ifmsg), ')'
!           end if

        end subroutine dmsg
        
        subroutine turbulence_read_error()

          print*, "ERROR: Turbulence model not recognised"
          STOP

        end subroutine turbulence_read_error

end module utils
