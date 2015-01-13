module global
    !-------------------------------------------------------------------
    ! The global module holds global variables
    !
    ! Global variables include string buffer lengths, file unit lengths,
    ! etc.
    !-------------------------------------------------------------------

    implicit none
    ! String buffer lengths
    integer, parameter :: FILE_NAME_LENGTH = 32
    integer, parameter :: SCHEME_NAME_LENGTH = 16
    integer, parameter :: STRING_BUFFER_LENGTH = 128
    integer, parameter :: ERROR_MESSAGE_LENGTH = 256
    ! File unit numbers
    integer, parameter :: CONFIG_FILE_UNIT = 1
    integer, parameter :: GRID_FILE_UNIT = 2
    integer, parameter :: STATE_FILE_UNIT = 10
    integer, parameter :: OUT_FILE_UNIT = 20

end module global
