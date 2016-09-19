!Variables for the sst module, should be included in main state.f90 file
!
    ! State variables
    !real, public, dimension(:, :, :, :), allocatable, target :: tqp !!!!!!not needed
    ! Infinity variables (free stream conditions)
    !real, public, dimension(:), allocatable, target :: tqp_inf !! nor needed added to n_var
    ! State variable component aliases
    integer, public :: sst_n_var = 2 !make it generalisse
    real, public, dimension(:, :, :), pointer :: tk
    real, public, dimension(:, :, :), pointer :: tw
    real, public, pointer :: tk_inf
    real, public, pointer :: tw_inf

