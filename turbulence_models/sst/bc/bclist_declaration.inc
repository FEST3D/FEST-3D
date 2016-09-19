! should be include in include file  "bclist_declaration" in main directory
! declaration do not use slect
    
    ! Fix Type Conditions
        integer, public :: sst_BC_FIX_tk
        integer, public :: sst_BC_FIX_tw

    ! Copy Type Conditions
        integer, public :: sst_BC_COPY_tk
        integer, public :: sst_BC_COPY_tw

        real :: mu, rho, T_l, T_r, T_face
        integer :: i,j,k
