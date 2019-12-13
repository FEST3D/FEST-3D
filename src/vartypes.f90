module vartypes
    implicit none

    type, public :: nodetype
        real :: x
        real :: y
        real :: z
    end type nodetype


    type, public :: extent
      integer :: imx
      integer :: jmx
      integer :: kmx
    end type extent


    type, public :: celltype
      real :: volume
      !< Store cell volume
      real :: centerx
      real :: centery
      real :: centerz
      !< Store Cell-center location 
    end type celltype


    type, public :: facetype
      real :: A
       !< Store magnitude of face area vector of direction faces
      real :: nx
      real :: ny
      real :: nz
       !< Store unit face normal vector for all faces 
    end type facetype


    type, public :: controltype
        real :: CFL=1.0
        !< Courant–Friedrichs–Lewy (CFL) (Read from input)
!        integer :: min_iter=1     
!        !< Minimum iteration value, starting iteration value
!        integer :: max_iters=1 
!        !< Maximum iteration value, stop after these many iteration
        integer :: start_from=0 
        !< Number of the folder (in time_directories) to load stored state from to restart computation
        integer                                           :: n_var=5
        ! Freestram variable used to read file before inf pointer are linked and allocated
    end type controltype
end module vartypes
