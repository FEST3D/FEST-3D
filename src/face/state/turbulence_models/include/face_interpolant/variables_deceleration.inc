! should be included in face_interpoant.f90 file 

!Variables for the sst module, should be included in main state.f90 file

! should be included in face_interpoant.f90 file 

    real, dimension(:, :, :), pointer :: x_tk_left, x_tk_right
    real, dimension(:, :, :), pointer :: y_tk_left, y_tk_right
    real, dimension(:, :, :), pointer :: z_tk_left, z_tk_right
    real, dimension(:, :, :), pointer :: x_tw_left, x_tw_right
    real, dimension(:, :, :), pointer :: y_tw_left, y_tw_right
    real, dimension(:, :, :), pointer :: z_tw_left, z_tw_right
    real, dimension(:, :, :), pointer :: x_tkl_left, x_tkl_right
    real, dimension(:, :, :), pointer :: y_tkl_left, y_tkl_right
    real, dimension(:, :, :), pointer :: z_tkl_left, z_tkl_right
    real, dimension(:, :, :), pointer :: x_tv_left, x_tv_right
    real, dimension(:, :, :), pointer :: y_tv_left, y_tv_right
    real, dimension(:, :, :), pointer :: z_tv_left, z_tv_right
    real, dimension(:, :, :), pointer :: x_te_left, x_te_right
    real, dimension(:, :, :), pointer :: y_te_left, y_te_right
    real, dimension(:, :, :), pointer :: z_te_left, z_te_right
    real, dimension(:, :, :), pointer :: x_tgm_left, x_tgm_right
    real, dimension(:, :, :), pointer :: y_tgm_left, y_tgm_right
    real, dimension(:, :, :), pointer :: z_tgm_left, z_tgm_right

    !public member
    public :: x_tk_left, x_tk_right
    public :: y_tk_left, y_tk_right
    public :: z_tk_left, z_tk_right
    public :: x_tw_left, x_tw_right
    public :: y_tw_left, y_tw_right
    public :: z_tw_left, z_tw_right
    public :: x_tkl_left, x_tkl_right
    public :: y_tkl_left, y_tkl_right
    public :: z_tkl_left, z_tkl_right
    public :: x_tv_left, x_tv_right
    public :: y_tv_left, y_tv_right
    public :: z_tv_left, z_tv_right
    public :: x_te_left, x_te_right
    public :: y_te_left, y_te_right
    public :: z_te_left, z_te_right
    public :: x_tgm_left, x_tgm_right
    public :: y_tgm_left, y_tgm_right
    public :: z_tgm_left, z_tgm_right
