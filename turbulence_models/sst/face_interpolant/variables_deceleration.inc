! should be included in face_interpoant.f90 file 

    real, dimension(:, :, :), pointer :: x_tk_left, x_tk_right
    real, dimension(:, :, :), pointer :: y_tk_left, y_tk_right
    real, dimension(:, :, :), pointer :: z_tk_left, z_tk_right
    real, dimension(:, :, :), pointer :: x_tw_left, x_tw_right
    real, dimension(:, :, :), pointer :: y_tw_left, y_tw_right
    real, dimension(:, :, :), pointer :: z_tw_left, z_tw_right

    !public member
    public :: x_tk_left, x_tk_right
    public :: y_tk_left, y_tk_right
    public :: z_tk_left, z_tk_right
    public :: x_tw_left, x_tw_right
    public :: y_tw_left, y_tw_right
    public :: z_tw_left, z_tw_right
