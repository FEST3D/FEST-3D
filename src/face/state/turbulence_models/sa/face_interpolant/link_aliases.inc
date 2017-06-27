!should be included in face_interpolant.f90 file
        

            ! Link xi faces left pointers
            x_tv_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)

            ! Link xi faces right pointers
            x_tv_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)

            ! Link eta faces left pointers
            y_tv_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)

            ! Link eta faces right pointers
            y_tv_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            
            ! Link zeta faces left pointers
            z_tv_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)

            ! Link zeta faces right pointers
            z_tv_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)

