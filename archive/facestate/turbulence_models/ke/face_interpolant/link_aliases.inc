!should be included in face_interpolant.f90 file
        

            ! Link xi faces left pointers
            x_tk_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 6)
            x_te_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 7)

            ! Link xi faces right pointers
            x_tk_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 6)
            x_te_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 7)

            ! Link eta faces left pointers
            y_tk_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 6)
            y_te_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 7)

            ! Link eta faces right pointers
            y_tk_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 6)
            y_te_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 7)
            
            ! Link zeta faces left pointers
            z_tk_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 6)
            z_te_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 7)

            ! Link zeta faces right pointers
            z_tk_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 6)
            z_te_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 7)

