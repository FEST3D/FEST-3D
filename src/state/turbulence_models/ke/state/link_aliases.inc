!should be included in state state, in subroutine with same name.

            call dmsg(1, 'turbulence_models/ke/state', 'link_aliases.inc')
            tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
            te(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
            tk_inf => qp_inf(6)
            te_inf => qp_inf(7)
