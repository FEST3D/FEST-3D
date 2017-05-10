!should be included in state state, in subroutine with same name.

            call dmsg(1, 'turbulence_models/sst/state', 'link_aliases.inc')
            tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 6)
            tw(-2:imx+2, -2:jmx+2, -2:kmx+2) => qp(:, :, :, 7)
            tk_inf => qp_inf(6)
            tw_inf => qp_inf(7)
