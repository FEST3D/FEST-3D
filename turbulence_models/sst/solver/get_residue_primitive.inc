!include in solver.f90
            
            dEdx(:, :, :, 6) = ( (-1 * tk(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( TKE_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 7) = ( (-1 * tw(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( omega_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
