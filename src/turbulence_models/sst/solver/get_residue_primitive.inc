!include in solver.f90
!            
!            dEdx(:, :, :, 6) = ( (-1 * tk(1:imx-1, 1:jmx-1, 1:kmx-1) / &
!                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
!                                     mass_residue) + &
!                             ( TKE_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
!            dEdx(:, :, :, 7) = ( (-1 * tw(1:imx-1, 1:jmx-1, 1:kmx-1) / &
!                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
!                                     mass_residue) + &
!                             ( omega_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
      beta = beta1*sst_F1(1:imx-1,1:jmx-1,1:kmx-1) + (1. - sst_F1(1:imx-1,1:jmx-1,1:kmx-1))*beta2

              dEdx(:, :, :, 5) = dEdx(:,:,:,5) + (-(gm-1)*TKE_residue)
              dEdx(:,:,:,6)=(-1.*tk(1:imx-1,1:jmx-1,1:kmx-1)*mass_residue/density(1:imx-1,1:jmx-1,1:kmx-1)) &
                       + (1./(1.+bstar*tw(1:imx-1,1:jmx-1,1:kmx-1)*delta_t(1:imx-1,1:jmx-1,1:kmx-1))&
                       * TKE_residue/density(1:imx-1,1:jmx-1,1:kmx-1))

              dEdx(:,:,:,7)=(-1.*tw(1:imx-1,1:jmx-1,1:kmx-1)*mass_residue/density(1:imx-1,1:jmx-1,1:kmx-1))&
                       +( 1./(1. + 2.*beta*tw(1:imx-1,1:jmx-1,1:kmx-1)*delta_t(1:imx-1,1:jmx-1,1:kmx-1))&
                       *omega_residue/density(1:imx-1,1:jmx-1,1:kmx-1))
