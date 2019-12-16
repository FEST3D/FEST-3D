! explicit residual smoothing
              residue(2:imx-2,2:jmx-2,2:kmx-2,:)=eps*(residue(1:imx-3,2:jmx-2,2:kmx-2,:)&
                                                    + residue(2:imx-2,1:jmx-3,2:kmx-2,:)&
                                                    + residue(2:imx-2,2:jmx-2,1:kmx-3,:)&
                                                    + residue(3:imx-1,2:jmx-2,2:kmx-2,:)&
                                                    + residue(2:imx-2,3:jmx-1,2:kmx-2,:)&
                                                    + residue(2:imx-2,2:jmx-2,3:kmx-1,:))&
                                           +(1-6*eps)*residue(2:imx-2,2:jmx-2,2:kmx-2,:) 
