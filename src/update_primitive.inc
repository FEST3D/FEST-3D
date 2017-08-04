 ! update primitive variable
  do k = 1,kmx-1
    do j = 1,jmx-1
      do i = 1,imx-1

        !check if user want to update from particular solution
        if(present(un))then
          u1(1:n_var) = un(i,j,k,1:n_var)
        else
          u1(1:n_var) = qp(i,j,k,1:n_var)
        end if

        ! finding primitive residue
        R(1) = mass_reidue(i,j,k)
        R(2) = -1*(u1(2)/u1(1))*mass_residue(i,j,k) + x_mom_residue(i,j,k)/u1(1)
        R(3) = -1*(u1(3)/u1(1))*mass_residue(i,j,k) + y_mom_residue(i,j,k)/u1(1)
        R(4) = -1*(u1(4)/u1(1))*mass_residue(i,j,k) + z_mom_residue(i,j,k)/u1(1)
        R(5) = 0.5*(gm-1.)*(sum(u1(2:4)**2)*mass_residue(i,j,k)) &
               -(gm-1.)*u1(2)*x_mom_residue(i,j,k)               &
               -(gm-1.)*u1(3)*y_mom_residue(i,j,k)               &
               -(gm-1.)*u1(4)*z_mom_residue(i,j,k)               &
               +(gm-1.)*energy_residue(i,j,k)

        select case(turbulence)
          case('none')
            !do nothing
            continue
          case('sst')
            beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
            R(5) = R(5) - (gm-1.)*TKE_residue(i,j,k)
            R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                   +(1./(1.+bstar*u1(6)*delta_t(i,j,k)))*TKE_residue(i,j,k)/u1(1)
            R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                   +(1./(1.+2.*beta*u1(6)*delta_t(i,j,k)))*omega_residue(i,j,k)/u1(1)
          case('kkl')
            eta  = u1(1)*dist(i,j,k)*(sqrt(0.3*u1(6))/(20*mu(i,j,k)))
            fphi = (1+cd1*eta)/(1+eta**4)
            R(5) = R(5) - (gm-1.)*TKE_residue(i,j,k)
            R(6) = -(u1(6)/u1(1))*mass_residue(i,j,k)&
                   + 1./(1.+(2.5*((cmu**0.75)*u1(1)*(u1(6)**1.5)/u1(7))&
                   -(2*mu(i,j,k)/dist(i,j,k))*delta_t(i,j,k)))*TKE_residue(i,j,k)/u1(1)
            R(7) = -(u1(7)/u1(1))*mass_residue(i,j,k)&
                   + 1./(1.-(6*mu(i,j,k)*fphi)*delta_t(i,j,k))*kl_residue(i,j,k)/u1(1)
          case DEFAULT
            Fatal_error
        end select

              
       !check if user want to store residue
        if(present(Rn)) then
          Rn(i,j,k,1:n_var) = Rn(i,j,k,1:n_var) + SF*R(1:n_var)
          if(TU) R(:) = Rn(i,j,k,:)
        end if
       

       !update
       u2(:) = u1(:) - R(:)*TF*delta_t(i,j,k)/volume(i,j,k)

       !check solution for non pyhysical results
       if(any(density < 0.) .or. any(pressure < 0.)then
         Fatal_error
       else !update
         qp(i,j,k,1:5) = u2(1:5)
         if(u2(6)>0.) qp(i,j,k,6) = u2(6)
         if(u2(7)>0.) qp(i,j,k,7) = u2(7)
       end if
      end do
    end do
  end do

