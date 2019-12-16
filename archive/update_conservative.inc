 ! update conservative variable
  do k = 1,kmx-1
    do j = 1,jmx-1
      do i = 1,imx-1

        ! getting conservative variable
        !check if user want to update from particular solution
        if(present(un))then
          u1(1)  = un(i,j,k,1)
          u1(2:) = un(i,j,k,2:)*u1(1)
          select case(turbulence)
            case('sst', 'kkl')
              KE = u1(6)
            case DEFAULT
              KE = 0.
          end select
          u1(5) = un(i,j,k,5)/(gm-1.) + (0.5*sum(u1(2:4)**2)/u1(1)) + KE
        else
          u1(1)  = qp(i,j,k,1)
          u1(2:) = qp(i,j,k,2:)*u1(1)
          select case(turbulence)
            case('sst', 'kkl')
              KE = u1(6)
            case DEFAULT
              KE = 0.
          end select
          u1(5) = qp(i,j,k,5)/(gm-1.) + (0.5*sum(u1(2:4)**2)/u1(1)) + KE
        end if

       ! get R
        R(:) = residue(i,j,k,:) 
       !check if user want to store residue
       if(present(Rn)) then
         Rn(i,j,k,1:n_var) = Rn(i,j,k,1:n_var) + SF*R(1:n_var)
         if(TU) R(:) = Rn(i,j,k,:)
       end if

       !update
       u2(:) = u1(:) - R(:)*TF*delta_t(i,j,k)/volume(i,j,k)

      ! getting primitve variable back variable
        u2(1)  = u2(1)
        u2(2:) = u2(2:)/u2(1)
        select case(turbulence)
          case('sst', 'kkl')
            KE = u2(6)
          case DEFAULT
            KE = 0.
        end select
        u2(5) = (gm-1.)*u2(5) - (0.5*u2(1)*sum(u2(2:4)**2)) - u2(1)*KE

       !check solution for non pyhysical results
       if((u2(1) < 0.) .or. (u2(5)) < 0.)then
         Fatal_error
       else !update
         qp(i,j,k,1:5) = u2(1:5)
         if(u2(6)>0.) qp(i,j,k,6) = u2(6)
         if(u2(7)>0.) qp(i,j,k,7) = u2(7)
       end if
      end do
    end do
  end do

