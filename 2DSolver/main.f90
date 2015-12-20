program main

    use solver
!   use ppm, only: output_data

    implicit none

    call setup_solver()

    do while (.not. converged())
        print *, 'Iteration ', iter
        call step()

        !if (mod(iter, 100) == 0) then
        !    write(filename, '(A,I5.5,A)') 'facenn', iter-1, '.data'
        !    call output_data(filename)
        !end if
        if (iter == max_iters) then
            exit
        end if
    end do

    call destroy_solver()

end program main
