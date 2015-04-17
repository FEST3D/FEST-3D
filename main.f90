program main

    use solver

    implicit none

    call setup_solver()

    do while (.not. converged())
        print *, 'Iteration ', iter
        call step()
        if (iter == max_iters) then
            exit
        end if
    end do

    call destroy_solver()

end program main
