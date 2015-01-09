program main

    use global, only: FILE_NAME_LENGTH
    use solver
    use state, only: writestate

    implicit none
    character(len=FILE_NAME_LENGTH) :: filename

    call setup_solver()

    do while (.not. converged())
        print *, 'Iteration ', iter
        if (mod(iter, 100) == 0) then
            print *, 'Iteration: ', iter, ' **********'
            write(filename, '(A,I5.5,A)') 'output', iter, '.fvtk'
            call writestate(filename)
        end if
        call step()
        if (iter == max_iters) then
            exit
        end if
    end do

    call destroy_solver()

end program main
