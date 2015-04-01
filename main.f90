program main

    use global, only: FILE_NAME_LENGTH
    use solver
    use state, only: writestate
    use ppm, only: output_data

    implicit none
    character(len=FILE_NAME_LENGTH) :: filename

    call setup_solver()

    do while (.not. converged())
        print *, 'Iteration ', iter
        if (mod(iter+1, 100) == 0) then
            print *, 'Iteration: ', iter, ' **********'
            write(filename, '(A,I5.5,A)') 'output', iter, '.fvtk'
            call writestate(filename)
        end if

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
