#!/bin/bash

ex='iitm_cfd_solver.out'
log='out'

alias compile_fortran="mpif90 -O3 -Wall -Wextra -Wconversion "\
"-Wno-compare-reals "\
"-fdefault-real-8 "\
"-Waliasing "\
"-Wsurprising "\
"-Wintrinsic-shadow "\
"-Werror "\
"-pedantic-errors "\
"-fbounds-check"
# Flags description: https://gcc.gnu.org/onlinedocs/gcc-4.5.4/gfortran/Error-and-Warning-Options.html

# Clear the log file before beginning the compile.
# Put in the timestamp as the header.
echo "Compile: `date +%Y/%m/%d-%H/%M/%S`" | tee $log

if [ -f $ex ]; then
    echo 'Deleting previous executable.' | tee -a $log
    rm -f $ex | tee -a $log
fi

# Compiling strategy:
# Remove the previously compiled object and module files.
# Compile the files one at a time.
# If the corresponding module file is missing, there is some issue.
# Stop compiling.

# Insert ppm after state. De comment lines in face_interpolant
# Insert ausm, ldfss0 and hlle adter van_leer. Decomment lines in
# scheme

filelist="
    global
    utils
    layout
    bitwise
    string
    grid
    geometry
    state
    ppm
    muscl
    face_interpolant
    van_leer
    viscous
    scheme
    parallel
    boundary_conditions
    solver
"

function process_modules {
    #espeak 'Processing module files.'
    for codefile in $filelist; do
        echo "Processing $codefile." | tee -a $log
        if [ -f "$codefile"".mod" ]; then
            rm "$codefile"".mod" | tee -a $log
        fi
        if [ -f "$codefile"".o" ]; then
            rm "$codefile"".o" | tee -a $log
        fi
        compile_fortran -c "$codefile"".f90" | tee -a $log
        if [ ! -f "$codefile"".mod" ]; then
            return -1
        fi
    done
    return 0
}

function create_executable {
    echo | tee -a $log
    echo 'Creating compiled file.' | tee -a $log

    # Create object file list
    obj_file_list=""
    for codefile in $filelist; do
        obj_file_list="$obj_file_list""$codefile"".o "
    done

    compile_fortran $obj_file_list main.f90 -o $ex | tee -a $log

    if [ -f $ex ]; then
        return 0
    else
        return -1
    fi
}

function output_todo {
    echo 'Check correctness of write_interface.py' | tee -a $log
    echo | tee -a $log
    echo 'TODO:' | tee -a $log
    grep -rn --exclude-dir={docs,results*} --exclude={run.sh,compile.sh,out} \
            'TODO' * | tee -a $log
}

function remove_compilation_files {
    echo | tee -a $log
    echo 'Removing generated files.' | tee -a $log
    #espeak 'Removing generated files.'
    for codefile in $filelist; do
        if [ -f "$codefile"".mod" ]; then
            rm "$codefile"".mod" | tee -a $log
        fi
        if [ -f "$codefile"".o" ]; then
            rm "$codefile"".o"  | tee -a $log
        fi
    done
}

if process_modules; then
    echo | tee -a $log
    echo 'All modules processed.' | tee -a $log

    if create_executable; then
        output_todo
    fi
fi

remove_compilation_files

unalias compile_fortran

echo | tee -a $log
echo 'End of compile script.' | tee -a $log

