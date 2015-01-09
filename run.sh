#!/bin/bash

ex='iitm_cfd_solver.out'
log='out'

alias compile_fortran="gfortran -Wall -Wextra -Wconversion "\
"-Wno-compare-reals "\
"-Waliasing "\
"-Wsurprising "\
"-Wintrinsic-shadow "\
"-Werror "\
"-pedantic-errors "\
"-fbounds-check"
# Flags description: https://gcc.gnu.org/onlinedocs/gcc-4.5.4/gfortran/Error-and-Warning-Options.html

echo "Run: `date +%Y/%m/%d-%H/%M/%S`" | tee $log

if [ -f $ex ]; then
    echo 'Deleting previous executable.' | tee -a $log
    rm -f $ex | tee -a $log
fi

# Compiling strategy:
# Remove the previously compiled object and module files.
# Compile the files one at a time.
# If the corresponding module file is missing, there is some issue.
# Stop compiling.

filelist="
    global
    utils
    string
    grid
    geometry
    state
    van_leer
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
    grep --exclude=run.sh 'TODO' * | tee -a $log
    #espeak 'Please finish the remaining tasks!'
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
    #espeak 'All module files processed.'

    if create_executable; then
        #espeak 'Compiled solver file prepped.'
        output_todo
    fi
fi

remove_compilation_files

if [ -f $ex ]; then
    echo 'Running solver.'
    echo >> $log
    echo 'Solver output:' >> $log
    echo >> $log
    espeak -a10 'Running solver.'
    ./$ex >> $log
fi

unalias compile_fortran
