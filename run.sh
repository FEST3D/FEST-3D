#!/bin/bash

ex='iitm_cfd_solver.out'
log='out'
output_dir='results'

echo "Run: `date +%Y/%m/%d-%H/%M/%S`" | tee -a $log

if ls output*.*vtk 1> /dev/null 2>&1; then
    # Output files exist
    # Back them up
    backup_dir="$output_dir`date +%Y%m%d%H%M%S`"
    mkdir $backup_dir
    mv output*.*vtk $backup_dir/
    echo "Previous output files moved to $backup_dir""." | tee -a $log
fi

if [ -f $ex ]; then
    echo 'Running solver.'
    echo >> $log
    echo 'Solver output:' >> $log
    echo >> $log
    espeak -a10 'Running solver.'
    ./$ex >> $log
fi

echo | tee -a $log
echo 'End of run script.' | tee -a $log
