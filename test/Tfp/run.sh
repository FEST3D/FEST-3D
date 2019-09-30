#!/bin/bash

ex='../../bin/FEST3D'
runlog='time_directories/aux/out'
total_process=2 #"$2" #parallel running - number of process

#echo "Run: `date +%Y/%m/%d-%H/%M/%S`" | tee -a $runlog

if [ -f $ex ]; then
    echo "Ran Test number 3  --->  Turbulent flow over a flat plate"
    echo "  __________Report__________   "
    mpiexec.hydra -np $total_process ./$ex >> $runlog
    cd pp
    python main.py
    cd ..
    rm -r time_directories/0000
    rm -r time_directories/0001
    rm -r time_directories/aux/out
    rm -r time_directories/aux/resnorm
    touch time_directories/aux/resnorm
fi
