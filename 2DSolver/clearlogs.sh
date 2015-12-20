#!/bin/bash
#----------------------------------------------------------------------
# This script does the following
# 1. Clears all the logs
# 2. Removes all the output files: .fvtk, .vtk and the output log
# 3. Stops the visdaemon if found to be running
#----------------------------------------------------------------------

output_dir='results'
runlog='out'
vislog='vislog'
resnorms='resnorms'
art_disspn='art_disspn'
visdaemon='visdaemon.py'

echo 'Clearing all logs' | tee -a $runlog

echo > $runlog
echo > $vislog
echo > $resnorms
echo > $art_disspn

if [ -f output00000.fvtk ]; then
    rm output*
fi

if [ -e .${visdaemon:0: -3}.pid ]; then
    python $visdaemon stop
    echo 'Deamon found to be running. Now stopped'
    rm visitlog.py
fi    
