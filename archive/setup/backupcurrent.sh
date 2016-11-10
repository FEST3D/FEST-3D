#!/bin/bash
# ---------------------------------------------------------------------
# backupcurrent.sh does the following:
# 1. Moves all vtk files in current directory from previous simulation
#    to a new directory named uniquely with date-and-time-stamp. Also,
#    get the latest fvtk file for restarting simulation
# 2. Stop visdaemon if found to be running
# 3. Copy all logs to the backup directory 
# 4. Copy all visit post processing scripts
# 5. Get images of residue norms as a function of number of iterations
#    and save to backup directory
# ---------------------------------------------------------------------

output_dir='results'
backup_dir="$output_dir`date +%Y%m%d%H%M%S`"
ex='iitm_cfd_solver.out'
runlog='out'
vislog='vislog'
resnorms='resnorms'
art_disspn='art_disspn'
visdaemon='visdaemon.py'
config='config.md'
plottinglist='visit_output_variables.md'

mkdir $backup_dir
mv output*.vtk $backup_dir/

if [ -e .${visdaemon:0: -3}.pid ]; then
    python $visdaemon stop
    echo 'Daemon found to be running. Now stopped' | tee -a $runlog
fi

cp $config $backup_dir/
cp $resnorms $backup_dir/
cp $runlog $backup_dir/
cp $plottinglist $backup_dir/
cp $art_disspn $backup_dir/
cp $vislog $backup_dir/

cp visitpostprocessing.py $backup_dir/
cp visit_display_last_iteration.py $backup_dir/
cp visit_open_database_gui.py $backup_dir/

