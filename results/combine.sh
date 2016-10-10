#!/bin/bash
# this script copies the output files from different process_x folders 
# into a new folder and creates a .visit master file to combine the blocks 
output_dir='net'
total_process=2 #"$1"
if [ -d "$output_dir" ]; then
rm -r "$output_dir"
fi
mkdir "$output_dir"
echo "!NBLOCKS $total_process" > net/net.visit
for (( i = 0; i< $total_process; i++ ))
do
folder="process_0$i"
cd "$folder"
tot_files=`ls output* | wc -l` 
for((j = 0; j<tot_files; j++))
do
if [ "$j" -lt 10 ]; then
k="0000" 
elif [ "$j" -lt 100 ]; then
k="000"
elif [ "$j" -lt 1000 ]; then
k="00"
fi
file="output$k$j.vtk"
cp "$file" ../"$output_dir"/"proc-$i-$k$j.vtk"
# write into .visit
if [ "$i" -eq 0 ]; then
for (( m = 0; m< $total_process; m++ ))
do
echo "proc-$m-$k$j.vtk" >> ../net/net.visit 
done
fi 
done
cd ..
done
