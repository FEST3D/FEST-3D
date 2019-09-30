#!/bin/bash
total_process="$1"
echo -n > grid_names
for((i =0; i<total_process;i++ ))
do
printf -v j "%02d" $i
echo "../gridfiles/grid_$j.txt" >> grid_names
done
