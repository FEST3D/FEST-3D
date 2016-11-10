#!/bin/bash
total_process="$1"
echo -n > grid_names
for((i =0; i<total_process;i++ ))
do
echo "grid_$i.txt" >> grid_names
done
