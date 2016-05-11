#!/bin/bash
cd results
for (( i = 0; i< $total_process; i++ ))
do
folder="process_0$i" 
if [ -d "$folder" ]; then
rm "$folder"/output*
else
mkdir "$folder"
fi
done

cd ..
