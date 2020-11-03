#/bin/bash

for dir in COR/*/; do
# get pair name from subdir
pair=`echo $dir | rev | cut -d'/' -f-2 | cut -d'/' -f2- | rev`
echo $pair.in

# list files in that subdur (SACpos and SACneg), write to a file for this pair
ls $dir*.SACpos > $pair.in
ls $dir*.SACneg >> $pair.in

done
