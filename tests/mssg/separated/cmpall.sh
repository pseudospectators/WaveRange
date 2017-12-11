#!/bin/bash

# Number of files
num=$(ls -l res.p_???? | grep -v ^d | grep -v ^t | wc -l)

# Loop for all files
for (( i=0; i<$num; i++ ))
do
./wrmssgenc res .enc 2 2 1 1e-6 $i
./wrmssgdec res .enc dec 2 2 1 $i
done

./compare.out

# Total encoded file size
find . -type f -name '*.enc' -exec du -ch {} + | grep total$
