#!/bin/bash

# Set tolerance vector
tolvec=( '1e-1' '1e-2' '1e-3' '1e-4' '1e-6' '1e-8' '1e-10' '1e-12' '1e-14' )

# Number of elements
num=${#tolvec[@]}

# Loop for all tolerances
for (( j=0; j<$num; j++ ))
do

# Encode
./wrmssgenc res .enc 1 2 1 ${tolvec[$j]} 0

# Move files
mkdir "encoded_${tolvec[$j]}"
mv *.enc "encoded_${tolvec[$j]}/"

done
