#!/bin/bash

# [Task 1] Print out the number of occurrences for each motif 
# [Task 2] output output a file called motif_count.txt 
# [Task 3] Create a fasta file for each motif (so 10 in total) which contains motif genes
# [Task 4] Each file should be named after the motif outputted to a new directory

###
# The idea of this script will be to first append the motifs into a list, then 
# Scan through the chromosome
###

motifs=()
printf "The following Motifs will be checked:\n"

#Internal Field Separator (IFS)
while IFS= read -r line
do
    printf "$line\n"
    motifs+=("$line")
    
done < motifs.txt

a=1
for i in "${motifs[@]}"
do 
   a+=1
   echo $a
done
