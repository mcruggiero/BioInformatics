#!/bin/bash

# [Task 1] Print out the number of occurrences for each motif 
# [Task 2] output output a file called motif_count.txt 
# [Task 3] Create a fasta file for each motif (so 10 in total) which contains motif genes
# [Task 4] Each file should be named after the motif outputted to a new directory

###
# The idea of this script will be to first append the motifs into a list, then 
# Scan through the chromosome
###


# Test to see if directory exists, delete if it does
File="/media/data/Documents/BioInformatics/001_RBIF-100-1DL/Week_2/bash_output"
if [ -d $File ]; then
    printf "$File exists, deleting file\n"
    rm -r bash_output
fi

#Make new directory for output
mkdir /media/data/Documents/BioInformatics/001_RBIF-100-1DL/Week_2/bash_output
mkdir /media/data/Documents/BioInformatics/001_RBIF-100-1DL/Week_2/bash_output/motifs

# I am going to make an array of the entire .fasta file 
# and then search the values by the motifs

chromosomes=()
chromo=()
while IFS= read -r line
do
    if [ ${line:0:4} == ">chr" ]; then
        chromosomes+=("$chromo")
        # Unsetting the chromo variable allows me to reset.
        unset chromo
        chromo+=$(echo "$line" | tr -cd "[:print:]\n")
    else
        chromo+=$(echo "$line" | tr -cd "[:print:]\n")
    fi
done < test1.fasta

#motifs will be our array where values are stored.
motifs=()
printf "The following Motifs will be checked:\n"

#Internal Field Separator (IFS)
while IFS= read -r line
do
    printf "$line\n"
    motifs+=("$line")
done < motifs.txt

printf "\n\n"

# This is the complete counting
for motif in "${motifs[@]}"
do 
    touch $File/motifs/"${motif^^}".txt
    printf "$motif," >> $File/motif_count.txt
    grep -o -i "$motif" test1.fasta | wc -l >> $File/motif_count.txt
done
