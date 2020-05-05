#!/usr/bin/env python

import os, shutil
import pandas as pd
from itertools import combinations, permutations

def copy_exomes():
    # Read in tab deliniated file
    df = pd.read_csv("clinical_data.txt", sep='\t', engine='python')
    
    # Filter out by needed terms
    mover = df[(20 <= df['Diamater (mm)']) & (df['Diamater (mm)'] <= 30)]
    print("The following entries are of interest:")
    print(mover)
    
    # Make Directory Structure
    if not os.path.exists("exomes"):
        os.mkdir("exomes")
    else:
        shutil.rmtree("exomes")
        os.mkdir("exomes")
    
    # Check to see if path exists, it it does, remove and create
    if not os.path.exists("exomes"):
        os.mkdir("exomes")
    else:
        shutil.rmtree("exomes")
        os.mkdir("exomes")

    # Copy files
    for i in mover["code_name"].tolist():
        print("copying {}.fasta file".format(i))
        shutil.copy2("exomes_given/{}.fasta".format(i), 
                     "exomes/{}.fasta".format(i))
