#!/usr/bin/env python

import os, shutil
import pandas as pd
import numpy as np
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
        
def create_crispy():
    animal_house = []
    
    #First build dataframe animals in exome dir
    for fasta in os.listdir("exomes"):
        with open('exomes/{}'.format(fasta)) as f:
            lines = [line.rstrip() for line in f]
        entry = {lines[2*i][1:]:lines[2*i + 1] for i in range(0,500)}
        entry["name"] = fasta[:-6]
        animal_house.append(entry)

    animal_house = pd.DataFrame(animal_house)
    animal_house.set_index('name', inplace=True)  
    
    #Iterate through motif_list
    with open("motif_list.txt") as f:
        motifs = [line.rstrip() for line in f]

    # Build counting dataframe for each animal
    animal_data = {}
    count_list = []
    top_per_animal = {}
    
    for animal in list(animal_house.index):
        animal_data[animal] = animal_house.loc[[animal]].T
        for motif in motifs:
            animal_data[animal][motif] = animal_data[animal][animal].apply(lambda a: a.count(motif))
    
        # Suming up columns and rows
        animal_data[animal].loc['{}'.format(animal),:]= animal_data[animal].sum(axis=0)
        
        # There is a quirk in Pandas where the Text column is added, will drop
        animal_data[animal].drop([animal], axis=1, inplace = True)
        
        # Casting to ints
        animal_data[animal] = animal_data[animal].astype(int)
        
        # Add the sum into a count_column list
        count_list.append(animal_data[animal].loc[['{}'.format(animal)]].T.dropna())
        
        # Populate top 3 lists
        top_per_animal[animal] = list(animal_data[animal].loc[animal].nlargest(3).index)


    # Join all of the columns into one dataframe
    count_values = pd.concat(count_list, axis=1)
    
    #Total sum per row: 
    count_values.loc[:,'sum'] = count_values.sum(axis=1)
    
    print("\nTop Ten Motifs")
    print(count_values["sum"].nlargest(10))
    print('\nTop Motifs Per Animal')
    for animal in top_per_animal:
        print(animal)
        print(top_per_animal[animal])
