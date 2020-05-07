#!/usr/bin/env python

import os, shutil
import pandas as pd
import numpy as np
from itertools import combinations, permutations
import datetime
import re

def copy_exomes():
 
    now = datetime.datetime.now()
    
    # Making Report as we go
    with open("report.txt", "a+") as f:
        f.write("\nAnalysis Run Started:\n" + now.strftime("%c") + "\n----------\n")
        f.close()
    
    
    # Read in tab deliniated file
    df = pd.read_csv("clinical_data.txt", sep='\t', engine='python')
    
    # Filter out by needed terms
    mover = df[(20 <= df['Diamater (mm)']) & (df['Diamater (mm)'] <= 30)]
    
    with open("report.txt", "a+") as f:
        f.write("The following entries are of interest:\n")
        f.write("{}".format(mover) + "\n")
        f.close()
    
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
        with open("report.txt", "a+") as f:
            f.write("copying {}.fasta file\n".format(i))
        shutil.copy2("exomes_given/{}.fasta".format(i), 
                     "exomes/{}.fasta".format(i))
        
def create_crispy():
    animal_apartment = []
    animal_list = {}
    gene_dict = {}

    #First build dataframe animals in exome dir
    for fasta in os.listdir("exomes"):
        with open('exomes/{}'.format(fasta)) as f:
            lines = [line.rstrip() for line in f]

        animal_list[fasta[:-6]] = lines
        entry = {lines[2*i]:lines[2*i + 1] for i in range(0,500)}
        entry["name"] = fasta[:-6]
        animal_apartment.append(entry)

    animal_house = pd.DataFrame(animal_apartment)
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

        # Look at only the populated data lists
        df = animal_data[animal][top_per_animal[animal]]
        df = df[(df.T != 0).any()]

        # Make animal dict to convert back to text file
        gene_dict[animal] = df.T.to_dict()

    # Join all of the columns into one dataframe
    count_values = pd.concat(count_list, axis=1)

    #Total sum per row: 
    count_values.loc[:,'sum'] = count_values.sum(axis=1)

    report_dict = {}

    # Build a report dictionary 
    for animal in gene_dict:
        report_dict[animal] = {}
        for gene in gene_dict[animal]:
            d = gene_dict[animal][gene]

            #Flatten Dictionary
            report_dict[animal][gene] = ' '.join("{!s}_{!r}".format(key,val) for (key,val) in d.items())

    # Build the final report with proper formating
    report = {}
    for animal in animal_list:
        text_list = []
        for line in animal_list[animal]:
            if line in report_dict[animal]:
                text_list.append(line + " " + report_dict[animal][line])
                text_list.append(animal_house.loc[animal][line])

        report["{}_topmotifs.fasta".format(animal)] = text_list

    # Send the report dictionary to .txt files, not sure where you want them
    # Sending to different directory "exomes_top"

    # Make Directory Structure
    if not os.path.exists("exomes_top"):
        os.mkdir("exomes_top")
    else:
        shutil.rmtree("exomes_top")
        os.mkdir("exomes_top")

    # Create Files
    for file in report:
        with open("exomes_top/{}".format(file), "w") as crispy:
            for line in report[file]:
                crispy.write("{}\n".format(line))

    similar = animal_house.T.to_dict()
    count_dict = {}
    gene_list = []

    # First make a list of all genes
    for animal in similar:
        for gene in similar[animal]:
            gene_list.append(similar[animal][gene])

    # Search for duplicates
    for animal in similar:
        count_dict[animal] = {}
        for gene in similar[animal]:
            count = gene_list.count(similar[animal][gene]) - 1
            count_dict[animal][gene] = count
            if count > 0: print("shared gene!")

    shared_genes = pd.DataFrame(count_dict)
    
    with open("report.txt", "a+") as f:
        f.write("\nShared Genes\n")
        f.write("{}".format(shared_genes[shared_genes != 0].dropna()))
        f.write("\n\nTop Ten Motifs\n")
        f.write("{}".format(count_values["sum"].nlargest(10)))
        f.write('\n\nTop Motifs Per Animal\n')
        for animal in top_per_animal:
            f.write("\n{}\n".format(animal))
            f.write("{}".format(top_per_animal[animal]))
        f.write("\n")
        f.close()
        
def locate():
    
    top_house = {}
    soggy_dict = {}

    # First build dictionary for candidates
    for fasta in os.listdir("exomes_top"):
        with open('exomes_top/{}'.format(fasta)) as f:
            lines = [line.rstrip() for line in f]

        entry = {lines[2*i]:lines[2*i + 1] for i in range(0,len(lines)//2)}
        top_house[fasta] = entry

    # Search through string for upsteam GG
    for animal in top_house:
        soggy_list = []
        have_gg = 0
        no_gg = 0

        # Reporter
        with open("report.txt", "a+") as f:
            f.write("\nSearching {} for NGG\n-------------\n".format(animal))
            f.close()

        for gene in top_house[animal]:
            candidate = top_house[animal][gene]
            search = re.finditer("GG", candidate)
            indices = [pos for pos in [m.start(0) for m in search] if pos > 19]

            if len(indices) > 0:
                soggy_list.append(gene)
                soggy_list.append(candidate)
                have_gg += 1

                # Reporter
                with open("report.txt", "a+") as f:
                    f.write("Found NGG at {}\n".format(gene.split(" ")[0]))
                    f.close()

            else:
                no_gg += 1

            # Reporter
        with open("report.txt", "a+") as f:
            f.write("On {}, found {} NGG genes out of {} candidates\n".format(animal,
                                                                              have_gg,
                                                                              have_gg + no_gg))
            f.write("-----------------\n")
            f.close()
        soggy_dict[animal] = soggy_list


    # Make Directory Structure
    if not os.path.exists("exomes_precrispr"):
        os.mkdir("exomes_precrispr")
    else:
        shutil.rmtree("exomes_precrispr")
        os.mkdir("exomes_precrispr")

    # Create Files
    for file in soggy_dict:
        listing = file.split("_")[0] + "_precrispr.fasta"
        with open("exomes_precrispr/{}".format(listing), "w") as crispy:
            for line in soggy_dict[file]:
                crispy.write("{}\n".format(line))
                
def edit_genome():
    big_house = {}
    crispy_dict = {}

    # First build dictionary for candidates
    for fasta in os.listdir("exomes_precrispr"):
        with open('exomes_precrispr/{}'.format(fasta)) as f:
            lines = [line.rstrip() for line in f]

        entry = {lines[2*i]:lines[2*i + 1] for i in range(0,len(lines)//2)}
        big_house[fasta] = entry

    # Search through string for upsteam GG
    for animal in big_house:
        crispy_list = []

        # Reporter
        with open("report.txt", "a+") as f:
            f.write("Inserting A at NGG sites on file {}\n-------------\n".format(animal))
            f.close()

        for gene in big_house[animal]:
            candidate = big_house[animal][gene]
            search = re.finditer("GG", candidate)
            indices = [pos for pos in [m.start(0) for m in search] if pos > 19]
            for ind in indices:
                candidate = candidate[:ind] + "A" + candidate[ind:]

            if len(indices) > 0:
                crispy_list.append(gene)
                crispy_list.append(candidate)

        crispy_dict[animal] = crispy_list

    # Make Directory Structure
    if not os.path.exists("exomes_postcrispr"):
        os.mkdir("exomes_postcrispr")
    else:
        shutil.rmtree("exomes_postcrispr")
        os.mkdir("exomes_postcrispr")

    # Create Files
    for file in crispy_dict:
        listing = file.split("_")[0] + "_postcrispr.fasta"
        with open("exomes_postcrispr/{}".format(listing), "w") as crispy:
            for line in crispy_dict[file]:
                crispy.write("{}\n".format(line))
