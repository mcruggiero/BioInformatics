#!/usr/bin/env python3

#####
# Michael Ruggiero August 18th, 2020
#####

#####
# Prompt
#####
# (2.1) Using NCBI or Unitprot to find the protein sequence of the gene PMS2. 
# Can you use BLAT and BLAST to search for its homologous genomic regions/genes, 
# using the whole or part of the PMS2 gene sequence (for example, the first 120 
# or last 120 amino acids)? Can you list the homologous genes you find and the 
# corresponding genomic coordinates as well as the PMS2 amino acids position 
# (e.g, position 1 to 120) if part of the whole sequence is used?
#
# (2.2) Find protein sequences of human, mouse,  rat and zebrafish sequences for 
# the protein FGFR4 in Unitprot, Align the 4 sequences with a MSA tool. How many 
# conserved regions can you identify? Can you use literature to verify your 
# results?
######


#Import all of these
import os, sys
import pandas as pd
import numpy as np
import argparse

# Calculate distance between strings
from Levenshtein import distance as l_dist

# Basic Libraries
from random import choices, randint
from collections import Counter, defaultdict
import re
from glob import glob

# Pretty Printer
import textwrap

# Online Libraries
import mygene
import ensembl_rest
import requests
import json

# Biopython scripts
from Bio import AlignIO, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna,generic_rna,generic_protein
from Bio.pairwise2 import format_alignment
from Bio.PDB import *

# This library can be a little bit of a pain to install. Please read online
# help
from Bio.Align.Applications import MafftCommandline

TF = textwrap.TextWrapper(width=80)
#####
# PART 1
######

# Change this value to explore other genes
ref_gene='PMS2'
ref_SPC='human'
ref_field='ensembl.gene'

# This is a little sloppy, but it works to query the ENSG id for the gene
mg        = mygene.MyGeneInfo()
out       = mg.querymany(ref_gene, scopes='symbol', fields=ref_field, species=ref_SPC)
gene      = out[0]['ensembl']['gene']
data      = ensembl_rest.sequence_id(gene) 
ensg      = data["query"]

# Now we need to grab the Fasta files from online
server    = "https://rest.ensembl.org"
ext       = "/xrefs/id/{}?".format(ensg)
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

# A little warning for us 
if not r.ok:
    r.raise_for_status()
    sys.exit()
 
# decodes json
decoded = r.json()

# Grabs the potential values
potential = []
for item in decoded:
    if item["dbname"] == "Uniprot_gn":
        potential.append(item["primary_id"])
        
link = 'https://www.uniprot.org/uniprot/'
urls = [link + protein + ".fasta" for protein in potential]
raw = [str(requests.get(x).text).split("\n") for x in urls]

fastas = {}
for protein in raw:
    fastas[protein[0].split("|")[1]] = "".join(protein[1:])

# Write part 1a to a fasta file
outF = open("part1a.fasta", "w")
for key, values in fastas.items():
    outF.write(">{}".format(key))
    outF.write("\n")
    outF.write(values)
    outF.write("\n")
outF.close()

# I kept getting timed out up on my BLAST api search, so I used the online api to get 
# the blast results instead as blast_results.fasta. (Prompt allowed)
# Giving up on this for now, maybe internet problem?
# from Bio.Blast import NCBIWWW
# result_handle = NCBIWWW.qblast("blastn", "nt", fastas["P54278"])

def aligner(fasta_file):
    with open(fasta_file) as f:
        # Set output file name
        aligned_output = "{}-aligned.fasta".format(fasta_file)

        # Apply MAFFT alignement
        mafft_cline = MafftCommandline(input= str(fasta_file))
        aligned, report = mafft_cline()
        
        # Uncomment to view MAFFT report
        # print(report)
        with open(aligned_output, "w") as handle:
            handle.write(aligned)

aligner("blast_results.fasta")

# Organize downloaded Fasta file
proteins = {}
for record in SeqIO.parse("blast_results.fasta", "fasta"):
    # Create an empty dictionary to populate
    proteins[record.id.split("|")[1]] = entry = {}
    
    # Set values to explore
    entry["description"] = record.description
    entry["species"] = entry["description"].split("=")[1][:-3]
    entry["organism_id"] = entry["description"].split("=")[2][:-3]
    entry["gene_name"] = entry["description"].split("=")[3][:-3]
    entry["protein"] = str(record.seq)

for record in SeqIO.parse("blast_results.fasta-aligned.fasta", "fasta"):
    key = record.id.split("|")[1]
    proteins[key]["aligned"] = str(record.seq)    

# Everything should be a nice little class organization, but I just don't have
# The time. Will pay back technical debt later. This class does a lot more than
# I can showcase here.

class Protein_Describer:
    def __init__(self, proteins_dict):
        """Creates useful information compiled from protein dictionary"""
    
        self.proteins_dict = proteins_dict
        self.proteins = list(proteins_dict.keys())
    
    def delta(self, protein_1, protein_2, place):
        "A simple Delta function for proteins"
        if protein_1[place] == protein_2[place]:
            return protein_1[place]
        else:
            return "-" 
    
    def shared(self, aligned_a, aligned_b):
        """A little shared protein creater from aligned sequences"""
        return "".join([self.delta(aligned_a, aligned_b, i) for i in range(len(aligned_a))])
    
    def proteins_align(self, protein_a, protein_b):
        """Aligns to proteins with BioPython"""
        # Set variables
        first = Seq(self.proteins_dict[protein_a]["protein"])
        second = Seq(self.proteins_dict[protein_b]["protein"])
        
        # Align proteins
        align = pairwise2.align.globalxx(first, second, one_alignment_only=True)
        aligned_a = align[0].seqA
        aligned_b = align[0].seqB
        
        # Calculate shared string
        shared = self.shared(aligned_a, aligned_b)

        # Returns dictionary of shared terms
        return {protein_a: aligned_a, 
                protein_b: aligned_b,
                "shared":  shared,
                "shared_count": Counter([x for x in shared.split("-") if x != ""]),
                "percent_simalarity": align[0].score / len(align[0].seqA),
                "score": align[0].score, 
                "levenshtein_distance": l_dist(str(first), str(second))}
    
    def conserved_regions(self):
        conserve_test_list = [self.proteins_dict[x]["aligned"] for x in self.proteins_dict]
        shared_regions = conserve_test_list[0]

        # Sloppy indexing, but I am in a rush
        for i in range(len(conserve_test_list)):
            for j in range(i + 1, len(conserve_test_list)):

                # return index of lists
                a, b = conserve_test_list[i], conserve_test_list[j]
                shared_regions = self.shared(shared_regions, self.shared(a, b))

        return shared_regions

# There is a lot of fun stuff to do here, but I just don't think I have enough
# time to cook up a distance matrix. I will throw it up on the discussion board

protein_comparisons = Protein_Describer(proteins)
conserved = protein_comparisons.conserved_regions()

# Write part 1 conclusion to file
outF = open("part1_conserved_for_all.txt", "w")
outF.write("Regions Conserved in all proteins\n")
outF.write("Prompt does not ask for limit in BLAST, thus fewer alignments found\n")
outF.write("\n")
pp = TF.fill(conserved)
outF.write(pp)
outF.write("\n")
outF.close()

#####
# PART 2
#####

# Change this value to explore other genes
ref_gene='FGFR4'

def unit_pro(ref_gene,ref_SPC):
    mg     = mygene.MyGeneInfo()
    out    = mg.querymany(ref_gene, scopes='symbol', fields='ensembl.gene', species=ref_SPC)
    gene   = out[0]['ensembl']['gene']
    server = "https://rest.ensembl.org"
    ext    = "/xrefs/id/" + gene + "?"
    r = requests.get(server + ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
        
    decoded = r.json()
    potential = []
    for item in decoded:
        if item["dbname"] == "Uniprot_gn":
            potential.append(item["primary_id"])
    return potential

# The species to explore from the prompt
species=['human', 'mouse','rat','zebrafish']

# Makes a dictionary with species and protein list
potential = {}
for each in species:
    potential[each] = unit_pro(ref_gene,each)

# Reverse key and values from dictionary 
search_items = {}
for i in potential:
    for j in potential[i]:
        search_items[j] = {}
        search_items[j]["species"] = i

# Build a fasta dictionary to explore
link = 'https://www.uniprot.org/uniprot/'
urls = [link + protein + ".fasta" for protein in search_items]
raw = [str(requests.get(x).text).split("\n") for x in urls]

fastas = {}
for protein in raw:
    key = protein[0].split("|")[1]
    fastas[key] = "".join(protein[1:])

# Add protein length to the dictionary as well
for key in search_items: 
    search_items[key]["protein"] = fastas[key]
    search_items[key]["length"]  = len(fastas[key])

# This is an important desision node here. I am going to take only the largest
# proteins from each species. The goal is to have a useful alignment with many
# overlaps. I also ran all proteins with no size restrictions, and I think
# the values were less useful, however it is included in the folder if you are
# interested in looking it over.

# Find largest proteins for each species
sizer = {k:0 for k in ['human', 'mouse','rat','zebrafish']}
for key in search_items:
    species = search_items[key]["species"]
    sizer[species] = max(sizer[species], search_items[key]["length"])

largest_searches = {}
for key in search_items:
    species = search_items[key]["species"]
    if search_items[key]["length"] == sizer[species]:
        largest_searches[key] = search_items[key]

# Build fasta list for largest proteins in 
outF = open("FGFR4_largest.fasta", "w")
for key in largest_searches:
    outF.write(">{} {}".format(key, largest_searches[key]["species"]))
    outF.write("\n")
    outF.write(largest_searches[key]["protein"])
    outF.write("\n")
outF.close()

# Use Aligner function from above
aligner("FGFR4_largest.fasta")

for record in SeqIO.parse("FGFR4_largest.fasta-aligned.fasta", "fasta"):
    # Create an empty dictionary to populate
    #key = record.id.split("|")[1]
    #proteins[key]["aligned"] = str(record.seq)    
    largest_searches[record.id]["aligned"] = str(record.seq)

largest_comparisons = Protein_Describer(largest_searches)
largest_conserved = largest_comparisons.conserved_regions()

# Write part 2 conclusion to file
outF = open("part1_conserved_for_all.txt", "w")
outF.write("Regions Conserved in all proteins\n")
outF.write("Prompt does not ask for limit in BLAST, thus fewer alignments found\n")
outF.write("\n")
pp = TF.fill(largest_conserved)
outF.write(pp)
outF.write("\n")
outF.close()