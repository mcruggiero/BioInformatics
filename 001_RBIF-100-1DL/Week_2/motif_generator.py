#!/usr/bin/env python

import random

###
# If you want a little python motif_generator to test out your bash scripts
###

class Motif_Generator:
    
    # The first thing this class will do is run the __init__
    def __init__(self, number, min_length, max_length):    
        # Make an empty list of motifs
        motifs = []
        
        # Use the motif_generator function to populate the motif list 
        for _ in range(number):
            motifs.append(self.motif_generator(min_length, max_length))
        
        # Write the motif list to file to a file
        with open('motifs.txt', 'w') as f:
            for motif in motifs:
                f.write("{}\n".format(motif))
        
    def motif_generator(self, min_length, max_length):
        # Set the size of the random item
        size = random.randint(min_length, max_length)
        
        # Make a blank list for the nucleotides
        nucleotides = []
        
        # Choose a nucleotide and append it to the list
        for _ in range(size):
            nucleotides.append(random.choice(["a", "c", "g", "t"]))

        # Turn that list into a string
        nucleotides = "".join(nucleotides)
        
        return nucleotides
    
if __name__ == '__main__':
    # When the class is called, it will make 10 motifs of sizes between 7 and 9
    Motif_Generator = Motif_Generator(10, 7, 9)
