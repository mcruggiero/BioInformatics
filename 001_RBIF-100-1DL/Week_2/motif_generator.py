#!/usr/bin/env python

# Numpy is a common library, I will use it here,
import numpy as np

###
# Until I get the motif list, I am going to cook myself something up
# I am not going to use numpy for the code below since I am going to translate
# the Python code to Bash.
###

class Motif_Generator:
    def __init__(self, number, min_length, max_length):
        motifs = []
        
        for _ in range(number):
            motifs.append(self.motif_generator(min_length, max_length))
        
        # Write the motifs to a file
        
        with open('motifs.txt', 'w') as f:
            for motif in motifs:
                f.write("{}\n".format(motif))
        
    def motif_generator(self, min_length, max_length):
        
        # Set the size of the random item
        size = np.random.randint(low = min_length, high = max_length)
        
        # Populate that list
        nucleotides = "".join(np.random.choice(["a", "c", "g", "t"], size))
        return nucleotides
    
if __name__ == '__main__':
    # Not really sure how large to make each motif, I suppose it doesn't matter.
    Motif_Generator = Motif_Generator(10, 7, 9)
