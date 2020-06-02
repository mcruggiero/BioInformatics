from Bio.Seq import translate
import re
import ensembl_rest

# This structure is overkill, but why not?
class Main:
    def __init__(self):
        fasta_dict, longest_chain = self.part_2()
        self.part_3()
        
    def part_2(self):
        fasta_dict = ensembl_rest.sequence_id(
             'ENSG00000258839',species='human',
             headers={'content-type': 'fasta'})

        dna_chain = fasta_dict["seq"]
        
        # Regex is standard
        chain_frame = re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',dna_chain)
        longest_chain = max(chain_frame, key=len)
        
        fast_a_fasta = open("MC1R.fasta","w") 
        fast_a_fasta.write(">{}\n".format(fasta_dict["desc"]))
        fast_a_fasta.write(fasta_dict["seq"])
        fast_a_fasta.write("\n+\n")
        fast_a_fasta.write(translate(longest_chain))
        fast_a_fasta.close()
        
        return fasta_dict, longest_chain
    
    def part_3(self):
        homology_dict = ensembl_rest.homology_ensemblgene('ENSG00000258839')
        
        # I am not quite sure what you mean by a "unique list", I am going to
        # make a csv with every line as a new target species
        holy_fasta = open("mc1r_homology_list.txt","w") 
        unique_species = []
        for species in homology_dict["data"][0]['homologies']:
            if species["target"]["species"] not in unique_species:
                holy_fasta.write("{},\n".format(species["target"]["species"]))
                unique_species.append(species["target"]["species"])
        holy_fasta.close()

if __name__== "__main__": 
    Main()
