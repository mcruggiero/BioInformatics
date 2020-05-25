from parseFastq import ParseFastQ
import re
import os
import shutil

class Main:
    def __init__(self):
        coding_dict = self.trimmer("harrington_clinical_data.txt", "hawkins_pooled_sequences.fastq")
        self.alignment()

    # Part I: Trim and report
    def trimmer(self, clinical_path, fast_path):
        #First we are going to make a dictionary with all of the clinical entries
        with open(clinical_path) as clinical:
            rows = (line.split('\t') for line in clinical)
            coding_dict = {row[1:][1].rstrip():{"name":row[0], "color":row[1]}  for row in rows}

        # Remove the first entry, we are not going to use a Pandas style index here
        coding_dict.pop('Barcode')

        fastqfile = ParseFastQ(fast_path)
        for i in fastqfile:
            # I don't think we need item 2 which is always a "+", so we will omit it

            coding_dict[i[1][0:5]][i[0]] = {"pretrimmed":i[1][5:], "prequality":i[3]}
            pretrim = coding_dict[i[1][0:5]][i[0]]["prequality"]

            #Search for consecutive D or F 
            trim_value = re.search(r'[DF][DF]', pretrim).start()
            coding_dict[i[1][0:5]][i[0]]["trim_value"] = trim_value
            coding_dict[i[1][0:5]][i[0]]["trimmed"] = coding_dict[i[1][0:5]][i[0]]["pretrimmed"][:trim_value]
            coding_dict[i[1][0:5]][i[0]]["quality"] = coding_dict[i[1][0:5]][i[0]]["prequality"][:trim_value]

        # Make a new directory unless it exists
        try:
            os.makedirs("fastqs")
            print("Made New Directory")
        except:
            print("Removing Old Directory")
            shutil.rmtree("fastqs")
            os.makedirs("fastqs")

        for listing in coding_dict:
            # First name file
            f= open("fastqs/{}_trimmed.fastq".format(coding_dict[listing]["name"]),"w+")
            for seq in coding_dict[listing]:
                # I kept the name and color in the dictionary, we can't add that to fastq file
                if seq not in {"name", "color"}:
                    f.write(seq + "\n")
                    f.write(coding_dict[listing][seq]["trimmed"] + "\n")
                    f.write("+" + "\n")
                    f.write(coding_dict[listing][seq]["quality"] + "\n")
            
        # Part 1 completed
        return coding_dict
    
    # Part 2&3 alignment
    def alignment(self):
        os.system("bwa index dgorgon_reference.fa")
        for name in os.listdir("fastqs"):
            split = name.split("_")[0]
            os.system("bwa mem dgorgon_reference.fa fastqs/{0}_trimmed.fastq > fastqs/{0}.sam".format(split))
            os.system("samtools view -bS fastqs/{0}.sam > fastqs/{0}.bam".format(split))
            os.system("samtools sort -m 100M -o fastqs/{0}.sorted.bam fastqs/{0}.bam".format(split))
            os.system("samtools index fastqs/{0}.sorted.bam".format(split))
            # Part 2&3 Completed
            
        
    
    
if __name__== "__main__":    
    Main()
