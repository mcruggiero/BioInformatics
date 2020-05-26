#!/usr/bin/env python

from parseFastq import ParseFastQ
import re
import os
import shutil
import pysam
import datetime

class Main:
    def __init__(self, har_clinic, hawkins_pooled):
        
        # Load reference
        with open("dgorgon_reference.fa") as f:
            dgogron = f.readlines()
        dgorgon_reference = [x.strip() for x in dgogron][1]
        
        # Create Report
        report = open("report.txt", "w")
        now = datetime.datetime.now()
        report.write("\nAnalysis Run Started:\n" + now.strftime("%c") + "\n---samples---\n")
        report.close()
        
        coding_dict = self.trimmer(har_clinic, hawkins_pooled)
        self.alignment()
        
        color_report = {}
        for i in coding_dict:
            file_path = "fastqs/{}.sorted.bam".format(coding_dict[i]["name"])
            coding_dict[i]["mutation"] = self.pileup(file_path)
            coding_mutation = coding_dict[i]["mutation"]
            expected = dgorgon_reference[coding_mutation["position"]]
            
            # This syntax is terrible, but I am running out of time
            # for big big sets, it actually might be faster, however
            mutation = coding_mutation["possible"] - {expected}
            mutation = list(mutation)[0]
            
            # Color Keeping: I think we should also have a count, but not asked for...
            if coding_dict[i]["color"] not in color_report:
                color_report[coding_dict[i]["color"]] = {"expected": expected,
                                                         "mutation": mutation,
                                                         "position": coding_mutation["position"]}
                
            
            with open("report.txt", "a") as myfile:
                report_text = "Sample {0} had a {1} mold, " \
                              "{2} reads, and had {3}% of the reads at " \
                              "position {4} had " \
                              "the mutation {5}. \n".format(coding_dict[i]["name"],
                                                            coding_dict[i]["color"],
                                                            coding_mutation["reads"],
                                                            round(coding_mutation[mutation],2),
                                                            coding_mutation["position"],
                                                            mutation)
                
                
                myfile.write(report_text)
        
        with open("report.txt", "a") as myfile:
            myfile.write("\n---mold-colors---\n")
            
        for color in color_report:
            with open("report.txt", "a") as myfile:
                report_text = "The {0} mold was caused by a mutation " \
                              "in position {1}. The wildtype base was {2} " \
                              "and the mutation was {3}. \n".format(color,
                                                                    color_report[color]["position"],
                                                                    color_report[color]["expected"],
                                                                    color_report[color]["mutation"])
            
                myfile.write(report_text)

        
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
    
    # Part 4 
    def pileup(self,path):
        samfile = pysam.AlignmentFile(path, "rb")
        count_dict = {}
        for item in samfile.pileup():
            total_count = item.get_num_aligned()
            summary = item.get_query_sequences()
            count_dict[item.pos] = {x:100*summary.count(x)/total_count for x in summary}

        mutation_dict = {}
        for i in count_dict.keys():
            if len(count_dict[i].keys()) > 1:
                mutation_dict = count_dict[i]
                mutation_dict["possible"] = set(count_dict[i].keys())
                mutation_dict["reads"] = total_count
                mutation_dict["position"] = i
                
        samfile.close()
        return mutation_dict

if __name__== "__main__":    
    
    ###
    # Important: changes these TXTs to match your file structure
    ###
    Main("harrington_clinical_data.txt", "hawkins_pooled_sequences.fastq")
