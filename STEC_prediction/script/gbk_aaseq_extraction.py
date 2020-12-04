import os
import sys
import re
import argparse 
from Bio import SeqIO


def path(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=path, help="input file path or file name, only 1 file - for many files, uses snakemake")
parser.add_argument("-pre", "--prefix", nargs = "?", default = "", help="output file prefix") 
args = parser.parse_args() 


def wrap(string, max_width):
        s=''
        for i in range(0,len(string),max_width):
            s+=string[i:i+max_width]
            s+='\n'
        return s

def cleanText(readData):
    text = re.sub('[-=+,#/\?:^$@*\"※~&%ㆍ!』\\‘|\(\)\[\]\<\>`\'…》]', '', readData)
    return text
    
def gbk_to_nt_fasta(file_path, outfile_name):
        if '.gbk' in os.path.basename(file_path):
            type = 'genbank'
            tag = re.sub('.gbk', '', os.path.basename(file_path))
                
        else:
            if 'gb' in record:
                type = 'genbank'
                tag = re.sub('.gb', '', os.path.basename(file_path))
            else:
                print("genbank file only, .gbk or .gb")
        
        if args.prefix == "":
            ofile = open("{}/{}_aa.fasta".format(os.path.dirname(os.path.abspath(file_path)), tag), "w")
        else:
            ofile = open("{}/{}.fasta".format(os.path.dirname(os.path.abspath(file_path)), outfile_name), "w")

        for rec in SeqIO.parse(file_path, type):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        #print(feature.location)
                        k = cleanText(str(feature.qualifiers["locus_tag"]))
                        ntseq = wrap(str(feature.location.extract(rec).seq.translate(table=11, to_stop=True)), 60)            
                        ofile.write(">" + k + "\n" + ntseq)
        ofile.close()  
        
import shutil
default_dir = os.path.dirname(os.path.abspath(args.input))
results_dir = default_dir + "/aa_fasta"


if not os.path.exists(results_dir):
    os.makedirs(results_dir) 
    

gbk_to_nt_fasta(args.input, args.prefix)

tag = re.sub('.gbk', '', os.path.basename(args.input))

shutil.move("{}/{}_aa.fasta".format(default_dir, tag), results_dir)
