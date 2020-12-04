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

args = parser.parse_args()

def cleanText(readData):
	text = re.sub('[-=+,#/\?:^$@*\"※~&%ㆍ!』\\‘|\(\)\[\]\<\>`\'…》]', '', readData)
	return text

def locustag(file_path):
    if '.gbk' in os.path.basename(file_path):
        type = 'genbank'
        tag = re.sub('.gbk', '', os.path.basename(file_path))

    else:
        if 'gb' in record:
            type = 'genbank'
            tag = re.sub('.gb', '', os.path.basename(file_path))
        else:
            print("genbank file only")

    ofile = open("{}/{}.tsv".format(os.path.dirname(os.path.abspath(file_path)), tag), "w")

    for rec in SeqIO.parse(file_path, type):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    #print(feature.location)
                    locus = cleanText(str(feature.qualifiers["locus_tag"]))
                    #seqlen = str(len(feature.location.extract(rec).seq))
                    #ofile.write(locus + "\t" + seqlen + "\n")
                    ofile.write(locus+"\n")
    ofile.close()   
    
    
    
import shutil
default_dir = os.path.dirname(os.path.abspath(args.input))
results_dir = default_dir + "/locus"


if not os.path.exists(results_dir):
    os.makedirs(results_dir) 
    
############################run####
locustag(args.input)

tag = re.sub('.gbk', '', os.path.basename(args.input))

shutil.move("{}/{}.tsv".format(default_dir, tag), results_dir)
