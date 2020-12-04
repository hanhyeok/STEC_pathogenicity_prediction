import os
import re
import numpy as np
import pandas as pd
import pickle
import argparse

def path(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--diamond_result", type=path, help="input file path or file name, only 1 file - for many files, uses snakemake1") 
parser.add_argument("-l", "--locus_tag", type=path, help="input file path or file name, only 1 file - for many files, uses snakemake3") 
args = parser.parse_args()

with open('/STEC_prediction/rep_seq_cluster.pickle', 'rb') as handle:
    rep_seq_cluster = pickle.load(handle)
with open('/STEC_prediction/threshold_cluster_pid.pickle', 'rb') as handle:
    threshold_cluster_pid = pickle.load(handle)
with open('/STEC_prediction/threshold_cluster.pickle', 'rb') as handle:
    threshold_cluster = pickle.load(handle)
with open('/STEC_prediction/pirate_cluster_list.pickle', 'rb') as handle:
    pirate_gene_list = pickle.load(handle)

def add_record(sseqid):
    cluster = rep_seq_cluster[sseqid][0]
    pid_threshold = threshold_cluster_pid[cluster]
    return cluster, pid_threshold


def pirate_gene_classifier(diamond_result, locus_tag):
  
    diamond = pd.read_csv(diamond_result, sep='\t', header=None)
    diamond.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    
    diamond_filtered = pd.DataFrame()
    diamond_filtered = diamond[diamond['pident'] >= 50].reset_index()
    
    diamond_filtered['cluster'] = np.nan
    diamond_filtered['pid'] = np.nan
    diamond_filtered['QC'] = np.nan
    
    diamond_filtered['cluster'] = diamond_filtered.apply(lambda x:add_record(x['sseqid'])[0], axis= 1)
    diamond_filtered['pid'] = diamond_filtered.apply(lambda x:add_record(x['sseqid'])[1], axis= 1)
    diamond_filtered['QC'] = diamond_filtered.apply(lambda x: 'PASS' if add_record(x['sseqid'])[1] <=  x['pident'] else 'FAIL', axis= 1)
    
    new = diamond_filtered[diamond_filtered['QC'] == 'PASS']
    
    locus = pd.read_csv(locus_tag, sep = '\t', header=None)
    
    assignment = {}
    
    
    for lc in locus[0]:
        if lc in list(set(new['qseqid'])):
            cluster = new[new['qseqid']==lc].iloc[0]['cluster']
            assignment[lc] = cluster
        else:
            assignment[lc] =  np.nan
            pass
        
    assignment_pd = pd.DataFrame.from_dict(assignment, orient='index')
    assignment_pd = assignment_pd.reset_index()
    result = {cluster: set(assignment['index'].to_list()) for cluster, assignment in assignment_pd.groupby(0)}
    
    tag = re.sub('.blastp', '', os.path.basename(diamond_result))
    output = pd.DataFrame(pirate_gene_list)
    output[tag] = np.nan
    output[tag] = output.apply(lambda x: list(result[x[0][0:x[0].find('__')]])[0] if x[0][0:x[0].find('__')] in result.keys() else 0, axis=1)  
       
    output.to_csv('{}_assigned_locus.tsv'.format(tag), sep='\t')
     



pirate_gene_classifier(args.diamond_result, args.locus_tag)
   
import shutil


results_dir =  "/data/classification"

if not os.path.exists(results_dir):
    os.makedirs(results_dir) 
    
tag = re.sub('.blastp', '', os.path.basename(args.diamond_result))

shutil.move("{}_assigned_locus.tsv".format(tag), results_dir)
