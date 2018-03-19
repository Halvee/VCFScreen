#!/usr/bin/env python

'''
Filename : VCFscreen_newlycomphet.py 
Author : Matt Halvorsen
Email : mhalvors1@gmail.com
Date created : 03/19/2018
Date last modified : 03/19/2018
Python version : 2.7
'''

import sys
import os
import argparse
import cyvcf2
from collections import defaultdict
import vcfscreen.screens as screens
import vcfscreen.misc as misc
import vcfscreen.cyvcf2_variant as cyvcf2_variant
from vcfscreen.samples import Samples
from vcfscreen.vcf_cnds import VcfCnds
from vcfscreen.annot import AnnotTxs
try:
    import pandas
except:
    print("ERROR : required module 'pandas' not installed.")
    sys.exit(1)

GT_VARNAMES=['genotypes', 'gt_alt_depths', 'gt_alt_freqs', 'gt_bases',
             'gt_depths', 'gt_phases', 
             'gt_phred_ll_homref','gt_phred_ll_het', 'gt_phred_ll_homalt',
             'gt_quals', 'gt_ref_depths', 'gt_types']

def main(ARGS = None):
    """
    workflow
    1. read fam file, tsv of transmitted het variants
    2. for each proband, cluster variants by max annot gene
    3. find combinations of het transmitted variants that are absent
       from both parents
    """
    if ARGS == None:
        ARGS = sys.argv[1:]
    args = parse_args(ARGS)

    """
    convert certain comma delim str args to lists
    """
    args.qual_impacts = set(misc.str_none_split(args.qual_impacts, ","))

    """
    read samples from fam file, send basic stats to stdout
    """
    samples_i = Samples(args.in_fam)
    samples_i.print_stats()
    n_samples = len(samples_i.samples)
    n_males = len(samples_i.males)
    n_females = len(samples_i.females)

    """
    read het transmitted table into pandas DataFrame
    """
    het_trans = pandas.read_csv(args.het_trans_tsv,
                                sep="\t", header=0)
    
    """
    iterate through all variants, performing cpht screen on each one
    """
    cpht_gts = {}
    prev_chrom = None
    linenum=0
    
    # first pass of transmitted het calls, get qualifying var calls,
    # carrying proband, and parent of origin
    for i in range(het_trans.shape[0]):
        gene=str(het_trans.loc[i,args.gene_col])
        chrom=str(het_trans.loc[i,"CHROM"])
        pos=str(het_trans.loc[i,"POS"]) 
        ref=str(het_trans.loc[i,"REF"]) 
        alt=str(het_trans.loc[i,"ALT"]) 
        par=str(het_trans.loc[i,"PAR_IID"])
        pro=str(het_trans.loc[i,"PRO_IID"])
        impact=str(het_trans.loc[i,"IMPACT_maximpact"])
        if impact not in args.qual_impacts: continue
        if gene not in cpht_gts: cpht_gts[gene]={}
        if pro not in cpht_gts[gene]: cpht_gts[gene][pro]={}
        var_id="-".join([chrom,pos,ref,alt])
        if par not in cpht_gts[gene][pro]:
            cpht_gts[gene][pro][par]=set()
        cpht_gts[gene][pro][par].add(var_id)

    # identify all gene-based instances of compound heterozygousity
    cpht_pro_gene=set()
    for gene in cpht_gts:
        for pro in cpht_gts[gene]:
            if len(cpht_gts[gene][pro])==2:
                cpht_pro_gene.add((pro, gene))
    
    # second pass, select all variants that are part of cpht genotypes
    # as variants to keep and to write to output tsv
    i_keep = []
    for i in range(het_trans.shape[0]):
        gene=str(het_trans.loc[i,args.gene_col])
        chrom=str(het_trans.loc[i,"CHROM"])
        pos=str(het_trans.loc[i,"POS"]) 
        ref=str(het_trans.loc[i,"REF"]) 
        alt=str(het_trans.loc[i,"ALT"]) 
        par=str(het_trans.loc[i,"PAR_IID"])
        pro=str(het_trans.loc[i,"PRO_IID"])
        impact=str(het_trans.loc[i,"IMPACT_maximpact"])
        var_id="-".join([chrom,pos,ref,alt])
        pro_gene = (pro, gene)
        if pro_gene in cpht_pro_gene:
            par_vars=cpht_gts[gene][pro]
            for par in par_vars:
                if var_id in par_vars[par]:
                    i_keep.append(i)

    # write vars that are part of cpht genotypes to out_tsv
    het_trans = het_trans.loc[i_keep, :]
    het_trans.to_csv(path_or_buf=args.out_tsv,
                     sep="\t", header=True, index=False)
    return

def parse_args(ARGS):
    args = argparse.ArgumentParser()
    args.add_argument("in_fam",                                               
                      action="store", type=str,                       
                      help="ped fam file.")
    args.add_argument("het_trans_tsv",	       	      
                      action="store", type=str, 
                      help="tab delim file that was output from " + \
                           "vcfscreen_parentaltransmitted.py.")
    args.add_argument("gene_col",
                      action="store", type=str,
                      help="name of column to get gene identifier from.")
    args.add_argument("impact_col",
                      action="store", type=str,
                      help="name of column to get max impact from")
    args.add_argument("qual_impacts",
       	       	      action="store", type=str,
       	       	      help="comma-delim set of max impacts to keep")
    args.add_argument("out_tsv",
                      action="store", type=str,
                      help="tab-delim output file with vars " + \
                           "that make up compound het genotypes.")
    return args.parse_args(args=ARGS)

if __name__ == "__main__":
    main()
