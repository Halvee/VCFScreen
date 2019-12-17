#!/usr/bin/env python

'''
Filename : VCFscreen_listcarriers.py
Author : Matt Halvorsen
Email : mhalvors1@gmail.com
Date created : 12/17/2019
Date last modified : 12/17/2019
Python version : 2.7
'''

import sys
import os
import argparse
import cyvcf2
from collections import defaultdict
import pandas as pd
import numpy as np

import vcfscreen.screens as screens
import vcfscreen.misc as misc
from vcfscreen.cyvcf2_variant import Cyvcf2Vcf,Cyvcf2Variant
from vcfscreen.samples import Samples
from vcfscreen.vcf_cnds import VcfCnds
from vcfscreen.annot import AnnotTxs

def main():
    
    """
    read user-provided args
    """
    args = parse_args()

    """
    read samples from fam file, send basic stats to stdout
    """
    samples_i = Samples(args.in_fam)
    n_samples = len(samples_i.samples)
    n_males = len(samples_i.males)
    n_females = len(samples_i.females)
    n_cases = len(samples_i.cases)
    n_ctrls = len(samples_i.ctrls)

    """
    load fam file into df. Coverage statistics will be stored here.
    """
    df = pd.read_table(args.in_fam, header=None, sep=" ")
    df = df.rename({0:"FID",1:"IID",2:"PID",3:"MID",4:"SEX",5:"PHENO"}, 
                   axis='columns')    
    df = df.set_index("IID", drop=False)

    """
    if varlist file defined, load into memory
    """
    if args.varlist != None:
        varlist = set([])
        varlist_fh = open(args.varlist, "r")
        for line in varlist_fh: varlist.add(line.rstrip())
        varlist_fh.close()
        args.varlist = varlist
    
    """
    init cyvcf2 VCF obj
    """
    vcf = cyvcf2.VCF(args.in_gt_vcf, strict_gt=True, gts012=True)

    """
    create sample idx
    """
    samples_i.get_vcf_idx(vcf.samples)

    """
    form dictionary of idxs
    """
    case_idxs = [samples_i.samples[x].idx for x in samples_i.cases]
    ctrl_idxs = [samples_i.samples[x].idx for x in samples_i.ctrls]
    idxs=dict()
    for x in samples_i.samples:
        idx = samples_i.samples[x].idx
        idxs[idx]=x

    """
    main loop for pulling out genotypes
    """
    nvar=0
    nbuff=0
    prev_chrom=None
    for vcf_variant in vcf:
        
        if args.varlist != None:
            if vcf_variant.ID not in args.varlist: continue

        # if nvar == 10000: break

        # get genotypes at variant site
        gt = vcf_variant.gt_types

        # get indeces for het or homalt genotypes
        if args.gts_include == "1":
            idxs_i = np.where(gt == 1)[0]
        elif args.gts_include == "2":
            idxs_i = np.where(gt == 2)[0]
        else:
            idxs_i = np.where((gt == 1) | (gt == 2))[0]

        
        # get iids to go along with idxs, append to output file
        for idx_i in idxs_i:
            iid_i = idxs[idx_i]
            out_str = "\t".join([vcf_variant.ID, iid_i])
            print(out_str)

        nvar += 1
                

    vcf.close()

    return

def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--varlist", action="store", type=str, default=None,
                      help="file with list of variants that qualify for " + \
                           "inclusion.")
    args.add_argument("--gts-include", action="store", type=str, 
                      default="12", choices=["1", "2", "12"],
                      help="which genotypes to include.")
    args.add_argument("in_fam",
                      help="input fam with with sample IDs.")
    args.add_argument("in_gt_vcf",
                      help="input genotype vcf")
    return args.parse_args()

if __name__ == "__main__":
    main()
