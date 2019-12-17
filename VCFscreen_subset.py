#!/usr/bin/env python

'''
Filename : VCFscreen_subset.py 
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

# if VCFscreen not in PYTHONPATH, make sure it is
#VCFSCREEN_DIR="/proj/sens2016011/mhalvors/src/VCFscreen/"
#sys.path.append(VCFSCREEN_DIR)

import vcfscreen.screens as screens
import vcfscreen.misc as misc
from vcfscreen.cyvcf2_variant import Cyvcf2Vcf,Cyvcf2Variant
from vcfscreen.samples import Samples
from vcfscreen.vcf_cnds import VcfCnds
from vcfscreen.annot import AnnotTxs

GT_VARNAMES=['genotypes', 'gt_alt_depths', 'gt_alt_freqs', 'gt_bases',
             'gt_depths', 'gt_phases', 'gt_phred_ll_het', 'gt_phred_ll_homalt',
             'gt_phred_ll_homref', 'gt_quals', 'gt_ref_depths', 'gt_types']

def main(ARGS = None):
    if ARGS == None:
        ARGS = sys.argv[1:]
    args = parse_args(ARGS)

    """
    convert certain comma delim str args to lists
    """
    args.qual_impacts = misc.str_none_split(args.qual_impacts, ",")
    args.max_impact_csqs = misc.str_none_split(args.max_impact_csqs, ",")
    args.max_csq_scores = misc.str_none_split(args.max_csq_scores, ",")
    args.min_csq_scores = misc.str_none_split(args.min_csq_scores, ",")

    """
    read cnds files
    """
    var_cnds = None
    if args.variant_cnds != None: var_cnds = VcfCnds(args.variant_cnds)

    """
    init cyvcf2 VCF obj, get info subfields, header for output
    """
    vcf = cyvcf2.VCF(args.in_vcf, strict_gt=True)
    cyvcf2_vcf = Cyvcf2Vcf(vcf)
    cyvcf2_vcf.get_info_subfields()
    if args.annotation_subfield == "ANN":
        cyvcf2_vcf.get_csq_keys(spliton="Functional annotations: ", delim="|",
                                chars_del=[" ", "'", '"'], 
                                ann_id=args.annotation_subfield)
    else:
        cyvcf2_vcf.get_csq_keys(spliton="Format: ", delim="|",
                                ann_id=args.annotation_subfield)
    vcf_header_str = cyvcf2_vcf.header_to_list(gt_varnames=GT_VARNAMES,
                                               max_impact=args.max_impact,
                                               max_impact_csqs=args.max_impact_csqs,
                                               max_csq_scores=args.max_csq_scores,
                                               min_csq_scores=args.min_csq_scores,
                                               delim="\t")
    

    """
    since we're writing to a VCF, if any new INFO items written, need to 
    add to header to reflect this.
    """
    if args.max_impact_csqs != None:
        for csq_name in args.max_impact_csqs:
            csq_name_ext = csq_name + "_maximpact"
            vcf.add_info_to_header({'ID': csq_name_ext, 
                                    'Description':'max '+csq_name+' to go along '+\
                                                  'with transcripts with max IMPACT',
                                    'Type':'Character',
                                    'Number':'1'})
    if args.max_csq_scores != None:
        for csq_name in args.max_csq_scores:
            csq_name_ext = csq_name + "_max"
            vcf.add_info_to_header({'ID': csq_name_ext, 
                                    'Description':'max value for '+csq_name + \
                                                  'along assessed transcripts '+\
                                                  'in CSQ field.',
                                    'Type':'Float',
                                    'Number':'1'})
    if args.min_csq_scores != None:
        for csq_name in args.min_csq_scores:
            csq_name_ext = csq_name + "_min"
            vcf.add_info_to_header({'ID': csq_name_ext, 
                                    'Description':'min value for '+csq_name + \
                                                  'along assessed transcripts '+\
                                                  'in CSQ field.',
                                    'Type':'Float',
                                    'Number':'1'})



    """
    init VCF writer object
    """
    w = cyvcf2.Writer(args.out_vcf, vcf)
    # to write variant record, for v in vcf: w.write_record(v)

    """
    iterate through all variants, performing de novo screen on each one
    """
    vargeno_counts = defaultdict(int)
    prev_chrom = None
    n_var=0
    n_var_keep=0

    """
    if intervals provided, make sure to parse over those, else whole vcf
    """
    if args.intervals != None:
        if os.path.isfile(args.intervals):
            intervals = open(args.intervals, "r").readlines()
            intervals = [x.rstrip() for x in intervals]
        else:
            intervals = [args.intervals]
    else:
        intervals = [""]

    """
    parse VCF file looking for de novo variant calls
    """
    for vcf_variant in cyvcf2_vcf.iterator(intervals):
        n_var+=1
        #if linenum == 1000000: break

        """
        create new Cyvcf2Variant instance
        """
        cyvcf2_variant=Cyvcf2Variant(vcf_variant)

        if vcf_variant.CHROM != prev_chrom:
            print("Extracting variants from chrom " + vcf_variant.CHROM)
            prev_chrom = vcf_variant.CHROM

        """
        assume single allele per site, exclude sites with call as '*'
        """
        alt = vcf_variant.ALT[0]
        if alt == '*': continue

        ## if no qualifying impact str found in CSQ, skip
        if args.qual_impacts != None:
            res = cyvcf2_variant.qual_impacts_screen(args.qual_impacts,
                                                     csq_subfield=args.annotation_subfield)
            if res == False: continue

        ## if desired, derive max impact annots from var, along with other
        ## user defined max or min scores in CSQ for variant
        csqs_maximpact_list = []
        max_csq_scores = []
        min_csq_scores = []
        if args.max_impact == True:
            cyvcf2_variant.get_annot_txs(cyvcf2_vcf.csq_keys, 
                                         csq_subfield=args.annotation_subfield)
            if args.annotation_subfield == "ANN":
                impact_subfield="Annotation_Impact"
            else:
                impact_subfield="IMPACT"
            res=cyvcf2_variant.maxmin_csqs(csq_subfield=args.annotation_subfield, 
                                           impact_subfield=impact_subfield,
                                           max_impact_csqs=args.max_impact_csqs,
                                           max_csq_scores=args.max_csq_scores, 
                                           min_csq_scores=args.min_csq_scores)
            (csqs_maximpact_list, max_csq_scores, min_csq_scores) = res

            """
            if corresponding values defined, add to vcf record
            """
            if args.max_impact_csqs!=None:
                for i in range(len(args.max_impact_csqs)):
                    max_impact_csq_name=args.max_impact_csqs[i] + "_maximpact"
                    max_impact_csq=csqs_maximpact_list[i]
                    vcf_variant.INFO[max_impact_csq_name] = max_impact_csq
            if args.min_csq_scores!=None:
                for i in range(len(args.min_csq_scores)):
                    min_csq_score_name=args.min_csq_scores[i] + "_min"
                    min_csq_score=float(min_csq_scores[i])
                    vcf_variant.INFO[min_csq_score_name] = min_csq_score
            if args.max_csq_scores!=None:
                for i in range(len(args.max_csq_scores)):
                    max_csq_score_name=args.max_csq_scores[i] + "_max"
                    max_csq_score=float(max_csq_scores[i])
                    vcf_variant.INFO[max_csq_score_name] = max_csq_score




        ## filter on variant cnds file provided
        if var_cnds.test_variant(vcf_variant) == False: continue

        ## if variant survives filters, retain record
        w.write_record(vcf_variant)
        n_var_keep += 1

    w.close()
    vcf.close()

    ## print basic stats on number of input variants, number of 
    ## variants to keep
    print("Number of variants in parent VCF : " + str(n_var))
    print("Number of variants retained post-filtration : " + str(n_var_keep))

    return

def parse_args(ARGS):
    args = argparse.ArgumentParser()
    import vcfscreen.arguments as a
    args = a.args_vars(args)
    args = a.args_annots(args)
    args.add_argument("--annotation-subfield", type=str, default="CSQ",
                      help="annotation subfield string to use")
    args.add_argument("variant_cnds",
                      help="input vcf file")
    args.add_argument("in_vcf",                                                 
                      help="input vcf file") 
    args.add_argument("out_vcf",
                      help="output vcf file.")
    return args.parse_args(args=ARGS)

if __name__ == "__main__":
    main()
