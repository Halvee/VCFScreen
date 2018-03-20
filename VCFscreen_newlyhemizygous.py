#!/usr/bin/env python

'''
Filename : VCFscreen_newlyhemizygous.py 
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
import vcfscreen.misc as misc
import vcfscreen.samples as samples
import vcfscreen.screens as screens
from vcfscreen.screens import variant_screen,hemi_screen
from vcfscreen.vcf_cnds import VcfCnds
from vcfscreen.misc import init_out_file,append_out_file,str_none_split
from vcfscreen.cyvcf2_variant import Cyvcf2Vcf,Cyvcf2Variant

X_CHROM_INTERVAL="X:1-155270560"
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
    args.qual_impacts = str_none_split(args.qual_impacts, ",")
    args.max_impact_csqs = str_none_split(args.max_impact_csqs, ",")
    args.max_csq_scores = str_none_split(args.max_csq_scores, ",")
    args.min_csq_scores = str_none_split(args.min_csq_scores, ",")

    """
    convert all comma delim args to list, or leave as none.
    """
    args_str_none=("filter_include","af_max_fields","vcf_info_flags_exclude")
    for arg in args_str_none:
        args.__dict__[arg] = misc.str_none_split(args.__dict__[arg], ",")

    """
    read samples from fam file, send basic stats to stdout
    """
    samples_i = samples.Samples(args.in_fam)
    samples_i.print_stats()

    """
    get all trios from samples with male proband
    """
    male_trios = {}
    for iid in samples_i.trios:
        if samples_i.samples[iid].gender == "M":
            male_trios[iid] = samples_i.samples[iid]
    if len(male_trios) == 0:
        print("ERROR" + \
              "Can't do hemizygous screen if no " + \
              "trios with male proband in cohort.")
        sys.exit(1)

    """
    read cnds files
    """
    var_cnds = None
    iid_cnds = None
    pid_cnds = None
    mid_cnds = None
    if args.variant_cnds != None: var_cnds = VcfCnds(args.variant_cnds)
    if args.iid_cnds != None: iid_cnds = VcfCnds(args.iid_cnds)
    if args.pid_cnds != None: pid_cnds = VcfCnds(args.pid_cnds)
    if args.mid_cnds != None: mid_cnds = VcfCnds(args.mid_cnds)

    """
    init cyvcf2 VCF obj, get info subfields, header for output
    """
    vcf = cyvcf2.VCF(args.in_vcf, strict_gt=True)
    cyvcf2_vcf = Cyvcf2Vcf(vcf)
    cyvcf2_vcf.get_info_subfields()
    cyvcf2_vcf.get_csq_keys(spliton="Format: ", delim="|")
    vcf_header_str = cyvcf2_vcf.header_to_list(gt_varnames=GT_VARNAMES,
                                               max_impact=args.max_impact,
                                               max_impact_csqs=args.max_impact_csqs,
                                               max_csq_scores=args.max_csq_scores,
                                               min_csq_scores=args.min_csq_scores,
                                               delim="\t")

    """
    create sample idx
    """
    samples_i.get_vcf_idx(vcf.samples)

    """
    iterate through all variants, performing newlyhemizygous screen on each one
    """
    hemi_counts = defaultdict(int)
    prev_chrom = None
    linenum=0

    """
    init output file
    """
    init_out_file(args.out_tsv, 
                  init_line = vcf_header_str + "\n")

    """
    only parse X chromosome, since this is only place where hemi vars can happen
    """
    for vcf_variant in vcf(args.x_chrom_interval):
        """
        assume single allele per site, exclude sites with call as '*'
        """
        alt = vcf_variant.ALT[0]
        if alt == '*': continue

        ## if no qualifying impact str found in CSQ, skip
        if args.qual_impacts != None:
            res = screens.qual_impacts_screen(vcf_variant, 
                                              args.qual_impacts,
                                              csq_subfield="CSQ")
            if res == False: continue

        ## if desired, derive max impact annots from var, along with other
        ## user defined max or min scores in CSQ for variant
        csqs_maximpact_list = []
        max_csq_scores = []
        min_csq_scores = []
        if args.max_impact == True:
            res = screens.get_maxmin_csqs(vcf_variant,
                                          csq_keys,
                                          max_impact_csqs=args.max_impact_csqs,
                                          max_csq_scores=args.max_csq_scores, 
                                          min_csq_scores=args.min_csq_scores,
                                          csq_subfield="CSQ",
                                          impact_subfield="IMPACT")
            (csqs_maximpact_list, max_csq_scores, min_csq_scores) = res

        ## variant cnds file provided, filter exclusively on that
        if var_cnds != None:
            if var_cnds.test_variant(vcf_variant) == False: continue
        ## otherwise, filter on user args
        else:
            var_pass = screens.variant_screen(vcf_variant, 
                                              min_qual=args.min_qual,
                                              filter_include=args.filter_include,
                                              vcf_info_flags_exclude=args.vcf_info_flags_exclude,
                                              internal_af_max=args.internal_af_max,
                                              af_max=args.af_max,
                                              af_max_fields=args.af_max_fields)
            if var_pass == False: continue

        if iid_cnds != None and pid_cnds != None and mid_cnds != None: 
            for iid in male_trios:
                iid_idx=samples_i.samples[iid].idx
                pid_idx=samples_i.samples[pid].idx
                mid_idx=samples_i.samples[mid].idx
                trio_idxs=(iid_idx,mid_idx)

                ## is iid hemi at site?
                trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs])
                if trio_gts[0] not in set([1,2]): continue

                ## is father hom ref at site?
                if trio_gts[1] not in set([0]): continue

                ## is mother het at site?
                if trio_gts[2] not in set([1]): continue

                ## test if proband and parents pass conditions in proband and parent
                ## cnds files, respectively
                if iid_cnds.test_gt(vcf_variant, iid_idx) == False: continue        
                if pid_cnds.test_gt(vcf_variant, pid_idx) == False: continue
                if mid_cnds.test_gt(vcf_variant, mid_idx) == False: continue
                hemi_carriers.add(iid)
                
        else:
            hemi_screen = screens.hemi_screen
            hemi_carriers = hemi_screen(vcf_variant, samples_i, male_trios,
                                        min_coverage = args.min_coverage,
                                        iid_min_perc_alt=args.iid_min_perc_alt,
                                        pid_max_perc_alt=args.pid_max_perc_alt,
                                        mid_min_perc_alt=args.mid_min_perc_alt,
                                        iid_het_phredmax=args.iid_het_phredmax,
                                        pid_hom_phredmax=args.pid_hom_phredmax,
                                        mid_het_phredmax=args.mid_het_phredmax,
                                        iid_hom_phredmin=args.iid_hom_phredmin,
                                        pid_het_phredmin=args.pid_het_phredmin,
                                        mid_hom_phredmin=args.mid_hom_phredmin)
        
        if len(hemi_carriers) > 0:
            iids = list(hemi_carriers)
            for iid in iids:
                samples_i.samples[iid].varcounts["newlyhemi"] += 1
            outs = cyvcf2_variant.variant_to_list(vcf_variant,
                                                  samples_i,
                                                  hemi_carriers,
                                                  info_subfields,
                                                  GT_VARNAMES,
                                                  csqs_maximpact=csqs_maximpact_list,
                                                  max_csq_scores=max_csq_scores,
                                                  min_csq_scores=min_csq_scores,
                                                  delim="\t")
            for out in outs:
                append_out_file(args.out_tsv, out)

    vcf.close()
    print("Screening of VCF for newly hemizygous variants complete.")
    samples_i.print_varcount_stats(var_types=["newlyhemi"])
    return

def parse_args(ARGS):
    args = argparse.ArgumentParser()
    import vcfscreen.arguments as a
    args = a.args_vars(args)
    args = a.args_annots(args)    
    args = a.args_phredthresh(args,
                              sampletypes=["iid","pid","mid"],
                              gt_types=["hom","het"],
                              threshes=["max","min"])
    args = a.args_perc_alt(args,
                           sampletypes=["iid","pid","mid"],
                           threshes=["max","min"])
    args = a.args_cnds(args, variant_cnds=True,
                        sampletypes=["iid","pid","mid"])
    args.add_argument("--x-chrom-interval", type=str,
                      action="store", default="X:1-155270560",
                      help="span of X chromosome to screen.")
    args = a.args_req(args)
    return args.parse_args(args=ARGS)

if __name__ == "__main__":
    main()
