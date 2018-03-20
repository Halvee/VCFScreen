#!/usr/bin/env python

'''
Filename : VCFscreen_denovo.py 
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
    read samples from fam file, send basic stats to stdout
    """
    samples_i = Samples(args.in_fam)
    samples_i.print_stats()
    n_samples = len(samples_i.samples)
    n_males = len(samples_i.males)
    n_females = len(samples_i.females)
    n_cases = len(samples_i.cases)
    n_ctrls = len(samples_i.ctrls)

    """
    get case, control idxs
    """
    case_idxs = [samples_i.samples[x].idx for x in samples_i.cases]
    ctrl_idxs = [samples_i.samples[x].idx for x in samples_i.ctrls]

    """
    read cnds files
    """
    var_cnds = None
    pro_cnds = None
    par_cnds = None
    if args.variant_cnds != None: var_cnds = VcfCnds(args.variant_cnds)
    if args.pro_cnds != None: pro_cnds = VcfCnds(args.pro_cnds)
    if args.par_cnds != None: par_cnds = VcfCnds(args.par_cnds)

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
    iterate through all variants, performing de novo screen on each one
    """
    dnm_counts = defaultdict(int)
    prev_chrom = None
    linenum=0

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
    init output file                                                         
    """                                                                         
    misc.init_out_file(args.out_tsv,                                                 
                       init_line = vcf_header_str)  

    """
    parse VCF file looking for de novo variant calls
    """
    for vcf_variant in cyvcf2_vcf.iterator(intervals):
        linenum+=1
        #if linenum == 1000000: break
        if vcf_variant.CHROM != prev_chrom:
            print("Extracting DNMs from chrom " + vcf_variant.CHROM)
            prev_chrom = vcf_variant.CHROM

        """
        create new Cyvcf2Variant instance 
        """
        cyvcf2_variant=Cyvcf2Variant(vcf_variant)

        """
        assume single allele per site, exclude sites with call as '*'
        """
        alt = cyvcf2_variant.cyvcf2_variant.ALT[0]
        if alt == '*': continue

        ## if no qualifying impact str found in CSQ, skip
        if args.qual_impacts != None:
            res = cyvcf2_variant.qual_impacts_screen(args.qual_impacts,
                                                     csq_subfield="CSQ")
            if res == False: continue

        ## if desired, derive max impact annots from var, along with other
        ## user defined max or min scores in CSQ for variant
        csqs_maximpact_list = []
        max_csq_scores = []
        min_csq_scores = []
        if args.max_impact == True:
            cyvcf2_variant.get_annot_txs(cyvcf2_vcf.csq_keys,
                                         csq_subfield="CSQ")
            res=cyvcf2_variant.maxmin_csqs(max_impact_csqs=args.max_impact_csqs,
                                           max_csq_scores=args.max_csq_scores,
                                           min_csq_scores=args.min_csq_scores,
                                           impact_subfield="IMPACT") 
            (csqs_maximpact_list, max_csq_scores, min_csq_scores) = res

        ## filter on internal ctrl-only maf, can't do in var cnds file
        if args.internal_ctrl_af_max != None:
            ctrl_af = cyvcf_variant.compute_maf(ctrl_idxs)
            if ctrl_af >= args.internal_ctrl_af_max:
                continue

        ## variant cnds file provided, filter exclusively on that
        if var_cnds != None:
            if var_cnds.test_variant(vcf_variant) == False: continue
        ## otherwise, filter on user args
        else:

            ## filter on FILTER column, default is PASS only, allow for 
            ## user defined FILTER classifs too
            if args.filter_include == None and vcf_variant.FILTER != None:
                continue
            elif args.filter_include != None and vcf_variant.FILTER != None:
                filter_include = set(",".split(args.filter_include))
                if vcf_variant not in filter_include: continue
           
            ## filter on user-defined VCF INFO flags
            if args.vcf_info_flags_exclude != None:
                vcf_info_flags_fail = False
                for vcf_info_flag in args.vcf_info_flags_exclude.split(","):
                    if vcf_variant.INFO.get(vcf_info_flag) == True:
                        vcf_info_flags_fail = True
                        break
                if vcf_info_flags_fail == True: continue

            ## filter on user-defined maximum internal MAF
            if args.internal_af_max != None:
                if vcf_variant.INFO.get("AF") > args.internal_af_max:
                    continue

            ## filter on user-defined maximum external MAF
            if args.af_max_fields != None and args.af_max != None:
                af_max_fields = args.af_max_fields.split(",")
                af_max_fail = False
                for af_max_field in af_max_fields:
                    af = vcf_variant.INFO.get(af_max_field)
                    if af > args.af_max: 
                        af_max_fail=True
                        break
                if af_max_fail == True: continue

        dnm_carriers = set()
        if pro_cnds != None and par_cnds != None:
            for iid in samples_i.trios:
                pid = samples_i.trios[iid].pid
                mid = samples_i.trios[iid].mid
                iid_idx=samples_i.samples[iid].idx
                pid_idx=samples_i.samples[pid].idx
                mid_idx=samples_i.samples[mid].idx
                trio_idxs=(iid_idx,pid_idx,mid_idx)

                ## is iid het or hom alt at site?
                trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs])
                if trio_gts[0] not in set([1,2]): continue

                ## are parent samples hom ref at site?
                if trio_gts[1] not in set([0]) or trio_gts[2] not in set([0]): 
                    continue
                
                ## test if proband and parents pass conditions in proband and parent
                ## cnds files, respectively
                if pro_cnds.test_gt(vcf_variant, iid_idx) == False: continue        
                if par_cnds.test_gt(vcf_variant, pid_idx) == False: continue
                if par_cnds.test_gt(vcf_variant, mid_idx) == False: continue
                dnm_carriers.add(iid)
                
        else:
            denovo_screen = screens.denovo_screen
            dnm_carriers = denovo_screen(vcf_variant, samples_i,
                                         min_coverage = args.min_coverage,
                                         pro_min_perc_alt=args.pro_min_perc_alt,
                                         par_max_perc_alt=args.par_max_perc_alt,
                                         pro_het_phredmax=args.pro_het_phredmax,
                                         par_hom_phredmax=args.par_hom_phredmax,
                                         pro_hom_phredmin=args.pro_hom_phredmin,
                                         par_het_phredmin=args.par_het_phredmin)

        if len(dnm_carriers) > 0:
            iids = list(dnm_carriers)
            for iid in iids: 
                samples_i.samples[iid].varcounts["dnm"] += 1
            outs = cyvcf2_variant.variant_to_list(samples_i,
                                                  dnm_carriers,
                                                  cyvcf2_vcf.info_subfields, 
                                                  GT_VARNAMES, 
                                                  csqs_maximpact=csqs_maximpact_list,
                                                  max_csq_scores=max_csq_scores,
                                                  min_csq_scores=min_csq_scores,
                                                  delim="\t")
            for out in outs:
                misc.append_out_file(args.out_tsv, out)

    vcf.close()
    samples_i.print_varcount_stats(var_types=["dnm"])
    return

def parse_args(ARGS):
    args = argparse.ArgumentParser()
    import vcfscreen.arguments as a
    args = a.args_vars(args)
    args = a.args_annots(args)
    args = a.args_phredthresh(args,
                              sampletypes=["pro","par"],
                              gt_types=["hom","het"],
                              threshes=["max","min"])
    args = a.args_perc_alt(args,
                           sampletypes=["pro","par"],
                           threshes=["max","min"])
    args = a.args_cnds(args, variant_cnds=True,
                        sampletypes=["pro","par"])
    args = a.args_req(args)
    return args.parse_args(args=ARGS)

if __name__ == "__main__":
    main()
