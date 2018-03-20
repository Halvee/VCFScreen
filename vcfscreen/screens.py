from collections import defaultdict
from annot import AnnotTxs

def variant_screen(vcf_variant, 
                   min_qual=0,
                   filter_include=None,
                   vcf_info_flags_exclude=None,
                   internal_af_max=None,
                   af_max=None,
                   af_max_fields=None):

    ## filter on FILTER field, default is PASS only, allow for 
    ## user defined FILTER classifs too
    if filter_include == None and vcf_variant.FILTER != None:
        return False
    elif filter_include != None and vcf_variant.FILTER != None:
        if vcf_variant not in filter_include: return False
   
    ## filter on QUAL field
    if vcf_variant.QUAL < min_qual: return False

    ## filter on user-defined VCF INFO flags
    if vcf_info_flags_exclude != None:
        vcf_info_flags_fail = False
        for vcf_info_flag in vcf_info_flags_exclude:
            if vcf_variant.INFO.get(vcf_info_flag) == True:
                vcf_info_flags_fail = True
                break
        if vcf_info_flags_fail == True: return False

    ## filter on user-defined maximum internal MAF
    if internal_af_max != None:
        if vcf_variant.INFO.get("AF") > internal_af_max:
            return False

    ## filter on user-defined maximum external MAF
    if af_max_fields != None and af_max != None:
        af_max_fail = False
        for af_max_field in af_max_fields:
            af = vcf_variant.INFO.get(af_max_field)
            if af > af_max: 
                af_max_fail=True
                break
        if af_max_fail == True: return False

    return True

def vargeno_screen(cyvcf2_variant_i, samples,           
                   qual_gts=set([1,2]),
                   min_coverage = 10,
                   min_perc_alt=0.2,
                   max_perc_alt=0.8,
                   homref_phredmin=None,
                   het_phredmax=None,   
                   het_phredmin=None,
                   homalt_phredmax=None,
                   homalt_phredmin=None): 
    vargeno_carriers=set()
    alt=cyvcf2_variant_i.ALT[0]
    if alt == '*': return vargeno_carriers
    variant_id="-".join([cyvcf2_variant_i.CHROM, 
                         str(cyvcf2_variant_i.POS),
                         cyvcf2_variant_i.REF, alt])
    for iid in samples.samples:
        iid_idx = samples.samples[iid].idx

        ## does sample have qualifying genotype?
        gt = cyvcf2_variant_i.gt_types[iid_idx]
        if gt not in qual_gts: continue

        ## does sample have >= min coverage at site?
        gt_cov = cyvcf2_variant_i.gt_depths[iid_idx]
        if gt_cov < min_coverage: continue

        ## if phred likelihood thresholds are provided, apply
        if homref_phredmin != None:
            if vcf_variant.gt_phred_ll_homref[iid_idx] < homref_phredmin:
                continue
        if het_phredmin != None:
            if vcf_variant.gt_phred_ll_het[iid_idx] < het_phredmin:
                continue
       	if het_phredmax != None:
            if vcf_variant.gt_phred_ll_het[iid_idx] > het_phredmax:
       	       	continue
       	if homalt_phredmin != None:
            if vcf_variant.gt_phred_ll_homalt[iid_idx] < homalt_phredmin:
       	       	continue
        if homalt_phredmax != None:
            if vcf_variant.gt_phred_ll_homalt[iid_idx] > homalt_phredmax:
                continue

        ## if min or max perc alt defined, apply
        alt_freq = cyvcf2_variant_i.gt_alt_freqs[iid_idx]
        if min_perc_alt != None: 
            if alt_freq < min_perc_alt: continue
        if max_perc_alt != None:
            if alt_freq > max_perc_alt: continue

        ## variant non-homref genotype detected
        vargeno_carriers.add(iid)

    return vargeno_carriers

def denovo_screen(vcf_variant, samples, iid_gts=set([1,2]),
                  par_gts=set([0]), min_coverage=10,
                  pro_min_perc_alt=None, par_max_perc_alt=None,
                  pro_het_phredmax=None, par_hom_phredmax=None,
                  pro_hom_phredmin=None, par_het_phredmin=None):
    dnm_carriers=set()
    alt=vcf_variant.ALT[0]
    if alt == '*': return dnm_carriers
    variant_id="-".join([vcf_variant.CHROM, str(vcf_variant.POS),
                         vcf_variant.REF, alt])
    for iid in samples.trios:
        pid = samples.trios[iid].pid
        mid = samples.trios[iid].mid
        iid_idx=samples.samples[iid].idx
        pid_idx=samples.samples[pid].idx
        mid_idx=samples.samples[mid].idx
        trio_idxs=(iid_idx,pid_idx,mid_idx)
    
        ## is iid het or hom alt at site?
        trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs])
        if trio_gts[0] not in iid_gts: continue

        ## are parent samples hom ref at site?
        if trio_gts[1] not in par_gts or trio_gts[2] not in par_gts: 
            continue

        ## do genotypes pass phred likelihood thresholds if defined?
        if pro_hom_phredmin != None:
            if vcf_variant.gt_phred_ll_homref[iid_idx] < pro_hom_phredmin:
                continue
        if pro_het_phredmax != None:
            if vcf_variant.gt_phred_ll_het[iid_idx] > pro_het_phredmax:
                continue
        if par_het_phredmin != None:
            het = min(vcf_variant.gt_phred_ll_het[pid_idx], 
                      vcf_variant.gt_phred_ll_het[mid_idx])
            homalt = min(vcf_variant.gt_phred_ll_homalt[pid_idx],
                         vcf_variant.gt_phred_ll_homalt[mid_idx])
            if het < par_het_phredmin or homalt < par_het_phredmin:
                continue
        if par_hom_phredmax != None:
            homref = max(vcf_variant.gt_phred_ll_homref[pid_idx], 
                         vcf_variant.gt_phred_ll_homref[pid_idx])
            if homref > par_hom_phredmax:
                continue

        ## only keep calls where min(iid_cov, pid_cov, mid_cov) > min_coverage:
        trio_covs = ([vcf_variant.gt_depths[id_x] for id_x in trio_idxs])
        if min(trio_covs) < min_coverage: continue

        ## only keep instances where var > pro_min_perc_alt in proband
        if pro_min_perc_alt != None:
            iid_gt_alt_depth = float(vcf_variant.gt_alt_depths[iid_idx])
            iid_perc_alt = iid_gt_alt_depth / trio_covs[0]
            if iid_perc_alt < pro_min_perc_alt: continue

        ## only keep instances where var < par_max_perc_alt in both parents
        if par_max_perc_alt != None:
            pid_gt_alt_depth = float(vcf_variant.gt_alt_depths[pid_idx])
            mid_gt_alt_depth = float(vcf_variant.gt_alt_depths[mid_idx])
            pid_perc_alt = pid_gt_alt_depth / trio_covs[1]
            mid_perc_alt = mid_gt_alt_depth / trio_covs[2]
            if max(pid_perc_alt, mid_perc_alt) > par_max_perc_alt: continue

        ## de novo mutation detected
        dnm_carriers.add(iid)

    return dnm_carriers

def mosaic_screen(vcf_variant, samples, iid_gts=set([1,2]),
                  transpar_gts=set([0,1]), ntranspar_gts=set([0]),
                  min_coverage=10,
                  pro_min_perc_alt=None, 
                  transpar_min_perc_alt=None, transpar_max_perc_alt=None,
                  ntranspar_max_perc_alt=None,
                  pro_homref_phredmin=None, pro_het_phredmax=None, 
                  pro_homalt_phredmin=None):
    mosaic_carriers=set()
    alt=vcf_variant.ALT[0]
    if alt == '*': return mosaic_carriers
    variant_id="-".join([vcf_variant.CHROM, str(vcf_variant.POS),
                         vcf_variant.REF, alt])
    for iid in samples.trios:
        pid = samples.trios[iid].pid
        mid = samples.trios[iid].mid
        iid_idx=samples.samples[iid].idx
        pid_idx=samples.samples[pid].idx
        mid_idx=samples.samples[mid].idx
        trio_idxs=(iid_idx,pid_idx,mid_idx)
    
        ## is iid het or hom alt at site?
        trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs])
        if trio_gts[0] not in iid_gts: continue

        ## identify more likely transmitting parent based on high
        ## percent alt reads at site
        pro_alt_freq = vcf_variant.gt_alt_freqs[iid_idx]
        pid_alt_freq = vcf_variant.gt_alt_freqs[pid_idx]                
        mid_alt_freq = vcf_variant.gt_alt_freqs[mid_idx]                
        if pid_alt_freq > mid_alt_freq:                                 
            transpar = pid                                              
            transpar_idx = pid_idx                                      
            transpar_gt = trio_gts[1]
            transpar_alt_freq = pid_alt_freq                                   
            ntranspar = mid                                             
            ntranspar_idx = mid_idx                                     
            ntranspar_gt = trio_gts[2]      
            ntranspar_alt_freq = mid_alt_freq                            
        else:                                                           
            transpar = mid                                              
            transpar_idx = mid_idx                                      
            transpar_gt = trio_gts[2]        
            transpar_alt_freq = mid_alt_freq                           
            ntranspar = pid                                             
            ntranspar_idx = pid_idx                                     
            ntranspar_gt = trio_gts[1]
            ntranspar_alt_freq = pid_alt_freq

        ## is more likely trasmitting parent homref or het at site, and is 
        ## less likely transmitting parent hom at site? If no then skip
        if transpar_gt not in transpar_gts or ntranspar_gt not in ntranspar_gts: 
            continue

        ## is proband var call found in greater than min percentage of reads
        if pro_alt_freq < pro_min_perc_alt: continue

        ## is var call between min and max allowed perc of reads in transmitting par?
        if  transpar_alt_freq < transpar_min_perc_alt: continue
        elif transpar_alt_freq > transpar_max_perc_alt: continue

        ## is var call less than or equal to max allowed perc of reads in ntrans par?
        if ntranspar_alt_freq < transpar_max_perc_alt: continue 
        
        ## subset on phred scores for proband genotype
        if vcf_variant.gt_phred_ll_homref[iid_idx] < pro_homref_phredmin:
            continue
       	if vcf_variant.gt_phred_ll_het[iid_idx] > pro_het_phredmax:
       	    continue
       	if vcf_variant.gt_phred_ll_homalt[iid_idx] < pro_homalt_phredmin:
       	    continue

        ## transmitted parental-mosaic variant detected
        mosaic_carriers.add(iid)

    return mosaic_carriers

def trans_screen(vcf_variant, samples, ntrans=False,
                 transpar_gts=set([0,1]), ntranspar_gts=set([0]),
                 min_coverage=10,
                 pro_min_perc_alt=None, pro_max_perc_alt=None,
                 transpar_min_perc_alt=0.2, transpar_max_perc_alt=0.8,
                 ntranspar_max_perc_alt=0.05,
                 pro_homref_phredmin=None, pro_het_phredmax=None, 
                 pro_homalt_phredmin=None):
    trans_carriers=set()
    alt=vcf_variant.ALT[0]
    if alt == '*': return mosaic_carriers
    variant_id="-".join([vcf_variant.CHROM, str(vcf_variant.POS),
                         vcf_variant.REF, alt])
    for iid in samples.trios:
        pid = samples.trios[iid].pid
        mid = samples.trios[iid].mid
        iid_idx=samples.samples[iid].idx
        pid_idx=samples.samples[pid].idx
        mid_idx=samples.samples[mid].idx
        trio_idxs=(iid_idx,pid_idx,mid_idx)
    
        ## get trio gts
        trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs]) 

        ## identify more likely transmitting parent based on high
        ## percent alt reads at site
        pro_alt_freq = vcf_variant.gt_alt_freqs[iid_idx]
        pid_alt_freq = vcf_variant.gt_alt_freqs[pid_idx]                
        mid_alt_freq = vcf_variant.gt_alt_freqs[mid_idx]                
        if pid_alt_freq > mid_alt_freq:                                 
            transpar = pid                                              
            transpar_idx = pid_idx                                      
            transpar_gt = trio_gts[1]
            transpar_alt_freq = pid_alt_freq                                   
            ntranspar = mid                                             
            ntranspar_idx = mid_idx                                     
            ntranspar_gt = trio_gts[2]      
            ntranspar_alt_freq = mid_alt_freq                            
        else:                                                           
            transpar = mid                                              
            transpar_idx = mid_idx                                      
            transpar_gt = trio_gts[2]        
            transpar_alt_freq = mid_alt_freq                           
            ntranspar = pid                                             
            ntranspar_idx = pid_idx                                     
            ntranspar_gt = trio_gts[1]
            ntranspar_alt_freq = pid_alt_freq

        ## is more likely trasmitting parent homref or het at site, and is 
        ## less likely transmitting parent hom at site? If no then skip
        if transpar_gt not in transpar_gts or ntranspar_gt not in ntranspar_gts: 
            continue
        
        ## is var call between min and max allowed perc of reads in transmitting par?
        if  transpar_alt_freq < transpar_min_perc_alt: continue
        elif transpar_alt_freq > transpar_max_perc_alt: continue

        ## is var call less than or equal to max allowed perc of reads in ntrans par?
        if ntranspar_alt_freq < transpar_max_perc_alt: continue 

        ## if defined, keep vars het in parent with high phred likelihood 
        if transpar_homref_phredmin != None:
            phred = vcf_variant.gt_phred_ll_homref[transpar_idx]
            if phred < transpar_homref_phredmin: continue
        if transpar_het_phredmax != None:
            phred = vcf_variant.gt_phred_ll_het[transpar_idx]
            if phred < transpar_het_phredmax: continue
        if transpar_homalt_phredmin != None:
            phred = vcf_variant.gt_phred_ll_homalt[transpar_idx]
            if phred < transpar_homalt_phredmin: continue

        ## if ntrans is false, then looking for transmitted het vars, 
        ## and proband must be '1', with perc_alt_min <= PERC_ALT < perc_alt_max. 
        ## if ntrans is true, looking for vars where proband gt is '0',
        ## with PERC_ALT < perc_alt_max.
        if ntrans == False:
            if trio_gts[0] not in set([1]): continue
            if pro_min_perc_alt !=None: 
                if pro_alt_freq < pro_min_perc_alt: continue
            if pro_max_perc_alt != None:
                if pro_alt_freq > pro_max_perc_alt: continue
        else:
            if trio_gts[0] not in set([0]): continue
            if pro_max_perc_alt != None:
                if pro_alt_freq > pro_max_perc_alt: continue

        ## parental transmitted (or nontransmitted) variant detected
        trans_carriers.add((transpar, iid))

    return trans_carriers


def hemi_screen(vcf_variant, samples, male_trios, iid_gts=set([1]),
                pid_gts=set([0]), mid_gts=set([1]), 
                min_coverage=10,
                iid_min_perc_alt=None, pid_max_perc_alt=None,
                mid_min_perc_alt=None,
                iid_het_phredmax=None, mid_het_phredmax=None,
                pid_hom_phredmax=None, pid_het_phredmin=None,
                iid_hom_phredmin=None, mid_hom_phredmin=None):
    hemi_carriers = set()
    alt = vcf_variant.ALT[0]
    for iid in male_trios:
        pid = samples.samples[iid].pid
        mid = samples.samples[iid].mid
        iid_idx=samples.samples[iid].idx
        pid_idx=samples.samples[pid].idx
        mid_idx=samples.samples[mid].idx
        trio_idxs=(iid_idx,pid_idx,mid_idx)

        ## get trio genotypes
        trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs])

        ## skip if proband not het at site
        if trio_gts[0] not in iid_gts: continue

        ## skip if father not hom ref at site
        if trio_gts[1] not in pid_gts: continue

        ## skip if mother not het at site
        if trio_gts[2] not in mid_gts: continue

        ## only keep calls where min(iid_cov, pid_cov, mid_cov) > min_coverage:
        trio_covs = ([vcf_variant.gt_depths[id_x] for id_x in trio_idxs])
        if min(trio_covs) < min_coverage: continue

        ## only keep instances where var > min_perc_alt in proband
        if iid_min_perc_alt != None:
            iid_gt_alt_depth = float(vcf_variant.gt_alt_depths[iid_idx])
            iid_perc_alt = iid_gt_alt_depth / trio_covs[0]
            if iid_perc_alt < pro_min_perc_alt: continue
        
        ## only keep instances where var > min_perc_alt in father
        if pid_max_perc_alt != None:
            pid_gt_alt_depth = float(vcf_variant.gt_alt_depths[pid_idx])
            pid_perc_alt = pid_gt_alt_depth / trio_covs[1]
            if pid_perc_alt > pid_max_perc_alt: continue
       
        ## only keep instances where var > min_perc_alt in mother
        if mid_min_perc_alt != None:
            mid_gt_alt_depth = float(vcf_variant.gt_alt_depths[mid_idx])
            mid_perc_alt = mid_gt_alt_depth / trio_covs[2]
            if mid_perc_alt < mid_min_perc_alt: continue

        ## do genotypes pass phred likelihood thresholds if defined?
        if iid_hom_phredmin != None:
            if vcf_variant.gt_phred_ll_homref[iid_idx] < pro_hom_phredmin:
                continue
        if iid_het_phredmax != None:
            if vcf_variant.gt_phred_ll_het[iid_idx] > pro_het_phredmax:
                continue
        
        ## do genotypes pass phred likelihood threshes for pid if def?
        if pid_hom_phredmax != None:
            if vcf_variant.gt_phred_ll_homref[pid_idx] > pid_hom_phredmax:
                continue
        if pid_het_phredmin != None:
            if vcf_variant.gt_phred_ll_het[pid_idx] < pid_hom_phredmin:
                continue

        ## do genotypes pass phred likelihood threshes for pid if def?
        if mid_hom_phredmin != None:
            if vcf_variant.gt_phred_ll_homref[mid_idx] < mid_hom_phredmin:
                continue
        if mid_het_phredmax != None:
            if vcf_variant.gt_phred_ll_homref[mid_idx] > mid_hom_phredmax:
                continue

        ## newly hemizygous variant detected
        hemi_carriers.add(iid)

    return hemi_carriers

def hom_screen(vcf_variant, samples, 
               pro_gts=set([2]), par_gts=set([1]),
               min_coverage=10,
               pro_min_perc_alt=None,
               par_min_perc_alt=None, par_max_perc_alt=None,
               pro_homref_phredmin=None, pro_het_phredmin=None,
               pro_homalt_phredmax=None, par_homref_phredmin=None,
               par_het_phredmax=None, par_homalt_phredmin=None):
    hom_carriers = set()
    alt = vcf_variant.ALT[0]
    for iid in samples.trios:
        pid = samples.samples[iid].pid
        mid = samples.samples[iid].mid
        iid_idx=samples.samples[iid].idx
        pid_idx=samples.samples[pid].idx
        mid_idx=samples.samples[mid].idx
        trio_idxs=(iid_idx,pid_idx,mid_idx)

        ## get trio genotypes
        trio_gts = ([vcf_variant.gt_types[id_x] for id_x in trio_idxs])

        ## skip if proband not hom alt at site
        if trio_gts[0] not in pro_gts: continue

        ## skip if father not het at site
        if trio_gts[1] not in par_gts: continue

        ## skip if mother not het at site
        if trio_gts[2] not in par_gts: continue
        
        ## only keep calls where min(iid_cov, pid_cov, mid_cov) > min_coverage:
        trio_covs = ([vcf_variant.gt_depths[id_x] for id_x in trio_idxs])
        if min(trio_covs) < min_coverage: continue

        ## only keep instances where var > min_perc_alt in proband
        if pro_min_perc_alt != None:
            iid_gt_alt_depth = float(vcf_variant.gt_alt_depths[iid_idx])
            iid_perc_alt = iid_gt_alt_depth / trio_covs[0]
            if iid_perc_alt < pro_min_perc_alt: continue

        ## only keep instances where var < max_perc_alt in parents
        if par_max_perc_alt != None:
            pid_gt_alt_depth = float(vcf_variant.gt_alt_depths[pid_idx])
            mid_gt_alt_depth = float(vcf_variant.gt_alt_depths[mid_idx])
            pid_perc_alt = pid_gt_alt_depth / trio_covs[1]
            mid_perc_alt = pid_gt_alt_depth / trio_covs[2]
            if pid_perc_alt > par_max_perc_alt: continue
            if mid_perc_alt > par_max_perc_alt: continue

        ## only keep instances where var > min_perc_alt in parents
        if par_min_perc_alt != None:
            pid_gt_alt_depth = float(vcf_variant.gt_alt_depths[pid_idx])
            mid_gt_alt_depth = float(vcf_variant.gt_alt_depths[mid_idx])
            pid_perc_alt = pid_gt_alt_depth / trio_covs[1]
            mid_perc_alt = pid_gt_alt_depth / trio_covs[2]
            if pid_perc_alt < par_min_perc_alt: continue
            if mid_perc_alt < par_min_perc_alt: continue
        
        ## do genotypes pass phred likelihood thresholds if defined?
        if pro_homref_phredmin != None:
            if vcf_variant.gt_phred_ll_homref[iid_idx] < pro_homref_phredmin:
                continue
        if pro_het_phredmin != None:
            if vcf_variant.gt_phred_ll_het[iid_idx] < pro_het_phredmin:
                continue
        if pro_homalt_phredmax != None:
            if vcf_variant.gt_phred_ll_homalt[iid_idx] > pro_homalt_phredmax:
                continue
        ## do genotypes pass phred likelihood threshes for parents?
        if par_homref_phredmin != None:
            if vcf_variant.gt_phred_ll_homref[pid_idx] < par_homref_phredmin:
                continue
            elif vcf_variant.gt_phred_ll_homref[mid_idx] < par_homref_phredmin:
                continue
        if par_het_phredmax != None:
            if vcf_variant.gt_phred_ll_het[pid_idx] > par_het_phredmax:
                continue
            elif vcf_variant.gt_phred_ll_het[mid_idx] > par_het_phredmax:
                continue
        if par_homalt_phredmin != None:
            if vcf_variant.gt_phred_ll_homalt[pid_idx] < par_homalt_phredmin:
                continue
            elif vcf_variant.gt_phred_ll_homalt[mid_idx] < par_homalt_phredmin:
                continue

        ## newly homzygous variant detected
        hom_carriers.add(iid)

    return hom_carriers


def qual_impacts_screen(cyvcf2_variant_i, 
                        qual_impacts,
                        csq_subfield="CSQ"):
    qual_impact_pass = False
    for qual_impact in qual_impacts:
        if cyvcf2_variant_i.INFO.get(csq_subfield) == None: break
        elif cyvcf2_variant_i.INFO.get(csq_subfield).find(qual_impact) != -1:
            qual_impact_pass = True
            break
    return qual_impact_pass

def get_maxmin_csqs(cyvcf2_variant_i,
                    csq_keys,
                    max_impact_csqs=None,
                    max_csq_scores=None, 
                    min_csq_scores=None,
                    csq_subfield="CSQ",
                    impact_subfield="IMPACT"):

    csqs_str = cyvcf2_variant_i.INFO.get(csq_subfield)
    if csqs_str == None: return [],[],[]

    annottxs_i = AnnotTxs(csq_keys, csqs_str)                                                

    csqs_maximpact_list = []
    if max_impact_csqs != None:
        csqs_max_impact = annottxs_i.max_csq(max_impact_csqs,                                      
                                             impact_subfield,                                      
                                             "max")     
        for csq_max_impact in csqs_max_impact:
            csqs_maximpact_list.append(csq_max_impact)
    
    max_csq_scores_list = []
    if max_csq_scores != None:
        for max_csq_score in max_csq_scores:
            csq_max = annottxs_i.max_csq([max_csq_score],                                      
                                         max_csq_score,                                      
                                         "max")
            max_csq_scores_list.append(csq_max[0])

    min_csq_scores_list = []
    if min_csq_scores != None:
        for min_csq_score in min_csq_scores:
            csq_min = annottxs_i.max_csq([min_csq_score],
                                         min_csq_score,
                                         "min")
            min_csq_scores_list.append(csq_min[0])
    
    return csqs_maximpact_list,max_csq_scores_list,min_csq_scores_list

