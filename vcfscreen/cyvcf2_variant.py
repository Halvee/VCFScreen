
'''
    Filename : cyvcf2_variant.py 
    Author : Matt Halvorsen
    Email : mhalvors1@gmail.com
    Date created : 03/19/2018
    Date last modified : 03/19/2018
'''

def iterator(intervals, cyvcf2_vcf_i, verbose=False):
    '''
    iterator function for going through every vcf in every provided
    interval in vcf
    '''
    for interval in intervals:
        if verbose == True:
            if interval == "": print("Parsing entire VCF.")
            else: print("Parsing " + interval + ".")
        for cyvcf2_variant_i in cyvcf2_vcf_i(interval):
            yield cyvcf2_variant_i

def get_info_subfields(cyvcf2_vcf_i):
    info_subfields = []
    for vcf_header_entry in cyvcf2_vcf_i.header_iter():
        field = str(vcf_header_entry['HeaderType'])
        if field == "INFO":
            info_subfields.append(str(vcf_header_entry['ID']))
    return info_subfields

def get_csq_keys(cyvcf2_vcf_i, spliton="Format: ", delim="|"): 
    csq_keys = None
    for vcf_header_entry in cyvcf2_vcf_i.header_iter(): 
        try:
            field = str(vcf_header_entry['HeaderType'])
            id = str(vcf_header_entry['ID'])
            if field == "INFO" and id == "CSQ":
                csq_str_full = vcf_header_entry['Description']
                csq_str = csq_str_full.split(spliton)[-1]
                csq_str = csq_str.replace('"','')
                csq_keys = csq_str.split(delim)
        except:
            pass
    return csq_keys

def header_to_list(cyvcf2_vcf_i, delim=None,
                   main_fields=["IID","CHROM","POS","ID","REF","ALT",
                                "QUAL","FILTER"],
                   max_impact=False,
                   max_impact_csqs=["IMPACT","Gene"],
                   max_csq_scores=None,
                   min_csq_scores=None,
                   gt_varnames=[]):
    '''
    from cyvcf2.VCF object, convert header info to output header list,
    or delimited string if delimiter is provided
    '''
    header_list=main_fields
    if max_impact == True:
        for max_impact_csq in max_impact_csqs:
            header_list.append(max_impact_csq + "_maximpact") 
    if max_csq_scores != None:
        for max_csq_score in max_csq_scores:
            header_list.append(max_csq_score + "_max")
    if min_csq_scores != None:
        for min_csq_score in min_csq_scores:
            header_list.append(min_csq_score + "_min") 
    info_subfields=get_info_subfields(cyvcf2_vcf_i)
    header_list.extend(info_subfields)
    header_list.extend(gt_varnames)

    if delim != None:
        return delim.join([str(x) for x in header_list])
    else:
        return header_list

def variant_to_list(cyvcf2_variant_i, samples, sample_carriers,
                    info_subfields, gt_varnames,
                    csqs_maximpact=[], 
                    max_csq_scores = [],
                    min_csq_scores = [],
                    delim=None):
    '''
    using cyvcf2.Variant instance, sample_idx, info subfields,
    gt varnames, return list of values for sample variant/genotype,
    or a delimited string if delimiter is provided
    '''
    outs = []
    alt = cyvcf2_variant_i.ALT[0]
    out_var=[cyvcf2_variant_i.CHROM, cyvcf2_variant_i.POS, cyvcf2_variant_i.ID, 
             cyvcf2_variant_i.REF, alt, cyvcf2_variant_i.QUAL,
             cyvcf2_variant_i.FILTER]
    for csq_maximpact in csqs_maximpact: out_var.append(csq_maximpact)
    for max_csq_score in max_csq_scores: out_var.append(max_csq_score)
    for min_csq_score in min_csq_scores: out_var.append(min_csq_score)
    for subfield in info_subfields:
        out_var.append(cyvcf2_variant_i.INFO.get(subfield))
    for iid in sample_carriers:
        if isinstance(iid, str): iid = [iid]
        else: iid = list(iid)
        out_i = []
        idx = samples.samples[iid[-1]].idx
        for gt_varname in gt_varnames:
            val = getattr(cyvcf2_variant_i, gt_varname)
            out_i.append(val[idx])
        out = iid + out_var + out_i
        outs.append(out)

    if delim!=None:
        for i in range(len(outs)):
            for j in range(len(outs[i])):
                #outs[i][j] = str(outs[i][j])
                try:
                    outs[i][j] = str(outs[i][j])
                except:
                    outs[i][j] = outs[i][j].encode("utf8")
            outs[i] = delim.join(outs[i])
    return outs

def get_sample_alts(cyvcf2_variant_i, sample_idx):
    '''
    for sample whose index is provided in input, get the two
    genotypes represented at callsite. Ref allele is represented
    by "None". A homozygous call might have appearance of (C, C)
    as return for example. A het call will usually look like
    (None, G) for example. In rare var studies, beware of genotypes 
    made of a mixture of two non-reference alleles. 
    '''
    ref = cyvcf2_varaint_i.REF
    alts = cyvcf2_variant_i.ALT
    sample_gts = cyvcf2_variant_i.genotypes[sample_idx]
    gt1_idx = sample_gts[0]
    gt2_idx = sample_gts[1]
    if gt1_idx == 0:
        gt1 = None
    else:
        gt1 = str(alts[gt1_idx-1])
    if gt2_idx == 0:
        gt2 = str(alts[gt2_idx-1])
    return gt1,gt2
