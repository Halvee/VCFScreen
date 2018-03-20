
'''
Filename : cyvcf2_variant.py 
Author : Matt Halvorsen
Email : mhalvors1@gmail.com
Date created : 03/19/2018
Date last modified : 03/20/2018
'''

import cyvcf2
from annot import AnnotTxs

class Cyvcf2Vcf(object):
    '''
    Wrapper class containing variables derived from cyvcf2.VCF instances, as 
    well as functions for accessing cyvcf2.VCF variables
    '''
    def __init__(self, cyvcf2_vcf_i):
        self.cyvcf2_vcf = cyvcf2_vcf_i
        self.info_subfields=[]
        self.csq_keys=None

    def get_info_subfields(self):
        '''
        get all info subfield names from cyvcf2.VCF, store as list in self
        '''
        for vcf_header_entry in self.cyvcf2_vcf.header_iter():
            field = str(vcf_header_entry['HeaderType'])
            if field == "INFO":
                self.info_subfields.append(str(vcf_header_entry['ID']))
        return self

    def get_csq_keys(self, spliton="Format: ", delim="|"): 
        '''
        get and store CSQ keys, store as list variable in self
        '''
        for vcf_header_entry in self.cyvcf2_vcf.header_iter(): 
            try:
                field = str(vcf_header_entry['HeaderType'])
                id = str(vcf_header_entry['ID'])
                if field == "INFO" and id == "CSQ":
                    csq_str_full = vcf_header_entry['Description']
                    csq_str = csq_str_full.split(spliton)[-1]
                    csq_str = csq_str.replace('"','')
                    self.csq_keys = csq_str.split(delim)
            except:
                pass
        return self

    def header_to_list(self, delim=None,
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
        self.get_info_subfields()
        header_list.extend(self.info_subfields)
        header_list.extend(gt_varnames)

        if delim != None:
            return delim.join([str(x) for x in header_list])
        else:
            return header_list

    def iterator(self, intervals, verbose=False):
        '''
        iterator function for going through every vcf in every provided
        interval in vcf
        '''
        for interval in intervals:
            if verbose == True:
                if interval == "": print("Parsing entire VCF.")
                else: print("Parsing " + interval + ".")
            for cyvcf2_variant_i in self.cyvcf2_vcf(interval):
                yield cyvcf2_variant_i

class Cyvcf2Variant(object):
    '''
    Wrapper class containing variables derived from cyvcf2.Variant instances, 
    as well as functions for accessing cyvcf2.Variant variables
    '''
    def __init__(self, cyvcf2_variant_i):
        self.cyvcf2_variant = cyvcf2_variant_i
        self.annot_txs = None
        self.metadata = {}

    def variant_to_list(self, samples, sample_carriers,
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
        alt = self.cyvcf2_variant.ALT[0]
        out_var=[self.cyvcf2_variant.CHROM, 
                 self.cyvcf2_variant.POS, 
                 self.cyvcf2_variant.ID, 
                 self.cyvcf2_variant.REF, 
                 alt, 
                 self.cyvcf2_variant.QUAL,
                 self.cyvcf2_variant.FILTER]
        for csq_maximpact in csqs_maximpact: out_var.append(csq_maximpact)
        for max_csq_score in max_csq_scores: out_var.append(max_csq_score)
        for min_csq_score in min_csq_scores: out_var.append(min_csq_score)
        for subfield in info_subfields:
            out_var.append(self.cyvcf2_variant.INFO.get(subfield))
        for iid in sample_carriers:
            if isinstance(iid, str): iid = [iid]
            else: iid = list(iid)
            out_i = []
            idx = samples.samples[iid[-1]].idx
            for gt_varname in gt_varnames:
                val = getattr(self.cyvcf2_variant, gt_varname)
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
    
    def get_sample_alts(self, sample_idx):
        '''
        for sample whose index is provided in input, get the two
        genotypes represented at callsite. Ref allele is represented
        by "None". A homozygous call might have appearance of (C, C)
        as return for example. A het call will usually look like
        (None, G) for example. In rare var studies, beware of genotypes 
        made of a mixture of two non-reference alleles. 
        '''
        ref = self.cyvcf2_variant.REF
        alts = self.cyvcf2_variant.ALT
        sample_gts = self.cyvcf2_variant.genotypes[sample_idx]
        gt1_idx = sample_gts[0]
        gt2_idx = sample_gts[1]
        if gt1_idx == 0:
            gt1 = None
        else:
            gt1 = str(alts[gt1_idx-1])
        if gt2_idx == 0:
            gt2 = str(alts[gt2_idx-1])
        return gt1,gt2

    def get_annot_txs(self, csq_keys, csq_subfield="CSQ"):
        """
        store annotation data from INFO field, consequence subfield
        """
        if self.annot_txs != None: return self
        csqs_str = self.cyvcf2_variant.INFO.get(csq_subfield)
        if csqs_str == None: return self
        self.annot_txs = AnnotTxs(csq_keys, csqs_str)
        return self

    def maxmin_csqs(self,
                    max_impact_csqs=None,
                    max_csq_scores=None, 
                    min_csq_scores=None,
                    csq_subfield="CSQ",
                    impact_subfield="IMPACT"):
        '''
        define consequence subfield as well as impact subfield within
        consequence, derive metadata to go along with gene/transcript with
        maximum impact classification, if provided by user also find max or 
        min defined score variable, return as 
        [max_impact_csqs,max_csq_scores,min_csq_scores]
        '''
        if self.annot_txs == None: return [],[],[]

        annottxs_i = self.annot_txs                                  

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

    def qual_impacts_screen(self, 
                            qual_impacts,
                            csq_subfield="CSQ"):
        '''
        screen consequence subfield in VCF INFO data to see if variant has 
        at least one transcript with one of the qualifying impact classifs
        '''
        qual_impact_pass = False
        for qual_impact in qual_impacts:
            csq_str = self.cyvcf2_variant.INFO.get(csq_subfield)
            if csq_str == None: break
            elif csq_str.find(qual_impact) != -1:
                qual_impact_pass = True
                break
        return qual_impact_pass

    def compute_maf(self, samples_include_idxs):
        gts = self.cyvcf2_variant.gt_types[samples_include_idxs]
        n_alt_alleles = 0
        n_samples = len(gts)
        for i in range(len(gts)):
            if gts[i] != 3: 
                n_alt_alleles += gts[i]
        maf = float(n_alt_alleles) / (2*n_samples) 
        return maf

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
