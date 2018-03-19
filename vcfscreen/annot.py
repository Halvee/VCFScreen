'''
    Filename : annot.py 
    Author : Matt Halvorsen
    Email : mhalvors1@gmail.com
    Date created : 03/19/2018
    Date last modified : 03/19/2018
'''

import sys

IMPACT_RANKINGS = {"HIGH":3,
                   "MODERATE":2,
                   "LOW":1,
                   "MODIFIER":0,
                   "NA":-1}
IMPACT_RANKINGS_INV = {3:"HIGH",
                       2:"MODERATE",
                       1:"LOW",
                       0:"MODIFIER",
                       -1:"NA"}
CSQ_KEYS_STR="Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature," + \
             "BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position," + \
             "Protein_position,Amino_acids,Codons,Existing_variation," + \
             "DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID," + \
             "CANONICAL,CCDS,GIVEN_REF,USED_REF,SIFT,PolyPhen," + \
             "DOMAINS,HGVS_OFFSET"

class AnnotTxs(object):
    def __init__(self, annot_keys, annot_line, 
                 delim=",", subdelim="|", format="ANN"):
        self.annot_keys = annot_keys
        self.annot_line = annot_line
        self.delim = delim
        self.subdelim = subdelim
        self.format = format
        self.max_eff = None
        self.load_annot_line()
    
    def load_annot_line(self):
        if self.annot_line == None: return
        self.annots = self.annot_line.split(self.delim)
        for i in range(len(self.annots)):
            tx = self.annots[i]
            annot_vals = tx.split(self.subdelim)
            self.annots[i] = AnnotTx(self.annot_keys, annot_vals)
        return self

    def max_csq(self, csqs_return, csq_key, csq_func):
        if csq_func != "max" and csq_func != "min":
            raise Exception(csq_func + " not applicable on " + \
                            csq_key + ", max or min only.")
        min_val = float("inf")
        max_val = - float("inf")
        csqs_return_out = []
        if self.annot_line == None:
            for csq_return in csqs_return:
                csqs_return_out.append("None")
            return csqs_return_out
        for csq in csqs_return:
            csqs_return_out.append([])
        for i in range(len(self.annots)):
            csq_i = self.annots[i].__dict__[csq_key]
            if csq_i == "": 
                csq_i = "NA"
                continue
            elif csq_key == "IMPACT":
                csq_i = IMPACT_RANKINGS[csq_i]
            else:
                csq_i = float(csq_i)
            if csq_i > max_val: max_val = csq_i
            if csq_i < min_val: min_val = csq_i
        i_keep={"max":[], "min":[]}
        for i in range(len(self.annots)):
            csq_i = self.annots[i].__dict__[csq_key]
            if csq_i == "": 
                csq_i = "NA"
                continue
            if csq_key == "IMPACT":
                csq_i = IMPACT_RANKINGS[csq_i]
            if float(csq_i) == max_val:
                i_keep["max"].append(i)
            if float(csq_i) == min_val:
                i_keep["min"].append(i)

        for i in i_keep[csq_func]:
            for j in range(len(csqs_return)):
                val = self.annots[i].__dict__[csqs_return[j]]
                if val == "": continue
                csqs_return_out[j].append(val)
        if len(csqs_return_out) == 0:
            csqs_return_out = "NA"
        for i in range(len(csqs_return_out)):
            if len(csqs_return_out[i]) == 0:
                csqs_return_out[i] = "NA"
            else:
                csqs_return_out[i] = ",".join(list(set(csqs_return_out[i])))
        return csqs_return_out

class AnnotTx(object):
    def __init__(self, keyvals, vals):
        for keyval in keyvals:
            self.__dict__[keyval] = None
        if len(keyvals) != len(vals):
            raise Exception("n(annot keys) != n(annot vals)")
        for i in range(len(vals)):
            key_i = keyvals[i]
            val_i = vals[i]
            self.__dict__[key_i] = val_i

def process_csq_header_desc(header_desc_str,                                    
                            replace_chars=["'",'"'," "],                        
                            main_split="Format: ",                              
                            subsplit="|"):                                      
    header_desc_list = header_desc_str.split(main_split)                        
    header_desc_str = header_desc_list[-1]                                      
    for replace_char in replace_chars:                                          
        header_desc_str = header_desc_str.replace(replace_char, "")             
    header_desc_list = header_desc_str.split(subsplit)                          
    return header_desc_list     
