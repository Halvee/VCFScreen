
'''
    Filename : VCFscreen_newlyhemizygous.py 
    Author : Matt Halvorsen
    Email : mhalvors1@gmail.com
    Date created : 03/19/2018
    Date last modified : 03/19/2018
'''

import os
import operator

class VcfCnds(object):
    def __init__(self, cnds_file):
        self.cnds_file = cnds_file
        self.cnds = []
        self.read_cnds_file()

    def read_cnds_file(self, print_to_stdout=False):
        if self.cnds_file == None: return self
        fh = open(self.cnds_file, "r")
        i = 0
        for line in fh:
            i += 1
            if line[0] == "#": continue
            elif line.rstrip() == "": continue
            try:
                data = line.rstrip().split()
                for i in range(3):
                    data[i] = data[i].replace('"','')
                [field, operand_str, value] = data[:3]
                value = format_value(value)
            except:
                raise Exception("Improper formatting in " + self.cnds_file + ": " + \
                                "line " + str(i) + " : " + line.rstrip())
            if operand_str in ["in","nin"]:
                operand = operand_str
                if os.path.isfile(value):
                    value = read_set_file(value)
                else:
                    value = set(value.split(","))
            elif operand_str in ["grep","grepv"]:
                operand = operand_str
                value = value
            else:
                try:
                    operand = operator.__dict__[operand_str]
                except:
                    raise Exception("Improper formatting in " + self.cnds_file + ": " + \
                                    "operand " + operand_str + " not defined.")
            cnd = [field, operand, value]
            self.cnds.append(cnd)
            if print_to_stdout == True:
                out_str = " ".join(["CONDITION",":",field,
                                    operand_str, value])
                print(out_str)
                
        fh.close()
        return self

    def test_gt(self, vcf_variant, idx):
        for cnd in self.cnds:
            [field, operand, value] = cnd
            vals = getattr(vcf_variant, field)
            vcf_val = vals[idx]
            if cmp_vals(vcf_val, value, operand) == False:
                return False
        return True

    def test_variant(self, vcf_variant):
        for cnd in self.cnds:
            [field, operand, value] = cnd
            if field.find("INFO:") == 0:
                x = field.split(":")
                subfield = x[1]
                vcf_val = vcf_variant.INFO.get(subfield)
            elif field.find("GT:") == 0:
                pass
            else:
                vcf_val = getattr(vcf_variant, field)
            if cmp_vals(vcf_val, value, operand) == False:
                return False

        return True

def read_set_file(set_filename):
    set_out = set()
    fh = misc.open_file(set_filename)
    for line in fh:
        set_out.add(line.rstrip())
    fh.close()
    return set_out

def cmp_vals(val, value, operand):
    if operand == "in":
        if val not in value:
            return False
    elif operand == "nin":
        if val in value:
            return False
    elif operand == "grep":
        if val.find(value) == -1:
            return False
    elif operand == "grepv":
        if val.find(value) != -1:
            return False
    else:
        if value != None and val == None: 
            val = 0
        if operand(val, value) == False:
            return False
    return True

def format_value(value):
    if value == "None":
        value = None
    elif os.path.isfile(value):
        value = read_set_file(value)
    elif value.find(",") != -1:
        value = set(value.split(","))
    else:
        try:
            value = float(value)
        except:
            value = str(value)
    return value
