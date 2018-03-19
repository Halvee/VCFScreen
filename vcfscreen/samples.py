
'''
    Filename : samples.py 
    Author : Matt Halvorsen
    Email : mhalvors1@gmail.com
    Date created : 03/19/2018
    Date last modified : 03/19/2018
'''

import sys
from collections import defaultdict

class Samples(object):
    """
    object for storing family and singleton sample data from a .fam file
    """
    def __init__(self, fam_file):
        self.fam_file = fam_file
        self.samples = {}
        self.trios = {}
        self.fam_trios = defaultdict(set)
        self.parentchild = defaultdict(set)
        self.cases = set()
        self.ctrls = set()
        self.males = set()
        self.females = set()
        self.case_fams = set()
        self.load_fam_file()
        self.get_males_females()
        self.get_cases_ctrls()
        self.get_trios()
        self.get_parentchild()

    def load_fam_file(self):
        """
        load family based sample info from .fam file
        """
        fam_fh = open(self.fam_file, "r")
        # store all data as cases/ctrls first
        for line in fam_fh:
            data = line.rstrip().split()
            [fid, iid, pid, mid, gender_num, pheno] = data[:6]
            sample_x = Sample(fid, iid, pid, mid, gender_num, pheno)
            self.samples[iid] = sample_x
        fam_fh.close()
        return self

    def get_cases_ctrls(self):
        """
        get dict of case iids, dict of control iids
        """
        for iid in self.samples:
            if self.samples[iid].pheno == "2":
                self.cases.add(iid)
            elif self.samples[iid].pheno == "1":
                self.ctrls.add(iid)
        return self

    def get_trios(self):
        """
        from case dict, get probands with pid + mid store to trios dict
        """
        for cohort in (self.cases, self.ctrls):
            for iid in self.samples:
                sample_i = self.samples[iid]
                if sample_i.pid!="0" and sample_i.mid!="0":
                    self.trios[iid] = self.samples[iid]
                    self.fam_trios[sample_i.fid].add(iid)
        return self

    def get_parentchild(self):
        """
        load dict of parent_x -> (child1, .., childN)
        """
        for iid in self.samples:
            for parentid in (self.samples[iid].pid, self.samples[iid].mid):
                if parentid != "0":
                    self.parentchild[parentid].add(iid)
        return self

    def get_males_females(self):
        for iid in self.samples:
            if self.samples[iid].gender == "M":
                self.males.add(iid)
            elif self.samples[iid].gender == "F":
                self.females.add(iid)
        return self

    def get_sonmother_pairs(self):
        sonmother={}
        for iid in self.samples:
            if self.samples[iid].gender =="M" and self.samples[iid].mid != "0":
                sonmother[iid] = self.samples[iid].mid
        return sonmother

    def print_stats(self):
        print("Summary stats for :")
        print(self.fam_file)
        print("")
        print("Males : " + str(len(self.males)))
        print("Females : " + str(len(self.females)))
        print("Cases : " + str(len(self.cases)))
        print("Controls : " + str(len(self.ctrls)))
        print("Trios : " + str(len(self.trios)))
        print("Parent/Child pairs : " + str(dict_set_len(self.parentchild)))
        print("")
        return 

    def print_varcount_stats(self, var_types=["dnm"], which="trios"):
        tot = 0
        if which not in self.__dict__:
            raise Exception("group " + which + " not in Samples obj.")
        iids = list(self.__dict__[which].keys())
        iids.sort()
        print("")
        for var_type in var_types:
            print("Variant count stats for samples : " + var_type)
            for iid in iids:
                try:
                    var_count = self.samples[iid].varcounts[var_type]
                except:
                    var_count = 0
                print(iid + " " + str(var_count))
                tot += var_count
        print("")
        print("Average : " + str(float(tot) / len(iids)))
        print("")
        return


    def get_vcf_idx(self, vcf_samples_list):
        idx = {}
        for i in range(len(vcf_samples_list)):
            idx[vcf_samples_list[i]] = i
        for iid in self.samples:
            if iid in idx: self.samples[iid].idx=idx[iid]
        return self

    def get_nmale_nfemale(self):
        nmale=0
        nfemale=0
        for iid in self.samples:
            if self.samples[iid].gender=="M":
                nmale += 1
            elif self.samples[iid].gender=="F":
                nfemale += 1
        return nmale,nfemale

class Sample(object):
    def __init__(self, fid, iid, pid, mid, gender_num, pheno):
        self.fid = fid
        self.iid = iid
        self.pid = pid
        self.mid = mid
        self.gender=gender_num_to_str(gender_num)
        self.pheno = pheno
        self.idx=None
        self.varcounts = defaultdict(int)

def gender_num_to_str(gender_num):
    gender_num = int(gender_num)
    if gender_num == 1:
        return "M"
    elif gender_num == 2:
        return "F"
    else:
        return "Unknown"

def dict_set_len(dict_set):
    k = 0
    for key_i in dict_set:
        for val_j in dict_set[key_i]:
            k += 1
    return k
