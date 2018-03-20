#!/usr/bin/env python

import os
import sys
import importlib
import textwrap

AUTHORSHIP="""
Filename : VCFscreen_newlyhomozygous.py 
Author : Matt Halvorsen
Email : mhalvors1@gmail.com
Date created : 03/19/2018
Date last modified : 03/19/2018
Python version : 2.7
"""

PROGS_AVAIL={"denovo":"VCFscreen_denovo",
             "newlyhom":"VCFscreen_newlyhomozygous",
             "newlyhem":"VCFscreen_newlyhemizygous",
             "newlycpht":"VCFscreen_newlycomphet",
             "parmosaic":"VCFscreen_parentalmosaic",
             "partrans":"VCFscreen_parentaltransmitted",
             "gt":"VCFscreen_vargeno"}
PROGS_DESC={"denovo":"Make de novo mutation calls on trio probands",
            "newlyhom":"Identify variants that are homozygous in a trio " + \
                       "proband and found in heterozygous from in each parent",
            "newlyhem":"Identify variant calls that are hemizygous in a " + \
                       "male trio proband, absent from the father, and " + \
                       "found in heterozygous form in the mother.",
            "newlycpht":"From tab-delimited file with parental transmitted " +\
                        "variant calls, identify variant transmissions that" +\
                        "form newly compound heterozygous genotypes absent "+\
                        "from each parent",
            "parmosaic":"Calls that are heterozygous in a trio proband, " +\
                        "absent from one parent, and found in a percentage " +\
                        "of reads in other parent inconsistent with " + \
                        "heterozygosity",
            "partrans":"Calls that are heterozygous in a trio parent, " + \
                       "absent in the other trio parent, and heterozygous " + \
                       "in the trio proband",
            "gt":"Single sample variant calls that are either heterozygous " +\
                 "or homozygous"}
            
def main():
    global PROGS_AVAIL
    try:
        ARGS = sys.argv[1:]
        prog = ARGS[0]
        PROGS_AVAIL[prog]
    except:
        usage()
    prog_import_str = PROGS_AVAIL[prog]
    func = importlib.import_module(prog_import_str)
    func.main(ARGS = ARGS[1:])
    return

def usage():
    global PROGS_AVAIL
    global PROGS_DESC
    global AUTHORSHIP
    print("")
    print("VCFscreen.py : extract genotypes for unrelated and family-based samples")
    print("")
    print("powered by   : cyvcf2")
    print("    Pedersen et al., Bioinformatics, Volume 33, Issue12, 15 June 2017")
    print("    https://github.com/brentp/cyvcf2")
    prog_names = list(PROGS_AVAIL.keys())
    print("")
    print("usage        : VCFscreen.py <function> [options]")
    prog_names.sort()
    print("")
    print("functions    : ")
    for prog_name in prog_names:
        print("")
        print("[ " + prog_name + " ]")
        for line in textwrap.wrap(PROGS_DESC[prog_name], 76): 
            print("    " + line)
    print("")
    sys.exit(1)

if __name__ == "__main__":
    main()
