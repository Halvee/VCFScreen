
'''
    Filename : VCFscreen_newlyhemizygous.py 
    Author : Matt Halvorsen
    Email : mhalvors1@gmail.com
    Date created : 03/19/2018
    Date last modified : 03/19/2018
'''

def args_vars(args):
    args.add_argument("--filter-include", type=str,
                      action="store", default=None,
                      help="comma delim set of FILTER classifs to allow, " + \
                           "in addition to PASS.")
    args.add_argument("--min-coverage", type=int,
                      action="store", default=10)
    args.add_argument("--min-qual", type=float,
                      action="store", default=0)
    args.add_argument("--af-max", type=float,
                      action="store", default=None,
                      help="maximum allowed allele freq")
    args.add_argument("--af-max-fields", type=str,
                      action="store", default=None,
                      help="comma delim set of allele freq fields in vcf INFO.")
    args.add_argument("--internal-af-max", type=float,
                      action="store", default=None,
                      help="maximum allowed internal leave-one-out allele freq")
    args.add_argument("--vcf-info-flags-exclude", type=str,
                      action="store", default=None,
                      help="comma-delimited flags in VCF info to skip vars " + \
                           "based on.")
    args.add_argument("--intervals", action="store", type=str, default=None,
                      help="intervals to screen for variants, " + \
                           "can be either a single interval or file.")
    return args

def args_annots(args):
    args.add_argument("--max-impact", action="store_true", default=False,
                      help="derive transcript(s) of maximum impact in CSQ subfield.")
    args.add_argument("--max-impact-csqs", action="store", 
                      default="IMPACT,Gene,SYMBOL,CCDS",
                      help="comma delim transcript lvl CSQ items for max " + \
                           "impact, to write to output tsv")
    args.add_argument("--qual-impacts", action="store", type=str,
                      default=None,
                      help="CSQ impact str requires one of these strings be " +\
                           "found, otherwise don't list variant (comma delim).")
    args.add_argument("--max-csq-scores", action="store", type=str,
                      default=None,
                      help="comma-delim list of CSQ keys to get max tx score on.")
    args.add_argument("--min-csq-scores", action="store", type=str,
                      default=None,
                      help="comma-delim	list of	CSQ keys to get	min tx score on.")
    return args

def args_phredthresh(args,
                     sampletypes=["iid","par","pid","mid"],
                     gt_types=["homref","het","homalt"],
                     threshes=["max","min"]):
    for sampletype in sampletypes:
        for gt_type in gt_types:
            for thresh in threshes:
                argname = "--" + sampletype + "-" + \
                          gt_type + "-phred" + thresh
                args.add_argument(argname, type=int,
                                  action="store", default=None)
    return args

def args_perc_alt(args,
                  sampletypes=["iid","par","pid","mid"],
                  threshes=["max","min"]):
    for sampletype in sampletypes:
        for thresh in threshes:
            argname = "--" + sampletype + "-" + thresh + "-perc-alt"
            args.add_argument(argname, action="store", type=float, default=None,
                              help="minimum percent alt allele in " + sampletype+".")
    return args

def args_cnds(args, variant_cnds=True,
              sampletypes=["iid","par","pid","mid"]):
    if variant_cnds == True:
        args.add_argument("--variant-cnds", action="store", type=str, default=None,
                          help=".cnds file for variant annotations.")
    for sampletype in sampletypes:
        argname = "--" + sampletype + "-cnds"
        args.add_argument(argname, action="store", type=str, default=None,
                          help=".cnds file for "+sampletype+" genotype annotations.")
    return args

def args_req(args):
    args.add_argument("in_fam",                                                 
                      help="input .fam file with sample info.")                 
    args.add_argument("in_vcf",                                                 
                      help="input vcf file") 
    args.add_argument("out_tsv",
                      help="output tsv file.")
    return args
