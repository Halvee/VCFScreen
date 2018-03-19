#!/usr/bin/env python

import os
import sys

PROGS_AVAIL=["denovo","newlyhom",
             "newlyhem", "newlycpht"]

def main():
    global PROGS_AVAIL
    try:
        ARGS = sys.argv[1:]
        prog = ARGS[0]
    except:
        usage()
    if prog == "denovo":
        import vcfscreen_denovo
        vcfscreen_denovo.main()
    else:
        usage()
    return

def usage():
    global PROGS_AVAIL
    progs_avail_str=" | ".join(PROGS_AVAIL)
    print("vcfscreen.py < " + \
          progs_avail_str + " >")
    sys.exit(1)

if __name__ == "__main__":
    main()
