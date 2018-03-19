
'''
    Filename : misc.py 
    Author : Matt Halvorsen
    Email : mhalvors1@gmail.com
    Date created : 03/19/2018
    Date last modified : 03/19/2018
'''

import os
import sys
import string

def str_none_split(string, split_on):
    try:
        str_list = string.split(split_on)
    except:
        return None
    return str_list

def init_out_file(out_filename, 
                  init_line = "",
                  check_existance=True):
    if os.path.exists(out_filename):
        yn=raw_input("File " + out_filename + " exists. " + \
                     "Overwrite? [Y/n] ")
        if yn != "Y":
            print("Not overwriting " + out_filename + ".")
            print("Exiting...")
            sys.exit(1)
    out_fh = open(out_filename, "w")
    out_fh.write(init_line + "\n")
    out_fh.close()
    return 

def append_out_file(out_filename, line):
    out_fh = open(out_filename, "a")
    out_fh.write(line + "\n")
    out_fh.close()
    return 

