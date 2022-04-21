#!/usr/bin/python

###     Author:     Xun Chen
###     Date:       2020/5/10
###     Contact:    xunchen85@gmail.com

import gzip
import re
import sys
import getopt

### define variables
try:
    opts,args = getopt.getopt(sys.argv[1:], '-h:-i:', ['help', 'input='])
except getopt.GetoptError:
    sys.exit()

for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        sys.exit()
    if opt_name in ('-i','--input'):
        inputFile = opt_value

### read input file
with open(inputFile, 'rt') as inputfile:
    for line in inputfile:
        line2 = re.split(r'\s+',line)
        if line2[0] == "sample_ID":
            continue
        medium = int(line2[1]) + int((int(line2[2])-int(line2[1]))/2)
        medium2 = medium + 1
        line2[3] = line2[3] + ":" + str(round((int(line2[2]) - int(line2[1])),0))
        line2[1] = medium
        line2[2] = medium2
        print("\t".join(map(str,line2)) + "\t+")
