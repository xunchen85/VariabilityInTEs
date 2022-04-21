### Author:     Xun Chen
### Date:       2020/6/10
### Contact:    xunchen85@gmail.com

#import gzip
import getopt
import re
import sys
#import numpy as np
#import pandas as pd
import glob, os


### example
# python Combined_TEenrichmentByTEfamily.py -i counts.anno.fileList -l hg19_rmsk_TE.bed.TEfamilyCounts >Peaks_public.counts
### define variables
try:
    opts,args = getopt.getopt(sys.argv[1:], '-h:-i:-l:', ['help', 'inputFile=', 'listTEFile='])
except getopt.GetoptError:
    print ('python Association_Genes_TEs.py -i <inputFile> -l <listTEFile>')
    sys.exit()
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
#        print("[*] Help info")
        sys.exit()
    if opt_name in ('-i','--inputFile'):
        input = opt_value
#        print('Input file: ' + input)
    if opt_name in ('-l','--listTEFile'):
        listTEf = opt_value
#        print('Bed file: ' + listTEf)

### list of files
fileList = []
with open(input,'rt') as fi:
    for line in fi:
        line2=re.split(r'\s+',line)
        fileList.append(line2[0])
#for file in glob.glob("*.anno"):

### TE list
listTE = {}
with open(listTEf, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        listTE[line2[0]] = line2
        listTE[line2[0]][2:(len(fileList)+2)] = [0] * len(fileList)

i = 2
listFiles = []
### table header
for file in fileList:
    listTE['TE_family'][i] = file.replace(".counts.anno","")
    #print(file)
    ### read anno file
    with open(file, 'rt') as inp:
        if os.stat(file).st_size == 0:
            listTE['TE_family'][i] = listTE['TE_family'][i] + ".notfinished"
            continue
        else:
            for line in inp:
                line2 = re.split(r'\s+',line)
                #print(line2)
                if line2[0] == "Binary":
                    listTE['TE_family'][i] = listTE['TE_family'][i] + ".notfinished"
                    continue
                elif len(line2) >= 4:
                    tmpName2 = line2[1] + ":" + line2[2] + ":" + line2[3]
                else:
                    listTE['TE_family'][i] = listTE['TE_family'][i] + ".notfinished"
                if tmpName2 in listTE:
                    listTE[tmpName2][i] = listTE[tmpName2][i] + int(line2[4])
    i = i+1

for TE in listTE:
    print("\t".join(map(str,listTE[TE])))
