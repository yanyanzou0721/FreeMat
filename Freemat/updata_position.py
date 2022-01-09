#!/usr/bin/env python
#
# The MIT License
#
# Copyright (c) 2022 yyzou.

### usage: python pairix_prepare.py in.bed out_prefix
### target: update chrom name & location
### input :  pairix_region_pairs.txt , bed.file
### output: updated_pairix_region_pairs.txt

import pandas as pd
from sys import argv

#### given series , change its position with 
### three params : row of pairix_region_pairs.txt  , chrom_info in bed.file  , col_number of chrom in pairs.txt 
def modify_chrom(s,mat,ind):
    for i in mat[mat[0]==s[ind]].index:
        if (s[ind+1]-mat[1].loc[i])*(s[ind+1]-mat[2].loc[i])<0:
            chrom = i
            break
    return chrom

pairs = argv[1]
in_bed = argv[2]
out_put = argv[3]

### set index 
new_chrom = pd.read_csv(in_bed,sep="\t",header=None)
new_chrom.index = new_chrom.apply(lambda x: str(x[0])+"-"+ str(x[1]) +"-"+ str(x[2]),axis=1) 

### pairix file readin
extract_data = pd.read_csv( pairs ,sep="\t",header=None)
extract_data["new_chrom1"] = extract_data.apply(lambda x: modify_chrom(x,new_chrom,1),axis=1)
extract_data["new_chrom2"] = extract_data.apply(lambda x: modify_chrom(x,new_chrom,3),axis=1)
extract_data["pos_1"] = extract_data.apply(lambda x: x[2]+1-new_chrom[1].loc[x["new_chrom1"]],axis=1)
extract_data["pos_2"] = extract_data.apply(lambda x: x[4]+1-new_chrom[1].loc[x["new_chrom2"]],axis=1)

 
extract_data[[0,"new_chrom1","pos_1","new_chrom2","pos_2",5,6,7]].to_csv( out_put , sep="\t" , header=False , index=False)
