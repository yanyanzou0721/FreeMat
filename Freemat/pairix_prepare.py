#!/usr/bin/env python
#
# The MIT License
#
# Copyright (c) 2022 yyzou.

### usage: python pairix_prepare.py in.bed out_prefix
### target: generate pairix region file
### input : bed file
### output: 
import pandas as pd
from sys import argv

in_bed = argv[1]
out_prefix = argv[2]

tmp = pd.read_csv(in_bed,header=None,sep="\t")
tmp[3] = tmp.apply(lambda x: str(x[0])+":"+ str(x[1]) +"-"+ str(x[2]),axis=1)  ### region for pairix
tmp[4] = tmp.apply(lambda x: str(x[0])+"-"+ str(x[1]) +"-"+ str(x[2]),axis=1)  ### new_chrom names
tmp[5] = tmp[2]-tmp[1]+1  #### new chrom length

res = []
for i in range(tmp.shape[0]):
    for j in range(i,tmp.shape[0]):
        res.append(tmp[3].loc[i]+"|"+tmp[3].loc[j])

out_pairix = out_prefix + "_pairix_region.txt"
with open(out_pairix,'w') as f:
    for i in res:
        f.write(i+"\n")
f.close()

out_chromsize = out_prefix + ".size"
tmp[[4,5]].to_csv(out_chromsize, header=False, index=False, sep='\t')
