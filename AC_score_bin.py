#!/usr/bin/env python3

import argparse, re
import pyBigWig
import pandas as pd
from scipy.stats import spearmanr
from scipy.stats import gmean
from statistics import mean
import pyranges as pr
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', default='', type=str, help='name of input file')
parser.add_argument('--out', dest='out', default='', type=str, help='output file of correction test')
parser.add_argument('--binsize', dest='binsize', default=5, type=int, help='size of bin')
##main
args = parser.parse_args()
input =args.input
out = args.out
binsize=args.binsize
if not out:
    input_prefix = os.path.splitext(input)[0]
    out = input_prefix + ".bin" + str(binsize) + ".txt"
df = pd.read_csv(input, header=None, sep="\t")
df = df.sort_values(by=[1])
gene_fc_list = list(df[2])
score_fc_list = list(df[1])

new_gene_fc =[]
new_score_fc = []
outlist = []
nbin = int(len(gene_fc_list)/binsize)+1
for i in range(0, nbin):
    start = binsize*i
    end = binsize * (i+1)
    tmpgenelist = gene_fc_list[start:end]
    tmpscorelist = score_fc_list[start:end]
    if tmpgenelist:
        gene_fc = mean(tmpgenelist)
        score_fc = mean(tmpscorelist)
        new_gene_fc.append(gene_fc)
        new_score_fc.append(score_fc)
        outlist.append(["bin"+str(i), score_fc, gene_fc])

corr, _ = spearmanr(new_gene_fc, new_score_fc)
print(corr,"\n")
df=pd.DataFrame(outlist)
df.to_csv(out,header=None,sep="\t",index=None)

outfh = open(out+".cor_test.txt","a")
outfh.write(input)
outfh.write("\t")
outfh.write(str(binsize))
outfh.write("\n"+str(corr))
outfh.write("\n")
outfh.close()
