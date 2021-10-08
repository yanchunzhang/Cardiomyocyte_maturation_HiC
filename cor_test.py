#!/usr/bin/env python3

import argparse
import pandas as pd
from scipy.stats import spearmanr
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', default='', type=str, help='name of input file')
parser.add_argument('--out', dest='out', default='', type=str, help='output file of correction test')

##main
args = parser.parse_args()
input = args.input
out = args.out

if not out:
    out = input + ".cor_test.txt"

df = pd.read_csv(input, header=None, sep="\t")
#gene_fc_list = list(np.log2(df[2]))
#score_fc_list = list(np.log2(df[1]))
gene_fc_list = list(df[2])
score_fc_list = list(df[1])
corr, _ = spearmanr(gene_fc_list, score_fc_list)
outfh = open(out,"a")
outfh.write(input)
outfh.write("\t")
outfh.write(str(corr))
outfh.write("\n")
outfh.close()
print(corr,"\n")
