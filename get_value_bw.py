#!/usr/bin/env python3

import argparse, re
import pyBigWig
import pandas as pd
from scipy.stats import spearmanr
import pyranges as pr


parser = argparse.ArgumentParser()
parser.add_argument('--bed', dest='bed', default='', type=str, help='name of bed file')
parser.add_argument('--bw', dest='bw', default='', type=str, help='filename of bigwig')
parser.add_argument('--out', dest='out', default='', type=str, help='output file of scores')

##main
args = parser.parse_args()

bwfile = args.bw
bedfile = args.bed
out = args.out

def read_bw(bwfile):
    bw = pyBigWig.open(bwfile)
    return bw

def get_value_bw(bw, interval):
    #bw is a readin bw from bwfile, interval is a list of [chr, start, end] and 0-based.
    mean_value = bw.stats(interval[0], interval[1], interval[2], exact=True)
    return mean_value[0]

def get_bed_value(bedfile, bw):
    print("read bed file")
    regions = {}
    bedlist = []
    with open(bedfile) as h:
        for line in h.readlines():
            line = line.strip().split()
            chr=line[0]
            start=int(line[1])
            end=int(line[2])
            strand=line[5]
            name=line[3]
            score=float(line[4])
            score=get_value_bw(bw, [chr,start,end])
            bedlist.append([chr, start,end,name,score,strand])
    return bedlist

def get_bed_value(bedfile, bw):
    print("read bed file")
    regions = {}
    bedlist = []
    with open(bedfile) as h:
        for line in h.readlines():
            line = line.strip().split()
            chr=line[0]
            start=int(line[1])
            end=int(line[2])
            strand=line[5]
            name=line[3]
            score=float(line[4])
            score=get_value_bw(bw, [chr,start,end])
            bedlist.append([chr, start,end,name,score,strand])
    return bedlist
    
bw = read_bw(bwfile)
bed = get_bed_value(bedfile, bw)
df = pd.DataFrame(bed, columns =['Chrom', 'Start', 'End', 'Name', 'Score', 'Strand'])
df.to_csv(out, sep='\t', index=False, header=False)
