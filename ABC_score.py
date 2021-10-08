#!/usr/bin/env python3

"""
Created on Aug 16 2021

##Calc ABC score from contact matrix and chip-seq bdg file
@author: yanchunzhang
##promtor can be any region around TSS location
"""

import argparse, re, time
import cooler as clr
import numpy as np
from scipy.stats import spearmanr
import pandas as pd
import os
import pickle
import pyranges as pr
from mpire import WorkerPool
import pyBigWig


parser = argparse.ArgumentParser()
#parser.add_argument('--matrix1', dest='matrix1', default='', type=str, help='name of dense matrix file1')
parser.add_argument('--matrixfile1', dest='matrixfile1', default='', type=str, help='name of matrix file1')
parser.add_argument('--matrixfile2', dest='matrixfile2', default='', type=str, help='name of matirx file2')
#parser.add_argument('--pickle1', dest='pickle1', default='', type=str, help='name of pickle1')
#parser.add_argument('--pickle2', dest='pickle2', default='', type=str, help='name of pickle2')
parser.add_argument('--promoter1', dest='promoter1', default='', type=str, help='name of file containing promoter information of sample1')
parser.add_argument('--promoter2', dest='promoter2', default='', type=str, help='name of file containing promoter information of sample2')
parser.add_argument('--peak1', dest='peak1', default='', type=str, help='name of bed file1 containing chip-seq peak signal')
parser.add_argument('--peak2', dest='peak2', default='', type=str, help='name of bed file2 containing chip-seq peak signal')
parser.add_argument('--bw1', dest='bw1', default='', type=str, help='name of bw file1 for chip-seq signal')
parser.add_argument('--bw2', dest='bw2', default='', type=str, help='name of bw file2 for chip-seq signal')
parser.add_argument('--chr', dest='chr', default='', type=str, help='')
parser.add_argument('--res', dest='res', default=5000, type=int, help='resolution of matrix')
parser.add_argument('--window', dest='window', default=5000000, type=int, help='size of window')
parser.add_argument('--remove_promoter', dest='remove_promoter', default=1, type=int, help='0 = not exclued other promoters in enhancer, 1 = yes')
parser.add_argument('--tss_intensity', dest='tss_intensity', default=0, type=int, help='0 = not count on tss intensity(score), 1 = count')
parser.add_argument('--tss_contact', dest='tss_contact', default=0, type=int, help='0 = not count on tss contact, 1 = count')
parser.add_argument('--gene_fc_file', dest='gene_fc_file', default='', type=str, help='file containing log2fc of gene expression')
parser.add_argument('--score_out1', dest='score_out1', default='', type=str, help='file name for output of scores of sample1')
parser.add_argument('--score_out2', dest='score_out2', default='', type=str, help='file name for output of scores of sample2')
parser.add_argument('--out', dest='out', default='', type=str, help='final output file containing score of sample1 and sample2, and foldchange of score and Gene_expression')
parser.add_argument('--nthreads', dest='nthreads', default=1, type=int, help='number of threads in certain steps')
parser.add_argument('--peak_signal_constant', dest='peak_signal_constant', default=0, type=int, help='if set peak_signal as 1 or not')
parser.add_argument('--contact_constant', dest='contact_constant', default=0, type=int, help='if set contact as 1 or not')

##main
args = parser.parse_args()

matrixfile1 = args.matrixfile1
matrixfile2 = args.matrixfile2

matfiletype = os.path.splitext(matrixfile1)[1]

promfile1 = args.promoter1
promfile2 = args.promoter2
peakfile1 = args.peak1
peakfile2 = args.peak2
bwfile1 = args.bw1
bwfile2 = args.bw2
res=args.res
chrom=args.chr
window=args.window
remove_promoter = args.remove_promoter
tss_intensity = args.tss_intensity
tss_contact = args.tss_contact

gene_fc_file = args.gene_fc_file
score_out1 = args.score_out1
score_out2 = args.score_out2

if not score_out1:
    score_out1 = peakfile1+".score_out.txt"
if not score_out2:
    score_out2 = peakfile2+".score_out.txt"

out=args.out
cor_test_out = out+".cor_test.txt"

nthreads = args.nthreads
contact_constant = args.contact_constant
peak_signal_constant = args.peak_signal_constant

if not score_out1:
    score_out1 = promfile1 + ".score.txt"
if not score_out2:
    score_out2 = promfile2 + ".score.txt"

def strlist(list):
    test_list = list.copy()
    test_list = [str(i) for i in test_list]
    return test_list

def overlap(region1, region2, cutoff=0):
    dist = 0
    min_start = min(int(region1[0]), int(region2[0]))
    max_start = max(int(region1[0]), int(region2[0]))
    min_end = min(int(region1[1]), int(region2[1]))
    max_end = max(int(region1[1]), int(region2[1]))

    if (max_start <= min_end + cutoff):
        return 1
    else:
        return 0

def read_matrix(matrixfile):
    matrix = {}
    with open(matrixfile) as h:
        for line in h.readlines():
            line = line.strip().split()
            start1=int(line[0])
            end1=int(line[1])
            value = float(line[2])
            bs1 = int(start1/res)
            be1 = int(end1/res)
            if bs1 not in matrix.keys():
                matrix[bs1] = {}
            matrix[bs1][be1] = value
    '''
    row = np.array([0, 0, 1, 2, 2, 2])
    col = np.array([0, 2, 2, 0, 1, 2])
    data = np.array([1, 2, 3, 4, 5, 6])
    mtx = sparse.coo_matrix((data, (row, col)), shape=(3, 3))
    '''
    return matrix

def read_pickle(file):
    pikd = open(file, 'rb')
    data = pickle.load(pikd)
    pikd.close()
    return data

def get_value_matrix(matrix, bin1, bin2):
    bstart = bin1
    bend = bin2
    if bin1 > bin2:
        bstart = bin2
        bend = bin1
    if bstart in matrix.keys() and bend in matrix[bstart].keys():
        value = matrix[bstart][bend]
    else:
        value = 0
    return value

def get_value_cooler(cooler, chr, bin1, bin2):
    start1 = bin1 * res
    end1 = (bin1+1) * res
    start2 = bin2 * res
    end2 = (bin2+1) * res
    region1 = str(chr) + ":"+str(start1)+"-"+str(end1)
    region2 = str(chr) + ":"+str(start2)+"-"+str(end2)
    value = float(cooler.matrix().fetch(region1, region2))
    if np.isnan(value):
        value = 0
    return value

def get_contact_value(chr,bin1,bin2):
    if matrixfile1:
        value = get_value_matrix(matrix1, bin1, bin2)
    if coolerfile:
        value = get_value_cooler(cooler, chr, bin1, bin2)
    return value

def calc_total_score(promoters, peaks, matrix, score_out, threads, gene_list = []):
    chrlist = promoters.keys()
    #cool = clr.Cooler(coolerfile)
    #scores = {}
    #f = open(score_out, "a")
    with WorkerPool(n_jobs=threads) as pool:
        results = pool.map(calc_score_per_chr, [(promoters, peaks, matrix[chr], chr, window, res, gene_list) for chr in chrlist])

    total_score_list = []
    for score_list in results:
        total_score_list = total_score_list + score_list

    total_score_df = pd.DataFrame(total_score_list, columns =['chrom', 'genename', 'score'])
    total_score_df = total_score_df.groupby('genename').sum()
    total_score_df['genename'] = total_score_df.index
    total_score_df.to_csv(score_out, header=True, index=False, sep = "\t")
    return total_score_df

def calc_score_per_chr(promoters, peaks, matrix, chr, window, res, gene_list):
    return_list = []
    for tss in promoters[chr]:
        start = tss[1]
        end = tss[2]
        genename = tss[3]
        if gene_list and genename not in gene_list:
            continue
        start_minus_window = int((start+end)/2)-window
        if start_minus_window<0:
            start_minus_window = 0
        tss_pr = pr.from_dict({"Chromosome": [chr], "Start": [start_minus_window], "End": [int((start+end)/2+window)]})
        peaks_pr = pr.PyRanges(pd.DataFrame(peaks[chr], columns=["Chromosome", "Start", "End", "name", "score"]))
        op_peaks = peaks_pr.overlap(tss_pr).as_df()

        bin1=int((start+end)/2/res)
        score_list = []
        tss_score = 1
        if tss_intensity:
            tss_score = tss[4]
        for i in range(0, len(op_peaks)):
            bin2 = int((op_peaks['Start'][i]+op_peaks['End'][i])/2/res)
            contact = 1
            peak_signal = 1
            if not contact_constant:
                if matfiletype == ".pickle":
                    contact = get_value_matrix(matrix, bin1, bin2)
            if not peak_signal_constant:
                peak_signal = op_peaks['score'][i]
            ACscore = contact * peak_signal * tss_score
            #print(tss, op_peaks.loc[i], matfiletype, contact, peak_signal, tss_score, ACscore)
            score_list.append(ACscore)
        #tss[4] = sum(score_list)
        return_list.append([chr, tss[3], sum(score_list)])
    print("Calc ", chr, " done", time.ctime())
    return return_list

def read_bw(bwfile):
    bw = pyBigWig.open(bwfile)
    return bw

def get_value_bw(bw, interval):
    #bw is a readin bw from bwfile, interval is a list of [chr, start, end] and 0-based.
    mean_value = bw.stats(interval[0], interval[1], interval[2], exact=True)
    return mean_value[0]

def read_promoter(promfile, bw):
    print("read promoter file", time.ctime())
    genes = {}
    promoters= {}
    with open(promfile) as h:
        for line in h.readlines():
            line = line.strip().split()
            chr=line[0]
            '''
            if chr!=chrom:
                continue
            '''
            start=int(line[1])
            end=int(line[2])
            genename=line[3]

            strand=line[5]
            score=float(line[4])
            if bw:
                score = get_value_bw(bw, [chr,start,end])
            if chr not in promoters.keys():
                promoters[chr] = []
            promoters[chr].append([chr,start,end,genename,score,strand])
            if genename not in genes.keys():
                genes[genename] = []
            genes[genename].append([chr,start,end,genename,score,strand])
    return genes, promoters

##read peakfile and remove enhancer overlap with promoter
def read_peak(peakfile, promoters, bw):
    print("read peak file", time.ctime())
    peaks = {}
    with open(peakfile) as h:
        for line in h.readlines():
            line = line.strip().split()
            chr=line[0]
            if chr not in promoters.keys():
                continue
            '''
            if chr!=chrom:
                continue
            '''
            start=int(line[1])
            end=int(line[2])
            name=line[3]
            score=float(line[4])
            if bw:
                score = get_value_bw(bw, [chr,start,end]) * (end-start)
            if chr not in peaks.keys():
                peaks[chr] = []
            tag = 0
            if remove_promoter:
                for i in promoters[chr]:
                    if i[0]!=chr:
                        continue
                    pstart = i[1]
                    pend = i[2]
                    genename = i[3]
                    tag += overlap([int(pstart),int(pend)],[start,end], 0)
            if tag==0:
                peaks[chr].append([chr,start,end,name,score])
    return peaks

#read promoter and peak intervals and scores
def prepare_tss_peak(promfile, peakfile, bwfile=""):
    bw = []
    if bwfile:
        print("read", bwfile)
        bw = read_bw(bwfile)
    genes, promoters = read_promoter(promfile, bw)
    peaks = read_peak(peakfile, promoters, bw)

    return [promoters, peaks]

gene_list = []
gene_fc_info = []
if gene_fc_file:
    gene_fc_info = pd.read_csv(gene_fc_file, header=None, sep = "\t")
    gene_list = list(gene_fc_info[0])

print("Start prepare promoter and peaks", time.ctime())
with WorkerPool(n_jobs=2) as pool:
    results = pool.map(prepare_tss_peak, [(promfile1, peakfile1, bwfile1), (promfile2, peakfile2, bwfile2)])
print("read genes and peaks done!", time.ctime())
promoters1 = results[0][0]
promoters2 = results[1][0]
peaks1 = results[0][1]
peaks2 = results[1][1]

matrix1=read_pickle(matrixfile1)
print("Start calc score for sample1", time.ctime())
scores_df1 = calc_total_score(promoters1, peaks1, matrix1, score_out1, nthreads, gene_list)
print("calc sample1 done", time.ctime())
matrix1.clear()
matrix2=read_pickle(matrixfile2)
print("Start calc score for sample2", time.ctime())
scores_df2 = calc_total_score(promoters2, peaks2, matrix2, score_out2, nthreads, gene_list)
print("calc score2 done", time.ctime())
matrix2.clear()

if gene_fc_file:
    #outfh = open(out, "a")
    cor_test_outfh = open(cor_test_out, "a")
    gene_fc_list = []
    score_fc_list = []
    fc_list = []
    for i in range(0, len(gene_fc_info)):
        gene = str(gene_fc_info[0][i])
        #if gene not in scores1.keys() or gene_fc[0][i] not in scores2.keys():
        if gene not in scores_df1.index or gene not in scores_df2.index:
            continue
        #score_fc = np.log2(scores_df1.loc[gene]['score'] / scores2[gene])
        #score_fc = np.log2(scores_df1.loc[gene]['score'] / scores_df2.loc[gene]['score'])
        if scores_df1.loc[gene]['score'] == 0 or scores_df2.loc[gene]['score'] == 0:
            continue
        score_fc = np.log2(scores_df1.loc[gene]['score'] / scores_df2.loc[gene]['score'])
        gene_fc = np.log2(gene_fc_info[1][i])
        score_fc_list.append(score_fc)
        gene_fc_list.append(gene_fc)
        fc_list.append([gene, score_fc, gene_fc])

    out_df = pd.DataFrame(fc_list)
    out_df.to_csv(out, header=False, index=False, sep ="\t")

    corr, _ = spearmanr(score_fc_list, gene_fc_list)
    print(corr,"\n")
    cor_test_outfh.write(str(corr)+"\n")

    #outfh.close()
    cor_test_outfh.close()
