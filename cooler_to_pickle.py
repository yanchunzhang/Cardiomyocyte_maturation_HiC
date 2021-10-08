#!/usr/bin/env python3

import argparse, re
#import pyBigWig
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
import cooler
import os
import pickle
import _pickle as cPickle
import bz2
import time
import mpire
from mpire import WorkerPool

parser = argparse.ArgumentParser()
parser.add_argument('--cooler', dest='cooler', default='', type=str, help='name of bed file')
parser.add_argument('--chrom', dest='chrom', default='', type=str, help='chrom')
parser.add_argument('--out', dest='out', default='', type=str, help='prefix of out pickle file')
parser.add_argument('--res', dest='res', default=5000, type=int, help='resolution of cooler file')
parser.add_argument('--threads', dest='threads', default=1, type=int, help='number of threads')

##main
args = parser.parse_args()

coolfile = args.cooler
chrom = args.chrom
out = args.out
res = args.res
nthreads = args.threads

def dump_intra_cooler_to_csr(coolfile, chrom, res):
    #matrix=[]
    rows=[]
    cols=[]
    data=[]
    cmd="cooler dump --join --balanced -r " + str(chrom) + " " + coolfile
    with os.popen(cmd) as pipe:
        for line in pipe:
            line = line.strip()
            line = line.split("\t")
            #print(line)
            if len(line)==8:
                #matrix.append([line[1],line[4],line[7]])
                rows.append(int(int(line[1])/res))
                cols.append(int(int(line[4])/res))
                value = float(line[7])
                data.append(value)
            else:
                continue
    matrix = csr_matrix((data, (rows, cols)))
    return matrix

def dump_intra_cooler_to_dict(coolfile, chrom, res):
    dict = {}
    cmd="cooler dump --join --balanced -r " + str(chrom) + " " + coolfile
    with os.popen(cmd) as pipe:
        for line in pipe:
            line = line.split()
            if len(line)==8:
                row = int(int(line[1])/res)
                col = int(int(line[4])/res)
                if row not in dict.keys():
                    dict[row] = {}
                dict[row][col] = float(line[7])
            else:
                continue
    return dict

def compressed_pickle(file, data):
    with bz2.BZ2File(file + '.pbz2', 'w') as f:
        cPickle.dump(data, f)

def decompress_pickle(file):
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data

def save_pickle(title, data):
    pikd = open(title + '.pickle', 'wb')
    pickle.dump(data, pikd)
    pikd.close()

def read_pickle(file):
    pikd = open(file, 'rb')
    data = pickle.load(pikd)
    pikd.close()
    return data

def dump_chrlist_cooler(coolfile, chrlist, res, nthreads, out):
    dicts = {}
    t = time.process_time()
    with WorkerPool(n_jobs=nthreads) as pool:
        results = pool.map(dump_intra_cooler_to_dict, [(coolfile, chr, res) for chr in chrlist])
    for i in range(0,len(chrlist)):
        dicts[chrlist[i]] = results[i]
    elapsed_time = time.process_time() - t
    print("dump_intra_cooler_to_dict", elapsed_time)
    out_dict = out
    t = time.process_time()
    try:
        save_pickle(out_dict, dicts)
    except:
        print("Failed to save mat_dict to pickle")
    elapsed_time = time.process_time() - t
    print("save dict to pickle", elapsed_time)
    return dicts

if __name__ == '__main__':
    if chrom:
        chrlist = chrom.split(",")
    else:
        chrlist = cooler.Cooler(coolfile).chromnames
    t = time.process_time()
    dump_chrlist_cooler(coolfile, chrlist, res=res, out=out, nthreads=nthreads)
    elapsed_time = time.process_time() - t
    print("Total elapsed time: ", elapsed_time)

'''
mats = {}
t = time.process_time()
for chr in chrlist:
    mats[chrom] = dump_intra_cooler_to_csr(coolfile, chr, res)
elapsed_time = time.process_time() - t
print("dump_intra_cooler_to_csr", elapsed_time)
t = time.process_time()
try:
    with open(out, 'wb') as fh:
        cPickle.dump(mats, out)
except:
    print("Failed to save csr_dict to pickle")
elapsed_time = time.process_time() - t
print("save csr to pickle", elapsed_time)
'''
