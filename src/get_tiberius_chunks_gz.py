import sys
import shutil
import logging

import numpy as np
import pandas as pd
import argparse
import os
import gzip
from glob import glob
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool
from Bio.Seq import Seq
import pickle
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GetChunks:
    def __init__(self, fasta_file='', pkl_file='', chunksize=199998, overlap=0, out_dir = "./"):
        self.fasta = fasta_file
        self.pkl = pkl_file
        self.chunksize = chunksize
        self.overlap = overlap
        self.out_dir = os.path.join(out_dir, "chunks_" + str(self.chunksize) + "_" + str(self.overlap))

        self.spec = self.pkl.split("/")[-3]
        basename = os.path.basename(self.pkl).split(".pkl")[-2].split("_")
        self.chrom = "_".join(basename[0:2])
        self.strand = basename[-1]

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        if self.fasta:
            self.read_fasta()
        else:
            print("-----Fasta file not found!-----")
        if self.pkl:
            self.read_pkl()
        else:
            print("-----Pickle file of one-hot labels not found!-----")
#        basename = os.path.basename(self.fasta_file).split(".txt")[0].split("_")
#        self.spec = "_".join(basename[0:2]) 
#        self.chrom = "_".join(basename[2:4])        


    def read_fasta(self):
        if self.fasta.endswith(".gz"):
            file = gzip.open(self.fasta, "rt")
        elif self.fasta.endswith(".txt"):
            file = open(self.fasta, "r")
        else:
            logger.error(f"Unsupported file format: {self.fasta}")

        seq = file.read()
        if self.strand == "backward":
            self.seq = Seq(seq).reverse_complement()
        else:
            self.seq = Seq(seq)
        file.close()

    def read_pkl(self):

        if self.pkl.endswith(".gz"):
            file = gzip.open(self.pkl, "rb")
        elif self.pkl.endswith(".pkl"):
            file = open(self.pkl, "rb")
        else:
            logger.error(f"Unsupported file format: {self.pkl}")
        label_matrix = pickle.load(file)
        label_matrix = label_matrix.toarray()
        if self.strand  == "backward":
            self.m = label_matrix[::-1]
        else:
            self.m = label_matrix
        file.close()

    def get_chunks(self):
        num_chunks = (len(self.seq) - self.overlap) // (self.chunksize - self.overlap) + 1
        for i in tqdm(range(num_chunks-1), desc = "Cut chr"):
            chunk = {}
            begin = i * (self.chunksize - self.overlap)
            end = i * (self.chunksize - self.overlap) + self.chunksize
            chunk_seq = self.seq[begin:end]    
            chunk_m = coo_matrix(self.m[begin:end,:])
#            sparse_matrix = coo_matrix(current_annotation['annotation'])
            chunk = {'Species': self.spec, 'chrom': self.chrom, 'seq': str(chunk_seq), \
                              'begin': begin, 'end': end, 'strand': self.strand, 'annotation': chunk_m}
            with open(os.path.join(f'{self.out_dir}',\
                 f"{self.spec}_{self.chrom}_{begin}-{end}_{self.strand}.pkl"), 'wb') as file:
                pickle.dump(chunk, file)

        print(f"{self.strand} of {self.fasta} has been cutted!")

"""        last_chunk = {}
        last_chunksize = (len(self.seq) - self.overlap)%(self.chunksize - self.overlap)
        last_seq = self.seq[-last_chunksize:]
        last_m = coo_matrix(self.m[-last_chunksize:])
        last_chunk = {'Species': self.spec, 'chrom': self.chrom, 'seq': str(last_seq), \
                'last_chunksize': last_chunksize, 'strand': self.strand, 'annotation': last_m}
        with open(os.path.join(f'{self.out_dir}',\
            f"{self.spec}_{self.chrom}_lastchunk_{last_chunksize}_{self.strand}.pkl"), 'wb') as file:
            pickle.dump(last_chunk, file)
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='fasta file')
    parser.add_argument('-pkl', '--pkl', help='pickle file of one-hot matrix of tiberius labels')
    parser.add_argument('-s', '--chunk_size', type=int, default=199998,
                        help='Chunksize to cut genome.')
    parser.add_argument('-over', '--overlap', type=int, default=0,
                        help='Overlap size to cut genome.')
    parser.add_argument('-o', '--out_dir', help='out path', default=r"D:\data\human\groups")
    
    args = parser.parse_args()
#    print(args)

#    out_path = os.path.join(args.out_dir, "chunks_" + str(args.chunk_size) + "_" + str(args.overlap))
    
#    files = os.listdir(out_path) if os.path.exists(out_path) else []
    files = []
    if len(files) > 100:
        print('---Skip cut sequences-------chunks exists!')
    else:
#        cpus = 10
#        p = Pool(cpus)
        chunks = GetChunks(args.fasta, args.pkl, args.chunk_size, args.overlap, args.out_dir)
        chunks.get_chunks()
#        p.close()
#        p.join()





