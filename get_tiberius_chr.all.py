import sys
import shutil

import numpy as np
import pandas as pd
import argparse
import os
from glob import glob
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool
from Bio.Seq import Seq
import pickle
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix



class GetChunks:
    def __init__(self, fasta_file='', pkl_file='', out_dir = "./"):
        self.fasta = fasta_file
        self.pkl = pkl_file
        # self.chunksize = chunksize
        # self.overlap = overlap
        self.out_dir = os.path.join(out_dir, "chunk_chr.all")

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
        with open(self.fasta, "r") as file:
            seq = file.read()
        if self.strand == "backward":
            self.seq = Seq(seq).reverse_complement()
        else:
            self.seq = Seq(seq)

    def read_pkl(self):
        with open(self.pkl, "rb") as f:
            label_matrix = pickle.load(f)
            label_matrix = label_matrix.toarray()
        if self.strand  == "backward":
            self.m = label_matrix[::-1]
        else:
            self.m = label_matrix

    def get_chr_all(self):
        file_name = f"{self.spec}_{self.chrom}_{self.strand}"
        with open(os.path.join(self.out_dir, file_name+".bin"), "wb") as seqf:
            seqf.write(str(self.seq).encode('ascii'))

        np.save(os.path.join(self.out_dir, file_name+".npy"), self.m)

#         num_chunks = (len(self.seq) - self.overlap) // (self.chunksize - self.overlap) + 1
#         for i in tqdm(range(num_chunks-1), desc = "Cut chr"):
#             chunk = {}
#             begin = i * (self.chunksize - self.overlap)
#             end = i * (self.chunksize - self.overlap) + self.chunksize
#             chunk_seq = self.seq[begin:end]    
#             chunk_m = coo_matrix(self.m[begin:end,:])
# #            sparse_matrix = coo_matrix(current_annotation['annotation'])
#             chunk = {'Species': self.spec, 'chrom': self.chrom, 'seq': str(chunk_seq), \
#                               'begin': begin, 'end': end, 'strand': self.strand, 'annotation': chunk_m}
#             with open(os.path.join(f'{self.out_dir}',\
#                  f"{self.spec}_{self.chrom}_{begin}-{end}_{self.strand}.pkl"), 'wb') as file:
#                 pickle.dump(chunk, file)

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
        chunks = GetChunks(args.fasta, args.pkl, args.out_dir)
        chunks.get_chr_all()
#        p.close()
#        p.join()





