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

import h5py

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GetChunks:
    def __init__(self, fasta_file='', pkl_file='', out_dir = "./", min_seq_len:int=1000000):
        self.fasta = fasta_file
        self.pkl = pkl_file
        self.min_seq_len = min_seq_len  # Minimum sequence length to save
        # self.chunksize = chunksize
        # self.overlap = overlap
        self.out_dir = os.path.join(out_dir, "chunk_chr.all_h5")

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
        if self.fasta.endswith(".txt.gz"):
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

        self.m = self.m.astype(np.bool_)

    def save_h5(self):
        file_name = f"{self.spec}_{self.chrom}_{self.strand}"
        h5_path = self.out_dir

        if len(self.seq) < self.min_seq_len:
            logger.info(f"seq len: {len(self.seq)} of {self.fasta} is too short, skip saving to h5.")

        if len(self.seq) != self.m.shape[0]:
            raise ValueError(f"Sequence length {len(self.seq)} does not match annotation matrix shape {self.m.shape[0]}.")

        # with open(os.path.join(h5_path, file_name+".bin"), "wb") as seqf:
        #     seqf.write(str(self.seq).encode('ascii'))

        
        with h5py.File(os.path.join(h5_path, file_name+".h5"), "w") as hf:
            hf.create_dataset('seq',
                              data=np.array(list(self.seq), dtype='S1'),
                              chunks=True,  # 启用分块存储
                                compression='gzip',
                                compression_opts=6,
                                shuffle=True
                              )
            hf.create_dataset('anno', 
                            data=self.m,
                            chunks=True,  # 启用分块存储
                            compression='gzip',
                            compression_opts=6,
                            shuffle=True)  # 启用字节洗牌提高压缩率
            
            # 存储元数据
            hf.attrs['spec'] = self.spec
            hf.attrs['chrom'] = self.chrom
            hf.attrs['strand'] = self.strand
        logger.info(f"{self.strand} of {self.fasta} has been cutted!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='fasta file')
    parser.add_argument('-pkl', '--pkl', help='pickle file of one-hot matrix of tiberius labels')
    parser.add_argument('-o', '--out_dir', help='out path', default=r"D:\data\human\groups")
    parser.add_argument('-min', '--min_seq_len', type=int, default=1000000,
                        help='Minimum sequence length to save. Default is 1000000.')
    
    args = parser.parse_args()
    print(args)

#    out_path = os.path.join(args.out_dir, "chunks_" + str(args.chunk_size) + "_" + str(args.overlap))
    
#    files = os.listdir(out_path) if os.path.exists(out_path) else []
    files = []
    if len(files) > 100:
        print('---Skip cut sequences-------chunks exists!')
    else:
#        cpus = 10
#        p = Pool(cpus)
        # chunks = GetChunks(args.fasta, args.pkl, args.out_dir, args.min_seq_len)
        # chunks.save_h5()
        pass
#        p.close()
#        p.join()





