import pandas as pd
import numpy as np
import statistics
import pickle
import argparse
import os
import gzip
from glob import glob
from tqdm import tqdm
from Bio import SeqIO


def stat_N_repeats(genome_file):

    with gzip.open(genome_file, 'rt') as f:
        genome_sequence = f.read()

    sequence = ''.join([line.strip() for line in genome_sequence.splitlines() if not line.startswith(">")])
    n_count = sequence.count('N')
    total_count = len(sequence)
    n_percentage = (n_count / total_count) * 100
    print(f"'N' 在基因组中的占比: {n_percentage:.2f}%")

    lowercase_atgc_count = sequence.count('a') + sequence.count('t') + sequence.count('g') + sequence.count('c')
    lowercase_atgc_percentage = (lowercase_atgc_count / total_count) * 100
    print(f"小写字母 'a'、't'、'g'、'c' 在基因组中的占比: {lowercase_atgc_percentage:.2f}%")

    sequence = sequence.lower()
    valid_count = sequence.count('a') + sequence.count('t') + sequence.count('g') + sequence.count('c') + sequence.count('n')
    non_valid_count = total_count - valid_count
    print(f"其他字符数目: {non_valid_count}")
    del sequence

#统计N和repeat在CDS中占的比例

def stat_N_repeats_cds(genome_file, gtf_file):
    genome_dict = {}
    with gzip.open(genome_file, 'rt') as f:
        for record in SeqIO.parse(f, "fasta"):
            genome_dict[record.id] = str(record.seq)

    gtf_data = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    cds_data = gtf_data[gtf_data[2] == 'CDS']
    cds_info = cds_data[[0, 3, 4]]
    cds_info.columns = ['chrom', 'start', 'end']
    cds_count = 0
    atgc_cds_count = 0
    N_cds_count = 0
    for index, row in cds_info.iterrows():
        chrom = row['chrom']
        start = row['start'] - 1 
        end = row['end']
        region_sequence = genome_dict[chrom][start:end]  #
        cds_count += len(region_sequence)

        atgc_cds_count += region_sequence.count('a')
        atgc_cds_count += region_sequence.count('t')
        atgc_cds_count += region_sequence.count('g')
        atgc_cds_count += region_sequence.count('c')
        N_cds_count += region_sequence.count('N')
    atgc_cds_percentage = (atgc_cds_count / cds_count) * 100
    N_cds_percentage = (N_cds_count / cds_count) * 100
    print(f"小写'atgc'在编码区（CDS）中的数目及占比: {atgc_cds_count} & {atgc_cds_percentage:.2f}%") #5%
    print(f"'N'在编码区（CDS）中的数目及占比: {N_cds_count} & {N_cds_percentage:.2f}%")


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-f', '--fasta', type=str, default='/home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna.gz')
    arg_parser.add_argument('-gtf', '--gtf', type=str, default='/home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/test/Homo_sapiens/Homo_sapiens.longest.tiberius.gtf')
    args = arg_parser.parse_args()
    stat_N_repeats(args.fasta)
    stat_N_repeats_cds(args.fasta, args.gtf)


