import argparse

import numpy as np
import pandas as pd
import os
import sys
import random
import gzip
import pickle
from tqdm import tqdm
from glob import glob
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
from Bio import SeqIO

class GenomeSequences:
    def __init__(self, fasta_file='', gtf_file = ''):
        """Initialize the GenomeSequences object.

        Arguments:
            fasta_file (str): Path to the FASTA file containing genome sequences.
            np_file (str): Path to the numpy file containing one-hot encoded sequences.
            chunksize (int): Size of each chunk for splitting sequences.
            overlap (int): Overlap size between consecutive chunks.
        """
        self.fasta_file = fasta_file
        self.gtf_file = gtf_file
        self.sequences = []
        self.sequence_names = []
        self.anno_seq_names = []
        self.chr_names = []
        if self.gtf_file:
            self.get_gtf_key()

        if self.fasta_file:
            self.read_fasta()
        # self.encode_sequences()

    def get_gtf_key(self):
        df = pd.read_csv(self.gtf_file, sep = "\t", comment = "#", header = None)
        self.anno_seq_names = df[0].unique().tolist()

    def read_fasta(self):
        """Read genome sequences from the specified FASTA file.
        """
        fasta_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(self.fasta_file, "rt"), "fasta"))
        chr_dict = {key: record for key, record in fasta_dict.items() \
                    if record.description.split(",")[0].split(" ")[-2].lower() == "chromosome"} 

        scaffold_dict = {key: record for key, record in fasta_dict.items() \
                         if "scaffold" in record.description.split(",")[0].split(" ")[-1].lower()}

#        chr_dict = {key: record for key, record in fasta_dict.items() if "chromosome" in record.description \
#                                                                        and "scaffold" not in record.description}
        #" ".join(fasta_dict["NC_010449.5"].description.split(",")[0].split(" ")[-2:])
#        scaffold_dict = {key: record for key, record in fasta_dict.items() if "scaffold" in record.description}
#        contig_dict = {key: record for key, record in fasta_dict.items() if "Contig" in record.description}
        seq_dict = chr_dict if chr_dict else scaffold_dict 
        anno_seq_dict = {key: value for key, value in seq_dict.items() if key in self.anno_seq_names}

        for key, record in anno_seq_dict.items():
            self.sequence_names.append(key)
            self.sequences.append(str(record.seq))
            self.chr_names.append(" ".join(record.description.split(",")[0].split(" ")[-4:]))

#    def filter_seq_by_gtf(self):
#        df = pd.read_csv(self.gtf_file, sep = "\t", comment = "#", header = None)
#        self.sequence_names = list(set(self.sequence_names) & set(df[0]))
    def save_fasta(self, args, spec):
        if not os.path.exists(rf'{args.output}/{spec}/fasta'):
            os.makedirs(rf'{args.output}/{spec}/fasta')

        if not os.path.exists(rf'{args.output}/{spec}/fasta_info'):
            os.makedirs(rf'{args.output}/{spec}/fasta_info')

        with open(fr'{args.output}/{spec}/fasta_info/fasta_info.txt', 'w') as chr_file:
            for chrom, chr_name in zip(self.sequence_names, self.chr_names):
                merged_output = chrom + ": " + chr_name
                chr_file.write(f"{merged_output}\n")

        for strand in ['forward', 'backward']:
            for chrom, seq in zip(self.sequence_names, self.sequences):
                if strand == 'backward':
                    continue
                with open(fr'{args.output}/{spec}/fasta/{spec}_{chrom}_{strand}.txt', 'w') as f:
                    f.write(seq)
       

# ==============================================================
# Authors: Lars Gabriel
#
# Class handling GTF information to generate trainings examples
# ==============================================================

import numpy as np
import sys


class GeneStructure:
    """Handles gene structure information from a gtf file,
    prepares one-hot encoded trainings examples"""

    def __init__(self, filename='', spec='', sequence_names=[]):
        """Initialize GeneStructure.

        Arguments:
            filename (str): Path to GTF file.
            np_file (str): Path to save/load numpy array.
            chunksize (int): Size of each chunk.
            overlap (int): Overlapping bases in chunks."""
        self.filename = filename
        self.spec = spec
        self.sequence_names = sequence_names
        self.gene_structures = []

        # one hot encoding for each sequence (intergenic [0], CDS [1], intron [2])
        self.one_hot = None
        # one hot encoding for phase of CDS (0 [0], 1 [1], 2 [2], non coding [3])
        self.one_hot_phase = None

        # chunks of one hot encoded numpy array fitted to genomic chunks
        self.chunks = None
        self.chunks_phase = None

        if self.filename:
            self.read_gtf()


    def read_gtf(self):
        """Read gene structure information from a GTF file.

        Arguments:
            filename (str): Path to GTF file."""
        with open(self.filename, 'rt') as f:
            for line in f:
                # Skip comments and header lines
                if line.startswith('#') or line.startswith('track'):
                    continue

                # Parse the GTF line
                line_parts = line.strip().split('\t')

                # Extract the gene structure information
                chromosome = line_parts[0]
                if chromosome not in self.sequence_names:
                    continue
                feature = line_parts[2]
                strand = line_parts[6]
                phase = line_parts[7]

                start = int(line_parts[3])
                end = int(line_parts[4])

                info = line_parts[8]

                # Store the gene structure information
                self.gene_structures.append((chromosome, feature, strand, phase, start, end, info))

        # Sort the gene structures by chromosome and end position
        self.gene_structures.sort(key=lambda x: (x[0], x[5]))



    def save_chr(self, args):
        if not os.path.exists(f'{args.output}/{self.spec}/anno_{args.type}'):
            os.makedirs(f'{args.output}/{self.spec}/anno_{args.type}')
        for strand in tqdm(self.one_hot, desc='save chr'):
            if strand == '+':
                strd = 'forward'
            if strand == '-':
                strd = 'backward'
            for chr in self.one_hot[strand]:
                sparse_matrix = coo_matrix(self.one_hot[strand][chr])
                with open(os.path.join(f'{args.output}/{self.spec}/anno_{args.type}', f'{chr}_{strd}.pkl'),
                          'wb') as file:
                    pickle.dump(sparse_matrix, file)


    def translate_to_one_hot_hmm_tiberius(self, sequence_lengths, args):
        """Translate gene structure one_hot_matrixormation to one-hot encoding.
            7 classes IR, I0, I1, I2, E0, E1, E2

        Arguments:
            sequence_names (list): Names of sequences.
            sequence_lengths (list): Lengths of sequences."""

        self.one_hot = {}

        numb_labels = 7
        if args.transition:
            numb_labels = 15

        # Initialize a numpy array to store the one-hot encoded positions
        for strand in ['+', '-']:
            self.one_hot[strand] = {seq: np.zeros((seq_l, numb_labels), dtype=np.int8) \
                                    for seq, seq_l in zip(self.sequence_names, sequence_lengths)}
            for seq in self.sequence_names:
                self.one_hot[strand][seq][:, 0] = 1

                # Set the one-hot encoded positions for each gene structure
        for chromosome, feature, strand, phase, start, end, info in self.gene_structures:
            if feature == 'CDS':
                exon_start = (3 - int(phase)) % 3
                self.one_hot[strand][chromosome][start - 1:end, 0] = 0

                one_help = (np.linspace(0, end - start, end - start + 1) + exon_start) % 3
                if strand == '-':
                    one_help = one_help[::-1]
                self.one_hot[strand][chromosome][start - 1:end, 4:7] = \
                    np.eye(3)[one_help.astype(int)]

        for chromosome, feature, strand, phase, start, end, info in self.gene_structures:
            if feature == 'intron':
                if strand == '+':
                    idx = start - 2
                    if args.transition:
                        idx = start - 3
                    exon_strand = np.argmax(self.one_hot[strand][chromosome][idx]) - 4
                else:
                    idx = end
                    if args.transition:
                        idx = end + 1
                    exon_strand = np.argmax(self.one_hot[strand][chromosome][idx]) - 4
                self.one_hot[strand][chromosome][start - 1:end, 1 + exon_strand] = 1
                self.one_hot[strand][chromosome][start - 1:end, 0] = 0

        if args.transition:
            def calculate_index(array, position, default, offset, condition):
                if condition:
                    return default
                else:
                    return np.argmax(array[position, :4]) + offset

            def update_one_hot(self, strand, chromosome, position, index):
                self.one_hot[strand][chromosome][position] = 0
                self.one_hot[strand][chromosome][position, index] = 1

            # states : Ir, I0, I1, I2, E0, E1, E2, START, EI0, EI1, EI2, IE0, IE1, IE2, STOP
            for chromosome, feature, strand, phase, start, end, info in self.gene_structures:
                if feature == 'CDS':
                    if strand == '+':
                        prev_condition = start - 2 < 0
                        prev = calculate_index(self.one_hot[strand][chromosome], start - 2, 7, 10, prev_condition)
                        prev = 7 if prev == 10 else prev

                        end_condition = end >= len(self.one_hot[strand][chromosome])
                        after = calculate_index(self.one_hot[strand][chromosome], end, 14, 7, end_condition)
                        after = 14 if after == 7 else after

                        update_one_hot(self, strand, chromosome, start - 1, prev)
                        update_one_hot(self, strand, chromosome, end - 1, after)
                    else:
                        end_condition = end >= len(self.one_hot[strand][chromosome])
                        prev = calculate_index(self.one_hot[strand][chromosome], end, 7, 10, end_condition)
                        prev = 7 if prev == 10 else prev

                        start_condition = start - 2 < 0
                        after = calculate_index(self.one_hot[strand][chromosome], start - 2, 14, 7, start_condition)
                        after = 14 if after == 7 else after

                        update_one_hot(self, strand, chromosome, end - 1, prev)
                        update_one_hot(self, strand, chromosome, start - 1, after)

        print(f'args.save_one_hot_matrix: {args.save_chr}')
        if args.save_chr:
            self.save_chr(args)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='fasta file',
                        default=r"/home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
    parser.add_argument('-a', '--anno', help='anno gtf file',
                        default=r"/home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz")
    parser.add_argument('-p', '--padding', help='padding', default=False)
    parser.add_argument('-I', '--filter_inframestop', action='store_true', default=True,
                        help='Filter out transcripts with in-frame stop codons.')
    parser.add_argument('-S', '--filter_short', type=int, default=90,
                        help='Filter out transcripts with in-frame stop codons.')
    parser.add_argument('-t', '--type', default='tiberius',
                        help='Same output as in Tiberius')
    parser.add_argument('-ts', '--transition', action='store_true',
                        help='transition')
    parser.add_argument('-sc', '--save_chr', action='store_true',
                        help='Save chromosome')
    parser.add_argument('-o', '--output', help='output', default=r"D:\data\human\groups")
    args = parser.parse_args()

    spec = args.fasta.replace('\\', '/').split('/')[-2]

    from get_longest_isoform.get_longest_isoform import main

    print('----------Start extracting longest isoform-----------')
    # if not os.path.exists(f'{args.output}/{spec}/anno'):
    #     os.makedirs(f'{args.output}/{spec}/anno')
    if args.type.lower() == 'tiberius':
        longest_anno_out = f'{args.output}/{spec}/{spec}.longest.tiberius.gtf'
        if os.path.exists(longest_anno_out):
            print('---Skip get longest isoform-------Longest isoform file exists!')
        else:
            main(args, longest_anno_out, filter_inframestop=args.filter_inframestop,
             filter_short=args.filter_short, quiet=False)
    else:
        longest_anno_out = f'{args.output}/{spec}/{spec}.longest.gtf'
        main(args, longest_anno_out, filter_inframestop=args.filter_inframestop,
             filter_short=args.filter_short, quiet=False)

    fasta_out = os.path.join(args.output, spec, "fasta")
    fasta_files = os.listdir(fasta_out) if os.path.exists(fasta_out) else []
    pkl_out = os.path.join(args.output, spec, "anno_tiberius")
    pkl_files = os.listdir(pkl_out) if os.path.exists(pkl_out) else []

    if len(fasta_files) > 15 and len(pkl_files) > 30:
        print('----SKIP--fasta and pkls exists!---------')
    elif len(fasta_files) > 15 and len(pkl_files) == 0:
        fasta = GenomeSequences(fasta_file=args.fasta, gtf_file=longest_anno_out)
        print('---Skip save sequence-------fasta exists!')
        seqs_length = [len(s) for s in fasta.sequences]
        seq_names = fasta.sequence_names
        del fasta
        ref_anno = GeneStructure(longest_anno_out, spec=spec, sequence_names=seq_names)
        ref_anno.translate_to_one_hot_hmm_tiberius(seqs_length, args)   
    else: 
        fasta = GenomeSequences(fasta_file=args.fasta, gtf_file=longest_anno_out)
        fasta.save_fasta(args, spec)
        seqs_length = [len(s) for s in fasta.sequences]
        seq_names = fasta.sequence_names
        del fasta
        ref_anno = GeneStructure(longest_anno_out, spec=spec, sequence_names=seq_names)
        ref_anno.translate_to_one_hot_hmm_tiberius(seqs_length, args)


