import argparse

import numpy as np
import os
import sys
import random
import gzip
import pickle
from tqdm import tqdm
from glob import glob
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix


class GenomeSequences:
    def __init__(self, fasta_file='', np_file='', chunksize=20000, overlap=1000):
        """Initialize the GenomeSequences object.

        Arguments:
            fasta_file (str): Path to the FASTA file containing genome sequences.
            np_file (str): Path to the numpy file containing one-hot encoded sequences.
            chunksize (int): Size of each chunk for splitting sequences.
            overlap (int): Overlap size between consecutive chunks.
        """
        self.fasta_file = fasta_file
        self.np_file = np_file
        self.chunksize = chunksize
        self.overlap = overlap

        self.sequences = []
        self.sequence_names = []
        self.one_hot_encoded = None
        self.chunks_one_hot = None
        self.chunks_seq = None
        if self.fasta_file:
            self.read_fasta()
        else:
            self.load_np_array(self.np_file)
        # self.encode_sequences()

    def read_fasta(self):
        """Read genome sequences from the specified FASTA file.
        """

        with gzip.open(self.fasta_file, "rb") as file:
            lines = file.readlines()
            current_sequence = []
            seq_name = None
            for line in tqdm(lines, desc='loading fasta'):
                line = line.decode()
                if line.startswith(">"):
                    new_seq_name = line[1:].strip().split(' ')[0]
                    if any(i in new_seq_name for i in ['NC', 'NW', 'X', 'Y']) or new_seq_name.isdigit():
                        self.sequence_names.append(new_seq_name)
                        if seq_name and current_sequence:
                            self.sequences.append(''.join(current_sequence))
                        current_sequence = []
                        seq_name = new_seq_name
                else:
                    if any(i in seq_name for i in ['NC', 'NW', 'X', 'Y']) or seq_name.isdigit():
                        current_sequence.append(line.strip())
            if any(i in seq_name for i in ['NC', 'NW', 'X', 'Y']) or seq_name.isdigit():
                self.sequences.append(''.join(current_sequence))

    def save_fasta(self, args, spec):
        if not os.path.exists(rf'{args.output}/{spec}/fasta'):
            os.makedirs(rf'{args.output}/{spec}/fasta')
        for strand in ['forward', 'backward']:
            for chrom, seq in zip(self.sequence_names, self.sequences):
                if strand == 'backward':
                    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                    seq = seq[::-1].translate(trantab)
                with open(fr'{args.output}/{spec}/fasta/{spec}_{chrom}_{strand}.txt', 'w') as f:
                    f.write(seq + '\n')

    def get_flat_chunks(self, sequence_name=None, strand='+', coords=False, pad=True):
        """Get flattened chunks of a specific sequence by name.

        Arguments:
            sequence_name (str): Name of the sequence to extract chunks from.
            strand (char): Strand direction ('+' for forward, '-' for reverse).

        Returns:
            chunks_one_hot (np.array): Flattened chunks of the specified sequence.
        """

        if not sequence_name:
            sequence_name = self.sequence_names
        sequences_i = [i for i in self.sequences]

        chunks = []
        chunk_coords = []
        for seq_name, sequence in zip(sequence_name, sequences_i):
            num_chunks = (len(sequence) - self.overlap) \
                         // (self.chunksize - self.overlap) + 1
            if num_chunks > 1:
                chunks += [sequence[i * (self.chunksize - self.overlap): \
                                    i * (self.chunksize - self.overlap) + self.chunksize] \
                           for i in range(num_chunks - 1)]
            if coords:
                num = num_chunks if pad else num_chunks - 1
                chunk_coords += [[
                    seq_name, strand,
                    i * (self.chunksize - self.overlap) + 1,
                    i * (self.chunksize - self.overlap) + self.chunksize] \
                    for i in range(num)]

            last_chunksize = (len(sequence) - self.overlap) % (self.chunksize - self.overlap)
            if pad and last_chunksize > 0:
                padding = np.zeros((self.chunksize, 6), dtype=np.uint8)
                padding[:, 4] = 1
                padding[0:last_chunksize] = sequence[-last_chunksize:]
                chunks.append(padding)

        # chunks = np.array(chunks)
        if strand == '-':
            for i in range(len(chunks)):
                trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                chunks[i] = chunks[i][::-1].translate(trantab)
            chunk_coords.reverse()
        if coords:
            return chunks, chunk_coords
        return chunks


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

    def __init__(self, filename='', np_file='', chunksize=20000, overlap=1000, spec=''):
        """Initialize GeneStructure.

        Arguments:
            filename (str): Path to GTF file.
            np_file (str): Path to save/load numpy array.
            chunksize (int): Size of each chunk.
            overlap (int): Overlapping bases in chunks."""
        self.filename = filename
        self.np_file = np_file
        self.chunksize = chunksize
        self.overlap = overlap
        self.spec = spec
        self.gene_structures = []

        # one hot encoding for each sequence (intergenic [0], CDS [1], intron [2])
        self.one_hot = None
        # one hot encoding for phase of CDS (0 [0], 1 [1], 2 [2], non coding [3])
        self.one_hot_phase = None

        # chunks of one hot encoded numpy array fitted to genomic chunks
        self.chunks = None
        self.chunks_phase = None

        if self.filename:
            self.read_gtf(self.filename)
        else:
            self.load_np_array(self.np_file)

    def read_gtf(self, filename):
        """Read gene structure information from a GTF file.

        Arguments:
            filename (str): Path to GTF file."""
        with open(filename, 'rt') as f:
            for line in f:
                # Skip comments and header lines
                if line.startswith('#') or line.startswith('track'):
                    continue

                # Parse the GTF line
                line_parts = line.strip().split('\t')

                # Extract the gene structure information
                chromosome = line_parts[0]
                if not any(i in chromosome for i in ['NC', 'NW', 'X', 'Y']) or chromosome.isdigit():
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

    def extract_exon_pkl(self, args, chrom_names, seqs, length=2000):
        if not os.path.exists(f'{args.output}/{spec}/exon_{length // 1000}k'):
            os.makedirs(f'{args.output}/{spec}/exon_{length // 1000}k')

        def split_by_length(start, end, length):
            if end - start <= length - 2:
                return [[start, end]]
            else:
                front = split_by_length(start, start + (end - start) // 2, length)
                last = split_by_length(start + (end - start) // 2 + 1, end, length)
                return front + last

        for chromosome, feature, strand, phase, start, end, info in tqdm(self.gene_structures,
                                                                         desc=f'extract exon {length // 1000}k'):
            start -= 1
            end -= 1
            seq_idx = chrom_names.index(chromosome)
            if feature == 'exon':
                locs = split_by_length(start, end, length)
                for loc in locs:
                    length_ = loc[1] - loc[0] + 1
                    if length_ > length:
                        continue
                    amend = (length - length_) // 2
                    front = loc[0] - amend
                    if front < 0:
                        front = 0
                    last = length + front
                    seq = seqs[seq_idx][front:last]
                    strd = 'forward'
                    if strand == "-":
                        strd = 'backward'
                        trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                        seq = seq[::-1].translate(trantab)
                    with open(f'{args.output}/{spec}/exon_{length // 1000}k/{chromosome}_{front}_{last - 1}_{strd}.txt',
                              'wt') as f:
                        f.write(seq + '\n')

    def save_to_file(self, filename):
        """Save one-hot encoding to a numpy file.

        Arguments:
            filename (str): Path to save numpy array."""
        self.np_file = filename
        np.save(filename, self.one_hot)

    def load_np_array(self):
        """Load one-hot encoding from a numpy file."""
        self.one_hot = np.load(self.np_file)

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

    def save_intergenic(self, args, sequence_names, sequences, num_samples=200000):
        def group_consecutive(test_list):
            if not test_list:
                return []

            result = []
            start = test_list[0]  # Start of the current sequence

            for i in range(1, len(test_list)):
                # Check if the current number is not consecutive
                if test_list[i] != test_list[i - 1] + 1:
                    # Append the current sequence to the result
                    result.append((start, test_list[i - 1]))
                    start = test_list[i]  # Update start to the current number

            # Append the last sequence
            result.append((start, test_list[-1]))

            return result

        if not os.path.exists(f'{args.output}/{spec}/intergenic_{args.type}'):
            os.makedirs(f'{args.output}/{spec}/intergenic_{args.type}')

        for strand in self.one_hot:
            for chrom in tqdm(self.one_hot[strand], desc='chrom intergenic'):
                idx_list = np.where(self.one_hot[strand][chrom][:, 0] == 1)[0].tolist()
                if args.type.lower() == '4':
                    channels = [1]
                else:
                    channels = [1, 2, 3]
                for channel in channels:
                    idx_list = np.concatenate((idx_list, np.where(self.one_hot[strand][chrom][:, channel] == 1)[0]))
                idx_list = np.sort(idx_list, axis=None).tolist()
                idx_list = group_consecutive(idx_list)
                for i in idx_list:
                    if i[1] - i[0] < 2000:
                        continue
                    seq_idx = sequence_names.index(chrom)
                    sequence = sequences[seq_idx][i[0]:i[1] + 1]
                    if strand == '-':
                        strd = 'backward'
                        trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                        sequence = sequence[::-1].translate(trantab)
                    else:
                        strd = 'forward'
                    with open(f'{args.output}/{spec}/intergenic_{args.type}/{chrom}_{i[0]}_{i[1]}_{strd}_full.txt',
                              'wt') as f:
                        f.write(sequence + '\n')

                    # if i[1] - i[0] >= 2000:
                    #     step = int((i[1] - i[0] - 2000) // 10) + 2000
                    #     for j in range(i[0], i[1] - 2000, step):
                    #         inter_list.append([j, j + 2000 - 1, strand, chrom])

        for intergentic in glob(f'{args.output}/{spec}/intergenic_{args.type}/*_full.txt'):
            chrom = intergentic.split('/')[-1].split('_')[0]
            strd = intergentic.split('/')[-1].split('_')[3]
            with open(intergentic, 'rt') as f:
                sequence = f.readlines()
                np.random.seed(42)
                size = len(sequence) // 2000
                print(f'intergenic size {size}')
                rand_idx_list = np.random.choice(len(sequence) - 2000, size, replace=False)
                for j in rand_idx_list:
                    sequence = sequence[j:j + 2000]
                    with open(f'{args.output}/{spec}/intergenic_{args.type}/{chrom}_{j}_{j + 2000 - 1}_{strd}.txt',
                              'wt') as f:
                        f.write(sequence + '\n')

    def translate_to_one_hot_4_labels(self, sequence_names, sequence_lengths, sequences, args):
        """Translate gene structure information to one-hot encoding.
            7 classes IR, I0, I1, I2, E0, E1, E2

        Arguments:
            sequence_names (list): Names of sequences.
            sequence_lengths (list): Lengths of sequences."""

        self.one_hot = {}

        numb_labels = 4

        # Initialize a numpy array to store the one-hot encoded positions
        for strand in ['+', '-']:
            self.one_hot[strand] = {seq: np.zeros((seq_l, numb_labels), dtype=np.int8) \
                                    for seq, seq_l in zip(sequence_names, sequence_lengths)}
            for seq in sequence_names:
                self.one_hot[strand][seq][:, 0] = 1

                # Set the one-hot encoded positions for each gene structure
        for chromosome, feature, strand, phase, start, end, info in self.gene_structures:
            start = start - 1
            end = end - 1
            if feature == 'intron':
                self.one_hot[strand][chromosome][start:end + 1, 1] = 1
            elif feature == 'CDS':
                self.one_hot[strand][chromosome][start:end + 1, 2] = 1
            elif feature == 'five_prime_utr':
                self.one_hot[strand][chromosome][start:end + 1, 3] = 1
            elif feature == 'three_prime_utr':
                self.one_hot[strand][chromosome][start:end + 1, 3] = 1
            self.one_hot[strand][chromosome][start:end + 1, 0] = 0

        if args.save_chr:
            self.save_chr(args)

        if args.save_intergenic:
            self.save_intergenic(args, sequence_names, sequences, num_samples=200000)

    def translate_to_one_hot_hmm(self, sequence_names, sequence_lengths, sequences, args):
        """Translate gene structure information to one-hot encoding.
            7 classes IR, I0, I1, I2, E0, E1, E2

        Arguments:
            sequence_names (list): Names of sequences.
            sequence_lengths (list): Lengths of sequences."""

        self.one_hot = {}

        numb_labels = 7
        if args.transition:
            # numb_labels = 17
            numb_labels = 15

        # Initialize a numpy array to store the one-hot encoded positions
        for strand in ['+', '-']:
            self.one_hot[strand] = {seq: np.zeros((seq_l, numb_labels), dtype=np.int8) \
                                    for seq, seq_l in zip(sequence_names, sequence_lengths)}
            for seq in sequence_names:
                self.one_hot[strand][seq][:, 0] = 1

                # Set the one-hot encoded positions for each gene structure
        for chromosome, feature, strand, phase, start, end, info in self.gene_structures:
            start = start - 1
            end = end - 1
            if feature == 'intron':
                if 'initial' in info:
                    self.one_hot[strand][chromosome][start:end + 1, 1] = 1
                if 'internal' in info:
                    self.one_hot[strand][chromosome][start:end + 1, 2] = 1
                if 'terminal' in info:
                    self.one_hot[strand][chromosome][start:end + 1, 3] = 1
            elif feature == 'exon':
                if 'initial' in info:
                    self.one_hot[strand][chromosome][start:end + 1, 4] = 1
                if 'internal' in info:
                    self.one_hot[strand][chromosome][start:end + 1, 5] = 1
                if 'terminal' in info:
                    self.one_hot[strand][chromosome][start:end + 1, 6] = 1
            self.one_hot[strand][chromosome][start:end + 1, 0] = 0

        if args.transition:
            # states : Ir, I0, I1, I2, E0, E1, E2, START, EI0, EI1, EI2, IE0, IE1, IE2, STOP
            for chromosome, feature, strand, phase, start, end, info in self.gene_structures:
                start = start - 1
                end = end - 1

                if feature == 'donor':
                    cur, total = info.split(';')[-2].split('"')[-2].split('_')
                    cur = int(cur)
                    total = int(total)
                    if strand == '+':
                        if int(total) >= 3:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 8] = 1
                            elif cur < total:
                                self.one_hot[strand][chromosome][start:end + 1, 9] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 10] = 1
                        else:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 8] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 10] = 1
                    else:
                        if int(total) >= 3:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 10] = 1
                            elif cur < total:
                                self.one_hot[strand][chromosome][start:end + 1, 9] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 8] = 1
                        else:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 10] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 8] = 1
                elif feature == 'acceptor':
                    cur, total = info.split(';')[-2].split('"')[-2].split('_')
                    cur = int(cur)
                    total = int(total)
                    if strand == '+':
                        if int(total) >= 3:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 11] = 1
                            elif cur < total:
                                self.one_hot[strand][chromosome][start:end + 1, 12] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 13] = 1
                        else:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 11] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 13] = 1
                    else:
                        if int(total) >= 3:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 13] = 1
                            elif cur < total:
                                self.one_hot[strand][chromosome][start:end + 1, 12] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 11] = 1
                        else:
                            if int(cur) == 1:
                                self.one_hot[strand][chromosome][start:end + 1, 13] = 1
                            else:
                                self.one_hot[strand][chromosome][start:end + 1, 11] = 1
                elif feature == 'start_codon':
                    self.one_hot[strand][chromosome][start:end + 1, 7] = 1
                elif feature == 'stop_codon':
                    self.one_hot[strand][chromosome][start:end + 1, 14] = 1
                # elif feature == 'five_prime_utr':
                #     self.one_hot[strand][chromosome][start:end + 1, 15] = 1
                # elif feature == 'three_prime_utr':
                #     self.one_hot[strand][chromosome][start:end + 1, 16] = 1

                self.one_hot[strand][chromosome][start:end + 1, 0] = 0

        if args.save_chr:
            self.save_chr(args)

        if args.save_intergenic:
            self.save_intergenic(args, seq_names, sequences, num_samples=200000)

    def translate_to_one_hot_hmm_tiberius(self, sequence_names, sequence_lengths, sequences, args):
        """Translate gene structure information to one-hot encoding.
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
                                    for seq, seq_l in zip(sequence_names, sequence_lengths)}
            for seq in sequence_names:
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

        print(f'args.save_chr: {args.save_chr}')
        if args.save_chr:
            self.save_chr(args)

        print(f'args.save_intergenic: {args.save_intergenic}')
        if args.save_intergenic:
            self.save_intergenic(args, sequence_names, sequences, 200000)

    def get_chunks_seq(self, seq_names, strand='+'):
        """Get all one-hot encoded chunks.

        Arguments:
            seq_names (list): Names of sequences to chunk.
            strand (str): Strand to process ('+' or '-').

        Returns:
            tuple: Tuple of chunks and phase chunks arrays.
        """
        chunks = []
        chunks_phase = []

        for seq_name in seq_names:
            num_chunks = len(self.one_hot[strand][seq_name]) // self.chunksize + 1

            chunks.extend([self.one_hot[strand][seq_name][i * self.chunksize: \
                                                          (i + 1) * self.chunksize, :] \
                           for i in range(num_chunks)])

            chunks_phase.extend([self.one_hot_phase[strand][seq_name][i * self.chunksize: \
                                                                      (i + 1) * self.chunksize, :] \
                                 for i in range(num_chunks)])

        if strand == '-':
            chunks = [c[::-1] for c in chunks]
            chunks.reverse()
            chunks_phase = [c[::-1] for c in chunks_phase]
            chunks_phase.reverse()

        return chunks, chunks_phase

    def get_flat_chunks_hmm_tiberius(self, seq_names, strand='+', coords=False):
        """Get one-hot encoded chunks, chunks smaller than chunksize are removed.

        Arguments:
            seq_names (list): Names of sequences to chunk.
            strand (str): Strand to process ('+' or '-').
            coords (bool): get coordinates of each chunk

        Returns:
            tuple: One hot encoded chunks of labels
        """
        self.chunks = []
        chunk_coords = []
        for seq_name in seq_names:
            num_chunks = (len(self.one_hot[strand][seq_name]) - self.overlap) \
                         // (self.chunksize - self.overlap) + 1

            if num_chunks - 1 == 0:
                continue
            if coords:
                chunk_coords += [[
                    seq_name, strand,
                    i * (self.chunksize - self.overlap) + 1,
                    i * (self.chunksize - self.overlap) + self.chunksize] \
                    for i in range(num_chunks - 1)]
            self.chunks += [self.one_hot[strand][seq_name][i * (self.chunksize - self.overlap): \
                                                           i * (self.chunksize - self.overlap) + self.chunksize, :] \
                            for i in range(num_chunks - 1)]

        self.chunks = np.array(self.chunks)
        if strand == '-':
            self.chunks = self.chunks[::-1, ::-1, :]
            chunk_coords.reverse()
        if coords:
            return self.chunks, chunk_coords
        return self.chunks


def save_chunks(chunks, coords, args, spec, strand):
    # chr; strand; start; end; seq;
    if not os.path.exists(rf'{args.output}/{spec}/chunks_{args.length}'):
        os.makedirs(rf'{args.output}/{spec}/chunks_{args.length}')
    with open(fr'{args.output}/{spec}/chunks_{args.length}/{spec}_{args.length}_{args.overlap}_{strand}.txt', 'w') as f:
        for chunk, coord in zip(chunks, coords):
            f.write(coord[0] + '\t' + str(coord[2]) + '\t' + str(coord[3]) + '\t' + chunk + '\n')


def save_chunks_hmm_tiberius(chunks, coords, args, spec, strand):
    # chr; strand; start; end; seq;
    if not os.path.exists(rf'{args.output}/{spec}/hmm_tiberius_{args.length}'):
        os.makedirs(rf'{args.output}/{spec}/hmm_tiberius_{args.length}')

    for chunk, coord in zip(chunks, coords):
        with open(fr'{args.output}/{spec}/hmm_tiberius_{args.length}/{spec}_{args.length}_{args.overlap}_{strand}.pkl',
                  'wb') as f:
            out = {'chr': coord[0], 'strand': coord[1], 'start': coord[2], 'end': coord[3], 'anno': chunk}
            pickle.dump(out, f)


def save_chunk_anno(args):
    anno_dict = {'+': {}, '-': {}}
    ann = f'anno_{args.type}'
    for anno_file in tqdm(glob(f'{args.output}/{spec}/{ann}/*.pkl'), desc='loading anno'):
        anno = pickle.load(open(anno_file, 'rb')).toarray()
        if 'forward' in anno_file:
            chr = os.path.basename(anno_file).split('_forward')[0]
            anno_dict['+'][chr] = anno
        else:
            chr = os.path.basename(anno_file).split('_backward')[0]
            anno_dict['-'][chr] = anno

    for chunk_file in glob(f'{args.output}/{spec}/chunks/*.txt'):
        if 'forward' in chunk_file:
            strand = 'forward'
        else:
            strand = 'backward'

        with open(chunk_file, 'r') as f:
            if not os.path.exists(f'{args.output}/{spec}/pickles_{args.type}'):
                os.makedirs(f'{args.output}/{spec}/pickles_{args.type}')
            for i, line in tqdm(enumerate(f), desc='create pickles'):
                chr, start, end, seq = line.strip().split('\t')
                start = int(start) - 1
                end = int(end) - 1
                pkl_dict = {'chr': chr, 'start': int(start), 'end': int(end), 'seq': seq}
                if 'forward' in chunk_file:
                    pkl_dict['anno'] = coo_matrix(anno_dict['+'][chr][int(start):int(end) + 1])
                else:
                    pkl_dict['anno'] = coo_matrix(anno_dict['-'][chr][int(start):int(end) + 1])
                with open(f'{args.output}/{spec}/pickles_{args.type}/{spec}_{chr}_{start}-{end}_{strand}.pkl',
                          'wb') as f:
                    pickle.dump(pkl_dict, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='fasta file',
                        default=r"D:\data\human\DNA_LLM\groups\Homo_sapiens\GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
    parser.add_argument('-a', '--anno', help='anno gtf file',
                        default=r"D:\data\human\DNA_LLM\groups\Homo_sapiens\GCF_000001405.40_GRCh38.p14_genomic.gtf.gz")
    parser.add_argument('-I', '--filter_inframestop', action='store_true', default=True,
                        help='Filter out transcripts with in-frame stop codons.')
    parser.add_argument('-S', '--filter_short', type=int, default=90,
                        help='Filter out transcripts with in-frame stop codons.')
    parser.add_argument('-t', '--type', default='tiberius',
                        help='Same output as in Tiberius')
    parser.add_argument('-o', '--output', help='output', default=r"D:\data\human\groups")
    args = parser.parse_args()
    print(args)

    spec = args.fasta.replace('\\', '/').split('/')[-2]

    from get_longest_isoform.get_longest_isoform import main

    print('----------Start extracting longest isoform-----------')

    if args.type.lower() == 'tiberius':
        longest_anno_out = f'{args.output}/{spec}/{spec}.longest.tiberius.gtf'
        main(args, longest_anno_out, filter_inframestop=args.filter_inframestop,
             filter_short=args.filter_short, quiet=False)
    else:
        longest_anno_out = f'{args.output}/{spec}/{spec}.longest.gtf'
        main(args, longest_anno_out, filter_inframestop=args.filter_inframestop,
             filter_short=args.filter_short, quiet=False)