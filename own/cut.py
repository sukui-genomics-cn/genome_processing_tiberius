import argparse
import logging
import os
import sys
import gzip
import pickle

from glob import glob
from tqdm import tqdm
import pandas as pd
import numpy as np

from hisat2_extract_splice_sites import extract_splice_sites
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
from multiprocessing import Pool
import concurrent.futures


class Fasta(object):
    def __init__(self, fasta_path, species):
        self.fasta_path = fasta_path
        self.species = species
        self.processed_fa = None

    def process_sequence(self, seq, chom):
        # seq_upper = seq.upper()  # 转换为大写
        # seq_filtered = ''.join([char for char in seq_upper if char in 'ATGCU'])  # 过滤非ATGCU字符
        # seq_processed = seq_filtered.replace('U', 'T')  # 将U替换为T
        # return seq_processed
        new_seq = []
        for i in tqdm(seq, desc=f'{self.species} chom: {chom} seq process', position=0, leave=True):
            # if i in 'ATGCU':
            #     if i == 'U':
            #         i = 'T'
            #     new_seq.append(i)
            new_seq.append(i)
        new_seq = ''.join(new_seq) + '\n'
        return new_seq

    def iterate_fasta(self, f):
        chom = None
        fa_chom = dict()
        lines = f.readlines()

        def is_valid_chr(chom):
            return chom[0].isdigit() or chom in ['X', 'Y']

        for i, line in tqdm(enumerate(lines), total=len(lines), position=0, leave=True):
            if not isinstance(line, str):
                line = line.decode()
            line = line.strip()
            if line.startswith('>'):  # 保留唯一标识符行
                new_chom = line.split(' ')[0][1:]
                if is_valid_chr(new_chom):
                    if chom and chr_seq:
                        fa_chom[chom].append(self.process_sequence(''.join(chr_seq), chom))
                    chom = new_chom
                    fa_chom[chom] = [line]
                    chr_seq = []
            else:
                if chom and is_valid_chr(chom):
                    chr_seq.append(line)
        if chom and chr_seq:
            fa_chom[chom].append(self.process_sequence(''.join(chr_seq), chom))
        return fa_chom

    def process(self):
        if self.fasta_path.endswith('.gz'):
            with gzip.open(self.fasta_path, 'r') as f:
                self.processed_fa = self.iterate_fasta(f)
        else:
            with open(self.fasta_path, 'r') as f:
                self.processed_fa = self.iterate_fasta(f)


class CutSeq(object):
    def __init__(self, cut_length, overlap):
        self.seq = None
        self.cut_length = cut_length
        if isinstance(overlap, float):
            self.overlap_length = int(self.cut_length * overlap)
        else:
            self.overlap_length = int(overlap)

    def cut_sequence(self, sequence, cut_length, overlap_length):
        for i in range(0, len(sequence), cut_length - overlap_length):
            if i + cut_length < len(sequence):
                yield (sequence[i:i + cut_length], i + 1, i + cut_length - 1 + 1)
            else:
                yield (sequence[i:], i + 1, len(sequence) - 1 + 1)
                continue

    def get_cut_seq_list(self, sequence):
        self.seq = sequence
        if len(sequence) < self.cut_length:
            return [(sequence, 0, len(sequence) - 1)]
        return list(self.cut_sequence(sequence, self.cut_length, self.overlap_length))


class Anno(object):
    def __init__(self, anno_file_list):
        self.anno_file_list = anno_file_list
        self.anno_dfs_dict = None
        self.anno_types = ['enhancer', 'promoter', 'TF_binding_site',
                           'exon', 'gene', 'mRNA', 'five_prime_UTR', 'CDS', 'three_prime_UTR', 'transcript',
                           'Intron', 'splicing_donor_site', 'splicing_acceptor_site']
        self.anno_types_map = {type: i for i, type in enumerate(self.anno_types)}

    def get_dfs(self):
        self.anno_dfs_dict = {}
        for annofile in self.anno_file_list:
            if not annofile: continue
            if 'regulatory' in annofile:
                self.anno_dfs_dict['regulatory'] = self.get_anno_df(annofile)
            else:
                self.anno_dfs_dict['gene'] = self.get_anno_df(annofile)
        return self.anno_dfs_dict

    def get_anno_df_info(self, anno_df, chrom, begin, end):
        select_df = anno_df[anno_df['chr'] == str(chrom)]
        select_df = select_df[(select_df['start'] >= begin) & (select_df['end'] <= end)]
        select_df = select_df[select_df['type'].isin(self.anno_types)]
        return select_df[['chr', 'type', 'start', 'end', 'strand']]

    def retrieve_ann_from_df(self, seq, selected_df, begin, end, seq_ann_info):
        if selected_df.empty:
            return
        selected_df['start'] -= begin
        selected_df['end'] -= begin
        for type, first, last, strand in selected_df[['type', 'start', 'end', 'strand']].values:
            # print(begin,end)
            if first < 0:
                first = 0
            if strand == '+':
                seq_ann_info[self.anno_types_map[type]][first:last + 1] = 1
            elif strand == '-':
                seq_ann_info[self.anno_types_map[type]][first:last + 1] = 2
            elif strand == '.':
                seq_ann_info[self.anno_types_map[type]][first:last + 1] = 3

        # add donor, acceptor
        for start, end, strand in selected_df[selected_df['type'] == 'Intron'][['start', 'end', 'strand']].values:
            donor = seq[start - 1:start - 1 + 2]
            acceptor = seq[end - 2:end]
            if strand == '-':
                trantab = str.maketrans('ACGT', 'TGCA')
                acceptor = donor[::-1].translate(trantab)
                donor = acceptor[::-1].translate(trantab)
            if (donor == 'GT' and acceptor == 'AG') or (donor == 'GC' and acceptor == 'AG') or (
                    donor == 'AT' and acceptor == 'AC'):
                # print('found!')
                if strand == '+':
                    seq_ann_info[self.anno_types_map['splicing_donor_site']][start:start + 2] = 1
                    seq_ann_info[self.anno_types_map['splicing_acceptor_site']][end - 2:end] = 1
                if strand == '-':
                    seq_ann_info[self.anno_types_map['splicing_donor_site']][end - 2:end] = 1
                    seq_ann_info[self.anno_types_map['splicing_acceptor_site']][start:start + 2] = 1

    # def add_intron_exon(self, seq, selected_df, begin, end, seq_ann_info):
    #     if seq_ann_info is None:
    #         seq_ann_info = [[0 for _ in range(len(seq))] for _ in range(len(self.anno_types))]
    #     selected_df['start'] -= begin
    #     selected_df['end'] -= begin
    #     for start, end, strand in selected_df[selected_df['type'] == 'Intron'][['start', 'end', 'strand']].values:
    #         donor = seq[start:start + 2]
    #         acceptor = seq[end - 2:end]
    #         if strand == '-':
    #             donor = donor[::-1]
    #             acceptor = acceptor[::-1]
    #         if donor == 'GT' and acceptor == 'AG':
    #             # print('found!')
    #             seq_ann_info[self.anno_types_map['splicing_donor_site']][start:start + 2] = 1
    #             seq_ann_info[self.anno_types_map['splicing_acceptor_site']][end - 2:end] = 1
    #         if donor == 'GC' and acceptor == 'AG':
    #             # print('found!')
    #             seq_ann_info[self.anno_types_map['splicing_donor_site']][start:start + 2] = 1
    #             seq_ann_info[self.anno_types_map['splicing_acceptor_site']][end - 2:end] = 1
    #         if donor == 'AT' and acceptor == 'AC':
    #             # print('found!')
    #             seq_ann_info[self.anno_types_map['splicing_donor_site']][start:start + 2] = 1
    #             seq_ann_info[self.anno_types_map['splicing_acceptor_site']][end - 2:end] = 1
    #         seq_ann_info[self.anno_types_map['Intron']][start:end + 1] = 1
    #
    #     for start, end, strand in selected_df[selected_df['type'] == 'exon'][['start', 'end', 'strand']].values:
    #         seq_ann_info[self.anno_types_map['exon']][start:end + 1] = 1
    #     return seq_ann_info

    @staticmethod
    def extract_intron(df):
        exons = {}
        exon_df = df[df['type'] == 'exon']
        for chr, start, end, strand, attributes in exon_df[['chr', 'start', 'end', 'strand', 'attributes']].values:
            if not (chr.isdigit() or chr in ['X', 'Y']):
                continue
            if chr not in exons:
                exons[chr] = {}
            for attr in attributes.split(';'):
                if attr.startswith('Parent='):
                    parent_id = attr.split(':')[1]
                    if parent_id not in exons[chr]:
                        exons[chr][parent_id] = [(start, end, strand)]
                    else:
                        exons[chr][parent_id].append((start, end, strand))

        for chr in exons:
            for parent_id in exons[chr]:
                exons[chr][parent_id].sort()

        introns = {}
        for chr in exons:
            introns[chr] = []
            for parent_id in exons[chr]:
                for i in range(1, len(exons[chr][parent_id])):
                    last = exons[chr][parent_id][i - 1]
                    cur = exons[chr][parent_id][i]
                    if cur[0] - last[1] < 4: continue
                    introns[chr].append((last[1], cur[0], last[2]))

        introns = [{'chr': chr, 'type': 'Intron', 'start': i[0], 'end': i[1], 'strand': i[2]} for chr in introns for i
                   in introns[chr]]
        introns_df = pd.DataFrame(introns)
        df = pd.concat([df, introns_df], ignore_index=True)
        return df

    def get_anno_df(self, anno=r"D:\data\human\Homo_sapiens.GRCh38.112.gff3.gz"):
        if anno.endswith('.gz'):
            df = pd.read_csv(anno,
                             compression='gzip',
                             sep='\t',
                             comment='#',
                             names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
                             dtype={'chr': str, 'source': str, 'type': str, 'start': int, 'end': int, 'score': str,
                                    'strand': str,
                                    'phase': str, 'attributes': str}
                             )
        else:
            df = pd.read_csv(anno,
                             sep='\t',
                             comment='#',
                             names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
                             dtype={'chr': str, 'source': str, 'type': str, 'start': int, 'end': int, 'score': str,
                                    'strand': str,
                                    'phase': str, 'attributes': str}
                             )
        df = df[df['type'].isin(self.anno_types)]
        df = self.extract_intron(df)
        df = df[['chr', 'type', 'start', 'end', 'strand']]
        df = df[(df['chr'].str.isdigit()) | (df['chr'].isin(['X', 'Y']))]
        return df

    def cut_df_by_seq(self, chrom, ann_gff3_dfs_dict, cut_seq_list, spec):
        sliced_df = {}
        for gff3_type in ann_gff3_dfs_dict:
            df = ann_gff3_dfs_dict[gff3_type]
            select_df = df[df['chr'] == str(chrom)]
            select_df = select_df[select_df['type'].isin(self.anno_types)]
            sliced_df[gff3_type] = {}
            for cut_seq, begin, end in tqdm(cut_seq_list, desc=f'{spec} {chrom} cutting df', position=0,
                                            leave=True):
                slice_df = select_df[(select_df['start'] <= end) & (select_df['end'] >= begin)]
                sliced_df[gff3_type][f'{begin}_{end}'] = slice_df
        return sliced_df


def write_files(current_annotation, args, spec, chrom, begin, end):
    with open(os.path.join(f'{args.species_output_path}/{spec}/{args.cut_length}_{args.overlap}',
                           f"{spec}_{chrom}_{begin}-{end}.pkl"),
              'wb') as file:
        pickle.dump(current_annotation, file)


def dump_pkl(species, args, ):
    # print(spec, chrom, cut_seq, begin, end, ann_gff3_dfs_dict, anno_gtf_file, species_output_path)
    try:
        # pool = concurrent.futures.ThreadPoolExecutor(max_workers=10)

        spec = species.replace('\\', '/').split('/')[-2]
        print(spec)
        if not os.path.exists(f'{args.species_output_path}/{spec}/{args.cut_length}_{args.overlap}'):
            os.makedirs(f'{args.species_output_path}/{spec}/{args.cut_length}_{args.overlap}')
        fasta_file = None
        for file in glob(f'{species}/*'):
            if file.endswith(('.fasta', '.fasta.gz', '.fna', '.fna.gz', '.fa', '.fa.gz')):
                fasta_file = file
        fa = Fasta(fasta_file, spec)
        fa.process()

        if fa.processed_fa is None:
            assert fa.processed_fa is None


        cutseq = CutSeq(args.cut_length, args.overlap)
        # cut_seq_ann_dict = dict()
        empty_cut = []
        for chrom in fa.processed_fa:
            seq = fa.processed_fa[chrom]

            chrom_length = int(seq[0].split(':')[-2])
            cut_seq_list = cutseq.get_cut_seq_list(seq[1])
            # pbar = tqdm(total=len(cut_seq_list))
            for cut_seq, begin, end in tqdm(cut_seq_list, desc=f'{spec} {chrom} cut seq processing', position=0,
                                            leave=True):
                current_annotation = {'Species': spec, 'chrom': chrom, 'chrom_length': chrom_length, 'seq': cut_seq,
                                      'begin': begin, 'end': end, 'strand': None}
                # for gff3_type in ann_gff3_dfs_dict:
                #     df = ann_gff3_dfs_dict[gff3_type]
                #     selected_df = ann_gff3.get_anno_df_info(df, chrom, begin, end)
                #     ann_gff3.retrieve_ann_from_df(cut_seq, selected_df, begin, end,
                #                                   current_annotation['annotation'])
                # ann_gff3.add_intron_exon(cut_seq, selected_df, begin, end, current_annotation['annotation'])

                with open(os.path.join(f'{args.species_output_path}/{spec}/{args.cut_length}_{args.overlap}',
                                       f"{spec}_{chrom}_{begin}-{end}.pkl"),
                          'wb') as file:
                    pickle.dump(current_annotation, file)

                # pool.submit(write_files, current_annotation, args, spec, chrom, begin, end)
            with open(os.path.join(f'{args.species_output_path}/{spec}/{args.cut_length}_{args.overlap}',
                                   'empty_cut.txt'), 'w') as file:
                for item in empty_cut:
                    file.write(f"{item}" + "\n")
        # pool.shutdown(wait=True)
    except Exception as e:
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)


def update(*a):
    pbar.update()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-f', '--fasta_path', type=str, default='./groups')
    arg_parser.add_argument('-o', '--species_output_path', type=str, default='../output')
    arg_parser.add_argument('-l', '--cut_length', type=int, default=8192)
    arg_parser.add_argument('-ol', '--overlap', default=0.1)
    arg_parser.add_argument('-c', '--cores', type=int, default=4)
    args = arg_parser.parse_args()
    # fasta = Fasta("./Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa.gz", './output')
    # fasta.process()
    pool = Pool(processes=args.cores)
    for species in glob(rf"{args.fasta_path}/*/"):
        pool.apply_async(dump_pkl, (species, args,))
    pool.close()
    pool.join()
