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


def fasta_stats(fasta_str):
    fasta_str = fasta_str
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0}
    for key in counts:
        counts[key] = fasta_str.count(key)
    return counts


def get_target_chr_fasta(path, chr_list):
    forward_out = {chr: [] for chr in chr_list}
    backward_out = {chr: [] for chr in chr_list}
    for chr in chr_list:
        p = glob(f'{path}/*_{chr}_forward.txt')
        forward_out[chr].extend(p)
        p = glob(f'{path}/*_{chr}_backward.txt')
        backward_out[chr].extend(p)
    return forward_out, backward_out


def merge_bwd_to_fwd(output, cds_file_names):
    if not os.path.exists(f'{output}/merged/'):
        os.makedirs(f'{output}/merged/')
    for name in tqdm(cds_file_names, desc='Merging'):
        fwd = [file for file in glob(f'{output}/{name}_forward_*.txt') if not file.endswith('cds_idx.txt')]
        bwd = [file for file in glob(f'{output}/{name}_backward_*.txt') if not file.endswith('cds_idx.txt')]
        fwd_idx = glob(f'{output}/{name}_forward_*_cds_idx.txt')
        bwd_idx = glob(f'{output}/{name}_backward_*_cds_idx.txt')

        if fwd_idx:
            fwd_idx_data = open(fwd_idx[0], 'rt').readline()
        else:
            fwd_idx_data = None

        if bwd_idx:
            bwd_idx_data = open(bwd_idx[0], 'rt').readline()
        else:
            bwd_idx_data = None

        if fwd_idx_data is None and bwd_idx_data is None:
            if fwd:
                shutil.copy(fwd[0], fr'{output}/merged/{name}_forward_0.txt')
            else:
                shutil.copy(bwd[0], fr'{output}/merged/{name}_forward_0.txt')
            continue

        if fwd_idx_data and bwd_idx_data:
            union = set(fwd_idx_data.split(',')).union(set(bwd_idx_data.split(',')))
        else:
            if fwd_idx_data:
                union = set(fwd_idx_data.split(','))
            if bwd_idx_data:
                union = set(bwd_idx_data.split(','))
        shutil.copy(fwd[0], f'{output}/merged/{name}_forward_{len(union)}.txt')
        with open(f'{output}/merged/{name}_forward_{len(union)}_cds_idx.txt', 'wt') as file:
            file.writelines(','.join(union))


def main(path, gtf, chrom, length, output, cpus=8, merge=True):
    p = Pool(cpus)
    if not os.path.exists(output):
        os.makedirs(output)
    species = os.path.basename(Path(path).parent)
    chr_list = chrom.split(',')
    forward_chr_file_path, backward_chr_file_path = get_target_chr_fasta(path, chr_list)

    gtf_df = pd.read_csv(gtf, sep='\t', header=None)
    cds_df_forward_list = {
        i: gtf_df[(gtf_df.iloc[:, 2] == 'CDS') & (gtf_df.iloc[:, 0] == i) & (gtf_df.iloc[:, 6] == '+')] for i in
        chr_list}
    cds_df_backward_list = {
        i: gtf_df[(gtf_df.iloc[:, 2] == 'CDS') & (gtf_df.iloc[:, 0] == i) & (gtf_df.iloc[:, 6] == '-')] for i in
        chr_list}

    for i in cds_df_forward_list:
        cds_df_forward_list[i] = cds_df_forward_list[i].copy()
        cds_df_forward_list[i].loc[:, 'cds_pos'] = cds_df_forward_list[i].iloc[:, 8].str.split(';').apply(
            lambda col: col[2].split(' ')[-1].replace('"', ''))
        cds_forward_range = cds_df_forward_list[i].iloc[:, [3, 4, 9]].values
        cds_df_forward_list[i] = cds_forward_range[cds_forward_range[:, 0].argsort()]

    for i in cds_df_backward_list:
        cds_df_backward_list[i] = cds_df_backward_list[i].copy()
        cds_df_backward_list[i].loc[:, 'cds_pos'] = cds_df_backward_list[i].iloc[:, 8].str.split(';').apply(
            lambda col: col[2].split(' ')[-1].replace('"', ''))
        cds_backward_range = cds_df_backward_list[i].iloc[:, [3, 4, 9]].values
        cds_df_backward_list[i] = cds_backward_range[cds_backward_range[:, 0].argsort()]

    results = []
    for chr_ in chr_list:
        result_fwd = p.apply_async(run,
                                   args=(
                                       forward_chr_file_path[chr_][0], cds_df_forward_list[chr_], chr_, length, output,
                                       species,))
        result_bwd = p.apply_async(run,
                                   args=(
                                       backward_chr_file_path[chr_][0], cds_df_backward_list[chr_], chr_, length,
                                       output, species,))
        results.append(result_fwd)
        results.append(result_bwd)
    p.close()
    p.join()

    for res in results:
        res.wait()

    cds_file_names = set()
    for r in results:
        cds_file_names = cds_file_names.union(r.get())
    if merge:
        merge_bwd_to_fwd(output, cds_file_names)


def run(chr_file_path, cds_range, chro, length, output, species):
    try:
        if 'forward' in chr_file_path:
            type_ = 'forward'
        else:
            type_ = 'backward'
        f = open(chr_file_path, 'rt')
        fasta = f.readlines()[0]
        fasta_counts = fasta_stats(fasta)
        print(chro, fasta_counts, flush=True)
        np.random.seed(42)
        sample_size = int(len(fasta) // length * 1.2)
        # samples = np.random.choice(len(fasta) - length, size=sample_size, replace=False)
        # samples = np.sort(samples)
        samples = np.array([i for i in range(3000, len(fasta), len(fasta) // sample_size)][:sample_size])
        # print(chr_path, samples)
        total = {0: 0, 1: 0}
        cds_file_names = set()
        for s in tqdm(samples, desc=f'{type_} {chr_file_path}', leave=True, file=sys.stdout):
            if samples.tolist().index(s) == len(samples) - 1:
                print(chro, type_, total[0], total[1], flush=True)

            if fasta[s:s + length].count('N') >= 1:
                continue

            nums = 0
            cds_idx = []
            cds_pos = {'initial': 0, 'internal': 0, 'terminal': 0}
            start_indices = np.searchsorted(cds_range[:, 0], np.arange(s, s + length), side='left')
            # print(chr_path, start_indices)
            for s_, idx in zip(range(s, s + length), start_indices):
                if np.any((cds_range[:idx, 0] <= s_) & (cds_range[:idx, 1] >= s_)):
                    nums += 1
                    cds_idx.append(str(s_ - s))
                    cds_index = ((cds_range[:idx, 0] <= s_) & (cds_range[:idx, 1] >= s_)).nonzero()[0]
                    cds_info = cds_range[cds_index, 2][0]
                    cds_pos[cds_info] += 1

            if nums == 0:
                total[0] += 1
            else:
                total[1] += 1

            with open(f'{output}/{species}_{chro}_{s}_{length}_{type_}_{nums}.txt', 'wt') as f:
                f.write(fasta[s:s + length])

            if nums > 0:
                cds_idx = ','.join(cds_idx)
                cds_pox = ','.join([key + ':' + str(cds_pos[key]) for key in cds_pos])
                with open(f'{output}/{species}_{chro}_{s}_{length}_{type_}_{nums}_cds_idx.txt', 'wt') as f:
                    f.write(cds_idx + '\t' + cds_pox)

            cds_file_names.add(f'{species}_{chro}_{s}_{length}')
            # break
        f.close()

        return cds_file_names

    except Exception as e:
        print(e, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', help='file path', default=r'D:\data\human\groups\Homo_sapiens\fasta')
    parser.add_argument('-g', '--gtf', help='longest gtf file',
                        default=r"D:\data\human\groups\Homo_sapiens\Homo_sapiens.longest.tiberius.gtf")
    parser.add_argument('-c', '--chr', help='chr list',
                        default="NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11")
    # NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10
    # NC_000020.11,NC_000021.9,NC_000022.11
    parser.add_argument('-l', '--length', help='cut length', type=int, default=2000)
    parser.add_argument('-t', '--threads', help='cpus to uses', type=int, default=10)
    parser.add_argument('-m', '--merge', help='merge forward and backward', default=True)
    parser.add_argument('-o', '--output', help='output', default=r"D:\data\human\groups")
    args = parser.parse_args()

    print(args)
    np.setbufsize(2 ** 23)
    main(args.path, args.gtf, args.chr, args.length, args.output, args.threads, args.merge)
