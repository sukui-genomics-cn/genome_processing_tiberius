#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# genome_anno.py: Handles the data structure for a genome annotation file
# ==============================================================

import os
import sys
import csv
import gzip
from Bio.Seq import Seq
import copy


class NotGtfFormat(Exception):
    pass


def assemble_transcript(exons, sequence, strand):
    parts = []
    exons.sort(reverse=strand == '-')
    for exon in exons:
        exon_seq = sequence.seq[exon[0] - 1:exon[1]]
        if strand == '-':
            exon_seq = exon_seq.reverse_complement()
        parts.append(str(exon_seq))  # Convert Seq object to string here

    coding_seq = Seq("".join(parts))
    if len(coding_seq) % 3 == 0:
        prot_seq = coding_seq.translate()
        return coding_seq, prot_seq
        # if prot_seq[-1] == '*':
        #     return coding_seq, prot_seq
    return None, None


# Check for in-frame stop codons
def check_in_frame_stop_codons(seq):
    return '*' in seq[:-1]


class Transcript:
    """
        Class handling the data structures and methods for a transcript
    """

    def __init__(self, id, gene_id, chr, source_anno, strand):
        """
            Args:
                id (str): Transcript ID
                gene_id (str): Gene ID
                chr (str): Chromosome/Sequence name where the transcript is located
                source_anno (str): Anno ID
                strand (str): Strand (+/-) on which the transctipt is located
        """
        self.id = id
        self.chr = chr
        self.gene_id = gene_id
        # self.transcript_lines[segment_type] = [lines of segment type]
        self.transcript_lines = {}
        self.gtf = []
        self.source_anno = source_anno
        self.start = -1
        self.end = -1
        self.cds_len = -1
        self.cds_coords = {}
        self.strand = strand
        self.source_method = ''

    def add_line(self, line):
        """
            Add a single line from the gtf file to the transcript data structure.

            Args:
                line (list): List of all elements of a line from a gtf file
        """
        if not (line[0] == self.chr or line[6] == self.strand):
            raise NotGtfFormat('File is not in gtf format. ' \
                               + 'Error in line {}\n'.format('\t'.join(map(str, line)))
                               + 'Transcript ID is not unique')

        if line[2] not in self.transcript_lines.keys():
            self.transcript_lines.update({line[2]: []})

        self.source_method = line[1]

        line[3] = int(line[3])
        line[4] = int(line[4])
        if self.start < 0 or line[3] < self.start:
            self.start = line[3]
        if self.end < 0 or line[4] > self.end:
            self.end = line[4]
        if self.gene_id == '' and not line[2] == 'transcript':
            self.gene_id = line[8].split('gene_id "')[1].split('";')[0]
        self.transcript_lines[line[2]].append(line)

    def set_gene_id(self, new_gene_id):
        self.gene_id = new_gene_id

    def get_cds_len(self):
        cds = self.get_type_coords('CDS', False)
        return sum([c[1] - c[0] + 1 for c in cds])

    def get_type_coords(self, type, frame=True):
        """
            Get the coordinates and reading frame of the coding regions
            Returns:
                (dict(list(list(int)))): Dictionary with list of type coords for
                                        each each frame phase (0,1,2)
        """
        # returns dict of cds_coords[phase] = [start_coord, end_coord] of all CDS
        if frame:
            coords = {'0': [], '1': [], '2': [], '.': []}
        else:
            coords = []
        if type == 'CDS' and type not in self.transcript_lines.keys():
            type = 'exon'
        if type not in self.transcript_lines.keys():
            return coords
        for line in self.transcript_lines[type]:
            if frame:
                coords[line[7]].append([line[3], line[4]])
            else:
                coords.append([line[3], line[4]])
        if frame:
            for k in coords.keys():
                coords[k].sort(key=lambda c: (c[0], c[1]))
            if type == 'CDS':
                coords['0'] += coords['.']
                del coords['.']
        else:
            coords.sort(key=lambda c: (c[0], c[1]))
        return coords

    def get_cds_coords(self):
        """
            Get the coordinates and reading frame of the coding regions

            Returns:
                (dict(list(list(int)))): Dictionary with list of CDS coords for
                                        each each frame phase (0,1,2)
        """
        # returns dict of cds_coords[phase] = [start_coord, end_coord] of all CDS
        if not self.cds_coords.keys():
            self.cds_coords = {'0': [], '1': [], '2': []}
            if 'CDS' in self.transcript_lines.keys():
                key = 'CDS'
            else:
                key = 'exon'
            for line in self.transcript_lines[key]:
                self.cds_coords[line[7]].append([line[3], line[4]])
            for k in self.cds_coords.keys():
                self.cds_coords[k].sort(key=lambda c: (c[0], c[1]))
        return self.cds_coords

    def add_missing_lines(self, genome, tiberius):
        """
            Add transcript, intron, CDS, exon coordinates if they were not
            included in the gtf file

            Returns:
                (boolean): FALSE if no cds were found for the tx, TRUE otherwise
        """
        # add intron lines
        self.find_introns(genome, tiberius)
        # check if tx has cds or exon
        if not self.check_cds_exons():
            return False
        # add transcript line
        self.find_transcript()
        # add start/stop codon line
        self.find_start_stop_codon(genome)
        self.merge_stop_codon_to_cds()
        # self.add_5p_3p_utr()
        self.cds_to_exon()
        if tiberius:
            self.modify_gene_transcript()
        return True

    def check_cds_exons(self):
        """
            Check if tx has CDS or exons.
        """
        if 'CDS' not in self.transcript_lines.keys() and 'exon' not in self.transcript_lines.keys():
            sys.stderr.write('Skipping transcript {}, no CDS nor exons in {}\n'.format(self.id, self.id))
            return False
        return True

    def modify_gene_transcript(self):
        self.start = self.transcript_lines['CDS'][0][3]
        self.transcript_lines['transcript'][0][3] = self.transcript_lines['CDS'][0][3]

        self.end = self.transcript_lines['CDS'][-1][4]
        self.transcript_lines['transcript'][0][4] = self.transcript_lines['CDS'][-1][4]

    def cds_to_exon(self):
        if 'CDS' not in self.transcript_lines:
            return
        self.transcript_lines['exon'] = copy.deepcopy(self.transcript_lines['CDS'])
        for i in range(len(self.transcript_lines['exon'])):
            self.transcript_lines['exon'][i][2] = 'exon'

    def add_5p_3p_utr(self):
        strand = self.transcript_lines['transcript'][0][6]

        trans_start = self.transcript_lines['transcript'][0][3]
        trans_end = self.transcript_lines['transcript'][0][4]
        if 'CDS' not in self.transcript_lines:
            return
        cds_start = self.transcript_lines['CDS'][0][3]
        cds_end = self.transcript_lines['CDS'][-1][4]

        # filter intron
        intro_left_list = []
        intro_right_list = []
        for intron in self.transcript_lines['intron']:
            if intron[4] < cds_start:
                intro_left_list.append(intron)
            if intron[3] > cds_end:
                intro_right_list.append(intron)

        self.transcript_lines['five_prime_utr'] = []
        self.transcript_lines['three_prime_utr'] = []

        if strand == '+':
            if len(intro_left_list) > 0:
                front = trans_start
                last = cds_start - 1

                for i in range(len(intro_left_list)):
                    intro = intro_left_list[i]
                    five_utr = self.transcript_lines['transcript'][0].copy()
                    five_utr[2] = 'five_prime_utr'
                    five_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                    five_utr[3] = front
                    if intro[4] < cds_start:
                        five_utr[4] = intro[3] - 1
                    else:
                        five_utr[4] = last
                    front = intro[4] + 1
                    self.transcript_lines['five_prime_utr'].append(five_utr)

                    if i == len(intro_left_list) - 1:
                        five_utr = five_utr.copy()
                        five_utr[3] = intro[4] + 1
                        five_utr[4] = last
                        self.transcript_lines['five_prime_utr'].append(five_utr)
            else:
                if trans_start >= cds_start - 1:
                    return
                five_utr = self.transcript_lines['transcript'][0].copy()
                five_utr[2] = 'five_prime_utr'
                five_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                five_utr[3] = trans_start
                five_utr[4] = cds_start - 1
                self.transcript_lines['five_prime_utr'].append(five_utr)

            if len(intro_right_list) > 0:
                begin = cds_end + 1
                end = trans_end
                for i in range(len(intro_right_list)):
                    intro = intro_right_list[i]
                    three_utr = self.transcript_lines['transcript'][0].copy()
                    three_utr[2] = 'three_prime_utr'
                    three_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                    three_utr[3] = begin
                    if intro[3] > cds_end:
                        three_utr[4] = intro[3] - 1
                    else:
                        three_utr[4] = end
                    begin = intro[4] + 1
                    self.transcript_lines['three_prime_utr'].append(three_utr)

                    if i == len(intro_right_list) - 1:
                        three_utr = three_utr.copy()
                        three_utr[3] = intro[4] + 1
                        three_utr[4] = end
                        self.transcript_lines['three_prime_utr'].append(three_utr)
            else:
                if trans_end <= cds_end:
                    return
                three_utr = self.transcript_lines['transcript'][0].copy()
                three_utr[2] = 'three_prime_utr'
                three_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                three_utr[3] = cds_end + 1
                three_utr[4] = trans_end
                self.transcript_lines['three_prime_utr'].append(three_utr)

        else:
            if len(intro_left_list) > 0:
                front = trans_start
                last = cds_start - 1

                for i in range(len(intro_left_list)):
                    intro = intro_left_list[i]
                    five_utr = self.transcript_lines['transcript'][0].copy()
                    five_utr[2] = 'three_prime_utr'
                    five_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                    five_utr[3] = front
                    if intro[4] < cds_start:
                        five_utr[4] = intro[3] - 1
                    else:
                        five_utr[4] = last
                    front = intro[4] + 1
                    self.transcript_lines['three_prime_utr'].append(five_utr)

                    if i == len(intro_left_list) - 1:
                        five_utr = five_utr.copy()
                        five_utr[3] = intro[4] + 1
                        five_utr[4] = last
                        self.transcript_lines['five_prime_utr'].append(five_utr)
            else:
                if trans_start >= cds_start - 1:
                    return
                five_utr = self.transcript_lines['transcript'][0].copy()
                five_utr[2] = 'three_prime_utr'
                five_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                five_utr[3] = trans_start
                five_utr[4] = cds_start - 1
                self.transcript_lines['three_prime_utr'].append(five_utr)

            if len(intro_right_list) > 0:
                begin = cds_end + 1
                end = trans_end
                for i in range(len(intro_right_list)):
                    intro = intro_right_list[i]
                    three_utr = self.transcript_lines['transcript'][0].copy()
                    three_utr[2] = 'five_prime_utr'
                    three_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                    three_utr[3] = begin
                    if intro[3] > cds_end:
                        three_utr[4] = intro[3] - 1
                    else:
                        three_utr[4] = end
                    begin = intro[4] + 1
                    self.transcript_lines['five_prime_utr'].append(three_utr)

                    if i == len(intro_right_list) - 1:
                        three_utr = three_utr.copy()
                        three_utr[3] = intro[4] + 1
                        three_utr[4] = end
                        self.transcript_lines['three_prime_utr'].append(three_utr)
            else:
                if trans_end <= cds_end:
                    return
                three_utr = self.transcript_lines['transcript'][0].copy()
                three_utr[2] = 'five_prime_utr'
                three_utr[-1] = f"gene_id \"{self.gene_id}\"; transcript_id \"{self.id}\";"
                three_utr[3] = cds_end + 1
                three_utr[4] = trans_end
                self.transcript_lines['five_prime_utr'].append(three_utr)

    def merge_stop_codon_to_cds(self):
        # merge stop_codon to last cds
        if 'CDS' in self.transcript_lines:
            cds_list = self.transcript_lines['CDS']
            if 'stop_codon' in self.transcript_lines:
                stop_list = self.transcript_lines['stop_codon']
                for stop in stop_list:
                    if stop[3] > cds_list[-1][4] + 1:
                        continue
                    if stop[6] == '+':
                        cds = cds_list[-1]
                        stop_end = stop[4]
                        cds[4] = stop_end
                        self.transcript_lines['CDS'][-1] = cds
                    else:
                        cds = cds_list[0]
                        stop_start = stop[3]
                        cds[3] = stop_start
                        self.transcript_lines['CDS'][0] = cds

    def find_introns(self, genome, tiberius):
        """
            Add intron lines.
        """
        if not 'intron' in self.transcript_lines.keys():
            self.transcript_lines.update({'intron': []})
            self.transcript_lines.update({'donor': []})
            self.transcript_lines.update({'acceptor': []})
            key = ''
            if 'exon' in self.transcript_lines.keys():
                key = 'exon'
            if tiberius and 'CDS' in self.transcript_lines.keys():
                key = 'CDS'
            if key:
                exon_lst = []
                for line in self.transcript_lines[key]:
                    exon_lst.append(line)
                exon_lst = sorted(exon_lst, key=lambda e: e[3])
                for i in range(1, len(exon_lst)):
                    if exon_lst[i][6] == '+':
                        intron = []
                        intron += exon_lst[i][0:2]
                        intron.append('intron')
                        if exon_lst[i - 1][4] + 1 >= exon_lst[i][3] - 1:
                            continue
                        intron.append(exon_lst[i - 1][4] + 1)
                        intron.append(exon_lst[i][3] - 1)
                        donor = genome[exon_lst[i][0]].seq[exon_lst[i - 1][4]:exon_lst[i - 1][4] + 2]
                        acceptor = genome[exon_lst[i][0]].seq[exon_lst[i][3] - 1 - 2:exon_lst[i][3] - 1]
                        intron += exon_lst[i][5:8]
                        intron.append("gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}_{}\";".format( \
                            self.gene_id, self.id, donor, acceptor))
                        self.transcript_lines['intron'].append(intron)

                        dss = []
                        dss += exon_lst[i][0:2]
                        dss.append('donor')
                        dss.append(exon_lst[i - 1][4] + 1)
                        dss.append(exon_lst[i - 1][4] + 1 + 1)
                        dss += exon_lst[i][5:8]
                        dss.append("gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}\";".format( \
                            self.gene_id, self.id, donor))
                        self.transcript_lines['donor'].append(dss)

                        accptr = []
                        accptr += exon_lst[i][0:2]
                        accptr.append('acceptor')
                        accptr.append(exon_lst[i][3] - 1 - 1)
                        accptr.append(exon_lst[i][3] - 1)
                        accptr += exon_lst[i][5:8]
                        accptr.append("gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}\";".format( \
                            self.gene_id, self.id, acceptor))
                        self.transcript_lines['acceptor'].append(accptr)
                    else:
                        intron = []
                        intron += exon_lst[i][0:2]
                        intron.append('intron')
                        if exon_lst[i - 1][4] + 1 >= exon_lst[i][3] - 1:
                            continue
                        intron.append(exon_lst[i - 1][4] + 1)
                        intron.append(exon_lst[i][3] - 1)
                        acceptor = genome[exon_lst[i][0]].seq[exon_lst[i - 1][4]:exon_lst[i - 1][4] + 2][::-1]
                        donor = genome[exon_lst[i][0]].seq[exon_lst[i][3] - 1 - 2:exon_lst[i][3] - 1][::-1]
                        trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                        acceptor = str(acceptor).translate(trantab)
                        donor = str(donor).translate(trantab)

                        intron += exon_lst[i][5:8]
                        intron.append(
                            "gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}_{}\";".format(self.gene_id, self.id,
                                                                                               donor, acceptor))
                        self.transcript_lines['intron'].append(intron)

                        accptr = []
                        accptr += exon_lst[i][0:2]
                        accptr.append('acceptor')
                        accptr.append(exon_lst[i - 1][4] + 1)
                        accptr.append(exon_lst[i - 1][4] + 1 + 1)
                        accptr += exon_lst[i][5:8]
                        accptr.append(
                            "gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}\";".format(self.gene_id, self.id,
                                                                                            acceptor))
                        self.transcript_lines['acceptor'].append(accptr)

                        dss = []
                        dss += exon_lst[i][0:2]
                        dss.append('donor')
                        dss.append(exon_lst[i][3] - 1 - 1)
                        dss.append(exon_lst[i][3] - 1)
                        dss += exon_lst[i][5:8]
                        dss.append(
                            "gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}\";".format(self.gene_id, self.id,
                                                                                            donor))
                        self.transcript_lines['donor'].append(dss)

    def find_transcript(self):
        """
            Add transcript lines.
        """
        if not 'transcript' in self.transcript_lines.keys():
            for k in self.transcript_lines.keys():
                for line in self.transcript_lines[k]:
                    if line[3] < self.start or self.start < 0:
                        self.start = line[3]
                    if line[4] > self.end:
                        self.end = line[4]
            tx_line = [self.chr, line[1], 'transcript', self.start, self.end, \
                       '.', line[6], '.', self.id]
            self.add_line(tx_line)

    def find_start_stop_codon(self, genome):
        """
            Add start/stop codon lines.
        """

        if not 'start_codon' in self.transcript_lines.keys():
            self.transcript_lines.update({'start_codon': []})
        else:
            # print(self.transcript_lines['start_codon'])
            start_seq = ''
            for idx in range(len(self.transcript_lines['start_codon'])):
                i = self.transcript_lines['start_codon'][idx]
                start_seq += genome[i[0]].seq[i[3] - 1:i[4] - 1 + 1]
                if i[6] == '-':
                    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                    start_seq = str(start_seq)[::-1].translate(trantab)
                self.transcript_lines['start_codon'][idx][-1] += f' site_seq "{start_seq}";'

        if not 'stop_codon' in self.transcript_lines.keys():
            self.transcript_lines.update({'stop_codon': []})
        else:
            # print(self.transcript_lines['stop_codon'])
            stop_seq = ''
            for idx in range(len(self.transcript_lines['stop_codon'])):
                i = self.transcript_lines['stop_codon'][idx]
                stop_seq += genome[i[0]].seq[i[3] - 1:i[4] - 1 + 1]
                if i[6] == '-':
                    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
                    stop_seq = str(stop_seq)[::-1].translate(trantab)
                self.transcript_lines['stop_codon'][idx][-1] += f' site_seq "{stop_seq}";'

        key = ''
        if 'CDS' in self.transcript_lines.keys():
            key = 'CDS'
        elif 'exon' in self.transcript_lines.keys():
            key = 'exon'

        if key:
            self.transcript_lines[key].sort(key=lambda x: x[3])
            tx = self.transcript_lines[key][0]
            line_1_seq = genome[tx[0]].seq[tx[3]:tx[3] + 3]
            line1 = [self.chr, tx[1], '', tx[3], tx[3] + 2, \
                     '.', self.strand, '0', "gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}\";".format( \
                    self.gene_id, self.id, line_1_seq)]
            tx = self.transcript_lines[key][-1]
            line_2_seq = genome[tx[0]].seq[tx[4] - 2:tx[4] + 1]
            line2 = [self.chr, tx[1], '', tx[4] - 2, tx[4], \
                     '.', self.strand, '0', "gene_id \"{}\"; transcript_id \"{}\"; site_seq \"{}\";".format( \
                    self.gene_id, self.id, line_2_seq)]

            fragmented_transcript = True
            if tx[6] == '+':
                line1[2] = 'start_codon'
                line2[2] = 'stop_codon'
                if self.transcript_lines[key][0][7] == 0:
                    fragmented_transcript = False
                start = line1
                stop = line2
            else:
                line1[2] = 'stop_codon'
                line2[2] = 'start_codon'
                if self.transcript_lines[key][-1][7] == 0:
                    fragmented_transcript = False
                stop = line1
                start = line2
            if not 'start_codon' in self.transcript_lines.keys() and not fragmented_transcript:
                if not fragmented_transcript:
                    self.add_line(start)
                else:
                    self.transcript_lines.update({'start_codon': []})
            if not 'stop_codon' in self.transcript_lines.keys():
                self.add_line(stop)

    def get_gtf(self, prefix=''):
        """
            Creates gtf output for the transcript.

            Returns:
                (list(list(str))): List of lines in gtf format as lists
        """
        gtf = []
        if prefix:
            prefix += '.'
        tx_line = []
        for k in self.transcript_lines.keys():
            total = len(self.transcript_lines[k])
            for i, g in enumerate(self.transcript_lines[k]):
                count = f'{i + 1}_{total}'
                site_seq = g[-1].split(';')[-2].split('"')[-2]
                if k == 'start_codon':
                    g[8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; site_seq \"{site_seq}";'
                elif k == 'stop_codon':
                    g[8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; site_seq \"{site_seq}";'
                elif k == 'exon':
                    cds_type = 'internal'
                    if len(self.transcript_lines[k]) == 1:
                        cds_type = 'initial'
                    elif (i == 0 and self.strand == '+') or (
                            i == len(self.transcript_lines[k]) - 1 and self.strand == '-'):
                        cds_type = 'initial'
                    elif (i == len(self.transcript_lines[k]) - 1 and self.strand == '+') or (
                            i == 0 and self.strand == '-'):
                        cds_type = 'terminal'
                    g[
                        8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; cds_type \"{cds_type}"; count \"{count}";'
                elif k == 'donor':
                    g[
                        8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; site_seq \"{site_seq}"; count \"{count}";'
                elif k == 'acceptor':
                    g[
                        8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; site_seq \"{site_seq}"; count \"{count}";'

                elif k == 'intron':
                    cds_type = 'internal'
                    if len(self.transcript_lines[k]) == 1:
                        cds_type = 'initial'
                    elif (i == 0 and self.strand == '+') or (
                            i == len(self.transcript_lines[k]) - 1 and self.strand == '-'):
                        cds_type = 'initial'
                    elif (i == len(self.transcript_lines[k]) - 1 and self.strand == '+') or (
                            i == 0 and self.strand == '-'):
                        cds_type = 'terminal'
                    g[
                        8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; cds_type \"{cds_type}"; count \"{count}";'
                elif k == 'transcript':
                    tx_line = g
                    tx_line[8] = prefix + self.id
                    continue
                elif k == 'CDS':
                    cds_type = 'internal'
                    if len(self.transcript_lines[k]) == 1:
                        cds_type = 'initial'
                    elif (i == 0 and self.strand == '+') or (
                            i == len(self.transcript_lines[k]) - 1 and self.strand == '-'):
                        cds_type = 'initial'
                    elif (i == len(self.transcript_lines[k]) - 1 and self.strand == '+') or (
                            i == 0 and self.strand == '-'):
                        cds_type = 'terminal'
                    g[
                        8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}"; cds_type \"{cds_type}"; count \"{count}";'
                else:
                    g[8] = f'transcript_id \"{prefix + self.id}\"; gene_id \"{self.gene_id}";'
                gtf.append(g)

        if not 'exon' in self.transcript_lines.keys():
            for g in self.transcript_lines['CDS']:
                gtf.append(g[:2] + ['exon'] + g[3:])

        gtf = sorted(gtf, key=lambda g: (g[3], g[4]))
        if tx_line:
            gtf = [tx_line] + gtf
        return gtf


class Anno:
    """
        Class handling the data structures and methods for a one genome annotation file
    """

    def __init__(self, path, id):
        """
            Args:
                path (str): Path to the annotation/gene prediction file in gtf format.
                id (str): Annotation ID
        """
        self.id = id
        self.genes = {'None': []}
        self.gene_gtf = {}
        self.transcripts = {}
        self.path = path
        self.translation_tab = []

    def addGtf(self):
        """
            Read a gtf file and create a dictionary of Transcript objects for
            all transcript in the file
        """
        open_fn = gzip.open if self.path.endswith(".gz") else open
        with open_fn(self.path, 'rt') as file:
            file_lines = csv.reader(file, delimiter='\t')
            for line in file_lines:
                line = [l.strip(' ') for l in line]
                if line[0][0] == '#':
                    continue
                line[3] = int(line[3])
                line[4] = int(line[4])
                if line[2] == 'gene':
                    gene_id = line[8].replace('"', '').split(';')[0].split(' ')[-1]
                    self.genes_update(gene_id)
                    if not gene_id in self.gene_gtf.keys():
                        self.gene_gtf.update({gene_id: line})
                    else:
                        sys.stderr.write('ERROR, gene_id not unique: {}\n'.format(gene_id))
                elif line[2] == 'transcript':
                    transcript_id = line[8].replace('"', '').split(';')[1].split(' ')[-1]
                    gene_id = ''
                    self.transcript_update(transcript_id, gene_id, line[0], line[6])
                    self.transcripts[transcript_id].add_line(line)
                else:
                    transcript_id = line[8].split('transcript_id "')
                    if len(transcript_id) > 1:
                        transcript_id = transcript_id[1].split('";')[0]
                    else:
                        raise NotGtfFormat('File: "{}" is not in gtf format. \n'.format( \
                            self.path) + 'Error in line {}\n'.format('\t'.join(map(str, line))))

                    gene_id = line[8].split('gene_id "')
                    if len(gene_id) > 1:
                        gene_id = gene_id[1].split('";')[0]
                    else:
                        gene_id = 'None'
                        for key, value in self.genes.items():
                            if value == transcript_id:
                                gene_id = key

                    self.transcript_update(transcript_id, gene_id, line[0], line[6])
                    self.genes_update(gene_id, transcript_id)
                    self.transcripts[transcript_id].add_line(line)

        for tx_id in self.genes['None']:
            gene_id = tx_id + '_g'
            self.genes_update(gene_id, tx_id)

    def norm_tx_format(self, genome, filter_inframestop, filter_short, tiberius):
        """
            Add to all Transcript objects transcript, intron, CDS, exon
            coordinates if they were not included in the gtf file.
            Delete all transripts that have no exons or CDS
        """
        for tx_id, tx in self.transcripts.copy().items():
            if 'CDS' not in tx.transcript_lines.keys():
                del self.transcripts[tx_id]
                continue
            exons = tx.get_type_coords('CDS', frame=False)
            # filter out tx with inframe stop codons
            if filter_inframestop:
                coding_seq, prot_seq = assemble_transcript(exons, genome[tx.chr], tx.strand)
                if not coding_seq or check_in_frame_stop_codons(prot_seq):
                    del self.transcripts[tx_id]
            # filter out transcripts with cds len shorter than args.filter_short
            if tx_id in self.transcripts and tx.get_cds_len() < filter_short:
                del self.transcripts[tx_id]

        tx_no_cds = []
        # add missing lines to all tx
        for k in self.transcripts.keys():
            # if len(self.transcripts[k]['transcript_lines']) <= 4:
            #     del self.transcripts[k]
            if not self.transcripts[k].add_missing_lines(genome, tiberius):
                tx_no_cds.append(k)
        for k in tx_no_cds:
            del self.transcripts[k]

    def genes_update(self, gene_id, transcript_id=''):
        """
            Update gene ID dict.
            Args:
                gene_id (str): Gene ID
                transcript_id (str): Transcript ID
        """
        # update gene ids
        if not gene_id in self.genes.keys():
            self.genes.update({gene_id: []})
        if transcript_id and transcript_id not in self.genes[gene_id]:
            self.genes[gene_id].append(transcript_id)
        if transcript_id in self.genes['None'] and not gene_id == 'None':
            self.genes['None'].remove(transcript_id)
            self.transcripts[transcript_id].gene_id = gene_id

    def transcript_update(self, t_id, g_id, chr, strand):
        """
            Update transcript ID dict.
            Args:
                t_id (str): Transcript ID
                g_id (str): Gene ID
                chr (str): Chromosome name
                strand (str): Strand (+/-)
        """
        if not t_id in self.transcripts.keys():
            self.transcripts.update({t_id: Transcript(t_id, g_id, chr, self.id, strand)})

    def find_genes(self):
        """
            Find all genes in the annotation and find the transcripts that
            belong to each gene. Also, cretae a dict with the gtf lines for each gene.
        """
        self.gene_gtf = {}
        self.genes = {}
        for tx in self.transcripts.values():
            if tx.gene_id in self.genes.keys():
                if not (tx.chr == self.gene_gtf[tx.gene_id][0] and \
                        tx.strand == self.gene_gtf[tx.gene_id][6]):
                    sys.stderr.write('ERROR, gene_id not unique: {}.'.format(tx.gene_id))
                    tx.gene_id = tx.gene_id + '.' + tx.chr + '.' + tx.strand
                    sys.stderr.write(' Adding new gene: {}\n'.format(tx.gene_id))
                else:
                    self.genes[tx.gene_id].append(tx.id)
                    self.gene_gtf[tx.gene_id][3] = min(self.gene_gtf[tx.gene_id][3], \
                                                       tx.start)
                    self.gene_gtf[tx.gene_id][4] = max(self.gene_gtf[tx.gene_id][4], \
                                                       tx.end)
                    continue
            self.genes.update({tx.gene_id: [tx.id]})
            self.gene_gtf.update({tx.gene_id: [tx.chr, tx.source_method, 'gene', \
                                               tx.start, tx.end, '.', tx.strand, '.', tx.gene_id]})

    def get_gtf(self):
        """
            Get annotaion file as gtf list.
            Returns:
                list(list(str)): Gtf file as list of lists
        """
        gtf = []
        gene_gtf = sorted(self.gene_gtf.values(), key=lambda g: (g[0], g[3], g[4]))
        for gene in gene_gtf:
            gtf.append(gene)
            for tx_id in self.genes[gene[8]]:
                gtf += self.transcripts[tx_id].get_gtf()
        return gtf

    def add_transcripts(self, txs, id_prefix=''):
        """
            Adds a dict of transcripts to the transcripts of the annotation.
            Args:
                dict(Transcript()): dictionary of Transcripts added to the annotation
        """
        if not id_prefix:
            self.transcripts.update({txs})
        else:
            for tx in txs.values():
                tx.id = id_prefix + tx.id
                self.transcripts.update({tx.id: tx})

    def get_subset(self, tx_list):
        """
            Get annotaion file for a subset of transcripts.
            Args:
                tx_list (list(str)): List of transcript IDs
            Returns:
                list(list(str)): Gtf file as list of lists
        """
        tx_subset = {}
        for tx in tx_list:
            tx_subset.update({tx: self.transcripts[tx]})
        return tx_subset

    def change_id(self, new_id):
        """
            Change annotation file ID.
        """
        self.id = new_id
        for k in self.transcripts.keys():
            self.transcripts[k].source_anno = self.id

    def get_transcript_list(self):
        """
            Returns:
                (List(Transcript)): List of all transcripts.
        """
        return list(self.transcripts.values())

    def rename_tx_ids(self, prefix=''):
        """
            Renames all tx and genes and returns translation table for old tx id to new tx id.
            Args:
                prefix (string): String added before each tx and gene ID.
            Returns:
                translation_tab (list(str, str)): Translation table for old tx id to new tx id.
        """
        self.translation_tab = []
        gene_numb = 1
        old_gene_gtf = sorted(self.gene_gtf.values(), key=lambda g: (g[0], g[3], g[4]))
        self.gene_gtf = {}
        old_genes = self.genes
        self.genes = {}
        old_txs = self.transcripts
        self.transcripts = {}
        if prefix:
            prefix += '_'
        for gene in old_gene_gtf:
            tx_numb = 1
            old_gene_id = gene[8]
            new_gene_id = "{}g{}".format(prefix, gene_numb)
            gene[8] = new_gene_id
            self.genes.update({new_gene_id: []})
            self.gene_gtf.update({new_gene_id: gene})
            for old_tx_id in old_genes[old_gene_id]:
                new_tx_id = "{}g{}.t{}".format(prefix, gene_numb, tx_numb)
                self.transcripts.update({new_tx_id: old_txs[old_tx_id]})
                self.transcripts[new_tx_id].id = new_tx_id
                self.transcripts[new_tx_id].gene_id = new_gene_id
                self.genes[new_gene_id].append(new_tx_id)
                tx_numb += 1
                self.translation_tab.append([new_tx_id, old_tx_id])
            gene_numb += 1
        return self.translation_tab

    def write_anno(self, out_path):
        """
            Write Annotation in gtf format to out_path.
            Args:
                (str) : path to the output file
        """
        if not os.path.exists(os.path.dirname(out_path)):
            os.makedirs(os.path.dirname(out_path))
        with open(out_path, 'w+') as file:
            out_writer = csv.writer(file, delimiter='\t', quotechar="|", lineterminator='\n')
            for line in self.get_gtf():
                out_writer.writerow(line)
