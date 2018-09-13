"""Utilities for translating ORF detection
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

import pysam
from tqdm import *

from .fasta import FastaReader
from .gtf import GTFReader
from .interval import Interval
from .common import is_read_uniq_mapping
from .common import merge_intervals
from .infer_protocol import infer_protocol


class PutativeORF:
    """Class for putative ORF."""

    def __init__(self, category, transcript_id, transcript_type, gene_id,
                 gene_name, gene_type, chrom, strand, intervals, seq, leader,
                 trailer):
        self.category = category
        self.tid = transcript_id
        self.ttype = transcript_type
        self.gid = gene_id
        self.gname = gene_name
        self.gtype = gene_type
        self.chrom = chrom
        self.strand = strand
        self.intervals = sorted(intervals, key=lambda x: x.start)
        start = intervals[0].start
        end = intervals[-1].end
        self.seq = seq
        self.oid = '{}_{}_{}_{}'.format(transcript_id, start, end, len(seq))
        self.leader = leader
        self.trailer = trailer

    @property
    def start_codon(self):
        if len(self.seq) < 3:
            return None
        return self.seq[:3]

    @classmethod
    def from_tracks(cls, tracks, category, seq='', leader='', trailer=''):
        """
        Parameters
        ----------
        tracks: list of GTFTrack
        """
        if not tracks:
            return None
        intervals = []
        tid = set()
        ttype = set()
        gid = set()
        gname = set()
        gtype = set()
        chrom = set()
        strand = set()
        for track in tracks:
            try:
                tid.add(track.transcript_id)
                ttype.add(track.transcript_type)
                gid.add(track.gene_id)
                gname.add(track.gene_name)
                gtype.add(track.gene_type)
                chrom.add(track.chrom)
                strand.add(track.strand)
                intervals.append(
                    Interval(track.chrom, track.start, track.end,
                             track.strand))
            except AttributeError:
                print('missing attribute {}:{}-{}'.format(
                    track.chrom, track.start, track.end))
                return None
        if (len(tid) != 1 or len(ttype) != 1 or len(gid) != 1
                or len(gname) != 1 or len(gtype) != 1 or len(chrom) != 1
                or len(strand) != 1):
            print('inconsistent tracks for one ORF')
            return None
        tid = list(tid)[0]
        ttype = list(ttype)[0]
        gid = list(gid)[0]
        gname = list(gname)[0]
        gtype = list(gtype)[0]
        chrom = list(chrom)[0]
        strand = list(strand)[0]
        return cls(category, tid, ttype, gid, gname, gtype, chrom, strand,
                   intervals, seq, leader, trailer)


def tracks_to_ivs(tracks):
    chrom = {track.chrom for track in tracks}
    strand = {track.strand for track in tracks}
    if len(chrom) != 1 or len(strand) != 1:
        print('fail to fetch seq: inconsistent chrom or strand')
        return None
    chrom = list(chrom)[0]
    strand = list(strand)[0]
    intervals = [
        Interval(chrom, track.start, track.end, strand) for track in tracks
    ]
    intervals = merge_intervals(intervals)
    return intervals


def transcript_to_genome_iv(start, end, intervals, reverse=False):
    total_len = sum(i.end - i.start + 1 for i in intervals)
    if reverse:
        start, end = total_len - end - 1, total_len - start - 1
    ivs = []
    start_genome = None
    end_genome = None

    ### find start in genome
    cur = 0
    for i in intervals:
        i_len = i.end - i.start + 1
        if cur + i_len > start:
            start_genome = i.start + start - cur
            break
        cur += i_len

    ### find end in genome
    cur = 0
    for i in intervals:
        i_len = i.end - i.start + 1
        if cur + i_len > end:
            end_genome = i.start + end - cur
            break
        cur += i_len

    ### find overlap with (start_genome, end_genome)
    for i in intervals:
        s = max(i.start, start_genome)
        e = min(i.end, end_genome)
        if s <= e:
            ivs.append(Interval(i.chrom, s, e, i.strand))
    return ivs


def fetch_seq(fasta, tracks):
    intervals = tracks_to_ivs(tracks)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    sequences = fasta.query(intervals)
    merged_seq = ''.join(sequences)
    strand = tracks[0].strand
    if strand == '-':
        return fasta.reverse_complement(merged_seq)
    return merged_seq


def search_orfs(fasta, intervals, min_len):
    if not intervals:
        return []

    orfs = []
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    intervals = merge_intervals(intervals)
    sequences = fasta.query(intervals)
    merged_seq = ''.join(sequences)
    reverse = False
    strand = intervals[0].strand
    if strand == '-':
        merged_seq = fasta.reverse_complement(merged_seq)
        reverse = True
    start_codons = ['ATG', 'CTG', 'GTG']
    stop_codons = ['TAG', 'TAA', 'TGA']
    for sc in start_codons:
        cur = 0
        while cur < len(merged_seq):
            start = merged_seq.find(sc, cur)
            if start == -1:
                break
            cur = start + 1
            for i in range(start, len(merged_seq), 3):
                if merged_seq[i:i + 3] in stop_codons:
                    ### found orf
                    size = i - start
                    if size < min_len:
                        break
                    ivs = transcript_to_genome_iv(start, i + 2, intervals,
                                                  reverse)
                    seq = merged_seq[start:i]
                    leader = merged_seq[:start][-99:]
                    trailer = merged_seq[i:][:100]
                    if ivs:
                        orfs.append((ivs, seq, leader, trailer))
                    break
    return orfs


def prepare_orfs(gtf, fasta, prefix, min_len=30):

    if not isinstance(gtf, GTFReader):
        gtf = GTFReader(gtf)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    print('preparing putative ORFs...')

    ### process CDS gtf
    print('searching cds...')
    cds_orfs = []
    for gid in tqdm(gtf.cds):
        for tid in gtf.cds[gid]:
            tracks = gtf.cds[gid][tid]
            # seq = fetch_seq(fasta, tracks)
            orf = PutativeORF.from_tracks(tracks, 'CDS')
            if orf:
                cds_orfs.append(orf)

    ### process UTR gtf
    utr5 = defaultdict(list)
    utr3 = defaultdict(list)
    for gid in gtf.utr:
        ### find first cds and last cds for gene
        gene_cds = []
        for tid in gtf.cds[gid]:
            gene_cds += gtf.cds[gid][tid]
        if not gene_cds:
            print('fail to find CDS for UTR')
            continue
        first_cds = gene_cds[0]
        for gc in gene_cds:
            if gc.start < first_cds.start:
                first_cds = gc
        last_cds = gene_cds[-1]
        for gc in gene_cds:
            if gc.end > last_cds.end:
                last_cds = gc

        for tid in gtf.utr[gid]:
            for track in gtf.utr[gid][tid]:
                if track.start < first_cds.start:
                    if track.end >= first_cds.start:
                        track.end = first_cds.start - 1
                    if track.strand == '+':
                        utr5[tid].append(track)
                    else:
                        utr3[tid].append(track)
                elif track.end > last_cds.end:
                    if track.start <= last_cds.end:
                        track.start = last_cds.end + 1
                    if track.strand == '+':
                        utr3[tid].append(track)
                    else:
                        utr5[tid].append(track)

    uorfs = []
    print('searching uORFs...')
    for tid in tqdm(utr5):
        tracks = utr5[tid]
        ttype = tracks[0].transcript_type
        gid = tracks[0].gene_id
        gname = tracks[0].gene_name
        gtype = tracks[0].gene_type
        chrom = tracks[0].chrom
        strand = tracks[0].strand

        ivs = tracks_to_ivs(tracks)
        orfs = search_orfs(fasta, ivs, min_len)
        for ivs, seq, leader, trailer in orfs:
            orf = PutativeORF('uORF', tid, ttype, gid, gname, gtype, chrom,
                              strand, ivs, seq, leader, trailer)
            uorfs.append(orf)

    dorfs = []
    print('searching dORFs...')
    for tid in tqdm(utr3):
        tracks = utr3[tid]
        ttype = tracks[0].transcript_type
        gid = tracks[0].gene_id
        gname = tracks[0].gene_name
        gtype = tracks[0].gene_type
        chrom = tracks[0].chrom
        strand = tracks[0].strand

        ivs = tracks_to_ivs(tracks)
        orfs = search_orfs(fasta, ivs, min_len)
        for ivs, seq, leader, trailer in orfs:
            orf = PutativeORF('dORF', tid, ttype, gid, gname, gtype, chrom,
                              strand, ivs, seq, leader, trailer)
            dorfs.append(orf)

    ### save to file
    print('saving putative ORFs file...')
    to_write = ('ORF_ID\tORF_type\ttranscript_id\ttranscript_type'
                '\tgene_id\tgene_name\tgene_type\tchrom'
                '\tstrand\tcoordinate\tseq\tleader\ttrailer\n')
    formatter = '{}\t' * 12 + '{}\n'
    for orf in tqdm(cds_orfs + uorfs + dorfs):
        coordinate = ','.join(
            ['{}-{}'.format(iv.start, iv.end) for iv in orf.intervals])
        to_write += formatter.format(orf.oid, orf.category, orf.tid, orf.ttype,
                                     orf.gid, orf.gname, orf.gtype, orf.chrom,
                                     orf.strand, coordinate, orf.seq,
                                     orf.leader, orf.trailer)

    with open('{}_putative_orfs.tsv'.format(prefix), 'w') as output:
        output.write(to_write)

    return (cds, utr5, utr3)


def split_bam(bam, protocol, prefix):
    """Split bam by read length and strand

    Parameters
    ----------
    bam : str
          Path to bam file
    protocol: str
          Experiment protocol [forward, reverse]
    prefix: str
            prefix for output files: {prefix}_xxnt_pos.wig and
            {prefix}__xxnt_neg.wig
    """
    coverages = defaultdict(lambda: defaultdict(Counter))
    qcfail = duplicate = secondary = unmapped = multi = valid = 0
    bam = pysam.AlignmentFile(bam, 'rb')
    total_count = bam.count()
    with tqdm(total=total_count) as pbar:
        for r in bam.fetch(until_eof=True):

            if r.is_qcfail:
                qcfail += 1
                continue
            if r.is_duplicate:
                duplicate += 1
                continue
            if r.is_secondary:
                secondary += 1
                continue
            if r.is_unmapped:
                unmapped += 1
                continue
            if not _is_read_uniq_mapping(r):
                multi += 1
                continue

            map_strand = '-' if r.is_reverse else '+'
            ref_positions = r.get_reference_positions()
            strand = None
            pos = None
            chrom = r.reference_name
            length = r.query_length
            if protocol == 'forward':
                if map_strand == '+':
                    strand = '+'
                    pos = ref_positions[0]
                else:
                    strand = '-'
                    pos = ref_positions[-1]
            elif protocol == 'reverse':
                if map_strand == '+':
                    strand = '-'
                    pos = ref_positions[-1]
                else:
                    strand = '+'
                    pos = ref_positions[0]
            coverages[length][strand][(chrom, pos)] += 1

            valid += 1

            pbar.update()

    summary = ('summary:\n\ttotal_reads: {}\n\tunique_mapped: {}\n'
              '\tqcfail: {}\n\tduplicate: {}\n\tsecondary: {}\n'
              '\tunmapped:{}\n\tmulti:{}\n\nlength dist:\n').format(total_count,
                       valid, qcfail, duplicate, secondary, unmapped, multi)

    for length in coverages:
        reads_of_length = 0
        for strand in coverages[length]:
            to_write = ''
            cur_chrom = ''
            for chrom, pos in sorted(coverages[length][strand]):
                if chrom != cur_chrom:
                    cur_chrom = chrom
                    to_write += 'variableStep chrom={}\n'.format(chrom)
                to_write += '{}\t{}\n'.format(
                    pos, coverages[length][strand][(chrom, pos)])
            fname = '{}_{}nt_{}.wig'.format(prefix, length, 'pos'
                                            if strand == '+' else 'neg')
            with open(fname, 'w') as output:
                output.write(to_write)
            reads_of_length += 1
        summary += '\t{}: {}\n'.format(length, reads_of_length)
    with open('{}_summary.txt'.format(prefix), 'w') as output:
        output.write(summary)

    return coverages


def align_coverages(coverages, base, saveto):
    """align coverages to determine the lag to the base

    Parameters
    ----------
    coverages: str
               Path to file which contains paths of all metagene
               from different lengths
               format:
               length (e.g. 28) path (e.g. metagene_28.tsv)
               length (e.g. 29) path (e.g. metagene_29.tsv)
    base: int
          The reference length to align against
    saveto: str
          Path to save the aligned offsets
    """
    base = int(base)
    with open(coverages) as f:
        cov_lens = f.readlines()
    cov_lens = {
        int(x.strip().split()[0]): x.strip().split()[1]
        for x in cov_lens
    }
    if base not in cov_lens:
        print('Failed to find base {} in coverages.'.format(base))
        return
    reference = pd.read_table(cov_lens[base])['count']
    to_write = 'relative lag to base: {}\n'.format(base)
    for length, path in cov_lens.items():
        cov = pd.read_table(path)['count']
        xcorr = np.correlate(reference, cov, 'full')
        origin = len(xcorr) // 2
        bound = min(base, length)
        xcorr = xcorr[origin - bound:origin + bound]
        lag = np.argmax(xcorr) - len(xcorr) // 2
        to_write += '\tlag of {}: {}\n'.format(length, lag)
    with open(saveto, 'w') as output:
        output.write(to_write)


def merge_wigs(wigs, offsets, strand, saveto):
    """merge wigs from different lengths into one with shift of offsets

    Parameters
    ----------
    wigs: str
          Path to file which contains paths of all wigs from differnt lengths
          format:
          length1 path1
          length2 path2
    offsets: str
             Path to file which contains offset for each length
             format:
             length1 offset1
             length2 offset2
    strand: str
            '+' for positive strand,
            '-' for negative strand
    saveto: str
            Path to save merged wig
    """
    coverages = defaultdict(int)
    with open(wigs) as wf:
        wigs = {
            int(x.strip().split()[0]): x.strip().split()[1]
            for x in wf.readlines()
        }
    with open(offsets) as of:
        offsets = {
            int(x.strip().split()[0]): int(x.strip().split()[1])
            for x in of.readlines()
        }
    for length, wig in wigs.items():
        with open(wig) as f:
            for line in f:
                if line.startswith('variableStep'):
                    line = line.strip()
                    chrom = line[line.index('=') + 1:]
                else:
                    pos, count = line.strip().split()
                    pos, count = int(pos), int(count)
                    if strand == '+':
                        pos_shifted = pos + offsets[length]
                    else:
                        pos_shifted = pos - offsets[length]
                    if pos_shifted >= 0:
                        coverages[(chrom, pos_shifted)] += count
    to_write = ''
    cur_chrom = ''
    for chrom, pos in sorted(coverages):
        if chrom != cur_chrom:
            cur_chrom = chrom
            to_write += 'variableStep chrom={}\n'.format(chrom)
        to_write += '{}\t{}\n'.format(pos, coverages[(chrom, pos)])
    with open(saveto, 'w') as output:
        output.write(to_write)


def parse_annotation(annotation):
    cds = []
    utr5 = []
    utr3 = []

    return (cds, utr5, utr3)


def detect_orfs(gtf, fasta, bam, prefix, annotation=None, protocol=None):

    if not isinstance(gtf, GTFReader):
        gtf = GTFReader(gtf)

    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    if annotation is None:
        cds, utr5, utr3 = prepare_orfs(gtf, fasta, prefix)
    else:
        cds, utr5, utr3 = parse_annotation(annotation)

    if protocol is None:
        protocol, _, _ = infer_protocol(bam, gtf, prefix)

    coverages = split_bam(bam, protocol, prefix)
