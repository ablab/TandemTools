import os
import re
import string
import sys
from joblib import Parallel, delayed
from bisect import bisect_left, insort
from collections import deque
from itertools import islice
from os.path import abspath, join, dirname, realpath

import numpy as np
from Bio import SeqIO

import config
from config import *

old_chars = "ACGT"
replace_chars = "TGCA"
try:
    maketrans = ''.maketrans
except AttributeError:
    # for Python 2
    from string import maketrans

tab = maketrans(old_chars, replace_chars)

cigar_pattern = re.compile(r'(\d+[M=XIDNSH])')


def get_ext_tools_dir():
    return abspath(join(dirname(dirname(realpath(__file__))), "ext_tools"))


def get_monomers_dict(monomers_fname):
    mm_records = list(SeqIO.parse(monomers_fname, "fasta"))
    monomers_names = [r.id for r in mm_records]
    symbols = list(string.ascii_uppercase+string.ascii_lowercase)
    monomers_dict = dict(zip(monomers_names, symbols))
    return monomers_dict


def is_python2():
    return sys.version_info[0] < 3


def compress_homopolymers(assemblies):
    for assembly in assemblies:
        real_coords = dict()
        with open(assembly.fname) as handle:
            with open(assembly.compressed_fname, "w") as out:
                for record in SeqIO.parse(handle, 'fasta'):
                    seq = record.seq
                compress_seq = ''
                prev_s = ''
                cur_i = 0
                for i,x in enumerate(seq):
                    if x == prev_s:
                        real_coords[i] = cur_i
                        continue
                    compress_seq+= x
                    real_coords[cur_i] = i
                    prev_s = x
                    cur_i+=1
                assembly.real_coords = real_coords
                out.write(">" + str(record.id) + "\n")
                out.write(compress_seq + "\n")


def run_parallel(_fn, fn_args, n_jobs=None, filter_results=False):
    n_jobs = n_jobs or config.MAX_THREADS
    parallel_args = {'n_jobs': n_jobs}
    try:
        new_style_parallel_args = {'backend': 'multiprocessing'}
        Parallel(**new_style_parallel_args)
        parallel_args.update(new_style_parallel_args)
    except TypeError:
        pass
    results_tuples = Parallel(**parallel_args)(delayed(_fn)(*args) for args in fn_args)
    results = []
    if results_tuples:
        if isinstance(results_tuples[0], list) or isinstance(results_tuples[0], tuple):
            results_cnt = len(results_tuples[0])
            results = [[result_list[i] for result_list in results_tuples] for i in range(results_cnt)]
        else:
            results = [result for result in results_tuples if result or not filter_results]
    return results


def add_suffix(fname, suffix):
    base, ext = os.path.splitext(fname)
    if ext in ['.gz', '.bz2', '.zip']:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base + (('.' + suffix) if suffix else '') + ext


def get_flye_cfg_fname():
    return abspath(join(dirname(dirname(realpath(__file__))), 'misc', 'asm_raw_reads.cfg'))


def rev_comp(seq):
    return str(seq).translate(tab)[::-1]


def get_fasta_len(fasta_fname):
    with open(fasta_fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            return len(record.seq)


def check_fasta_files(fasta_fnames):
    for fasta_fname in fasta_fnames:
        num_seqs = 0
        with open(fasta_fname) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                num_seqs += 1
                if num_seqs > 1:
                    print("ERROR! TandemQUAST is currently designed for the assemblies consisting of one contig. "
                          "If you still want to assess the quality of this assembly, you can concatenate contigs into one file (separating them with N-s) and re-run tandemQUAST.")
                    return False
    return True


def get_non_canonical_kmers(fasta_fname, kmers):
    non_canonical_kmers = set()
    with open(fasta_fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            assembly_seq = str(record.seq)
            assembly_len = len(record.seq)
            for i in range(assembly_len - KMER_SIZE+1):
                kmer = assembly_seq[i:i + KMER_SIZE]
                if kmer in kmers or rev_comp(kmer) in kmers:
                    non_canonical_kmers.add(kmer)
    return non_canonical_kmers


def get_kmers_positions(fasta_fname, unique_kmers):
    ref_kmers_pos = dict()
    kmer_by_pos = []
    with open(fasta_fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            assembly_seq = str(record.seq)
            assembly_len = len(record.seq)
            kmer_by_pos = [0] * len(record.seq)
            for i in range(assembly_len - KMER_SIZE+1):
                kmer = assembly_seq[i:i + KMER_SIZE]
                if kmer in unique_kmers:
                    ref_kmers_pos[kmer] = i
                    kmer_by_pos[i] = kmer
                elif rev_comp(kmer) in unique_kmers:
                    ref_kmers_pos[rev_comp(kmer)] = i
                    kmer_by_pos[i] = rev_comp(kmer)
    return ref_kmers_pos, kmer_by_pos


def get_kmers(kmers_fname):
    unique_kmers = set()
    with open(kmers_fname) as f:
        for line in f:
            kmer = line.strip().split()[0]
            unique_kmers.add(kmer)
    return unique_kmers


def get_canon_kmer(kmer):
    return min(kmer, rev_comp(kmer))


def mean_ignore_zeros(arr):
    arr = [x for x in arr if x!=0]
    if not arr:
        return 0
    return np.mean(arr)


def running_median(seq, N):
    seq = iter(seq)
    arr = [item for item in islice(seq, N)]
    d = deque(arr)
    arr.sort()
    medians = [mean_ignore_zeros(arr)]
    for item in seq:
        old = d.popleft()
        d.append(item)
        del arr[bisect_left(arr, old)]
        insort(arr, item)
        medians.append(mean_ignore_zeros(arr))
    return medians
