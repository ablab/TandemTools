import os
import re
import string
import subprocess
import sys
from bisect import bisect_left, insort
from collections import deque
from itertools import islice
from os.path import abspath, join, dirname, realpath, basename, exists, isdir

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

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


def run_parallel(_fn, fn_args, n_jobs=None, filter_results=False):
    n_jobs = n_jobs or config.MAX_THREADS
    parallel_args = {'n_jobs': n_jobs}
    import joblib
    from joblib import Parallel, delayed
    try:
        # starting from joblib 0.10 the default backend has changed to 'loky' which causes QUAST crashes,
        # since it uses default values in qconfig module. So, we explicitly require 'multiprocessing' here.
        # Note that Parallel doesn't have 'require' argument in joblib 0.9 and earlier.
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


def mask_files(assemblies, out_dir, no_reuse=False):
    repeatmasker_bin = "RepeatMasker"
    if not exists(repeatmasker_bin):
        print("Warning: RepeatMasker is not found! TEs will not be masked that can affect k-mer based metrics")
        for assembly in assemblies:
            assembly.fname = assembly.raw_fname
        return
    rm_dir = join(out_dir, "rm_output")
    if not isdir(rm_dir):
        os.makedirs(rm_dir)
    for assembly in assemblies:
        masked_fname = join(out_dir, "%s.masked.fa" % assembly.name)
        if exists(masked_fname) and not no_reuse:
            assembly.fname = masked_fname
            continue
        cmd = [repeatmasker_bin, "-pa", str(MAX_THREADS), assembly.raw_fname, "-dir", rm_dir]
        subprocess.call(cmd)
        #28728    0.5  0.7  0.0  T2T_CENX_57-61.1k_590-3860k_221691-2944894  2464731 2470723  (252481) + L1HS
        te_coords = []
        if exists(join(rm_dir, basename(assembly.raw_fname) + ".out")):
            with open(join(rm_dir, basename(assembly.raw_fname) + ".out")) as f:
                for line in f:
                    fs = line.split()
                    if not fs:
                        continue
                    if 'LINE' in fs[10] and int(fs[6]) - int(fs[5]) >= 500 and float(fs[1]) < 10:  #TODO: CHANGE
                        te_coords.append((int(fs[5]), int(fs[6])))
                    if 'Alu' in fs[10] and int(fs[6]) - int(fs[5]) >= 200 and float(fs[1]) < 10:
                        te_coords.append((int(fs[5]), int(fs[6])))

        with open(assembly.raw_fname) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seq = str(record.seq).upper()
                for coord in te_coords:
                    seq = seq[:coord[0]] + 'N' * (coord[1] - coord[0]) + seq[coord[1]:]
                record.seq = Seq(seq)
                with open(masked_fname, "w") as handle:
                    SeqIO.write(record, handle, "fasta")
                break
        assembly.fname = masked_fname


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
