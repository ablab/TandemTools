import os
import subprocess
import sys
from collections import defaultdict
from os.path import abspath, join, exists, dirname

from Bio import SeqIO
from slugify import slugify

from config import *
from ext_tools.Flye.flye.polishing.polish import ASSEMBLY_BIN
from scripts.utils import get_fasta_len, get_flye_cfg_fname, rev_comp, get_kmers, run_parallel, get_kmers_positions


def make_flye():
    if exists(ASSEMBLY_BIN):
        return
    flye_bin_dir = dirname(ASSEMBLY_BIN)
    flye_dir = dirname(flye_bin_dir)
    if not exists(flye_bin_dir):
        os.makedirs(flye_bin_dir)

    print('Compiling Flye...')
    return_code = subprocess.call(['make', '-C', flye_dir], stdout=open(join(flye_dir, "make.log"),"w"), stderr=open(join(flye_dir, "make.err"),"w"))
    if return_code != 0:
        raise
    print('Flye is compiled successful!')


def run_flye(assembly, reads_fname, out_dir, threads):
    try:
        make_flye()
    except:
        print('Failed to compile Flye! Please try to compile it manually: create %s folder and run "make" in %s'
              % (dirname(ASSEMBLY_BIN), dirname(dirname(ASSEMBLY_BIN))))
        sys.exit(2)
    cmd = [ASSEMBLY_BIN, '--reads', reads_fname,
           '--asm', assembly.compressed_fname or assembly.fname,
           '--kmers', abspath(assembly.kmers_fname),
           '--out-file', abspath(assembly.chains_fname),
           '--out-asm', 'draft_assembly.fasta', '--max-diff', str(assembly.max_aln_diff),
           '--genome-size', str(get_fasta_len(assembly.fname)), '--config', abspath(get_flye_cfg_fname()),
           '--log', join(out_dir, 'mapping.log'),
           '--threads', str(threads), '--min-ovlp', str(MIN_CHAIN_LEN), '--kmer', str(KMER_SIZE)]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))


def postprocess_chains(assembly):
    rare_kmers = get_kmers(assembly.kmers_fname)
    unique_kmers = get_kmers(assembly.solid_kmers_fname)

    rare_kmers_by_pos = []
    unique_kmers_by_pos = []
    assembly_seq = ""
    with open(assembly.fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            assembly_len = len(record.seq)
            assembly_seq = str(record.seq)
            rare_kmers_by_pos = [0] * assembly_len
            unique_kmers_by_pos = [0] * assembly_len
            for i in range(len(assembly_seq) - KMER_SIZE + 1):
                kmer = assembly_seq[i:i + KMER_SIZE]
                if kmer in rare_kmers or rev_comp(kmer) in rare_kmers:
                    rare_kmers_by_pos[i] = 1
                if kmer in unique_kmers or rev_comp(kmer) in unique_kmers:
                    unique_kmers_by_pos[i] = 1

    read_alignments = defaultdict(list)
    read_lengths = dict()
    read_seeds = defaultdict(lambda: defaultdict(list))
    with open(assembly.chains_fname) as f:
        #Aln     -a3aa51bb-4cZ 4 23715 24082 +cenX 1894557 1918322 3053970 -1894573 1135315 140 0.27543
        for line in f:
            fs = line.split()
            if "Aln" in line and len(fs) >= 8:
                read_name, align_start, align_end, read_len, \
                ref_name, ref_start, ref_end, ref_len = fs[1], fs[2], fs[3], fs[4], fs[5], fs[6], fs[7], fs[8]
                align_start, align_end, read_len, ref_start, ref_end, ref_len = map(int, (align_start, align_end, read_len, ref_start, ref_end, ref_len))
                read_lengths[read_name[1:]] = read_len
                if read_name.startswith('-'):
                    align_start, align_end = read_len - align_end - 1, read_len - align_start - 1
                if ref_name.startswith('-'):
                    ref_start, ref_end = ref_len - ref_end - 1, ref_len - ref_start - 1
                read_alignments[read_name[1:]].append((ref_start, ref_end, align_start, align_end))
            elif len(fs) >= 2:
                read_pos, ref_pos = int(fs[0]), int(fs[1])
                if read_name.startswith('-'):
                    read_pos = read_len - read_pos - KMER_SIZE
                if ref_name.startswith('-'):
                    ref_pos = ref_len - ref_pos - KMER_SIZE
                read_seeds[read_name[1:]][(ref_start, ref_end, align_start, align_end)].append((read_pos, ref_pos))

    num_alignments = 0
    with open(assembly.bed_fname, "w") as f:
        for read_name, aligns in read_alignments.items():
            max_kmers = 0
            max_len = 0
            selected_chain = []
            for align in aligns:
                seeds = read_seeds[read_name][align]
                seeds.sort(key=lambda x: x[1])
                best_chain = None
                best_kmers = 0
                best_len = 0
                if len(seeds) >= MIN_CHAIN_KMERS:
                    prev_pos = 0
                    unique_seeds = []
                    for seed in seeds:
                        read_pos, ref_pos = seed
                        if ref_pos - prev_pos >= KMER_SIZE or not unique_seeds:
                            unique_seeds.append((read_pos, ref_pos))
                            prev_pos = ref_pos
                    unique_seeds.append(seeds[-1])
                    seeds = unique_seeds
                    new_chains = []
                    breakpoints = []
                    for i in range(1, len(seeds)):
                        if seeds[i][1] - seeds[i-1][1] >= MAX_REF_GAP:
                            gap_s, gap_e = seeds[i-1][1]+KMER_SIZE, seeds[i][1]
                            if sum(rare_kmers_by_pos[gap_s:gap_e]) > MAX_MISSED_KMERS:
                                breakpoints.append(i-1)
                                #print(read_name,gap_s,gap_e)
                    if breakpoints:
                        chain_start1, chain_end1, chain_start2, chain_end2 = seeds[0][1], seeds[-1][1], seeds[0][0], seeds[-1][0]
                        start_n = 0
                        for p in breakpoints:
                            chain_end1, chain_end2 = seeds[p][1], seeds[p][0]
                            new_chains.append((chain_start1, chain_end1, chain_start2, chain_end2, p-start_n+1))
                            if p < len(seeds):
                                chain_start1, chain_start2, start_n = seeds[p+1][1], seeds[p+1][0], p+1

                        chain_end1, chain_end2 = seeds[-1][1], seeds[-1][0]
                        if chain_end1 > chain_start1:
                            new_chains.append((chain_start1, chain_end1, chain_start2, chain_end2, len(seeds) - p))
                        chains = []
                        total_kmers = 0
                        total_len = 0
                        for c in new_chains:
                            chain_start1, chain_end1, chain_start2, chain_end2, chain_kmers = c
                            if chain_kmers > MIN_CHAIN_KMERS and chain_end1 - chain_start1 >= MIN_CHAIN_LEN:
                                chains.append((chain_start1, chain_end1, chain_start2, chain_end2))
                                total_kmers += chain_kmers
                                total_len+= chain_end1-chain_start1
                        if total_kmers > best_kmers:
                            best_kmers = total_kmers
                            best_len = total_len
                            best_chain = chains
                    else:
                        best_kmers = len(seeds) if len(seeds) > MIN_CHAIN_KMERS/2 else 0
                        best_len = seeds[-1][1] - seeds[0][1]
                        best_chain = [[seeds[0][1], seeds[-1][1], seeds[0][0], seeds[-1][0]]]

                    if best_len < 100:
                        continue
                    if best_kmers > max_kmers:
                        max_kmers = best_kmers
                        max_len = best_len
                        selected_chain = best_chain

            for c in selected_chain:
                ref_start, ref_end, align_start, align_end = c
                if assembly.real_coords:
                    ref_start, ref_end = assembly.real_coords[ref_start], assembly.real_coords[ref_end]
                if (ref_end-ref_start) < MIN_CHAIN_LEN:
                    continue
                num_alignments += 1
                f.write("seq\t%d\t%d\t%s\t%d\t%d\t%d\n" % (
                ref_start, ref_end, slugify(read_name), align_start, align_end, read_lengths[read_name]))

    print("  Total %d alignments" % num_alignments)
    print("  Longest chains saved to %s" % assembly.bed_fname)


def do(assemblies, reads_fname, out_dir, threads, no_reuse):
    print("")
    print("*********************************")
    print("Read mapping started...")
    run_parallel(run_flye, [(assembly, reads_fname, out_dir, max(1, threads // len(assemblies)))
                            for assembly in assemblies if not exists(assembly.bed_fname) or no_reuse], n_jobs=threads)
    for assembly in assemblies:
        postprocess_chains(assembly)
    print("Read mapping finished")

