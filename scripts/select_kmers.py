import random
import subprocess
import sys
from collections import defaultdict
from os.path import join, exists, basename

from Bio import SeqIO
import numpy as np

from config import *
from scripts.reporting import make_plot
from scripts.utils import rev_comp, get_ext_tools_dir, get_kmers_positions

kmc_bin = join(get_ext_tools_dir(), 'KMC', sys.platform, 'kmc')
kmctools_bin = join(get_ext_tools_dir(), 'KMC', sys.platform, 'kmc_tools')


def draw_plot(assembly, kmers, plot_fname, kmers_type):
    unique_kmer_pos = []
    assembly_len = 0
    with open(assembly.fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            assembly_seq = str(record.seq)
            assembly_len = len(assembly_seq)
            unique_kmer_pos = [0] * assembly_len
            for i in range(len(assembly_seq) - KMER_SIZE + 1):
                kmer = assembly_seq[i:i + KMER_SIZE]
                if kmer in kmers or rev_comp(kmer) in kmers:
                    unique_kmer_pos[i] = 1
    unique_bins = [sum(unique_kmer_pos[i:i + KMER_WINDOW_SIZE]) for i in range(0, assembly_len, KMER_WINDOW_SIZE)]

    make_plot(plot_fname, "Distribution of %s k-mers" % kmers_type, assembly.label, xlabel="Position", ylabel="unique $\it{k}$--mer counts",
              bar_values=unique_bins, plot_color="blue", max_x=assembly_len)


def run_kmc(seq, out_dir):
    kmers = set()
    tmp_db_fname = join(out_dir, "tmp_kmc")
    tmp_out_fname = join(out_dir, "tmp_kmers.txt")
    cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-b", "-fm", "-ci1", "-cx20", "-cs10000", seq, tmp_db_fname, out_dir]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    cmd = [kmctools_bin, "transform", tmp_db_fname, "dump", tmp_out_fname]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    with open(tmp_out_fname) as f:
        for line in f:
            kmer, freq = line.strip().split()
            kmers.add(kmer)
    return kmers


def find_rare_kmers(assembly, assembly_kmc_db, tmp_dir):
    potential_kmers = run_kmc(assembly.fname, tmp_dir)
    num_kmers = 0
    with open(assembly.all_kmers_fname, "w") as f:
        for i, kmer in enumerate(potential_kmers):
            f.write(">%d\n" % i)
            f.write("%s\n" % kmer)
            num_kmers += 1
    cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-cx1", "-ci1", assembly.all_kmers_fname, assembly_kmc_db, tmp_dir]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    print("  %d unique k-mers without filtration" % num_kmers)


def filter_kmers(reads_fname, kmers):
    bad_kmers = set()
    with open(reads_fname) as handle:
        for record in SeqIO.parse(handle, reads_fname.split('.')[-1]):
            read_seq = str(record.seq)
            read_len = len(read_seq)
            read_counter = defaultdict(int)
            for i in range(read_len - KMER_SIZE + 1):
                kmer = str(read_seq[i:i + KMER_SIZE])
                if kmer in kmers or rev_comp(kmer) in kmers:
                    read_counter[kmer] +=1
            for kmer, c in read_counter.items():
                if c > 1:
                    bad_kmers.add(min(kmer, rev_comp(kmer)))

    solid_kmers = kmers-bad_kmers
    return solid_kmers


def do(assemblies, reads_fname, hifi_reads_fname, out_dir, tmp_dir, no_reuse, only_polish=False):
    print("")
    print("*********************************")
    print("K-mers selection started...")

    for assembly in assemblies:
        print("")
        if exists(assembly.solid_kmers_fname) and not no_reuse:
            print("Reusing previously selected k-mers...")
            continue
        print("Analyzing %s assembly" % assembly.label)
        readkmers_kmc_db = join(tmp_dir, basename(reads_fname))
        readkmers_fname = join(out_dir, basename(reads_fname)+".kmers.txt")
        if not exists(readkmers_fname) or no_reuse:
            cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-ci1", "-cs10000", reads_fname, readkmers_kmc_db, tmp_dir]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
            cmd = [kmctools_bin, "transform", readkmers_kmc_db, "dump", readkmers_fname]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        tmp_all_kmc_db = join(out_dir, assembly.name+".all_kmers")
        tmp_all_kmers_fname = join(out_dir, assembly.name+".all_kmers.txt")
        cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-b", "-fm", "-ci1", assembly.fname, tmp_all_kmc_db, tmp_dir]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "transform", tmp_all_kmc_db, "dump", tmp_all_kmers_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        kmer_freq = defaultdict()
        with open(tmp_all_kmers_fname) as f:
            for line in f:
                kmer, freq = line.split()
                kmer_freq[min(kmer,rev_comp(kmer))] = int(freq)

        assembly_kmc_db = join(tmp_dir, assembly.name)
        find_rare_kmers(assembly, assembly_kmc_db, tmp_dir)

        tmp_kmers_fname = join(tmp_dir, assembly.name+".kmers.txt")
        missed_kmers_db = join(tmp_dir, assembly.name+".missed")
        missed_kmers_fname = join(out_dir, assembly.name+".missed.kmers.txt")
        intersect_kmc_db = join(tmp_dir, assembly.name+".intersect")
        kmers_hist_fname = join(out_dir, assembly.name+".hist.kmers.txt")
        cmd = [kmctools_bin, "simple", assembly_kmc_db, readkmers_kmc_db, "intersect", intersect_kmc_db, "-cs%d" % MAX_OCC_IN_READS, "-ocright"]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "simple", assembly_kmc_db, readkmers_kmc_db, "kmers_subtract", missed_kmers_db]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "transform", missed_kmers_db, "dump", missed_kmers_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "transform", intersect_kmc_db, "-ci1", "-cx%d" % MAX_OCCURRENCES, "dump", tmp_kmers_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "transform", intersect_kmc_db, "histogram", kmers_hist_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        hist_values = [0] * (MAX_OCCURRENCES+1)
        missed_kmers = set()
        with open(missed_kmers_fname) as f:
            for line in f:
                hist_values[0] += 1
                missed_kmers.add(line.split()[0])
        print("  %d k-mers do not occur in reads!" % hist_values[0])
        plot_fname = join(out_dir, assembly.name + "_missedkmers_hist.png")
        draw_plot(assembly, missed_kmers, plot_fname, "missed")

        with open(kmers_hist_fname) as f:
            for line in f:
                freq, num_kmers = map(int, line.strip().split())
                hist_values[min(freq, MAX_OCCURRENCES)] += num_kmers

        plot_fname = join(out_dir, assembly.name + "_kmers_hist.png")
        make_plot(plot_fname, "Histogram of k-mer occurrences", assembly.label, xlabel="Occurrences in reads", ylabel="$\it{k}$-mer counts",
                  bar_values=hist_values, plot_color="green")

        total_kmers = sum(hist_values[1:-1])
        cum_occ = 0

        real_min_occ = 0
        real_max_occ = MAX_OCCURRENCES

        for i, freq in enumerate(hist_values[1:]):
            cum_occ += freq
            if not real_min_occ and cum_occ >= total_kmers*MIN_FREQ:
                real_min_occ = i + 1
                continue
            if cum_occ >= total_kmers*MAX_FREQ:
                real_max_occ = i + 1
                break

        # print("Max k-mer occurrences in read-set: %d" % real_max_occ)

        hifi_kmers = defaultdict(int)
        if hifi_reads_fname:
            hifi_kmers_kmc_db = join(tmp_dir, basename(hifi_reads_fname))
            hifi_kmers_fname = join(out_dir, basename(hifi_reads_fname)+".kmers.txt")
            if not exists(hifi_kmers_fname) or no_reuse:
                cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-ci1", "-cs10000", hifi_reads_fname, hifi_kmers_kmc_db, tmp_dir]
                subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
                cmd = [kmctools_bin, "transform", hifi_kmers_kmc_db, "dump", hifi_kmers_fname]
                subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
            with open(hifi_kmers_fname) as f:
                for line in f:
                    k, freq = line.split()
                    hifi_kmers[k] = int(freq)
            # print("Total %d k-mers in HiFi read-set" % len(hifi_kmers))

        selected_kmers = set()
        solid_kmers = set()

        with open(tmp_kmers_fname) as f:
            for line in f:
                kmer, occ = line.strip().split()
                if real_min_occ <= int(occ) <= real_max_occ and (not hifi_kmers or hifi_kmers[kmer] >= 1):
                    selected_kmers.add(kmer)

        ref_kmers_pos, kmer_by_pos = get_kmers_positions(assembly.fname, selected_kmers)
        kmer_markers = [1 if i else 0 for i in kmer_by_pos]
        kmer_density = [sum(kmer_markers[i:i+KMER_SELECT_WINDOW_SIZE]) for i in range(0, len(kmer_by_pos), KMER_SELECT_WINDOW_SIZE/2)]
        kmer_density = [k for k in kmer_density if k > 1]
        max_density = int(np.mean(kmer_density) + 2*np.std(kmer_density))
        filtered_kmers = []
        for i in range(0, len(kmer_by_pos), KMER_SELECT_WINDOW_SIZE/2):
            cur_kmers = [k for k in kmer_by_pos[i:i+KMER_SELECT_WINDOW_SIZE] if k]
            if not cur_kmers:
                continue
            if sum(kmer_markers[i:i+KMER_SELECT_WINDOW_SIZE]) > max_density:
                filtered_kmers.extend(random.sample(cur_kmers, max_density))
            else:
                filtered_kmers.extend(cur_kmers)
        selected_kmers = set(filtered_kmers)
        print("  %d locally unique k-mers were selected" % len(selected_kmers))

        plot_fname = join(out_dir, assembly.name + "_selected_kmers.png")
        draw_plot(assembly, selected_kmers, plot_fname, "selected")

        #print("  Filtering unique solid k-mers...")
        with open(assembly.kmers_fname, "w") as out_f:
            for kmer in selected_kmers:
                out_f.write("%s\n" % kmer)
                if kmer_freq[min(kmer,rev_comp(kmer))] == 1:
                    solid_kmers.add(kmer)

        if not only_polish:
            solid_kmers = filter_kmers(reads_fname,solid_kmers)
        with open(assembly.solid_kmers_fname, "w") as f:
            for kmer in solid_kmers:
                f.write("%s\n" % kmer)
        print("  %d unique solid k-mers were selected" % len(solid_kmers))

