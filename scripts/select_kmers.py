import random
import re
import subprocess
import sys
from collections import defaultdict
from os.path import join, exists, basename

from Bio import SeqIO
import numpy as np

from config import *
from scripts.reporting import make_plot
from scripts.utils import rev_comp, get_kmers_positions, get_fasta_len


def draw_plot(assembly_fname, label, kmers, plot_fname, kmers_type):
    unique_kmer_pos = []
    assembly_len = 0
    with open(assembly_fname) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            assembly_seq = str(record.seq)
            assembly_len = len(assembly_seq)
            unique_kmer_pos = [0] * assembly_len
            for i in range(len(assembly_seq) - KMER_SIZE + 1):
                kmer = assembly_seq[i:i + KMER_SIZE]
                if kmer in kmers or rev_comp(kmer) in kmers:
                    unique_kmer_pos[i] = 1
    unique_bins = [sum(unique_kmer_pos[i:i + KMER_WINDOW_SIZE]) for i in range(0, assembly_len, KMER_WINDOW_SIZE)]

    make_plot(plot_fname, "Distribution of %s k-mers" % kmers_type, label, xlabel="Position", ylabel="$\it{k}$--mer counts",
              bar_values=unique_bins, plot_color="blue", max_x=assembly_len)


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


def calculate_thresholds(kmers, ref_kmers_pos, reads_fname):
    all_reads_kmer_pos = defaultdict(dict)
    with open(reads_fname) as handle:
        for record in SeqIO.parse(handle, reads_fname.split('.')[-1]):
            read_name = record.id
            read_seq = str(record.seq)
            read_len = len(read_seq)
            read_counter = defaultdict(int)
            read_pos = dict()
            for i in range(read_len - KMER_SIZE + 1):
                kmer = str(read_seq[i:i + KMER_SIZE])
                if kmer in kmers:
                    read_pos[kmer] = i
                    read_counter[kmer] +=1
                elif rev_comp(kmer) in kmers:
                    read_pos[rev_comp(kmer)] = i
                    read_counter[rev_comp(kmer)] +=1
            for kmer, c in read_counter.items():
                if c == 1:
                    all_reads_kmer_pos[read_name][kmer] = read_pos[kmer]
    diff_distances = []
    positioned_kmers = sorted([(kmer, ref_kmers_pos[kmer]) for kmer in kmers], key=lambda x: x[1])

    prev_kmer, prev_pos = positioned_kmers[0]
    for i in range(1, len(kmers)):
        kmer, pos = positioned_kmers[i]
        if pos - prev_pos < 19:
            continue
        for read_kmer_pos in all_reads_kmer_pos.values():
            if kmer in read_kmer_pos and prev_kmer in read_kmer_pos:
                dist1, dist2 = abs(pos-prev_pos), abs(read_kmer_pos[kmer]-read_kmer_pos[prev_kmer])
                dist_diff = abs(dist1-dist2)*1.0/min(dist1,dist2)
                if dist_diff <= 1:
                    diff_distances.append(dist_diff)
        prev_kmer, prev_pos = kmer, pos
    iqr_diff = np.percentile(diff_distances, 75) - np.percentile(diff_distances, 25)
    max_diff = np.median(diff_distances) + iqr_diff
    return max_diff


def check_jf_version():
    p = subprocess.Popen(["jellyfish", '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    version_pattern = re.compile('(?P<major_version>\d+)\.(?P<minor_version>\d+)')
    v = version_pattern.search(str(stdout))
    if not v.group('major_version') or not v.group('minor_version'):
        return False
    version, minor_version = 2, 2
    if int(v.group('major_version')) == version and int(v.group('minor_version')) >= minor_version:
        return True


def do(assemblies, raw_reads_fname, reads_fname, hifi_reads_fname, out_dir, tmp_dir, no_reuse, only_polish=False):
    print("")
    print("*********************************")
    print("K-mers selection started...")

    jellyfish_bin = "jellyfish"
    try:
        subprocess.call(["jellyfish", "--version"], stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    except:
        print('Jellyfish is not found in your PATH! Please install Jellyfish and rerun TandemTools. Jellyfish can be installed through conda install jellyfish')
        sys.exit(2)
    if not check_jf_version():
        print('Jellyfish is obsolete. Please install new version of Jellyfish and rerun TandemTools. '
              'Jellyfish can be installed via "conda install jellyfish"')
        sys.exit(2)

    for assembly in assemblies:
        print("")
        if exists(assembly.solid_kmers_fname) and not no_reuse:
            print("Reusing previously selected k-mers...")
            continue
        print("Analyzing %s assembly" % assembly.label)
        assembly_fname = assembly.compressed_fname or assembly.fname
        readkmers_db = join(tmp_dir, basename(reads_fname) + ".jf")
        readkmers_fname = join(out_dir, basename(reads_fname)+".kmers.txt")
        if not exists(readkmers_fname) or no_reuse:
            cmd = [jellyfish_bin, "count", "-m%d" % KMER_SIZE, "-s500M", "-C", "-t%d" % MAX_THREADS, "-o", readkmers_db, reads_fname]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
            cmd = [jellyfish_bin, "dump", readkmers_db, "-o", readkmers_fname]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        assembly_db = join(tmp_dir, assembly.name + ".jf")
        cmd = [jellyfish_bin, "count", "-m%d" % KMER_SIZE, "-s100M", "-o", assembly_db, assembly_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [jellyfish_bin, "dump", assembly_db, "-o", assembly.all_kmers_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        read_counts_fname = join(out_dir, assembly.name + ".read_counts.txt")
        cmd = [jellyfish_bin, "query", "-s", assembly.all_kmers_fname, readkmers_db, "-o", read_counts_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        unique_kmers = set()
        freq_in_assembly = defaultdict()
        with open(assembly.all_kmers_fname) as f:
            for record in SeqIO.parse(f, 'fasta'):
                kmer, freq = str(record.seq), int(record.id)
                freq_in_assembly[min(kmer,rev_comp(kmer))] = freq
                if freq == 1:
                    unique_kmers.add(kmer)

        freq_in_reads = defaultdict()
        max_occ_in_assembly = max(1, get_fasta_len(assembly_fname) // 100000)
        hist_values = [0] * (MAX_OCCURRENCES+1)
        rare_hist_values = [0] * (MAX_OCCURRENCES+1)
        missed_kmers = set()
        with open(read_counts_fname) as f:
            for line in f:
                kmer, freq = line.split()
                freq = int(freq)
                if freq_in_assembly[kmer] <= max_occ_in_assembly:
                    rare_hist_values[min(freq, MAX_OCCURRENCES)] += 1
                hist_values[min(freq, MAX_OCCURRENCES)] += 1
                if freq > 0:
                    freq_in_reads[kmer] = freq
                else:
                    missed_kmers.add(kmer)

        print("  %d k-mers do not occur in reads!" % hist_values[0])
        plot_fname = join(out_dir, assembly.name + "_missedkmers_hist.png")
        draw_plot(assembly_fname, assembly.label, missed_kmers, plot_fname, "missed")

        plot_fname = join(out_dir, assembly.name + "_kmers_hist.png")
        make_plot(plot_fname, "Histogram of k-mer occurrences", assembly.label, xlabel="Occurrences in reads", ylabel="$\it{k}$-mer counts",
                  bar_values=hist_values, plot_color="green")

        total_kmers = sum(hist_values[1:-1])
        cum_occ = 0

        real_min_occ = 0
        real_max_occ = MAX_OCCURRENCES

        for i, freq in enumerate(rare_hist_values[1:]):
            cum_occ += freq
            if not real_min_occ and cum_occ >= total_kmers*MIN_FREQ:
                real_min_occ = i + 1
                continue
            if cum_occ >= total_kmers*MAX_FREQ:
                real_max_occ = i + 1
                break

        #print("Max k-mer occurrences in read-set: %d" % real_max_occ)

        hifi_kmers = defaultdict(int)
        hifi_read_counts_fname = join(out_dir, assembly.name + ".hifi_read_counts.txt")
        if hifi_reads_fname:
            hifi_kmers_db = join(tmp_dir, basename(hifi_reads_fname) + ".jf")
            cmd = [jellyfish_bin, "count", "-m%d" % KMER_SIZE, "-s500M", "-C", "-t%d" % MAX_THREADS, "-o", hifi_kmers_db, hifi_reads_fname]
            subprocess.call(cmd)
            cmd = [jellyfish_bin, "query", "-s", assembly.all_kmers_fname, hifi_kmers_db, "-o", hifi_read_counts_fname]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
            with open(hifi_read_counts_fname) as f:
                for line in f:
                    k, freq = line.split()
                    hifi_kmers[k] = int(freq)
            # print("Total %d k-mers in HiFi read-set" % len(hifi_kmers))

        selected_kmers = set()
        solid_kmers = set()

        for kmer, occ in freq_in_reads.items():
            if freq_in_assembly[kmer] <= max_occ_in_assembly and real_min_occ <= int(occ) <= real_max_occ and (not hifi_kmers or hifi_kmers[kmer] >= 1):
                selected_kmers.add(kmer)

        ref_kmers_pos, kmer_by_pos = get_kmers_positions(assembly_fname, selected_kmers)
        kmer_markers = [1 if i else 0 for i in kmer_by_pos]
        kmer_density = [sum(kmer_markers[i:i+KMER_SELECT_WINDOW_SIZE]) for i in range(0, len(kmer_by_pos), KMER_SELECT_WINDOW_SIZE)]
        kmer_density = [k for k in kmer_density if k > 1]
        max_density = int(np.median(kmer_density) + 2*np.std(kmer_density))
        filtered_kmers = []
        for i in range(0, len(kmer_by_pos), KMER_SELECT_WINDOW_SIZE):
            cur_kmers = [k for k in kmer_by_pos[i:i+KMER_SELECT_WINDOW_SIZE] if k]
            if not cur_kmers:
                continue
            if sum(kmer_markers[i:i+KMER_SELECT_WINDOW_SIZE]) > max_density:
                filtered_kmers.extend(random.sample(cur_kmers, max_density))
            else:
                filtered_kmers.extend(cur_kmers)
        selected_kmers = set(filtered_kmers)
        print("  %d rare k-mers were selected" % len(selected_kmers))
        max_diff = calculate_thresholds(selected_kmers, ref_kmers_pos, reads_fname)
        assembly.max_aln_diff = max(0.1, min(max_diff, 0.25))

        plot_fname = join(out_dir, assembly.name + "_selected_kmers.png")
        draw_plot(assembly_fname, assembly.label, selected_kmers, plot_fname, "selected")

        #print("  Filtering unique solid k-mers...")
        with open(assembly.kmers_fname, "w") as out_f:
            for kmer in selected_kmers:
                out_f.write("%s\n" % kmer)
                if freq_in_assembly[kmer] == 1:
                    solid_kmers.add(kmer)

        if platform == "pacbio" and not only_polish:
            assembly_db = join(tmp_dir, assembly.name + "_uncompressed.jf")
            cmd = [jellyfish_bin, "count", "-m%d" % KMER_SIZE, "-s100M", "-o", assembly_db, assembly.fname]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
            cmd = [jellyfish_bin, "dump", assembly_db, "-o", assembly.all_kmers_fname]
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

            unique_kmers = set()
            with open(assembly.all_kmers_fname) as f:
                for record in SeqIO.parse(f, 'fasta'):
                    kmer, freq = str(record.seq), int(record.id)
                    if freq == 1:
                        unique_kmers.add(kmer)
            solid_kmers = filter_kmers(raw_reads_fname, unique_kmers)
        elif not only_polish:
            solid_kmers = filter_kmers(reads_fname, solid_kmers)
        with open(assembly.solid_kmers_fname, "w") as f:
            for kmer in solid_kmers:
                f.write("%s\n" % kmer)
        print("  %d unique solid k-mers were selected" % len(solid_kmers))

