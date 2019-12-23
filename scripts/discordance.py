import itertools
from os.path import join

import matplotlib.pyplot as plt

from slugify import slugify

from scripts.coverage_test import calculate_coverage
from scripts.kmer_analysis import get_kmers_read_pos
from scripts.utils import get_kmers, get_kmers_positions, get_fasta_len
from config import *


def draw_discordance_coverage(assembly_len, assembly_name1, assembly_name2, bed_file1, bed_file2,
                              read_names1, read_names2, plot_title, out_dir, plot_filename):
    plt.figure(figsize=(20,10))
    plt.title(plot_title, fontsize=40)
    plt.xlabel('Position', fontsize=32)
    plt.ylabel('Coverage', fontsize=32)

    plt.xticks(fontsize=32)
    plt.yticks(fontsize=32)

    coverage1 = calculate_coverage(assembly_len, bed_file1, read_names1)
    p1 = plt.plot(range(len(coverage1)), coverage1)
    coverage2 = calculate_coverage(assembly_len, bed_file2, read_names2)
    p2 = plt.plot(range(len(coverage2)), coverage2, color="orange")
    plt.legend((p1[0], p2[0]), ("Preference for %s" % assembly_name1, "Preference for %s" % assembly_name2), prop={'size': 32})
    plt.ylim(0, max([max(coverage1)+2, max(coverage2)+2, 10]))
    plot_file = join(out_dir, plot_filename)
    plt.savefig(plot_file)
    print("Discordance coverage plot saved to %s" % plot_file)


def do(assemblies, reads_file, out_dir):
    if len(assemblies) < 2:
        return
    print("")
    print("*********************************")
    print("Discordance test started...")
    for (assembly1, assembly2) in itertools.combinations(assemblies, 2):
        kmers1 = get_kmers(assembly1.good_kmers_fname)
        kmers2 = get_kmers(assembly2.good_kmers_fname)
        shared_kmers = kmers1.intersection(kmers2)
        with open(join(out_dir, "shared_kmers.txt"), "w") as f:
            for k in shared_kmers:
                f.write(k + "\n")

        ref_kmers_pos1, kmer_by_pos1 = get_kmers_positions(assembly1.fname, shared_kmers)
        ref_kmers_pos2, kmer_by_pos2 = get_kmers_positions(assembly2.fname, shared_kmers)

        score1 = 0
        voting_reads1 = 0
        score2 = 0
        voting_reads2 = 0

        read_kmer_pos1, reads_coords1 = get_kmers_read_pos(assembly1, reads_file, shared_kmers)
        read_kmer_pos2, reads_coords2 = get_kmers_read_pos(assembly2, reads_file, shared_kmers)

        selected_reads1 = []
        selected_reads2 = []
        for read_name in reads_coords1.keys():
            read_name = slugify(read_name)
            read_score = 0
            for kmer in shared_kmers:
                if kmer in ref_kmers_pos1 and read_name in read_kmer_pos1 and kmer in read_kmer_pos1[read_name] and \
                        read_kmer_pos1[read_name][kmer] in reads_coords1[read_name]:
                    if abs(ref_kmers_pos1[kmer] - reads_coords1[read_name][read_kmer_pos1[read_name][kmer]]) <= 1000 or \
                            abs(ref_kmers_pos2[kmer] - reads_coords1[read_name][
                                read_kmer_pos1[read_name][kmer]]) <= 1000:
                        score1 += 1
                        read_score += 1
                if kmer in ref_kmers_pos2 and read_name in read_kmer_pos2 and kmer in read_kmer_pos2[read_name] and \
                        read_kmer_pos2[read_name][kmer] in reads_coords2[read_name]:
                    if abs(ref_kmers_pos2[kmer] - reads_coords2[read_name][read_kmer_pos2[read_name][kmer]]) <= 1000:
                        score2 += 1
                        read_score -= 1
            if read_score > KMER_SIZE:
                voting_reads1 += 1
                selected_reads1.append(read_name)
            elif read_score < -KMER_SIZE:
                voting_reads2 += 1
                selected_reads2.append(read_name)
        total_discordance = score1 - score2
        print("Discordance between %s and %s: %d. There are %d (%d) discordant reads voting for %s (%s)." %
              (assembly1.name, assembly2.name, total_discordance, voting_reads1, voting_reads2, assembly1.name,
               assembly2.name))
        plot_fname = join(out_dir, "report", "discordance_%s_vs_%s.png" % (assembly1.name, assembly2.name))
        draw_discordance_coverage(max(get_fasta_len(assembly1.fname), get_fasta_len(assembly2.fname)), assembly1.name,
                                  assembly2.name,
                                  assembly1.bed_fname, assembly2.bed_fname, selected_reads1, selected_reads2,
                                  "Discordance coverage", out_dir, plot_fname)
    print("Discordance test finished.")

