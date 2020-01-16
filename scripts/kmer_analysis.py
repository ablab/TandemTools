from collections import defaultdict
from os.path import exists, join

from Bio import SeqIO
from slugify import slugify

from config import *
from scripts.reporting import make_plot, draw_report_table
from scripts.utils import rev_comp, get_kmers, get_fasta_len, get_kmers_positions, cigar_pattern, \
    get_non_canonical_kmers


def get_clusters(arr):
    arr.sort()
    clusters = []
    cluster = [arr[0]]
    for i in range(1, len(arr)):
        if arr[i] - arr[i - 1] < MAX_CLUMP_DIST:
            cluster.append(arr[i])
        else:
            if len(cluster) >= MIN_CLUMP_SIZE:
                clusters.append(cluster)
            cluster = [arr[i]]
    if len(cluster) >= MIN_CLUMP_SIZE:
        clusters.append(cluster)
    return clusters


def get_kmers_read_pos(assembly, reads_fname, solid_kmers=None):
    reads_coords = defaultdict(dict)
    read_lengths = dict()
    with open(assembly.bed_fname) as f:
        for i, line in enumerate(f):
            fs = line.split()
            if len(fs) < 7:
                continue
            ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs[:7]
            read_lengths[read_name] = int(read_len)

    strand_in_ref = defaultdict()
    with open(assembly.sam_fname) as f:  ## TODO: parse in parallel
        for i, line in enumerate(f):
            fs = line.split()
            if len(fs) < 7:
                continue
            #ref     16      63e6e431        167     60      37M4I78M1I17M1I7M2I7M1I20M1D18M2I16M1D39M1I12M1I74M1I26M1I30M2I81M2D16M1I32M1I22M4I15M2D
            read_name, ref_s, cigar = slugify(fs[0]), int(fs[3]), fs[5]
            if read_name not in read_lengths:
                continue
            strand_direction = -1 if fs[1] == '16' else 1
            strand_in_ref[read_name] = strand_direction
            align_start = 0 #if strand_direction == 1 else read_lengths[read_name]
            operations = cigar_pattern.findall(cigar)
            align_len = 0
            ref_len = 0
            for op in operations:
                n_bases, operation = int(op[:-1]), op[-1]
                if operation == 'S' or operation == 'H':
                    align_start += n_bases*strand_direction
                elif operation == 'M':
                    for i in range(n_bases):
                        reads_coords[read_name][align_start+(align_len*strand_direction)] = ref_s + ref_len
                        ref_len += 1
                        align_len += 1
                elif operation == 'D':
                    ref_len += n_bases
                elif operation == 'I':
                    for i in range(n_bases):
                        reads_coords[read_name][align_start+(align_len*strand_direction)] = ref_s + ref_len
                        align_len += 1

    read_kmer_pos = defaultdict()
    if solid_kmers:
        non_canonical_kmers = get_non_canonical_kmers(assembly.fname, solid_kmers)
        with open(reads_fname) as handle:
            for record in SeqIO.parse(handle, reads_fname.split('.')[-1]):
                read_name = slugify(record.name)
                if read_name not in strand_in_ref:
                    continue
                read_seq = str(record.seq)
                read_len = len(read_seq)
                read_kmer_pos[read_name] = dict()
                for i in range(read_len - KMER_SIZE+1):
                    kmer = str(read_seq[i:i + KMER_SIZE])
                    if kmer in solid_kmers:
                        read_kmer_pos[read_name][kmer] = i
                    elif rev_comp(kmer) in solid_kmers:
                        read_kmer_pos[read_name][rev_comp(kmer)] = i
    return read_kmer_pos, reads_coords


def do(assemblies, reads_fname, out_dir, no_reuse=False):
    print("")
    print("*********************************")
    print("K-mer analysis started...")

    kmer_stats_table = [['Assembly'] + [assembly.name for assembly in assemblies]]
    kmer_stats_table.append(["K-mers forming single clump"] + ["-" for assembly in assemblies])
    kmer_stats_table.append(["K-mers forming multiple clumps"] + ["-" for assembly in assemblies])
    kmer_stats_table.append(["K-mers forming no clumps"] + ["-" for assembly in assemblies])
    for i,assembly in enumerate(assemblies):
        if exists(assembly.good_kmers_fname) and exists(join(out_dir, "report", assembly.name + "_kmer_stats.txt")) and not no_reuse:
            print("Reusing latest results...")
            with open(join(out_dir, "report", assembly.name + "_kmer_stats.txt")) as f:
                line = f.readline()
                kmer_stats_table[1][i+1] = line.split("\t")[1]
                line = f.readline()
                kmer_stats_table[2][i+1] = line.split("\t")[1]
                line = f.readline()
                kmer_stats_table[3][i+1] = line.split("\t")[1]
            continue

        solid_kmers = get_kmers(assembly.solid_kmers_fname)
        assembly_len = get_fasta_len(assembly.fname)
        ref_kmers_pos, kmer_by_pos = get_kmers_positions(assembly.fname, solid_kmers)
        read_kmer_pos, reads_coords = get_kmers_read_pos(assembly, reads_fname, solid_kmers)

        no_clumps = []
        one_clump = []
        good_kmers = []
        multi_clumps = []
        bad_kmers1 = []
        bad_kmers2 = []
        multi_clump_pos = []
        no_clump_pos = []
        for kmer, pos in ref_kmers_pos.items():
            pos_in_ref = []
            read_pos = []
            reads = []
            for read_name, kmers_pos in read_kmer_pos.items():
                if kmer in kmers_pos and kmers_pos[kmer] in reads_coords[read_name]:
                    k_pos = kmers_pos[kmer]
                    read_pos.append(k_pos)
                    pos_in_ref.append(reads_coords[read_name][k_pos])
                    reads.append(read_name)
            if read_pos and len(read_pos)>= MIN_CLUMP_SIZE:
                clusters = get_clusters(pos_in_ref)
                if not clusters:
                    no_clumps.append(pos)
                    bad_kmers1.append(kmer)
                    multi_clump_pos.extend(pos_in_ref)
                elif len(clusters) == 1:
                    one_clump.append(pos)
                    good_kmers.append(kmer)
                else:
                    multi_clumps.append(pos)
                    bad_kmers2.append(kmer)
                    no_clump_pos.extend(pos_in_ref)

        all_kmers = len(one_clump) + len(multi_clumps) + len(no_clumps)
        with open(assembly.good_kmers_fname, "w") as f:
            for kmer in good_kmers:
                f.write("%s\n" % kmer)
        with open(join(out_dir, "report", assembly.name + "_kmer_stats.txt"), "w") as f:
            f.write("Single clump\t%.2f (%d)\n" % (len(one_clump)*100.0/all_kmers,len(one_clump)))
            f.write("Multi clump\t%.2f (%d)\n" % (len(multi_clumps)*100.0/all_kmers,len(multi_clumps)))
            f.write("No clumps\t%.2f (%d)\n" % (len(no_clumps)*100.0/all_kmers,len(no_clumps)))
        kmer_stats_table[1][i+1] = "%.2f (%d)" % (len(one_clump)*100.0/all_kmers,len(one_clump))
        kmer_stats_table[2][i+1] = "%.2f (%d)" % (len(multi_clumps)*100.0/all_kmers,len(multi_clumps))
        kmer_stats_table[3][i+1] = "%.2f (%d)" % (len(no_clumps)*100.0/all_kmers,len(no_clumps))

        one_clump_dist = [0] * assembly_len
        no_clump_dist = [0] * assembly_len
        multi_clump_dist = [0] * assembly_len
        for p in one_clump:
            one_clump_dist[p] = 1
        for p in no_clumps:
            no_clump_dist[p] = 1
        for p in multi_clumps:
            multi_clump_dist[p] = 1

        one_clump_vals = [sum(one_clump_dist[i:i+KMER_WINDOW_SIZE]) for i in range(0, assembly_len, KMER_WINDOW_SIZE)]
        multi_clump_vals = [sum(multi_clump_dist[i:i+KMER_WINDOW_SIZE]) for i in range(0, assembly_len, KMER_WINDOW_SIZE)]
        no_clump_vals = [sum(no_clump_dist[i:i+KMER_WINDOW_SIZE]) for i in range(0, assembly_len, KMER_WINDOW_SIZE)]

        plot_fname = join(out_dir, "report", assembly.name + "_kmer_analysis.png")
        make_plot(plot_fname, "K-mer analysis", assembly.label, xlabel="Position", ylabel="$\it{k}$-mer counts", list_vals=[one_clump_vals, multi_clump_vals, no_clump_vals],
                  legend=("Single clump", "Multiple clumps", "No clumps"), max_x=assembly_len)
    #draw_report_table("K-mer statistics", "", kmer_stats_table)
    print("K-mer analysis finished.")
