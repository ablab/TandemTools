import subprocess
import sys
from collections import defaultdict
from os.path import join, exists, basename
from slugify import slugify

from config import *
from scripts.kmer_analysis import get_kmers_read_pos
from scripts.reporting import make_plot, make_plotly_html
from scripts.unit_analysis import analyze_unit_structure
from scripts.utils import get_monomers_dict, cigar_pattern, get_fasta_len, get_ext_tools_dir

decomposer_bin = join(get_ext_tools_dir(), "string_decomposer", "run_decomposer.py")


def approx_binary_search(arr, arr_i, l, r, x):
    while l <= r and  l + (r - l) / 2 < len(arr):
        mid = l + (r - l) // 2
        if abs(x - arr[mid][arr_i]) <= 50:
            return mid
        elif arr[mid][arr_i] < x:
            l = mid + 1
        else:
            r = mid - 1
    return -1


def run_decomposer(fasta_fname, out_fname, monomers_fname):
    print("Running StringDecomposer, it can take a while...")
    cmd = [decomposer_bin, "-o", out_fname, "-t", str(MAX_THREADS), "-i", "70", fasta_fname, monomers_fname]
    return_code = subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    if return_code != 0 or not exists(out_fname):
        print("Error! Cannot decompose %s to monomers" % fasta_fname)
        sys.exit(1)


def get_reads_monomer_structure(reads_fname, monomers_fname, out_dir):
    decomposition_fname = join(out_dir, "decomposition_%s.tsv" % basename(reads_fname).split('.')[0])
    if not exists(decomposition_fname):
        run_decomposer(reads_fname, decomposition_fname, monomers_fname)

    reads_mm_structure = defaultdict(list)
    monomers_dict = get_monomers_dict(monomers_fname)
    with open(decomposition_fname) as handle:
        for line in handle:
            read_name, monomer_name, s, e = line.split('\t')[:4]
            monomer_name = monomer_name.replace("'", "")
            s, e = min(int(s), int(e)), max(int(s), int(e))
            reads_mm_structure[slugify(read_name)].append((monomers_dict[monomer_name], s, e))
    return reads_mm_structure


def get_ref_monomers(assembly, monomers_fname, out_dir):
    decomposition_fname = join(out_dir, "decomposition_%s.tsv" % assembly.name)
    if not exists(decomposition_fname):
        run_decomposer(assembly.raw_fname, decomposition_fname, monomers_fname)

    monomers = []
    monomer_stats = defaultdict(list)
    monomers_dict = get_monomers_dict(monomers_fname)
    with open(decomposition_fname) as handle:
        for line in handle:
            # polished_repeat_4       X02418.1:1870-2055'     54726   54912   99.00   +
            read_name, monomer_name, s, e, idy = line.split('\t')[:5]
            monomer_name = monomer_name.replace("'", "")
            monomer_name = monomer_name.split()[0]
            s, e = map(int, (s,e))
            monomers.append((monomers_dict[monomer_name], s, e))
            monomer_stats[monomers_dict[monomer_name]].append((s, e-s+1, float(idy)))
    return monomers, monomer_stats


def parse_read_alignments(assembly, nucl_align_dir):
    reads_coords = defaultdict(dict)
    with open(assembly.bed_fname) as f:
        for line in f:
            fs = line.split()
            if len(fs) < 7:
                continue
            ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs[:7]
            ref_s, ref_e, align_start, align_end, read_len = map(int, (ref_s, ref_e, align_start, align_end, read_len))

            # if read_len < MIN_LONGREAD_LEN:
            #    continue
            nucl_align_fname = join(nucl_align_dir, read_name + ".out")
            if not exists(nucl_align_fname):
                continue
            with open(nucl_align_fname) as nucl_align_f:
                for line in nucl_align_f:
                    fs = line.split()
                    if len(fs) < 8:
                        continue
                    read_name, read_len, s_nucl, e_nucl, strand, ref_s_nucl, ref_e_nucl, cigar = \
                        fs[0], fs[1], fs[2], fs[3], fs[4], fs[7], fs[8], fs[-1]

                    s_nucl, e_nucl, ref_s_nucl, ref_e_nucl = map(int, (s_nucl, e_nucl, ref_s_nucl, ref_e_nucl))
                    align_start += s_nucl
                    align_end = align_start + (e_nucl - s_nucl)
                    ref_s += ref_s_nucl
                    ref_e = ref_s + (ref_e_nucl - ref_s_nucl)
                    cigar = cigar.split(':')[-1]

                    strand_direction = 1
                    if strand == '-':
                        s_nucl, e_nucl = e_nucl, s_nucl
                        align_start, align_end = align_end, align_start
                        strand_direction = -1

                    operations = cigar_pattern.findall(cigar)
                    align_len = 0
                    ref_len = 0
                    for op in operations:
                        n_bases, operation = int(op[:-1]), op[-1]
                        if operation == 'S' or operation == 'H':
                            align_start += n_bases
                        elif operation == 'M' or operation == '=' or operation == 'X':
                            for i in range(n_bases):
                                reads_coords[read_name][align_start + (align_len * strand_direction)] = ref_s + ref_len
                                ref_len += 1
                                align_len += 1
                        elif operation == 'D':
                            for i in range(n_bases):
                                reads_coords[read_name][align_start + (align_len * strand_direction)] = ref_s + ref_len
                                ref_len += 1
                            # ref_len += n_bases
                        elif operation == 'I':
                            # align_len += n_bases
                            for i in range(n_bases):
                                reads_coords[read_name][align_start + (align_len * strand_direction)] = ref_s + ref_len
                                align_len += 1
                    break
        return reads_coords


def do(assemblies, reads_fname, monomers_fname, out_dir):
    print("")
    print("*********************************")
    print("Monomer analysis started...")

    reads_mm_structure = get_reads_monomer_structure(reads_fname, monomers_fname, out_dir)

    for assembly in assemblies:
        print("")
        print("Processing %s assembly..." % assembly.label)
        ref_mm_structure, ref_stats = get_ref_monomers(assembly, monomers_fname, out_dir)
        ref_mm_structure.sort(key=lambda x:x[1])

        make_plotly_html(assembly, ref_stats, out_dir)

        _, reads_coords = get_kmers_read_pos(assembly, reads_fname)

        assembly_len = get_fasta_len(assembly.fname)
        reads_monomers = [[] for i in range(len(ref_mm_structure))]
        coverage = [0] * assembly_len
        for read_name, coord_dict in reads_coords.items():
            for (mm_name, mm_start, mm_end) in reads_mm_structure[read_name]:
                ref_i = -1
                ref_i2 = -1
                mm_len = mm_end-mm_start,
                if mm_start in reads_coords[read_name] and mm_end in reads_coords[read_name]:
                    mm_start = reads_coords[read_name][mm_start]
                    mm_end = reads_coords[read_name][mm_end]
                    mm_start, mm_end = min(mm_start, mm_end), max(mm_start, mm_end)
                    ref_i = approx_binary_search(ref_mm_structure, 1, 0, len(ref_mm_structure), mm_start)
                    if mm_end - mm_start > 50:
                        ref_i2 = approx_binary_search(ref_mm_structure, 2, 0, len(ref_mm_structure), mm_end)
                    else:
                        ref_i2 = ref_i
                if ref_i > -1 and ref_i2 == ref_i:
                    coverage[ref_i] += 1
                    reads_monomers[ref_i].append((mm_name, mm_len, mm_start))

        read_support = []
        for i in range(len(ref_mm_structure)):
            if len(reads_monomers[i]) >= MIN_COV:
                read_support.append(sum([1 for m in reads_monomers[i] if m[0] == ref_mm_structure[i][0]])*1.0/coverage[i])
            else:
                read_support.append(1)

        plot_fname = join(out_dir, "report", assembly.name + "_monomer_analysis.png")
        make_plot(plot_fname, "Monomer analysis", assembly.label, xlabel="Position", ylabel="MonomersRatio", plot_values=read_support, plot_color="blue",
                  ymax=1.05, max_x=assembly_len)

        #### UNIT ANALYSIS
        ref_unit_structure, monomers_pattern = analyze_unit_structure(ref_mm_structure)

        unit_occ = defaultdict(int)
        units_fname = join(out_dir, "%s_units.txt") % assembly.name
        with open(units_fname, "w") as f:
            f.write("\t".join(["Unit", "Start", "End", "Monomer sequence\n"]))
            for i in range(len(ref_unit_structure)):
                if not ref_unit_structure[i]:
                    continue
                unit_str = ref_unit_structure[i][2]
                unit_occ[unit_str] += 1
                f.write("\t".join(str(s) for s in [i+1, ref_unit_structure[i][0], ref_unit_structure[i][1], unit_str])+"\n")

        print("Total units: %d, units sequences saved to %s" % (len(ref_unit_structure), units_fname))