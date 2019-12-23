from os.path import join

from scripts.reporting import make_plot
from scripts.utils import get_fasta_len


def calculate_coverage(assembly_len, bed_fname, read_names=None):
    coverage = [0] * assembly_len
    starts = [0] * assembly_len
    ends = [0] * assembly_len
    with open(bed_fname) as f:
        for line in f:
            fs = line.split()
            ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs
            if read_names is not None and read_name not in read_names:
                continue
            ref_s, ref_e, align_start, align_end, read_len = map(int, (ref_s, ref_e, align_start, align_end, read_len))
            starts[ref_s] += 1
            ends[ref_e - 1] += 1
    cur_cov = 0
    for i in range(assembly_len):
        cur_cov += starts[i]
        cur_cov -= ends[i]
        coverage[i] = cur_cov
    return coverage


def do(assemblies, out_dir):
    print("")
    print("*********************************")
    print("Coverage test started...")

    for assembly in assemblies:
        coverage = calculate_coverage(get_fasta_len(assembly.fname), assembly.bed_fname)
        plot_fname = join(out_dir, "report", assembly.name + "_coverage.png")
        make_plot(plot_fname, "Coverage", assembly.label, xlabel="Position", ylabel="Coverage", fill_values=coverage, fill_color="blue")
    print("Coverage test finished.")