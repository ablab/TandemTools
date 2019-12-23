from os.path import join

import numpy as np

from config import *
from scripts.reporting import make_plot
from scripts.utils import get_fasta_len


def do(assemblies, out_dir):
    print("")
    print("*********************************")
    print("Breakpoint analysis started...")

    for assembly in assemblies:
        assembly_len = get_fasta_len(assembly.fname)
        starts = [0] * assembly_len
        ends = [0] * assembly_len
        ideal_starts = [0] * assembly_len
        ideal_ends = [0] * assembly_len
        coverage = [0] * assembly_len
        ideal_coverage = [0] * assembly_len
        tips = [0] * assembly_len
        with open(assembly.bed_fname) as f:
            for line in f:
                fs = line.split()
                ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs
                ref_s, ref_e, align_start, align_end, read_len = map(int, (ref_s, ref_e, align_start, align_end, read_len))

                align_start, align_end = min(align_start, align_end), max(align_start, align_end)
                tips[ref_s] += 1
                tips[ref_e-1] += 1
                starts[ref_s] += 1
                ends[ref_e - 1] += 1
                ideal_starts[max(0,ref_s-align_start)] += 1
                ideal_ends[min(assembly_len-1,ref_e+read_len - align_end-1)] += 1

        cur_cov = 0
        ideal_cur_cov = 0
        uncovered_regions = []
        prev_s, prev_e = -1,-1
        for i in range(assembly_len):
            cur_cov += starts[i]
            cur_cov -= ends[i]
            ideal_cur_cov += ideal_starts[i]
            ideal_cur_cov -= ideal_ends[i]
            ideal_coverage[i] = ideal_cur_cov
            coverage[i] = cur_cov
            if cur_cov < MIN_BP_COV:
                if prev_s != -1:
                    prev_e = i
                else: prev_s = i
            elif prev_s != -1 and prev_e != -1:
                uncovered_regions.append((prev_s,prev_e))
                prev_s,prev_e=-1,-1
            else:
                prev_s,prev_e=-1,-1

        factor = 2
        step = BP_WINDOW_SIZE//factor
        real_bp_ratio = [sum(tips[i:i+BP_WINDOW_SIZE])*1.0/max(1,coverage[i]+sum(starts[i+1:i+BP_WINDOW_SIZE]))
                        if max(coverage[i:i+BP_WINDOW_SIZE]) >= MIN_BP_COV else 0
                        for i in range(0, len(coverage), step)]
        ideal_bp_ratio = [(sum(ideal_starts[i:i+BP_WINDOW_SIZE])+sum(ideal_ends[i:i+BP_WINDOW_SIZE]))*1.1/max(1,ideal_coverage[i]+sum(ideal_starts[i+1:i+BP_WINDOW_SIZE]))
                        if max(ideal_coverage[i:i+BP_WINDOW_SIZE]) >= MIN_BP_COV else 0
                        for i in range(0, len(coverage), step)]

        def running_mean(data):
            cumsum = np.cumsum(np.insert(data, 0, 0))
            return (cumsum[MOVING_AVG_WINDOW_SIZE:] - cumsum[:-MOVING_AVG_WINDOW_SIZE]) / 10

        for i in range(2):
            ideal_bp_ratio[i], ideal_bp_ratio[-i-1] = 0, 0
            real_bp_ratio[i], real_bp_ratio[-i-1] = 0, 0
        real_vals = [min(1,v) for v in running_mean(real_bp_ratio)]
        ideal_vals = [min(1,v) for v in running_mean(ideal_bp_ratio)]

        plot_fname = join(out_dir, "report", assembly.name + "_bp_analysis.png")
        uncovered_bars = [(r[0] / step, r[1] / step) for r in uncovered_regions if (r[1] / step - r[0] / step) > 10]
        make_plot(plot_fname, "Breakpoint", assembly.label, xlabel="position", ylabel="breakpointRatio", plot_values=real_vals, plot_color="red",
                  fill_values=ideal_vals, fill_color="gray", ymax=1, max_x=assembly_len, bg_bars=uncovered_bars)
    print("Breakpoint analysis finished.")