import itertools
from os.path import join

import matplotlib
import matplotlib.pyplot as plt


def draw_bimapping_plot(assembly1, assembly2, out_dir):
    alignments_1 = dict()
    with open(assembly1.bed_fname) as f:
        for line in f:
            fs = line.split()
            ref, ref_s, ref_e, read_name = fs[:4]
            ref_s, ref_e = map(int, (ref_s, ref_e))
            alignments_1[read_name] = (int(ref_s),int(ref_e))
    alignments_2 = dict()
    with open(assembly2.bed_fname) as f:
        for line in f:
            fs = line.split()
            ref, ref_s, ref_e, read_name, align_start, align_end, read_len = fs
            ref_s, ref_e, align_start, align_end, read_len = map(int, (ref_s, ref_e, align_start, align_end, read_len))
            alignments_2[read_name] = (int(ref_s),int(ref_e))

    coords_x = []
    coords_y = []
    diff_n = 0
    missed_n = 0
    alignments = []
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.title("Dot plot", fontsize=40)
    plt.xlabel(assembly1.label, fontsize=32)
    plt.ylabel(assembly2.label, fontsize=32)

    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    for read_name in alignments_1:
        if read_name in alignments_2:
            coords_x.append(alignments_1[read_name][0])
            coords_y.append(alignments_2[read_name][0])
            alignments.append((alignments_1[read_name][0], alignments_2[read_name][0]))
            if abs(alignments_1[read_name][0] - alignments_2[read_name][0]) >= 2000:
                diff_n += 1
        else:
            missed_n += 1

    print("Reads mapped to %s: %d, mapped to %s: %d." %
          (assembly1.name, len(alignments_1), assembly2.name, len(alignments_2)))

    plt.xlim(0,max(max(coords_x), max(coords_y))+1000)
    plt.ylim(0,max(max(coords_x), max(coords_y))+1000)
    ax.scatter(coords_x, coords_y, s=12, color="black")
    line = matplotlib.lines.Line2D([0, 1], [0, 1], color='red')
    transform = ax.transAxes
    line.set_transform(transform)
    ax.add_line(line)

    #print("Common mappings: %d, new mappings in %s: %d, new mappings in %s: %d" % (len(coords_x) - diff_n, missed_n, len(alignments_2) - len(coords_x)))
    plot_fname = join(out_dir, "report", "%s_vs_%s.png" % (assembly1.name, assembly2.name))
    plt.savefig(plot_fname)
    print("Bi-mapping plot for %s and %s is saved to %s" % (assembly1.name, assembly2.name, plot_fname))


def do(assemblies, out_dir):
    if len(assemblies) < 2:
        return
    print("")
    print("*********************************")
    print("Pairwise comparison started...")
    for (assembly1, assembly2) in itertools.combinations(assemblies, 2):
        draw_bimapping_plot(assembly1, assembly2, out_dir)
    print("Pairwise comparison finished.")

