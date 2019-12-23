import itertools
import subprocess
from collections import defaultdict
from os.path import join, exists, basename

from Bio import SeqIO

from config import *
from scripts.assembly import Assembly
from scripts.reporting import make_plot
from scripts.utils import rev_comp, get_kmers_positions, get_fasta_len, get_kmers

kmc_bin = "/Users/alla/git/quast/external_tools/KMC/macosx/kmc"
kmctools_bin = "/Users/alla/git/quast/external_tools/KMC/macosx/kmc_tools"


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

    make_plot(plot_fname, "Distribution of %s k-mers" % kmers_type, assembly.label, xlabel="Position", ylabel="$\it{k}$-mer counts",
              bar_values=unique_bins, plot_color="blue", max_x=assembly_len, ymax=600)


def run_kmc(seq, out_dir, unique_kmers):
    tmp_seq_fname = join(out_dir, "seq.tmp")
    tmp_fname = join(out_dir, "kmc.tmp")
    with open(tmp_seq_fname, "w") as f:
        f.write(">1\n")
        f.write(seq)
    cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-cx1", "-ci1", "-cs10000", tmp_seq_fname, "kmers", out_dir]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    cmd = [kmctools_bin, "transform", "kmers", "dump", tmp_fname]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
    with open(tmp_fname) as f:
        for line in f:
            kmer = line.strip().split()[0]
            unique_kmers[min(kmer,rev_comp(kmer))]+=1


def get_locally_unique_kmers(assembly, assembly_kmc_db, out_dir, tmp_dir):
    seq = ''
    with open(assembly.fname) as f:
        for record in SeqIO.parse(f, 'fasta'):
            seq = str(record.seq)

    unique_kmers = defaultdict(int)
    for i in range(0, len(seq), CHUNK_SIZE):
        run_kmc(seq[i:i+CHUNK_SIZE], tmp_dir, unique_kmers)

    with open(assembly.all_kmers_fname, "w") as f:
        for i, (k, freq) in enumerate(unique_kmers.items()):
            f.write(">%d_%d\n" % (i, freq))
            f.write("%s\n" % k)
    cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-cx1", "-ci1", assembly.all_kmers_fname, assembly_kmc_db, tmp_dir]
    subprocess.call(cmd)
    return unique_kmers


def filter_kmers(reads_fname, kmers):
    bad_kmers = set()
    j=0
    with open(reads_fname) as handle:
        for record in SeqIO.parse(handle, reads_fname.split('.')[-1]):
            if j%100==0:
                print(j)
            read_seq = str(record.seq)
            read_len = len(read_seq)
            read_counter = defaultdict(int)
            for i in range(read_len - KMER_SIZE):
                kmer = str(read_seq[i:i + KMER_SIZE])
                if kmer in kmers:
                    read_counter[kmer] +=1
                elif rev_comp(kmer) in kmers:
                    read_counter[rev_comp(kmer)] +=1
            for kmer, c in read_counter.items():
                if c > 1:
                    bad_kmers.add(kmer)
            j+=1
    solid_kmers = kmers-bad_kmers
    return solid_kmers


def do(assemblies, reads_fname, out_dir, tmp_dir, no_reuse):
    print("")
    print("*********************************")
    print("K-mers selection started...")

    no_reuse=1
    for assembly in assemblies:
        readkmers_kmc_db = join(tmp_dir, basename(reads_fname))
        readkmers_fname = join(out_dir, basename(reads_fname)+".kmers.txt")
        if not exists(readkmers_fname):
            cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-ci1", "-cs10000", reads_fname, readkmers_kmc_db, tmp_dir]
            print(' '.join(cmd))
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
            cmd = [kmctools_bin, "transform", readkmers_kmc_db, "dump", readkmers_fname]
            print(' '.join(cmd))
            subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        assembly_kmc_db = join(tmp_dir, assembly.name)
        tmp_allkmers_fname = join(tmp_dir, assembly.name+".allkmers.txt")
        tmp_kmers_fname = join(tmp_dir, assembly.name+".kmers.txt")
        if not exists(assembly_kmc_db) or no_reuse:
            cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-ci1", "-cs10000", assembly.fname, assembly_kmc_db, tmp_dir]
            print(' '.join(cmd))
            subprocess.call(cmd)
            cmd = [kmctools_bin, "transform", assembly_kmc_db, "dump", tmp_allkmers_fname]
            print(' '.join(cmd))
            subprocess.call(cmd)

        missed_kmers_db = join(tmp_dir, assembly.name+".missed")
        missed_kmers_fname = join(out_dir, assembly.name+".missed.kmers.txt")
        intersect_kmc_db = join(tmp_dir, assembly.name+".intersect")
        kmers_hist_fname = join(out_dir, assembly.name+".hist.kmers.txt")
        cmd = [kmctools_bin, "simple", assembly_kmc_db, readkmers_kmc_db, "intersect", intersect_kmc_db, "-cs%d" % MAX_OCC_IN_READS, "-ocright"]
        subprocess.call(cmd)
        cmd = [kmctools_bin, "simple", assembly_kmc_db, readkmers_kmc_db, "kmers_subtract", missed_kmers_db]
        print(' '.join(cmd))
        subprocess.call(cmd,)
        cmd = [kmctools_bin, "transform", missed_kmers_db, "dump", missed_kmers_fname]
        print(' '.join(cmd))
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "transform", intersect_kmc_db, "-ci%d" % MIN_OCCURRENCES, "-cx%d" % MAX_OCCURRENCES, "dump", tmp_kmers_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
        cmd = [kmctools_bin, "transform", intersect_kmc_db, "histogram", kmers_hist_fname]
        subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

        hist_values = [0] * (MAX_OCCURRENCES+1)
        missed_kmers = set()
        with open(missed_kmers_fname) as f:
            for line in f:
                #hist_values[0] += 1
                missed_kmers.add(line.split()[0])

        kmers_pos = dict()
        k =19
        positions=[]
        with open(assembly.fname) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                assembly_len = len(record.seq)
                assembly_seq = str(record.seq)
                for i in range(len(record.seq) - k):
                    kmer = assembly_seq[i:i + k]
                    if kmer in missed_kmers:
                        kmers_pos[kmer] = i
                        #print(kmer,i)
                        positions.append(i)
                    elif rev_comp(kmer) in missed_kmers:
                        kmers_pos[rev_comp(kmer)] = i
                        #print(kmer,i)
                        positions.append(i)
        print(positions)
        print("%d kmers do not occur in reads!" % len(missed_kmers))
        plot_fname = join(out_dir, assembly.name + "_missedkmers_hist.png")
        draw_plot(assembly, missed_kmers, plot_fname, "missed")
        continue
        with open(kmers_hist_fname) as f:
            for line in f:
                freq, num_kmers = map(int, line.strip().split())
                hist_values[min(freq, MAX_OCCURRENCES)] += num_kmers
        total_kmers = sum(hist_values[:-1])
        cum_occ = 0

        real_max_occ = MAX_OCCURRENCES
        for i, freq in enumerate(hist_values):
            cum_occ += freq
            if cum_occ >= total_kmers*0.85:
                real_max_occ = i
                break

        print(real_max_occ)
        selected_kmers = set()
        solid_kmers = set()
        with open(tmp_kmers_fname) as f:
            for line in f:
                kmer, occ = line.strip().split()
                if int(occ) <= real_max_occ:
                    selected_kmers.add(kmer)
        with open(assembly.kmers_fname, "w") as out_f:
            for kmer in selected_kmers:
                out_f.write("%s\n" % kmer)
                selected_kmers.add(kmer)
                if local_kmers_counter[kmer] == 1 or local_kmers_counter[rev_comp(kmer)] == 1:
                    solid_kmers.add(kmer)

        solid_kmers=filter_kmers(reads_fname,solid_kmers)
        with open(assembly.solid_kmers_fname, "w") as f:
            for kmer in solid_kmers:
                f.write("%s\n" % kmer)
        print("%d k-mers (%d solid) were selected" % (len(selected_kmers), len(solid_kmers)))

        plot_fname = join(out_dir, assembly.name + "_selected_kmers.png")
        draw_plot(assembly, selected_kmers, plot_fname, "selected")

        plot_fname = join(out_dir, assembly.name + "_kmers_hist.png")
        make_plot(plot_fname, "Histogram of k-mer occurrences", assembly.label, xlabel="Occurrences in reads", ylabel="K-mer counts",
                  bar_values=hist_values, plot_color="green")

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

if __name__ == '__main__':
    '''
    reads_seq = dict()
    freq_kmers = []
    uniq_kmers = []
    with open("rel3_new2/rel3-cenx_solid_kmers.txt") as f:
        for line in f:
            if not line.startswith(">"):
                kmer= line.split()[0]
                uniq_kmers.append(min(kmer,rev_comp(kmer)))
    with open("rel3_new2/freq_kmers.txt") as f:
        for line in f:
            kmer= line.split()[0]
            freq_kmers.append(min(kmer,rev_comp(kmer)))
    uniq_kmers.sort()
    freq_kmers.sort()
    bad_kmers = set()
    import editdistance
    for n,i in enumerate(uniq_kmers):
        if n%100 == 0:
            print(n, len(bad_kmers))
        for j in uniq_kmers[n+1:]:
            #dist = hamming_distance(i,j)
            if editdistance.eval(i,j) < 2:
            #if hamming_distance(i,j) < 2:
                bad_kmers.add(i)
                #print(i,j,dist)
                break
    with open("rel3_new2/cleared_kmers.fasta", "w") as f:
        for k in uniq_kmers:
            if k not in bad_kmers:
                f.write("%s\n" % k)
    for n,i in enumerate(uniq_kmers):
        if i in bad_kmers:
            continue
        if n%100 == 0:
            print(n, len(bad_kmers))
        for j in freq_kmers:
            if editdistance.eval(i,j) < 2:
                bad_kmers.add(i)
                break
    print(len(bad_kmers))
    with open("rel3_new2/good_solid_kmers.fasta", "w") as f:
        for k in uniq_kmers:
            if k not in bad_kmers:
                f.write("%s\n" % k)

    with open("centromeric_reads/centromeric_reads.fasta") as handle:
        from Bio import SeqIO
        for record in SeqIO.parse(handle, 'fasta'):
            from slugify import slugify
            reads_seq[slugify(record.name)] = str(record.seq)
    read_names =set()
    with open("cenX_filtered_reads.txt") as f:
        for line in f:
            read_names.add(line.strip())
    with open("cenX_filtered_reads.fasta", "w") as f:
        for name, seq in reads_seq.items():
            if name in read_names:
                f.write(">%s\n" % name)
                f.write("%s\n" % seq)
    exit()
    bad_kmers = set()
    hifi_kmers = {}
    with open("rel3_new2/rel3-cenx_bad_kmers.txt") as f:
        for line in f:
            bad_kmers.add(line.strip())
    with open("hifi/centromeric_reads_hifi.fasta.kmers.txt") as f:
        for line in f:
            k, freq = line.split()
            hifi_kmers[k]=int(freq)
    a =[]
    for k in bad_kmers:
        if k in hifi_kmers:
            a.append(hifi_kmers[k])
    print(a.count(1), a.count(2), min(a), max(a))

    tmp_db = join(tmp_dir, "freq_kmers")
    tmp_file = join(out_dir, "freq_kmers.txt")
    cmd = [kmc_bin, "-k%d" % KMER_SIZE, "-fm", "-ci1000", "-cs10000", assembly.fname, tmp_db, tmp_dir]
    subprocess.call(cmd)
    cmd = [kmctools_bin, "transform", tmp_db, "dump", tmp_file]
    subprocess.call(cmd, stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))'''
    assembly1 = Assembly(out_dir="rel3_new2", name="T2T7",raw_fname="rel3-cenx",fname="rel3_new2/rel3-cenx.masked.fa")
    assembly1 = Assembly(out_dir="polish_t2t4_50reads_new", name="T2T4_polish",raw_fname="polished-1",fname="polish_t2t4_50reads_new/polished-1.masked.fa")
    assembly1 = Assembly(out_dir="t2t7_05dec", name="T2T7",raw_fname="cenx-t2t7",fname="t2t7_05dec/cenx-t2t7.masked.fa")
    assembly1 = Assembly(out_dir="chm13_fix", name="T2T7",raw_fname="chm13-cenx",fname="chm13_fix/chm13-cenx.masked.fa")
    #assembly1 = Assembly(out_dir="chm13_fix", name="T2T7",raw_fname="chm13-cenx",fname="t2t7_rel3_1.fasta")
    assembly1 = Assembly(out_dir="res_centroflye",raw_fname="cenx-v0-8-3",fname="res_centroflye/cenx-v0-8-3.masked.fa", name="centroFlye")
    assembly1 = Assembly(out_dir="chm13_fix", name="centroFlye",raw_fname="chm13-cenx",fname="cenx_new_polish.fa")
    #assembly1 = Assembly(out_dir="chm13_fix", name="T2T4",raw_fname="chm13-cenx",fname="t2t4_new_polish.fa")

    #assembly2 = Assembly(out_dir="polish_v083_rel2_local_unique_15_new",raw_fname="polished-1",fname="polish_v083_rel2_local_unique_15_new/polished-1.masked.fa", name="centroFlye_polish")
    #assembly1 = Assembly(out_dir="t2t4_new", name="T2T4",raw_fname="cenx-v04", fname="t2t4_new/cenx-v04.masked.fa")
    #assembly2 = Assembly(out_dir="polish_t2t4_50reads", raw_fname="cenx-v04",name="centroFlye_polish")
    #assembly1 = Assembly(out_dir="polish_v083_rel2_and_hifi",raw_fname="polished-1",fname="polish_v083_rel2_and_hifi/polished-1.masked.fa", name="centroFlye_polish")
    #assembly2 = Assembly(out_dir="rel3_new",name="rel3-cenx")
    #assembly1 = Assembly(out_dir="t2t6_new", name="T2T6",raw_fname="cenx-v06",fname="cenX_v06.fasta")

    do([assembly1], "hifi/centromeric_reads_hifi.fasta", "hifi", "hifi/tmp", True)
