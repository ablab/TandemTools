import os
import subprocess
import sys
from os.path import exists, join, dirname

import config
from ext_tools.Flye.flye.polishing.polish import polish
from scripts import select_kmers
from scripts.make_alignments import make_flye
from scripts.utils import get_fasta_len, get_flye_cfg_fname, get_ext_tools_dir

POLISH_BIN = join(get_ext_tools_dir(), "Flye", "bin", "flye-polish")


def run_flye_polish(assembly_fname, reads_fname, out_dir, kmers_fname):
    platform = "nano"
    polish(assembly_fname, reads_fname, out_dir, kmers_fname, get_fasta_len(assembly_fname),
           2, config.MAX_THREADS, platform, get_flye_cfg_fname())


def do(assemblies, reads_fname, hifi_reads_fname, out_dir, tmp_dir, no_reuse):
    print("")
    print("*********************************")
    print("Running polishing module...")
    out_dir = join(out_dir, "polished")
    if not exists(out_dir):
        os.makedirs(out_dir)
    try:
        make_flye()
    except:
        print('Failed to compile Flye! Please try to compile it manually: create %s folder and run "make" in %s'
              % (dirname(POLISH_BIN), dirname(dirname(POLISH_BIN))))
        sys.exit(2)
    for i in range(4):
        select_kmers.do(assemblies, reads_fname, hifi_reads_fname, out_dir, tmp_dir, no_reuse=True, only_polish=True)
        for assembly in assemblies:
            print("Polishing genome (%d/%d)" % (i+1, 4))
            assembly.raw_fname = polish(assembly.raw_fname, reads_fname, out_dir, assembly.kmers_fname,
                                        get_fasta_len(assembly.raw_fname), config.MAX_THREADS, "nano", get_flye_cfg_fname(), i)
    print("Polished assemblies saved to %s" % out_dir)