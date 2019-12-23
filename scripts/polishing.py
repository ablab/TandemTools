import os
import subprocess
from os.path import exists, join

import config
from ext_tools.Flye.flye.polishing.polish import polish
from scripts.utils import get_fasta_len, get_flye_cfg_fname


def run_flye_polish(assembly_fname, reads_fname, out_dir, kmers_fname):
    platform = "nano"
    polish(assembly_fname, reads_fname, out_dir, kmers_fname, get_fasta_len(assembly_fname),
           2, config.MAX_THREADS, platform, get_flye_cfg_fname(), output_progress=True)


def do(assemblies, reads_fname, out_dir, no_reuse):
    print("")
    print("*********************************")
    print("Running polishing module...")
    out_dir = join(out_dir, "polished")
    if not exists(out_dir):
        os.makedirs(out_dir)
    for assembly in assemblies:
        run_flye_polish(assembly.fname, reads_fname, out_dir, assembly.kmers_fname)
    print("Polished assemblies saved to %s" % out_dir)