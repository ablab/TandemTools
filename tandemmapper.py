#!/usr/bin/env python

import datetime
import os
import sys
from os.path import abspath, isdir, join, basename
import click

import config
from scripts.utils import compress_homopolymers


@click.command()
@click.argument('assembly_fnames', type=click.Path(exists=True), nargs=-1)
@click.option('--nano', 'nano_reads_fname', type=click.Path(exists=True), help='File with ONT reads')
@click.option('--pacbio', 'pacbio_reads_fname', type=click.Path(exists=True), help='File with PacBio CLR reads')
@click.option('-o', 'out_dir', type=click.Path(), required=True, help='Output folder')
@click.option('-t', 'threads', type=click.INT, help='Threads')
@click.option('-l', 'labels', help='Comma separated list of assembly labels')
@click.option('--hifi', 'hifi_reads_fname',  type=click.Path(), help='File with PacBio HiFi reads')
@click.option('-f', '--no-reuse', 'no_reuse', is_flag=True, help='Do not reuse old files')
@click.option('--no-nucl-align', 'no_nucl_alignment', is_flag=True, help='Do not perform nucleotide alignment '
                                                                         '(use with caution)')

def main(assembly_fnames, nano_reads_fname, pacbio_reads_fname, hifi_reads_fname, out_dir, labels, threads,
         no_reuse, no_nucl_alignment):
    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("%s TandemMapper started" % date)

    if nano_reads_fname and pacbio_reads_fname or (not nano_reads_fname and not pacbio_reads_fname):
        print("ERROR! You should specify ONE path to a file with reads (ONT or Pacbio CLR reads)")
        sys.exit(2)

    if nano_reads_fname:
        raw_reads_fname = nano_reads_fname
        config.platform = "nano"
    else:
        raw_reads_fname = pacbio_reads_fname
        config.platform = "pacbio"

    if not assembly_fnames:
        print("ERROR! You should specify at least one assembly file.")
        sys.exit(2)

    from scripts import select_kmers, make_alignments
    from scripts.assembly import Assembly
    out_dir = abspath(out_dir)
    if not isdir(out_dir):
        os.makedirs(out_dir)
    tmp_dir = join(out_dir, "tmp")
    if not isdir(tmp_dir):
        os.makedirs(tmp_dir)

    list_labels = [None] * len(assembly_fnames)
    if labels:
        list_labels = labels.replace('"', '').split(',')
        if len(list_labels) != len(assembly_fnames):
            print("ERROR! Number of labels must correspond to the number of analyzed assemblies")
            sys.exit(2)
    assemblies = [Assembly(assembly_fnames[i], name=list_labels[i], out_dir=out_dir) for i in range(len(assembly_fnames))]

    reads_real_coords = dict()
    if config.platform == "pacbio":
        for assembly in assemblies:
            assembly.real_coords = compress_homopolymers(assembly.fname, assembly.compressed_fname)
        reads_fname = join(out_dir, "compress_" + basename(raw_reads_fname).replace('fastq', 'fasta'))
        reads_real_coords = compress_homopolymers(raw_reads_fname, reads_fname)
    else:
        reads_fname = raw_reads_fname

    select_kmers.do(assemblies, raw_reads_fname, reads_fname, hifi_reads_fname, out_dir, tmp_dir, no_reuse)
    make_alignments.do(assemblies, reads_fname, reads_real_coords, out_dir, threads, no_reuse, no_nucl_alignment)

    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("")
    print("%s TandemMapper finished" % date)


if __name__ == '__main__':
    main()