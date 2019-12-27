#!/usr/bin/env python

import datetime
import os
import sys
from os.path import abspath, isdir, join
import click


@click.command()
@click.argument('assembly_fnames', type=click.Path(exists=True), nargs=-1)
@click.option('-l', 'labels', help='Comma separated list of assembly labels')
@click.option('-r', 'reads_fname', type=click.Path(exists=True), required=True, help='File with reads')
@click.option('--hi-fi', 'hifi_reads_fname',  type=click.Path(), help='File with PacBio HiFi reads (optional)')
@click.option('-o', 'out_dir',  type=click.Path(), help='Output folder')
@click.option('-t', 'threads', type=click.INT, help='Threads')
@click.option('-f', '--no-reuse', 'no_reuse', is_flag=True, help='Do not reuse old files')
def main(assembly_fnames, reads_fname, hifi_reads_fname, out_dir, labels, threads, no_reuse):
    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("%s TandemMapper started" % date)

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
    assemblies = [Assembly(assembly_fnames[i], out_dir, assembly_fnames[i], name=list_labels[i]) for i in range(len(assembly_fnames))]
    select_kmers.do(assemblies, reads_fname, hifi_reads_fname, out_dir, tmp_dir, no_reuse)
    make_alignments.do(assemblies, reads_fname, out_dir, threads, no_reuse)

    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("")
    print("%s TandemMapper finished" % date)


if __name__ == '__main__':
    main()