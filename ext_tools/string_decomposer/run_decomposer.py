#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import os
from os import listdir
from os.path import join, isfile
import sys
import argparse

import subprocess
from subprocess import check_output

import re
import edlib

def edist(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = str(ar[0]), str(ar[1])
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    if ed == -1:
        return 0
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai*100


def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            records[r] = records[r].upper() 
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def add_rc_monomers(monomers):
    res = []
    for m in monomers:
        res.append(m)
        res.append(make_record(m.seq.reverse_complement(), m.name + "'", m.id + "'"))
    return res

def convert_read(decomposition, read, monomers):
    res = []
    for d in decomposition:
        monomer, start, end = d["m"], d["start"], d["end"]
        scores = {}
        for m in monomers:
            score = aai([read.seq[start:end + 1], m.seq])
            scores[m.name] = score
        
        if monomer == None:
            for s in scores:
                if monomer == None or scores[s] > scores[monomer]:
                    monomer = s
        secondbest, secondbest_score = None, -1 
        for m in scores:
            if m != monomer and abs(scores[m] - scores[monomer]) < 5:
                if not secondbest or secondbest_score < scores[m]:
                    secondbest, secondbest_score = m, scores[m]
        res.append({"m": monomer, "start": str(d["start"]), "end": str(d["end"]), "score": scores[monomer], \
                                "second_best": str(secondbest), "second_best_score": secondbest_score, "alt": scores, "q": "+"})

    window = 2
    for i in range(len(res)):
        sm, cnt = 0, 0
        for j in range(i - window, i + window + 1):
            if j >= 0 and j < len(res):
                sm += res[j]["score"]
                cnt += 1
        if sm/cnt < 80:
            res[i]["q"] = "?"

    return res

def print_read(fout, fout_alt, dec, read, monomers, identity_th):
    dec = convert_read(dec, read, monomers)
    for d in dec:
        if d["score"] >= identity_th:
            fout.write("\t".join([read.name, d["m"], d["start"], d["end"], "{:.2f}".format(d["score"]), \
                                                    d["second_best"], "{:.2f}".format(d["second_best_score"]), d["q"]]) + "\n")
            for a in d["alt"]:
                star = "-"
                if a == d["m"]:
                    star = "*"
                fout_alt.write("\t".join([read.name, a, d["start"], d["end"], "{:.2f}".format(d["alt"][a]), star]) + "\n")

def convert_tsv(decomposition, reads, monomers, outfile, identity_th):
    with open(outfile[:-len(".tsv")] + "_alt.tsv", "w") as fout_alt:
        with open(outfile, "w") as fout:
            cur_dec = []
            prev_read = None
            for ln in decomposition.split("\n")[:-1]:
                read, monomer, start, end = ln.split("\t")[:4]
                if read != prev_read and prev_read != None:
                    print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers, identity_th)
                    cur_dec = []
                prev_read = read
                start, end = int(start), int(end)
                cur_dec.append({"m": monomer, "start": start, "end": end})
            if len(cur_dec) > 0:
                print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers, identity_th)

def run(sequences, monomers, num_threads, scoring):
    ins, dels, mm, match = scoring.split(",")
    p = os.path.abspath(__file__)
    sd_exec_file = p[:-len("run_decomposer.py")] + "/src/dp"
    #print("Run", sd_exec_file, " with parameters ", sequences, monomers, num_threads, scoring, file=sys.stderr)
    out = check_output([sd_exec_file, sequences, monomers, num_threads, "5000", ins, dels, mm, match])
    #print("Calculated raw decomposition", file=sys.stderr)
    return out.decode("utf-8")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Decomposes string into blocks alphabet')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('-t', '--threads',  help='number of threads (by default 1)', default="1", required=False)
    parser.add_argument('-o', '--out-file',  help='output tsv-file (by default final_decomposition.tsv)', default="final_decomposition.tsv", required=False)
    parser.add_argument('-i', '--min-identity',  \
                         help='only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0)', type=int, default=0, required=False)
    parser.add_argument('-s', '--scoring', \
                         help='set scoring scheme for SD in the format "insertion,deletion,match,mismatch" (by default "-1,-1,-1,1")', default="-1,-1,-1,1", required=False)
    parser.add_argument('-r', '--raw',  help='save initial monomer decomposition to raw_decomposition.tsv (by default False)', action="store_true")

    args = parser.parse_args()
    raw_decomposition = run(args.sequences, args.monomers, args.threads, args.scoring)
    if args.raw:
        with open("raw_decomposition.tsv", "w") as fout:
            fout.write(raw_decomposition)

    reads = load_fasta(args.sequences, "map")
    monomers = load_fasta(args.monomers)
    monomers = add_rc_monomers(monomers)
    convert_tsv(raw_decomposition, reads, monomers, args.out_file, int(args.min_identity))