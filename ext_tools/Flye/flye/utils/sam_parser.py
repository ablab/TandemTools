#(c) 2019 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Provides multithreaded parser for SAM files
"""

from __future__ import absolute_import
from __future__ import division

import os
import re
import sys
from collections import namedtuple, defaultdict
import subprocess
import logging
import multiprocessing
import ctypes

#In Python2, everything is bytes (=str)
#In Python3, we are doing IO in bytes, but everywhere else strngs = unicode
from ext_tools.Flye.flye.six import iteritems
from ext_tools.Flye.flye.six.moves import range

if sys.version_info < (3, 0):
    from string import maketrans
    _STR = lambda x: x
    _BYTES = lambda x: x
else:
    maketrans = bytes.maketrans
    _STR = bytes.decode
    _BYTES = str.encode


import ext_tools.Flye.flye.utils.fasta_parser as fp

logger = logging.getLogger()

Alignment = namedtuple("Alignment", ["qry_id", "trg_id", "qry_start", "qry_end",
                                     "qry_sign", "qry_len", "trg_start",
                                     "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate",
                                     "is_secondary"])


class AlignmentException(Exception):
    pass


class PafHit(object):
    """
    Stores paf alignment
    """
    __slots__ = ("query", "query_length", "query_start", "query_end",
                 "target", "target_length", "target_start", "target_end")
    def __init__(self, raw_hit):
        hit = raw_hit.split()

        self.query = hit[0]
        self.query_length = int(hit[1])
        self.query_start = int(hit[2])
        self.query_end = int(hit[3])

        self.target = hit[5]
        self.target_length = int(hit[6])
        self.target_start = int(hit[7])
        self.target_end = int(hit[8])

    def query_mapping_length(self):
        return self.query_end - self.query_start + 1

    def target_mapping_length(self):
        return self.target_end - self.target_start + 1

    def query_left_overhang(self):
        return self.query_start

    def query_right_overhang(self):
        return self.query_length - self.query_end + 1

    def target_left_overhang(self):
        return self.target_start

    def target_right_overhang(self):
        return self.target_length - self.target_end + 1


def read_paf(filename):
    """
    Streams out paf alignments
    """
    with open(filename, "rb") as f:
        for raw_hit in f:
            yield PafHit(_STR(raw_hit))


def read_paf_grouped(filename):
    """
    Outputs chunks of alignments for each (query, target)pair.
    Assumes that PAF alignment is already sorted by query.
    """
    prev_hit = None
    target_hits = defaultdict(list)
    for hit in read_paf(filename):
        if prev_hit is not None and hit.query != prev_hit.query:
            for trg in sorted(target_hits):
                yield target_hits[trg]
            target_hits = defaultdict(list)

        target_hits[hit.target].append(hit)
        prev_hit = hit

    if len(target_hits):
        for trg in sorted(target_hits):
            yield target_hits[trg]


class SynchronizedSamReader(object):
    """
    Parses SAM file in multiple threads.
    """
    def __init__(self, sam_alignment, reference_fasta,
                 max_coverage=None, use_secondary=False):
        #will not be changed during exceution, each process has its own copy
        self.aln_path = sam_alignment
        self.aln_file = None
        self.ref_fasta = {_BYTES(h) : _BYTES(s)
                          for (h, s) in iteritems(reference_fasta)}
        self.change_strand = True
        self.max_coverage = max_coverage
        self.seq_lengths = {}
        self.use_secondary = use_secondary
        self.cigar_parser = None
        self.processed_contigs = None

        #reading SAM header
        if not os.path.exists(self.aln_path):
            raise AlignmentException("Can't open {0}".format(self.aln_path))

        with open(self.aln_path, "rb") as f:
            for line in f:
                if not line or not _is_sam_header(line):
                    break
                if line.startswith(b"@SQ"):
                    seq_name = None
                    seq_len = None
                    for tag in line.split():
                        if tag.startswith(b"SN"):
                            seq_name = tag[3:]
                        if tag.startswith(b"LN"):
                            seq_len = int(tag[3:])
                    if seq_name and seq_len:
                        self.seq_lengths[_STR(seq_name)] = seq_len

        #will be shared between processes
        self.lock = multiprocessing.Lock()
        self.eof = multiprocessing.Value(ctypes.c_bool, False)
        self.position = multiprocessing.Value(ctypes.c_longlong, 0)

    def init_reading(self):
        """
        Call from the reading process, initializing local variables
        """
        self.aln_file = open(self.aln_path, "rb")
        self.processed_contigs = set()
        self.cigar_parser = re.compile(b"[0-9]+[MIDNSHP=X]")

    def stop_reading(self):
        """
        Call when the reading is done
        """
        self.aln_file.close()

    def is_eof(self):
        return self.eof.value

    def parse_cigar(self, cigar_str, read_str, ctg_name, ctg_pos):
        ctg_str = self.ref_fasta[ctg_name]
        trg_seq = []
        qry_seq = []
        trg_start = ctg_pos - 1
        trg_pos = ctg_pos - 1
        qry_start = 0
        qry_pos = 0

        left_hard = True
        left_soft = True
        hard_clipped_left = 0
        hard_clipped_right = 0
        soft_clipped_left = 0
        soft_clipped_right = 0
        for token in self.cigar_parser.findall(cigar_str):
            size, op = int(token[:-1]), token[-1:]
            if op == b"H":
                if left_hard:
                    qry_start += size
                    hard_clipped_left += size
                else:
                    hard_clipped_right += size
            elif op == b"S":
                qry_pos += size
                if left_soft:
                    soft_clipped_left += size
                else:
                    soft_clipped_right += size
            elif op == b"M":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                qry_pos += size
                trg_pos += size
            elif op == b"I":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(b"-" * size)
                qry_pos += size
            elif op == b"D":
                qry_seq.append(b"-" * size)
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                trg_pos += size
            else:
                raise AlignmentException("Unsupported CIGAR operation: " + str(op))
            left_hard = False
            if op != b"H":
                left_soft = False

        trg_seq = b"".join(trg_seq)
        qry_seq = b"".join(qry_seq)
        matches = 0
        for i in range(len(trg_seq)):
            if trg_seq[i] == qry_seq[i]:
                matches += 1
        err_rate = 1 - float(matches) / len(trg_seq)

        trg_end = trg_pos
        qry_end = qry_pos + hard_clipped_left
        qry_len = qry_end + hard_clipped_right
        qry_start += soft_clipped_left
        qry_end -= soft_clipped_right

        return (trg_start, trg_end, len(ctg_str), trg_seq,
                qry_start, qry_end, qry_len, qry_seq, err_rate)

    def get_chunk(self):
        """
        Alignment file is expected to be sorted!
        """

        chunk_buffer = []
        parsed_contig = None

        with self.lock:
            self.aln_file.seek(self.position.value)
            if self.eof.value:
                return None, []

            current_contig = None
            while True:
                self.position.value = self.aln_file.tell()
                line = self.aln_file.readline()
                if not line: break
                if _is_sam_header(line): continue

                tokens = line.strip().split()
                if len(tokens) < 11:
                    continue
                    #raise AlignmentException("Error reading SAM file")

                read_contig = tokens[2]
                flags = int(tokens[1])
                is_unmapped = flags & 0x4
                is_secondary = flags & 0x100
                is_supplementary = flags & 0x800    #allow supplementary

                #if is_unmapped or is_secondary: continue
                if is_unmapped: continue
                if is_secondary and not self.use_secondary: continue
                if read_contig in self.processed_contigs:
                    raise AlignmentException("Alignment file is not sorted")

                if read_contig != current_contig:
                    prev_contig = current_contig
                    current_contig = read_contig

                    if prev_contig is not None:
                        self.processed_contigs.add(prev_contig)
                        parsed_contig = prev_contig
                        break
                    else:
                        chunk_buffer = [tokens]
                else:
                    chunk_buffer.append(tokens)

            if not parsed_contig:
                self.eof.value = True
                parsed_contig = current_contig
        #end with

        sequence_length = 0
        alignments = []
        for tokens in chunk_buffer:
            read_id = tokens[0]
            read_contig = tokens[2]
            cigar_str = tokens[5]
            qry_seq = tokens[10]
            trg_seq = tokens[9]
            ctg_pos = int(tokens[3])
            flags = int(tokens[1])
            is_reversed = flags & 0x16
            is_secondary = flags & 0x100

            if qry_seq == b"*":
                raise Exception("Error parsing SAM: record without read sequence")

            trg_start = ctg_pos
            qry_start = 0
            qry_end = len(qry_seq) - qry_seq.count(b'-')
            trg_end = trg_start + len(trg_seq) - trg_seq.count(b'-')
            qry_len= len(qry_seq) - qry_seq.count(b'-')
            trg_len = len(self.ref_fasta[read_contig])

            #(trg_start, trg_end, trg_len, trg_seq,
            #qry_start, qry_end, qry_len, qry_seq, err_rate) = \
            #        self.parse_cigar(cigar_str, read_str, read_contig, ctg_pos)

            #OVERHANG = cfg.vals["read_aln_overhang"]
            #if (float(qry_end - qry_start) / qry_len > self.min_aln_rate or
            #        trg_start < OVERHANG or trg_len - trg_end < OVERHANG):
            matches = 0
            for i in range(len(trg_seq)):
                if trg_seq[i] == qry_seq[i]:
                    matches += 1
            err_rate = 1 - float(matches) / len(trg_seq)
            aln = Alignment(_STR(read_id), _STR(read_contig),
                            qry_start, qry_end, "-" if is_reversed else "+", qry_len,
                            trg_start, trg_end, "+", trg_len,
                            _STR(qry_seq), _STR(trg_seq),
                            err_rate, is_secondary)
            alignments.append(aln)

            sequence_length += qry_end - qry_start
            #In rare cases minimap2 does not output SQ tag, so need to check
            if _STR(parsed_contig) in self.seq_lengths:
                contig_length = self.seq_lengths[_STR(parsed_contig)]
                if sequence_length // contig_length > self.max_coverage:
                    break

        if parsed_contig is None:
            return None, []
        return _STR(parsed_contig), alignments


def preprocess_sam(sam_file, work_dir):
    """
    Proprocesses minimap2 output by adding SEQ
    to secondary alignments, removing
    unaligned reads and then sorting
    file by reference sequence id
    """
    expanded_sam = sam_file + "_expanded"
    merged_file = sam_file + "_merged"
    sorted_file = sam_file + "_sorted"

    #puting SAM headers to the final postprocessed file first
    with open(sam_file, "rb") as hdr_in, open(merged_file, "wb") as fout:
        for line in hdr_in:
            if not _is_sam_header(line):
                break
            fout.write(line)

    #adding SEQ fields to secondary alignments
    with open(sam_file, "rb") as fin, open(expanded_sam, "wb") as fout:
        prev_id = None
        prev_seq = None
        primary_reversed = None
        for line in fin:
            if _is_sam_header(line):
                continue

            tokens = line.strip().split()
            flags = int(tokens[1])
            is_unmapped = flags & 0x4
            is_secondary = flags & 0x100
            is_supplementary = flags & 0x800
            is_reversed = flags & 0x16

            if is_unmapped:
                continue

            read_id, cigar_str, read_seq = tokens[0], tokens[5], tokens[9]
            has_hard_clipped = b"H" in cigar_str

            #Checking format assumptions
            if has_hard_clipped:
                if is_secondary:
                    raise Exception("Secondary alignment with hard-clipped bases")
                if not is_supplementary:
                    raise Exception("Primary alignment with hard-clipped bases")
            if not is_secondary and read_seq == b"*":
                raise Exception("Missing SEQ for non-secondary alignment")

            if read_seq == b"*":
                if read_id != prev_id:
                    raise Exception("SAM file is not sorted by read names")
                if is_reversed == primary_reversed:
                    tokens[9] = prev_seq
                else:
                    tokens[9] = fp.reverse_complement_bytes(prev_seq)

            #Assuming that the first read alignmnent in SAM is primary
            elif prev_id != read_id:
                if has_hard_clipped:
                    raise Exception("Hard clipped bases in the primamry read")
                prev_id = read_id
                prev_seq = read_seq
                primary_reversed = is_reversed

            fout.write(b"\t".join(tokens) + b"\n")

    #don't need the original SAM anymore, cleaning up space
    os.remove(sam_file)

    #logger.debug("Sorting alignment file")
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    subprocess.check_call(["sort", "-k", "3,3", "-T", work_dir, expanded_sam],
                          stdout=open(sorted_file, "wb"), env=env)

    #don't need the expanded file anymore
    os.remove(expanded_sam)

    #appending to the final file, that already contains headers
    with open(sorted_file, "rb") as sort_in, open(merged_file, "ab") as fout:
        for line in sort_in:
            if not _is_sam_header(line):
                fout.write(line)

    os.remove(sorted_file)
    os.rename(merged_file, sam_file)


def _is_sam_header(line):
    return line[:3] in [b"@PG", b"@HD", b"@SQ", b"@RG", b"@CO"]
