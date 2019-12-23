#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs Minimap2 and parses its output
"""

from __future__ import absolute_import
from __future__ import division
import os
from collections import namedtuple
import subprocess
import logging

import ext_tools.Flye.flye.utils.fasta_parser as fp
from ext_tools.Flye.flye.utils.utils import which
from ext_tools.Flye.flye.utils.sam_parser import AlignmentException, preprocess_sam
from ext_tools.Flye.flye.six import iteritems
from ext_tools.Flye.flye.six.moves import range


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"


ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])


def check_binaries():
    if not which(MINIMAP_BIN):
        raise AlignmentException("Minimap2 is not installed")
    if not which("sort"):
        raise AlignmentException("UNIX sort utility is not available")


def make_alignment(reference_file, reads_file, num_proc,
                   work_dir, platform, out_alignment, reference_mode,
                   sam_output):
    """
    Runs minimap2 and sorts its output
    """
    minimap_ref_mode = {False: "ava", True: "map"}
    minimap_reads_mode = {"nano": "ont", "pacbio": "pb"}
    mode = minimap_ref_mode[reference_mode] + "-" + minimap_reads_mode[platform]

    _run_minimap(reference_file, reads_file, num_proc, mode,
                 out_alignment, sam_output)

    if sam_output:
        preprocess_sam(out_alignment, work_dir)


def get_contigs_info(contigs_file):
    contigs_info = {}
    contigs_fasta = fp.read_sequence_dict(contigs_file)
    for ctg_id, ctg_seq in iteritems(contigs_fasta):
        contig_type = ctg_id.split("_")[0]
        contigs_info[ctg_id] = ContigInfo(ctg_id, len(ctg_seq),
                                          contig_type)

    return contigs_info


def shift_gaps(seq_trg, seq_qry):
    """
    Shifts all ambigious query gaps to the right
    """
    lst_trg, lst_qry = list("$" + seq_trg + "$"), list("$" + seq_qry + "$")
    is_gap = False
    gap_start = 0
    for i in range(len(lst_trg)):
        if is_gap and lst_qry[i] != "-":
            is_gap = False
            swap_left = gap_start - 1
            swap_right = i - 1

            while (swap_left > 0 and swap_right >= gap_start and
                   lst_qry[swap_left] == lst_trg[swap_right]):
                lst_qry[swap_left], lst_qry[swap_right] = \
                            lst_qry[swap_right], lst_qry[swap_left]
                swap_left -= 1
                swap_right -= 1

        if not is_gap and lst_qry[i] == "-":
            is_gap = True
            gap_start = i

    return "".join(lst_qry[1 : -1])


def get_uniform_alignments(alignments, seq_len):
    """
    Leaves top alignments for each position within contig
    assuming uniform coverage distribution
    """
    def _get_median(lst):
        if not lst:
            raise ValueError("_get_median() arg is an empty sequence")
        sorted_list = sorted(lst)
        if len(lst) % 2 == 1:
            return sorted_list[len(lst) // 2]
        else:
            mid1 = sorted_list[(len(lst) // 2) - 1]
            mid2 = sorted_list[(len(lst) // 2)]
            return (mid1 + mid2) / 2

    WINDOW = 100
    MIN_COV = 10
    COV_RATE = 1.25

    #split contig into windows, get median read coverage over all windows and
    #determine the quality threshold cutoffs for each window
    wnd_primary_cov = [0 for _ in range(seq_len // WINDOW + 1)]
    wnd_aln_quality = [[] for _ in range(seq_len // WINDOW + 1)]
    wnd_qual_thresholds = [1.0 for _ in range(seq_len // WINDOW + 1)]
    for aln in alignments:
        for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW):
            if not aln.is_secondary:
                wnd_primary_cov[i] += 1
            wnd_aln_quality[i].append(aln.err_rate)

    #for each window, select top X alignmetns, where X is the median read coverage
    cov_threshold = max(int(COV_RATE * _get_median(wnd_primary_cov)), MIN_COV)
    for i in range(len(wnd_aln_quality)):
        if len(wnd_aln_quality[i]) > cov_threshold:
            wnd_qual_thresholds[i] = sorted(wnd_aln_quality[i])[cov_threshold]

    #for each alignment, count in how many windows it passes the threshold
    filtered_alignments = []
    total_sequence = 0
    filtered_sequence = 0
    for aln in alignments:
        good_windows = 0
        total_windows = aln.trg_end // WINDOW - aln.trg_start // WINDOW
        total_sequence += aln.trg_end - aln.trg_start
        for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW):
            if aln.err_rate <= wnd_qual_thresholds[i]:
                good_windows += 1

        if good_windows > total_windows // 2:
            filtered_alignments.append(aln)
            filtered_sequence += aln.trg_end - aln.trg_start

    #filtered_reads_rate = 1 - float(len(filtered_alignments)) / len(alignments)
    #filtered_seq_rate = 1 - float(filtered_sequence) / total_sequence
    #logger.debug("Filtered {0:7.2f}% reads, {1:7.2f}% sequence"
    #                .format(filtered_reads_rate * 100, filtered_seq_rate * 100))

    return filtered_alignments


def split_into_chunks(fasta_in, chunk_size):
    out_dict = {}
    for header, seq in iteritems(fasta_in):
        #print len(seq)
        for i in range(0, max(len(seq) // chunk_size, 1)):
            chunk_hdr = "{0}$chunk_{1}".format(header, i)
            start = i * chunk_size
            end = (i + 1) * chunk_size
            if len(seq) - end < chunk_size:
                end = len(seq)

            #print(start, end)
            out_dict[chunk_hdr] = seq[start : end]

    return out_dict


def merge_chunks(fasta_in, fold_function=lambda l: "".join(l)):
    """
    Merges sequence chunks. Chunk names are in format `orig_name$chunk_id`.
    Each chunk is as dictionary entry. Value type is arbitrary and
    one can supply a custom fold function
    """
    def name_split(h):
        orig_hdr, chunk_id = h, 0
        return orig_hdr, 0

    out_dict = {}
    cur_seq = []
    cur_contig = None
    for hdr in sorted(fasta_in, key=name_split):
        orig_name, dummy_chunk_id = name_split(hdr)
        if orig_name != cur_contig:
            if cur_contig != None:
                out_dict[cur_contig] = fold_function(cur_seq)
            cur_seq = []
            cur_contig = orig_name
        cur_seq.append(fasta_in[hdr])

    if cur_seq:
        out_dict[cur_contig] = fold_function(cur_seq)

    return out_dict


def _run_minimap(reference_file, reads_files, num_proc, mode, out_file,
                 sam_output):
    cmdline = [MINIMAP_BIN, reference_file]
    cmdline.extend(reads_files)
    cmdline.extend(["-x", mode, "-t", str(num_proc)])
    if sam_output:
        #a = SAM output, p = min primary-to-seconday score
        #N = max secondary alignments
        cmdline.extend(["-a", "-p", "0.5", "-N", "10"])

    try:
        devnull = open(os.devnull, "wb")
        #logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline, stderr=devnull,
                              stdout=open(out_file, "wb"))
    except (subprocess.CalledProcessError, OSError) as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AlignmentException(str(e))
