#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Configuration file for the Python part of the pipeline
"""

import os

vals = {
        "pkg_root" : os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "pipeline_version" : 2,

        #additional configuration files for binary modules
        "bin_cfg" : {
            "raw" : "config/bin_cfg/asm_raw_reads.cfg",
            "corrected" : "config/bin_cfg/asm_corrected_reads.cfg",
            "subasm" : "config/bin_cfg/asm_subasm.cfg"
        },

        #assembly parameters
        "big_genome_kmer" : 29000000,
        "kmer_size" : {
            "raw" : [15, 17],
            "corrected" : [17, 17],
            "subasm" : [31, 31]
        },
        "min_overlap_range" : {
            "raw" : [1000, 5000],
            "corrected" : [1000, 5000],
            "subasm" : [1000, 1000]
        },
        #"reduced_asm_cov"  : 40,

        #polishing
        "simple_kmer_length" : 4,
        "solid_kmer_length" : 10,
        "max_bubble_length" : 500,
        "max_bubble_branches" : 50,
        "max_read_coverage" : 1000,
        "min_polish_aln_len" : 500,

        #final coverage filtering
        "relative_minimum_coverage" : 5,
        "hard_minimum_coverage" : 3,

        "err_modes" : {
            "pacbio" : {
                "subs_matrix" : "config/bin_cfg/pacbio_substitutions.mat",
                "hopo_matrix" : "config/bin_cfg/pacbio_homopolymers.mat",
                "solid_missmatch" : 0.2,
                "solid_indel" : 0.2,
                "max_aln_error" : 0.25
            },
            "nano" : {
                "subs_matrix" : "config/bin_cfg/nano_r94_substitutions.mat",
                "hopo_matrix" : "config/bin_cfg/nano_r94_homopolymers.mat",
                "solid_missmatch" : 0.2,
                "solid_indel" : 0.2,
                "max_aln_error" : 0.25
            },
        },

        #scaffolding
        "scaffold_gap" : 100
    }
