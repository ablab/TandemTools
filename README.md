## Quick start
```shell
    tandemquast.py --nano test_data/simulated_reads.fasta test_data/simulated_polished.fa -o simulated_res
```

## Introduction

**TandemTools** package includes **TandemQUAST** tool for evaluating and improving assemblies of extra-long tandem repeats (ETR) and **TandemMapper** tool for mapping long error-prone reads to ETRs.

Note: TandemTools is designed specifically for ETR (range in length from hundreds of thousands to millions of nucleotides). It is strongly not recommended to run TandemTools on shorter TRs.


## Installation

Requirements are listed in ```requirements.txt``` and can be installed through Conda as ```conda install --file requirements.txt```

## Usage

```shell
tandemquast.py [options] --nano/--pacbio <reads_file> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    --nano      PATH                    File with Oxford Nanopore reads used for ETR assembly 
      or
    --pacbio    PATH                    File with PacBio CLR reads used for ETR assembly
    -o          PATH                    Folder to store all result files

Optional arguments:    
    -t          INT                     Maximum number of threads [default: 4]
    -l          \"label,label,...\"     Human-readable names of assemblies to use in reports, comma-separated. If contain spaces, use quotes 
    -m          PATH                    FASTA file with monomer sequences
    --hifi     PATH                    File with accurate PacBio HiFi reads
    --polish                            Run tandemQUAST polishing module (metrics will not be calculated).
```

```shell
tandemmapper.py [options] -r <reads_file> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    --nano      PATH                    File with Oxford Nanopore reads used for ETR assembly 
      or
    --pacbio    PATH                    File with PacBio CLR reads used for ETR assembly
    -o          PATH                    Folder to store alignment results

Optional arguments:    
    -t          INT                     Maximum number of threads [default: 4]
```

## Output files

The following files are contained in <output_dir> directory (specified by `-o`) and include results 
for all input assemblies. See TandemQUAST paper for the detailed explanation of different metrics.

### TandemMapper output

`<output_dir>/*_alignment.bed` - TandemMapper alignments in BED format
`<output_dir>/*_alignment.sam` - TandemMapper alignments in SAM format

### TandemQUAST output

####Basic metrics:

`<output_dir>/report/*_coverage.png` - plot of read coverage.

`<output_dir>/report/*_bp_analysis.png` - plot of breakpoint ratio values
 (see TandemQUAST paper for the details). The red peaks in this plot may correspond
 to large-scale assembly errors.

`<output_dir>/report/*_kmer_analysis.png` - distribution of different types of 
unique k-mers along the assembly. Each bar 
shows the number of different types of k-mers in a bin of length 20 kb. The blue bars
represent single-clump k-mers. The high number of single-clump k-mers suggests the good 
base-level quality of the assembly. The orange (multiple-clumps) and green (no-clumps) bars
suggest a low base-level quality in the region, caused by lack of polishing or an assembly error.

`<output_dir>/report/*_kmer_stats.txt` - distribution of different types of unique k-mers in TXT format.

####Pairwise comparison:

`<output_dir>/report/*_vs_*.png` - a dot plot comparing mappings for two assemblies

`<output_dir>/report/discordance_*_vs_*.png` - a discordance plot showing coverage of two assemblies by discordant reads.
The peaks of coverage for one assembly suggest that this assembly is more "supported" by reads in this region than other assembly. 

####Centromeric metrics:

`<output_dir>/report/*_monomer_lengths.html` - an interactive HTMl-page showing 
monomer length distribution along the assembly.

`<output_dir>/report/*_units.txt` - file with a list of HOR units.

## Citation

Alla Mikheenko, Andrey V. Bzikadze, Alexey Gurevich, Karen H. Miga, and Pavel A. Pevzner. TandemMapper and TandemQUAST: mapping long reads and assessing/improving assembly quality in extra-long tandem repeats, 2019, bioRxiv


## Contacts

Please report any problems to the [issue tracker](https://github.com/ablab/tandemQUAST/issues). Alternatively, you can write directly to [a.mikheenko@spbu.ru](mailto:a.mikheenko@spbu.ru).