## Quick start
```shell
    tandemquast.py -r test_data/simulated_reads.fasta test_data/simulated_polished.fa -o simulated_res
```

## Introduction

**TandemTools** package includes **TandemQUAST** tool for evaluating and improving assemblies of extra-long tandem repeats (ETR) and **TandemMapper** tool for mapping long error-prone reads to ETRs.

Note: TandemTools is designed specifically for ETR (range in length from hundreds of thousands to millions of nucleotides). It is strongly not recommended to run TandemTools on shorter TRs.


## Installation

Requirements are listed in ```requirements.txt``` and can be installed through Conda as ```conda install --file requirements.txt```

## Usage

```shell
tandemquast.py [options] -r <reads_file> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    -r          PATH                    File with reads used for ETR assembly (only PacBio CLR or Oxford Nanopore reads)
    -o          PATH                    Folder to store all result files

Optional arguments:    
    -t          INT                     Maximum number of threads [default: 4]
    -l          \"label,label,...\"     Human-readable names of assemblies to use in reports, comma-separated. If contain spaces, use quotes 
    -m          PATH                    File with monomer sequences
    --hi-fi     PATH                    File with accurate PacBio HiFi reads
    --polish                            Run tandemQUAST polishing module (metrics will not be calculated).
```

```shell
tandemmapper.py [options] -r <reads_file> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    -r          PATH                    File with reads used for ETR assembly (only PacBio CLR or Oxford Nanopore reads)
    -o          PATH                    Folder to store alignment results

Optional arguments:    
    -t          INT                     Maximum number of threads [default: 4]
```
## Citation

Alla Mikheenko, Andrey V. Bzikadze, Alexey Gurevich, Karen H. Miga, and Pavel A. Pevzner. TandemMapper and TandemQUAST: mapping long reads and assessing/improving assembly quality in extra-long tandem repeats, 2019, bioRxiv


## Contacts

Please report any problems to the [issue tracker](https://github.com/ablab/tandemQUAST/issues). Alternatively, you can write directly to [a.mikheenko@spbu.ru](mailto:a.mikheenko@spbu.ru).