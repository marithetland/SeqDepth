# SeqDepth

Calculate average coverage (sequencing depth) and standard deviation for Illumina short-read paired-end sequences.

## Table of Contents

[Requirements](#Requirements)  
[Usage](#Usage)  
[Example command](#example-command)  
[Output](#Output)  

## Requirements

* Linux or MacOS
* BWA v0.7.17-r1188
* SAMtools v1.7
* Picard tools

## Usage

For each sample, paired-end reads (preferably trimmed) and their assembled fasta-file must be provided. Please ensure all three files have the same prefix, e.g: "sample_1.fastq.gz" "sample_2.fastq.gz" and "sample.fasta".

``` 
SeqDepth: Calculate the sequencing depth (overall coverage) and standard deviation of trimmed FASTQ-files. 

Usage: seqDepth.py [-h] [-v] -r READS [READS ...] -a ASSEMBLIES
                   [ASSEMBLIES ...] [-k] [-o OUTPUT]

Input (required):
  -r READS [READS ...], --reads READS [READS ...]
                        Provide the full file path to the (trimmed) reads
                        (ending in e.g. *fastq.gz, fq.gz)
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        Provide the full file path to the assembly-files
                        (ending in *fasta)

Optional arguments:
  -k, --keep            Use this flag if you want to keep the BAM files used
                        in the calculation.
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Output options:
  -o OUTPUT, --output OUTPUT
                        Optional outfile name. Default="seqdepth_results.csv"

```

For each successful calculation, a "sample.Successful" file is created. If the program is interrupted or you want to add files to your outfile, please leave these files as they are when you (re)start the program (otherwise you can just delete them). Also remember to use the same outfile with the -o option (if you specified).


## Example command

``` python
seqdepth.py -a ./assemblies/*fasta -r ./reads/*fastq.gz -k 
```

## Output

* The output is stored in "seqdepth_results.txt" if no outfile has been specified.
* If you use the --keep option, all BAM-files will be kept.
