# SingleCellTranscriptEvidence

## Introduction
The analysis of single cell RNA-Seq data involves the use of a disparate set of software packages to identify the genes that are expressed within the sample. Although sophisticated software packages are now just a click away from being employed for these analyses, they expect the user to string the outputs from one step of the analyses to the next piece of software. We address this challenge by creating a pipeline that brings together a set of well-tested tools and combines them into a single script that is simple and easy-to-use. The pipeline described here allows users to analyze data from the Sequence Read Archive (SRA) by providing just the BioProject accession or specific run identifiers in combination with HISAT2 index files for the genome of interest. Upon execution, this pipeline identifies sequencing reads that correspond to novel genes that are yet to be annotated.

## Workflow
![Workflow](https://github.com/NCBI-Hackathons/SingleCellTranscriptEvidence/blob/master/workflow.jpg)


## Dependencies
The pipeline expects the following packages are installed and available in the PATH:
* Python (version 3.5 or above)
* HISAT2 (version 2.1.0)
* StringTie (version 1.3.3b)
* Samtools (version 1.3.1)

## Usage
```
usage: main.py [-h] [--sra_acc SRA_ACC] [--file FILE] [-r REFERENCE]
               [-p PROCESSES] [-o OUTDIR] [-b BAM]
               [-nso NOVEL_SPLICESITE_OUTFILE] [-sf STRINGTIE_FILE]
               [-a ABUNDANCE] [-m MULTI_MAP_FRAC]
               genome

positional arguments:
  genome                path of genome file

optional arguments:
  -h, --help            show this help message and exit
  --sra_acc SRA_ACC     SRR or PRJNA number (default: )
  --file FILE           path to a file with newline separated SRR numbers
                        (default: )
  -r REFERENCE, --reference REFERENCE
                        path to reference .gtf file (default: )
  -p PROCESSES, --processes PROCESSES
                        number of cores to use in run (default: 4)
  -o OUTDIR, --outdir OUTDIR
                        name of directory to save everything to (default: )
  -b BAM, --bam BAM     name of hisat2 output bam file (default:
                        hisat.sorted.bam)
  -nso NOVEL_SPLICESITE_OUTFILE, --novel_splicesite_outfile NOVEL_SPLICESITE_OUTFILE
                        Set stringtie novel_splicesite_outfile parameter
                        (default: splicesite.tab)
  -sf STRINGTIE_FILE, --stringtie_file STRINGTIE_FILE
                        name of stringtie output gtf file (default:
                        stringtie_file.gtf)
  -a ABUNDANCE, --abundance ABUNDANCE
                        Set stringtie -A parameter (default: abundance.tab)
  -m MULTI_MAP_FRAC, --multi_map_frac MULTI_MAP_FRAC
                        Set stringtie -M parameter (default: .95)
```
## Examples
main.py --sra_acc <SRR_accession> /path/to/HISAT2_index_files

## Output explanation

## Limitations/To-do

## People/Team
* Ashis Saha <ashis@jhu.edu>
* Michael Chambers <greenkidneybean@gmail.com>
* Allissa Dillman <allissa.dillman@gmail.com>
* Jessime Kirk <jessime@email.unc.edu>
* Sara Kimiko Suzuki <sksuzuki@ad.unc.edu>
* Wes Crouse <wcrouse@email.unc.edu>
* Vamsi Kodali <vkkodali@gmail.com>
