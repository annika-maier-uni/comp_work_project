# Simple RNA-Seq Pipeline using Nextflow

This pipeline is a RNA-Seq analysis workflow implemented using Nextflow. The workflow includes quality 
control, trimming, and alignment steps using popular tools like FastQC, TrimGalore, and HISAT2. 

## Workflow Overview
The pipeline performs the following steps:

1. Read Samplesheet: The input samplesheet (CSV) contains information about the samples, such as the sample name, 
strandedness, and the paths to paired-end FASTQ files. FASTQ Quality Control (FastQC): Performs quality control 
checks on raw FASTQ files using FastQC.

2. Trimming (TrimGalore): Trims adapters and low-quality regions from the FASTQ files.

3. Reference Genome Indexing (HISAT2): Builds the reference genome index required for alignment using HISAT2.

4. Alignment (HISAT2): Aligns the trimmed FASTQ files to the reference genome.


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please get familiar to these tools. Make sure to install nextflow version 23.10, which only runs in mac-os, linux, or WSL.

First, prepare a samplesheet with your input data that should have the same structure, file name, and column names as this example samplesheet:

**samplesheet.csv**:

```csv
sample,fastq_1,fastq_2,strandedness
SRR6357070,'./data/SRR6357070_1.fastq.gz','./data/SRR6357070_2.fastq.gz',auto
```

Each row represents a pair of fastq files (paired end). This tool does not allow single end data.
The strandedness refers to the library preparation and will be automatically inferred if set to `auto`.

## Executing the Pipeline

### Parameters
-- samplesheet: The path to the samplesheet CSV file containing sample information. The default value is ./modules/trimgalore/samplesheet_new.csv.
-- reference:   The path to the reference genome file (FASTA format). The default value is ./modules/hisat2/align/genome.fa. You can use the genome 
                data from here: https://github.com/nf-core/test-datasets/tree/rnaseq/reference

### Now, you can run the pipeline using:

```bash
nextflow run main.nf 
```

```bash
nextflow run main.nf --samplesheet <SAMPLESHEET>  
```

```bash
nextflow run main.nf --samplesheet <SAMPLESHEET> --reference <REFERENCE>
```

## Notes
This pipeline uses Nextflow's view() function to print intermediate channel contents, which is helpful for debugging.
Make sure the input files (samplesheet and reference genome) are correctly formatted and accessible before running the workflow.


## License
This example shows a simplified version. Once the file is pushed to GitHub, the README will be automatically displayed when someone 
visits the repository.