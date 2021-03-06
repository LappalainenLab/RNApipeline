# RNApipeline
#### A pipeline to automate RNA sequence aligning and quality control workflows via list-based batch submission and parallel processing
---

## What is RNApipeline for?

RNApipeline is a series of scripts, or handlers, to automate and speed up DNA sequence alignment and quality control.
Currently, RNApipeline is designed to work with Illumina paired- and single-end bulk RNAseq data data.
The workflow is designed to process samples in batches and in parallel.
It is also intended to be easy for users to configure.
The handlers use a list of sequences, with full sequence paths as input.
RNApipeline uses [GNU Parallel](http://www.gnu.org/software/parallel/) to speed up analysis.
The handlers are designed to run as jobs submitted to a job scheduler.

## Configuration File

The included configuration file, `proj.conf`, provides information needed to run each of the handlers within it.
No other information is needed as RNApipeline pulls all necessary information from `proj.conf`.
Variables that are used by more than one handler are located at the top of `proj.conf`, followed by handler-specific variables, ending with software definitions.
Please read `proj.conf` <!--or visit [the wiki page](https://github.com/MorrellLAB/sequence_handling/wiki/Configuration)--> for more usage information.

`proj.conf` is broken up into several sections.
The first section, at the top of `Config` contains variables that are used by more than one handler.
Each section below is headed by a block of hash (`#`) marks and contains variables for one specific handler only.

For example, the section headed by
```shell
############################################
##########    Adapter_Trimming    ##########
############################################
```
contains variables for Adapter_Trimming only.
These variables are completely ignored by other handlers.

Please note, some of the variables are pre-defined in `proj.conf`.
These have been set for using the entirety of RNApipeline, and follows naming conventions used by all of the handlers.
If you choose to not use some of the handlers in your analyses (See *Do I have to use the entire workflow as is?* below), please modify variables as needed.

## Why use list-based batch submission?

Piping one sample alone through this workflow can take over 12 hours to run to completion.
Most sequence handling jobs are not dealing with one sample, so the amount of time to run this workflow increases drastically.
List-based batch submission simplifies the amount of typing that one has to do, and enables parallel processing to decrease time spent waiting for samples to finish.

An example list is shown below
>/home/path\_to\_sample/sample\_001\_R1.fastq.gz
>/home/path\_to\_sample/sample\_001\_R2.fastq.gz
>/home/path\_to\_sample/sample\_003_R1.fastq.gz
>/home/path\_to\_sample/sample\_003\_R2.fastq.gz

## Why use parallel processing?

Parallel processing decreases the amount of time by running multiple jobs at once and keeping track of which are done, which are running, and which have yet to be run.
This workflow, with the list-based batch submissions and parallel processing, both simplifies and quickens the process of sequence handling.

## Do I have to use the entire workflow as is?

No; no two handlers are entirely dependent on one another.
While all these handlers are designed to easily use the output from one to the next, these handlers are not required to achieve the end result of RNApipeline.
If you prefer tools other than the ones used within this workflow, you can modify or replace any or all of the handlers offered in RNApipeline.
This creates a pseudo-modularity for the entire workflow that allows for customization for each user.

## Dependencies

Due to the pseudo-modularity of this workflow, dependencies for each individual handler are listed below.
 The pipeline as a whole depends on BASH and a compute cluster that uses the [Grid Engine](http://www.univa.com/) as a scheduler

<!-- Please note that this is not a complete list of dependencies. Check the [dependencies wiki page](https://github.com/MorrellLab/sequence_handling/wiki/Dependencies) for  dependencies for each handler. -->

___

## Basic Usage

To run RNApipeline, use the following command, assuming you are in the RNApipeline directory:
```shell
./main.sh <handler> proj.conf
```
Where `<handler>` is one of the handlers listed below and `proj.conf` is the full file path to the configuration file.
A brief usage message can be viewed by passing no arguments to RNApipeline:
```shell
./main.sh
```

## Handlers

### Quality\_Assessment

The Quality_Assessment handler runs [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on a series of samples organized in a project directory for quality control.
The Quality_Assessment handler depends on:

 - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### Sequence\_Trimming

The Sequence_Trimming handler runs Trimmomatic on a series of samples, removing adapters and trimming based on quality.
This handler supports both paired-end and single-ended samples.
A list of all trimmed samples will be output at the end of all runs.
The Sequence_Trimming handler depends on:

 - [Java](https://www.java.com/en/)
 - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

### Read\_Mapping

The Read_Mapping handler maps reads to a reference genome using STAR.
This handler supports both paired-end and single-ended samples.
A list of all mapped samples will be output at the end of all runs.
The Read_Mapping handler depends on:

 - [STAR](https://github.com/alexdobin/STAR)

### SAM\_Processing

The SAM_Processing handler converts the SAM files from read mapping to the BAM format using [SAMTools](http://www.htslib.org/) and Picard.
In the conversion process, it will sort the reads and mark duplicates for the finished BAM file.
The SAM_Processing handler depends on:

- [SAMTools](http://www.htslib.org/)
- [Java](https://www.java.com/en/)
- [Picard](http://broadinstitute.github.io/picard/)

### Quantify\_Summarize

The Quanitfy_Summarize handler depends on:

 - [featureCounts](http://subread.sourceforge.net/)
 - [SAMTools](http://www.htslib.org/)
 - [Bedtools](https://bedtools.readthedocs.io/en/latest/)
 - [bc](https://www.gnu.org/software/bc/)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)
 - [GNU datamash](https://www.gnu.org/software/datamash/)
 - [Rscript](https://cran.r-project.org/)
 - [dstat](https://github.com/mojaveazure/ma_tools)

## Acknowledgements

RNApipeline was inspired by the [`sequence_handling`](https://github.com/MorrellLAB/sequence_handling) pipeline written by the [Morrell Lab](https://morrelllab.github.io/)
