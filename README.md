# `sequence_handling`
#### A pipeline to automate RNA sequence aligning and quality control workflows via list-based batch submission and parallel processing
---

## What is `sequence_handling` for?

`sequence_handling` is a series of scripts, or handlers, to automate and speed up DNA sequence alignment and quality control. Currently, `sequence_handling` is designed to work with Illumina paired-end whole-genome and exome capture data. Work is underway to expand `sequence_handling` to accept single-end and GBS data.
The workflow is designed to process samples in batches and in parallel. It is also intended to be easy for users to configure. The handlers use a list of sequences, with full sequence paths as input. `sequence_handling` uses [GNU Parallel](http://www.gnu.org/software/parallel/) to speed up  analysis.The handlers are designed to run as jobs submitted to a job scheduler, specifically the [Portable Batch System](http://www.pbsworks.com/).
> **NOTE:** This workflow is designed to use the [Portable Batch System](http://www.pbsworks.com/) and run on the [Minnesota Supercomputing Institute](https://www.msi.umn.edu). Heavy modifications will need to be made if not using these systems.
## Configuration File
The included configuration file, `proj.conf`, provides information needed to run each of the handlers within it. No other information is needed as `sequence_handling` pulls all necessary information from `proj.conf`. Variables that are used by more than one handler are located at the top of `proj.conf`, followed by handler-specific variables, ending with software definitions. Please read `proj.conf` <!--or visit [the wiki page](https://github.com/MorrellLAB/sequence_handling/wiki/Configuration)--> for more usage information.
`proj.conf` is broken up into several sections. The first section, at the top of `Config` contains variables that are used by more than one handler. Each section below is headed by a block of hash (`#`) marks and contains variables for one specific handler only. For example, the section headed by
```shell
############################################
##########    Adapter_Trimming    ##########
############################################
```
contains variables for Adapter_Trimming only. These variables are completely ignored by other handlers.
Please note, some of the variables are pre-defined in `proj.conf`. These have been set for using the entirety of `sequence_handling`, and follows naming conventions used by all of the handlers. If you choose to not use some of the handlers in your analyses (See *Do I have to use the entire workflow as is?* below), please modify variables as needed.

## Why use list-based batch submission?

Piping one sample alone through this workflow can take over 12 hours to run to completion. Most sequence handling jobs are not dealing with one sample, so the amount of time to run this workflow increases drastically. List-based batch submission simplifies the amount of typing that one has to do, and enables parallel processing to decrease time spent waiting for samples to finish. An example list is shown below.
>/home/path\_to\_sample/sample\_001\_R1.fastq.gz
>/home/path\_to\_sample/sample\_001\_R2.fastq.gz
>/home/path\_to\_sample/sample\_003_R1.fastq.gz
>/home/path\_to\_sample/sample\_003\_R2.fastq.gz

## Why use parallel processing?

Parallel processing decreases the amount of time by running multiple jobs at once and keeping track of which are done, which are running, and which have yet to be run. This workflow, with the list-based batch submissions and parallel processing, both simplifies and quickens the process of sequence handling.

## Do I have to use the entire workflow as is?

No. No two handlers are entirely dependent on one another. While all these handlers are designed to easily use the output from one to the next, these handlers are not required to achieve the end result of `sequence_handling`. If you prefer tools other than the ones used within this workflow, you can modify or replace any or all of the handlers offered in `sequence_handling`. This creates a pseudo-modularity for the entire workflow that allows for customization for each user.

<!--
## Dependencies

Due to the pseudo-modularity of this workflow, dependencies for each individual handler are listed below. Some general dependencies for the workflow as a whole are also listed here:
 - [GNU Parallel](http://www.gnu.org/software/parallel/)
 - A quality control mechanism, such as [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - A quality trimmer, such as the [Seqqs](https://github.com/morrelllab.seqqs)/[Sickle](https://github.com/najoshi/sickle) combo
 - An adapter trimmer, such as [Scythe](https://github.com/vsbuffalo/scythe)
 - A read mapper, such as [The Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/) (BWA)
 - Tools for plotting results, such as [R](http://cran.r-project.org/)
 - SAM file processing utilities, such as [SAMTools](http://www.htslib.org/) and/or [Picard](http://broadinstitute.github.io/picard/)
Please note that this is not a complete list of dependencies. Check the [dependencies wiki page](https://github.com/MorrellLab/sequence_handling/wiki/Dependencies) for  dependencies for each handler.
-->

___

## Basic Usage

To run `sequence_handling`, use the following command, assuming you are in the `sequence_handling` directory:
```shell
./main.sh <handler> proj.conf
```
Where `<handler>` is one of the handlers listed below and `proj.conf` is the full file path to the configuration file.
A brief usage message can be viewed by passing no arguments to `sequence_handling`:
```shell
./main.sh
```
