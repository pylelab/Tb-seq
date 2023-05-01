# Tb-seq data processing

This processing pipeline is adapted from Sexton AN, Wang PY, Rutenberg-Schoenberg M, Simon MD. Interpreting Reverse Transcriptase Termination and Mutation Events for Greater Insight into the Chemical Probing of RNA. Biochemistry. 2017 Sep 5;56(35):4713-4721. doi: 10.1021/acs.biochem.7b00323. Epub 2017 Aug 18. PMID: 28820243; PMCID: PMC5648349. 

This workflow is designed to calculate reverse transcription events from aligned next generation sequencing data. 

Requirement:

Python 2.7
 with modules: sys, os, re, traceback, time

RStudio 1.2.5001
with packages: ggplot2, tidyr, dplyr, lazyeval

Files needed:
1. Config file (.csv)
	Configuration file; file located in working directory/configs
2. Reference file (.fa, with indexes built by Hisat2 or Bowtie2 as required)
	A FASTA file containing reference sequence of target; file located in working directory/genomes
3. Primer list (if applicable)
	Primer file containing sequences of oligo nucleotides that pair with RNA; file located in working_directory/primers.dir
4. RTEventsCounter.py
	Script used to calculate RT termination events; file located in working directory/scripts
5. conf.py (for RTEventsCounter.py)
	Configuration file for RTEventsCounter.py; file located in working directory/scripts
6. Terbium.sh
	Full pipeline to process fastq files; located in working directory/scripts
7. Terbium_submission.pbs
	Configuration file for Terbium.sh; located in working directory/scripts

Output of the RTEventsCounter.py is a .CSV file containing a table of counts. This file is further processed using the following scripts:

8. 1910_datawrangling_Terbium.R
This script will sort and filter the data into a relevant data frame; output will be a .l.csv file
9. Terbium_data_analysis.R
This script will produce final normalized reactivity values and allow for visualization of the data; output will be a .csv file

