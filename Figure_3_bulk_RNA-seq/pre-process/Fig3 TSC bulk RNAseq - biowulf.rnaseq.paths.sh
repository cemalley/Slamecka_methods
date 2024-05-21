#!/bin/bash


# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley

# use "source biowulf.rnaseq.paths.sh" in the beginning of each analysis script to use these variables


#-------------------------------------------------------
# define sequencing run, software and reference versions
#-------------------------------------------------------


# sequencing run
SEQRUN="PIPE.RNA-seq"

# FASTQ file suffix
SUFR1="_R1_001.fastq.gz"
SUFR2="_R2_001.fastq.gz"

# software versions
BBTV="38.86" # BBTools, not used, just for reference
SUBRV="2.0.1" # Subread version
STARV="2.7.5a" # STAR version
SALV="0.14.2" # salmon version

# GTF file version
GTFV="100"


#-------------------------------------------------------
# define paths to directories
#-------------------------------------------------------


# define base directory for the analysis depending on the computer and drive used
BASEPATH="/data/username/${SEQRUN}/"

# define paths to directories
# FASTQ directory
FQDIR="${BASEPATH}FASTQ/"
# path the alignment BAM files
BAMPATH="${BASEPATH}BAM/"
# path to REFERENCE folder
REFPATH="/data/username/REFERENCE/"

# define directory containing programs
PROGPATH="/home/username/PROGRAMS/"
gffread="${PROGPATH}cufflinks-2.2.1.Linux_x86_64/gffread" # the version has not changed since 2014
salmon="${PROGPATH}salmon-${SALV}_linux_x86_64/bin/salmon"
fastqc="${PROGPATH}FastQC/fastqc"
bbmap="${PROGPATH}bbmap/bbmap.sh"
featureCounts="${PROGPATH}subread-${SUBRV}-Linux-x86_64/bin/featureCounts"
star="${PROGPATH}STAR-${STARV}/bin/Linux_x86_64/STAR"

# define REFERENCE file names and INDEX folders
GEN="${REFPATH}Homo_sapiens.GRCh38.dna.primary_assembly.fasta" # genome
GTF="${REFPATH}Homo_sapiens.GRCh38.${GTFV}.gtf" # GTF file
TR="${REFPATH}Homo_sapiens.GRCh38.${GTFV}.transcripts.fa" # transcriptome file
# path to salmon index
SALIDX="${REFPATH}salmon-index-${SALV}-GTF${GTFV}"
