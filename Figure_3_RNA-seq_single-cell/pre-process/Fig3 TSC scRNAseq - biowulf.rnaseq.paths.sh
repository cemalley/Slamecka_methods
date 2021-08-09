#!/bin/bash


# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley

# this script serves to set all paths to files needed for the entire analysis pipeline
# use "source rnaseq.paths.sh" in the beginning of each analysis script to use these variables


# sequencing run
SEQRUN="2020-7-TE.exp.sc"

# define base directory for the analysis depending on the computer and drive used
BASEPATH=/data/username/${SEQRUN}/

# FASTQ directory
FQDIR="${BASEPATH}FASTQ/"

# BCL directory
BCLDIR="${BASEPATH}BCL/200703_A00898_0008_AHGGVLDRXX/"