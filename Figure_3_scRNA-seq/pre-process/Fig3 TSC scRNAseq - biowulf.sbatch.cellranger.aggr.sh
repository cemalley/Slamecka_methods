#!/bin/bash


# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


source ../INPUTS/biowulf.rnaseq.paths.sh

module load cellranger

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& cellranger aggr \
--id=H9TE_p24 \
--csv=../INPUTS/samplesheet.aggr.csv \
--normalize=mapped

