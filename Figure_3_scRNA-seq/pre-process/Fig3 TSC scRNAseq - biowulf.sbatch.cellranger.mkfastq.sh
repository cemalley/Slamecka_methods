#!/bin/bash


# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


source ../INPUTS/biowulf.rnaseq.paths.sh

module load cellranger

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& cellranger mkfastq \
--id=IS026 \
--run=${BCLDIR} \
--csv=../INPUTS/samplesheet.mkfastq.csv \
--localcores=$SLURM_CPUS_PER_TASK


