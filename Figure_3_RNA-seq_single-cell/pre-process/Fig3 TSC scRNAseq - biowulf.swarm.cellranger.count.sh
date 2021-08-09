
export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/username/2020-7-TE.exp.sc/INPUTS/biowulf.rnaseq.paths.sh \
&& cellranger count \
--id=H9TE_p24_r1 \
--transcriptome=/fdb/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=${FQDIR} \
--sample=H9TE_p24_r1 \
--localcores=$SLURM_CPUS_PER_TASK \
--expect-cells=5000


export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/username/2020-7-TE.exp.sc/INPUTS/biowulf.rnaseq.paths.sh \
&& cellranger count \
--id=H9TE_p24_r2 \
--transcriptome=/fdb/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=${FQDIR} \
--sample=H9TE_p24_r2 \
--localcores=$SLURM_CPUS_PER_TASK \
--expect-cells=5000

