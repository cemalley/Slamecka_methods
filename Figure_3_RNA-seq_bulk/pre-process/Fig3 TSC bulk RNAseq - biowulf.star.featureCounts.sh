source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_hPSC_1_R1_001.fastq.gz ${FQDIR}H7_hPSC_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_hPSC_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_hPSC_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_hPSC_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_hPSC_1-star-counts.txt ${FQDIR}H7_hPSC_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_hPSC_2_R1_001.fastq.gz ${FQDIR}H7_hPSC_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_hPSC_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_hPSC_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_hPSC_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_hPSC_2-star-counts.txt ${FQDIR}H7_hPSC_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_hPSC_3_R1_001.fastq.gz ${FQDIR}H7_hPSC_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_hPSC_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_hPSC_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_hPSC_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_hPSC_3-star-counts.txt ${FQDIR}H7_hPSC_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_hPSC_4_R1_001.fastq.gz ${FQDIR}H7_hPSC_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_hPSC_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_hPSC_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_hPSC_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_hPSC_4-star-counts.txt ${FQDIR}H7_hPSC_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_hPSC_5_R1_001.fastq.gz ${FQDIR}H7_hPSC_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_hPSC_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_hPSC_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_hPSC_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_hPSC_5-star-counts.txt ${FQDIR}H7_hPSC_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_TE_exp_early_1_R1_001.fastq.gz ${FQDIR}H7_TE_exp_early_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_TE_exp_early_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_TE_exp_early_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_TE_exp_early_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_TE_exp_early_1-star-counts.txt ${FQDIR}H7_TE_exp_early_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_TE_exp_early_2_R1_001.fastq.gz ${FQDIR}H7_TE_exp_early_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_TE_exp_early_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_TE_exp_early_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_TE_exp_early_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_TE_exp_early_2-star-counts.txt ${FQDIR}H7_TE_exp_early_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_TE_exp_early_3_R1_001.fastq.gz ${FQDIR}H7_TE_exp_early_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_TE_exp_early_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_TE_exp_early_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_TE_exp_early_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_TE_exp_early_3-star-counts.txt ${FQDIR}H7_TE_exp_early_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_TE_p0_1_R1_001.fastq.gz ${FQDIR}H7_TE_p0_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_TE_p0_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_TE_p0_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_TE_p0_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_TE_p0_1-star-counts.txt ${FQDIR}H7_TE_p0_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_TE_p0_2_R1_001.fastq.gz ${FQDIR}H7_TE_p0_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_TE_p0_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_TE_p0_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_TE_p0_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_TE_p0_2-star-counts.txt ${FQDIR}H7_TE_p0_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H7_TE_p0_3_R1_001.fastq.gz ${FQDIR}H7_TE_p0_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H7_TE_p0_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H7_TE_p0_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H7_TE_p0_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H7_TE_p0_3-star-counts.txt ${FQDIR}H7_TE_p0_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_hPSC_1_R1_001.fastq.gz ${FQDIR}H9_hPSC_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_hPSC_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_hPSC_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_hPSC_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_hPSC_1-star-counts.txt ${FQDIR}H9_hPSC_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_hPSC_2_R1_001.fastq.gz ${FQDIR}H9_hPSC_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_hPSC_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_hPSC_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_hPSC_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_hPSC_2-star-counts.txt ${FQDIR}H9_hPSC_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_hPSC_3_R1_001.fastq.gz ${FQDIR}H9_hPSC_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_hPSC_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_hPSC_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_hPSC_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_hPSC_3-star-counts.txt ${FQDIR}H9_hPSC_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_hPSC_4_R1_001.fastq.gz ${FQDIR}H9_hPSC_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_hPSC_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_hPSC_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_hPSC_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_hPSC_4-star-counts.txt ${FQDIR}H9_hPSC_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_hPSC_5_R1_001.fastq.gz ${FQDIR}H9_hPSC_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_hPSC_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_hPSC_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_hPSC_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_hPSC_5-star-counts.txt ${FQDIR}H9_hPSC_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_early_1_R1_001.fastq.gz ${FQDIR}H9_TE_exp_early_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_early_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_early_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_early_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_early_1-star-counts.txt ${FQDIR}H9_TE_exp_early_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_early_2_R1_001.fastq.gz ${FQDIR}H9_TE_exp_early_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_early_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_early_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_early_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_early_2-star-counts.txt ${FQDIR}H9_TE_exp_early_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_early_3_R1_001.fastq.gz ${FQDIR}H9_TE_exp_early_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_early_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_early_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_early_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_early_3-star-counts.txt ${FQDIR}H9_TE_exp_early_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_hypoxia_1_R1_001.fastq.gz ${FQDIR}H9_TE_exp_hypoxia_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_hypoxia_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_hypoxia_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_hypoxia_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_hypoxia_1-star-counts.txt ${FQDIR}H9_TE_exp_hypoxia_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_hypoxia_2_R1_001.fastq.gz ${FQDIR}H9_TE_exp_hypoxia_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_hypoxia_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_hypoxia_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_hypoxia_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_hypoxia_2-star-counts.txt ${FQDIR}H9_TE_exp_hypoxia_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_hypoxia_3_R1_001.fastq.gz ${FQDIR}H9_TE_exp_hypoxia_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_hypoxia_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_hypoxia_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_hypoxia_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_hypoxia_3-star-counts.txt ${FQDIR}H9_TE_exp_hypoxia_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_hypoxia_4_R1_001.fastq.gz ${FQDIR}H9_TE_exp_hypoxia_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_hypoxia_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_hypoxia_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_hypoxia_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_hypoxia_4-star-counts.txt ${FQDIR}H9_TE_exp_hypoxia_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_hypoxia_5_R1_001.fastq.gz ${FQDIR}H9_TE_exp_hypoxia_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_hypoxia_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_hypoxia_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_hypoxia_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_hypoxia_5-star-counts.txt ${FQDIR}H9_TE_exp_hypoxia_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_hypoxia_6_R1_001.fastq.gz ${FQDIR}H9_TE_exp_hypoxia_6_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_hypoxia_6-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_hypoxia_6-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_hypoxia_6-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_hypoxia_6-star-counts.txt ${FQDIR}H9_TE_exp_hypoxia_6-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_late_1_R1_001.fastq.gz ${FQDIR}H9_TE_exp_late_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_late_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_late_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_late_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_late_1-star-counts.txt ${FQDIR}H9_TE_exp_late_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_late_2_R1_001.fastq.gz ${FQDIR}H9_TE_exp_late_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_late_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_late_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_late_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_late_2-star-counts.txt ${FQDIR}H9_TE_exp_late_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_exp_late_3_R1_001.fastq.gz ${FQDIR}H9_TE_exp_late_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_exp_late_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_exp_late_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_exp_late_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_exp_late_3-star-counts.txt ${FQDIR}H9_TE_exp_late_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_p0_1_R1_001.fastq.gz ${FQDIR}H9_TE_p0_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_p0_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_p0_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_p0_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_p0_1-star-counts.txt ${FQDIR}H9_TE_p0_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_p0_2_R1_001.fastq.gz ${FQDIR}H9_TE_p0_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_p0_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_p0_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_p0_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_p0_2-star-counts.txt ${FQDIR}H9_TE_p0_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_p0_3_R1_001.fastq.gz ${FQDIR}H9_TE_p0_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_p0_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_p0_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_p0_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_p0_3-star-counts.txt ${FQDIR}H9_TE_p0_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}H9_TE_p0_4_R1_001.fastq.gz ${FQDIR}H9_TE_p0_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}H9_TE_p0_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}H9_TE_p0_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}H9_TE_p0_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}H9_TE_p0_4-star-counts.txt ${FQDIR}H9_TE_p0_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_hPSC_1_R1_001.fastq.gz ${FQDIR}JHU191i_hPSC_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_hPSC_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_hPSC_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_hPSC_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_hPSC_1-star-counts.txt ${FQDIR}JHU191i_hPSC_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_hPSC_2_R1_001.fastq.gz ${FQDIR}JHU191i_hPSC_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_hPSC_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_hPSC_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_hPSC_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_hPSC_2-star-counts.txt ${FQDIR}JHU191i_hPSC_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_hPSC_3_R1_001.fastq.gz ${FQDIR}JHU191i_hPSC_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_hPSC_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_hPSC_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_hPSC_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_hPSC_3-star-counts.txt ${FQDIR}JHU191i_hPSC_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_early_1_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_early_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_early_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_early_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_early_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_early_1-star-counts.txt ${FQDIR}JHU191i_TE_exp_early_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_early_2_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_early_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_early_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_early_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_early_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_early_2-star-counts.txt ${FQDIR}JHU191i_TE_exp_early_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_early_3_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_early_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_early_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_early_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_early_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_early_3-star-counts.txt ${FQDIR}JHU191i_TE_exp_early_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_early_4_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_early_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_early_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_early_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_early_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_early_4-star-counts.txt ${FQDIR}JHU191i_TE_exp_early_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_late_1_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_late_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_late_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_late_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_late_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_late_1-star-counts.txt ${FQDIR}JHU191i_TE_exp_late_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_late_2_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_late_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_late_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_late_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_late_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_late_2-star-counts.txt ${FQDIR}JHU191i_TE_exp_late_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_late_3_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_late_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_late_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_late_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_late_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_late_3-star-counts.txt ${FQDIR}JHU191i_TE_exp_late_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_late_4_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_late_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_late_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_late_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_late_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_late_4-star-counts.txt ${FQDIR}JHU191i_TE_exp_late_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_late_5_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_late_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_late_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_late_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_late_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_late_5-star-counts.txt ${FQDIR}JHU191i_TE_exp_late_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_exp_late_6_R1_001.fastq.gz ${FQDIR}JHU191i_TE_exp_late_6_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_exp_late_6-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_exp_late_6-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_exp_late_6-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_exp_late_6-star-counts.txt ${FQDIR}JHU191i_TE_exp_late_6-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_p0_1_R1_001.fastq.gz ${FQDIR}JHU191i_TE_p0_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_p0_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_p0_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_p0_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_p0_1-star-counts.txt ${FQDIR}JHU191i_TE_p0_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_p0_2_R1_001.fastq.gz ${FQDIR}JHU191i_TE_p0_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_p0_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_p0_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_p0_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_p0_2-star-counts.txt ${FQDIR}JHU191i_TE_p0_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU191i_TE_p0_3_R1_001.fastq.gz ${FQDIR}JHU191i_TE_p0_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU191i_TE_p0_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU191i_TE_p0_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU191i_TE_p0_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU191i_TE_p0_3-star-counts.txt ${FQDIR}JHU191i_TE_p0_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_hPSC_1_R1_001.fastq.gz ${FQDIR}JHU198i_hPSC_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_hPSC_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_hPSC_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_hPSC_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_hPSC_1-star-counts.txt ${FQDIR}JHU198i_hPSC_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_hPSC_2_R1_001.fastq.gz ${FQDIR}JHU198i_hPSC_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_hPSC_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_hPSC_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_hPSC_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_hPSC_2-star-counts.txt ${FQDIR}JHU198i_hPSC_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_hPSC_3_R1_001.fastq.gz ${FQDIR}JHU198i_hPSC_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_hPSC_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_hPSC_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_hPSC_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_hPSC_3-star-counts.txt ${FQDIR}JHU198i_hPSC_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_hPSC_4_R1_001.fastq.gz ${FQDIR}JHU198i_hPSC_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_hPSC_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_hPSC_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_hPSC_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_hPSC_4-star-counts.txt ${FQDIR}JHU198i_hPSC_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_early_1_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_early_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_early_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_early_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_early_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_early_1-star-counts.txt ${FQDIR}JHU198i_TE_exp_early_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_early_2_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_early_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_early_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_early_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_early_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_early_2-star-counts.txt ${FQDIR}JHU198i_TE_exp_early_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_early_3_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_early_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_early_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_early_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_early_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_early_3-star-counts.txt ${FQDIR}JHU198i_TE_exp_early_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_early_4_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_early_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_early_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_early_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_early_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_early_4-star-counts.txt ${FQDIR}JHU198i_TE_exp_early_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_hypoxia_1_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_hypoxia_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_hypoxia_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_hypoxia_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_hypoxia_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_hypoxia_1-star-counts.txt ${FQDIR}JHU198i_TE_exp_hypoxia_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_hypoxia_2_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_hypoxia_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_hypoxia_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_hypoxia_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_hypoxia_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_hypoxia_2-star-counts.txt ${FQDIR}JHU198i_TE_exp_hypoxia_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_hypoxia_3_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_hypoxia_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_hypoxia_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_hypoxia_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_hypoxia_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_hypoxia_3-star-counts.txt ${FQDIR}JHU198i_TE_exp_hypoxia_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_late_1_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_late_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_late_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_late_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_late_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_late_1-star-counts.txt ${FQDIR}JHU198i_TE_exp_late_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_late_2_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_late_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_late_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_late_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_late_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_late_2-star-counts.txt ${FQDIR}JHU198i_TE_exp_late_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_exp_late_3_R1_001.fastq.gz ${FQDIR}JHU198i_TE_exp_late_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_exp_late_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_exp_late_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_exp_late_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_exp_late_3-star-counts.txt ${FQDIR}JHU198i_TE_exp_late_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_p0_1_R1_001.fastq.gz ${FQDIR}JHU198i_TE_p0_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_p0_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_p0_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_p0_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_p0_1-star-counts.txt ${FQDIR}JHU198i_TE_p0_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_p0_2_R1_001.fastq.gz ${FQDIR}JHU198i_TE_p0_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_p0_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_p0_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_p0_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_p0_2-star-counts.txt ${FQDIR}JHU198i_TE_p0_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_p0_3_R1_001.fastq.gz ${FQDIR}JHU198i_TE_p0_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_p0_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_p0_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_p0_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_p0_3-star-counts.txt ${FQDIR}JHU198i_TE_p0_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_p0_4_R1_001.fastq.gz ${FQDIR}JHU198i_TE_p0_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_p0_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_p0_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_p0_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_p0_4-star-counts.txt ${FQDIR}JHU198i_TE_p0_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_p0_5_R1_001.fastq.gz ${FQDIR}JHU198i_TE_p0_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_p0_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_p0_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_p0_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_p0_5-star-counts.txt ${FQDIR}JHU198i_TE_p0_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JHU198i_TE_p0_6_R1_001.fastq.gz ${FQDIR}JHU198i_TE_p0_6_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JHU198i_TE_p0_6-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JHU198i_TE_p0_6-star-Aligned.sortedByCoord.out.bam ${FQDIR}JHU198i_TE_p0_6-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JHU198i_TE_p0_6-star-counts.txt ${FQDIR}JHU198i_TE_p0_6-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_hPSC_1_R1_001.fastq.gz ${FQDIR}MCW032i_hPSC_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_hPSC_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_hPSC_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_hPSC_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_hPSC_1-star-counts.txt ${FQDIR}MCW032i_hPSC_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_hPSC_2_R1_001.fastq.gz ${FQDIR}MCW032i_hPSC_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_hPSC_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_hPSC_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_hPSC_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_hPSC_2-star-counts.txt ${FQDIR}MCW032i_hPSC_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_hPSC_3_R1_001.fastq.gz ${FQDIR}MCW032i_hPSC_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_hPSC_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_hPSC_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_hPSC_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_hPSC_3-star-counts.txt ${FQDIR}MCW032i_hPSC_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_hPSC_4_R1_001.fastq.gz ${FQDIR}MCW032i_hPSC_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_hPSC_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_hPSC_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_hPSC_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_hPSC_4-star-counts.txt ${FQDIR}MCW032i_hPSC_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_hPSC_5_R1_001.fastq.gz ${FQDIR}MCW032i_hPSC_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_hPSC_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_hPSC_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_hPSC_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_hPSC_5-star-counts.txt ${FQDIR}MCW032i_hPSC_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_exp_early_1_R1_001.fastq.gz ${FQDIR}MCW032i_TE_exp_early_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_exp_early_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_exp_early_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_exp_early_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_exp_early_1-star-counts.txt ${FQDIR}MCW032i_TE_exp_early_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_exp_early_2_R1_001.fastq.gz ${FQDIR}MCW032i_TE_exp_early_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_exp_early_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_exp_early_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_exp_early_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_exp_early_2-star-counts.txt ${FQDIR}MCW032i_TE_exp_early_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_exp_early_3_R1_001.fastq.gz ${FQDIR}MCW032i_TE_exp_early_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_exp_early_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_exp_early_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_exp_early_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_exp_early_3-star-counts.txt ${FQDIR}MCW032i_TE_exp_early_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_exp_early_4_R1_001.fastq.gz ${FQDIR}MCW032i_TE_exp_early_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_exp_early_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_exp_early_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_exp_early_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_exp_early_4-star-counts.txt ${FQDIR}MCW032i_TE_exp_early_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_p0_1_R1_001.fastq.gz ${FQDIR}MCW032i_TE_p0_1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_p0_1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_p0_1-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_p0_1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_p0_1-star-counts.txt ${FQDIR}MCW032i_TE_p0_1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_p0_2_R1_001.fastq.gz ${FQDIR}MCW032i_TE_p0_2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_p0_2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_p0_2-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_p0_2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_p0_2-star-counts.txt ${FQDIR}MCW032i_TE_p0_2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_p0_3_R1_001.fastq.gz ${FQDIR}MCW032i_TE_p0_3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_p0_3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_p0_3-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_p0_3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_p0_3-star-counts.txt ${FQDIR}MCW032i_TE_p0_3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_p0_4_R1_001.fastq.gz ${FQDIR}MCW032i_TE_p0_4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_p0_4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_p0_4-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_p0_4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_p0_4-star-counts.txt ${FQDIR}MCW032i_TE_p0_4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_p0_5_R1_001.fastq.gz ${FQDIR}MCW032i_TE_p0_5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_p0_5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_p0_5-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_p0_5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_p0_5-star-counts.txt ${FQDIR}MCW032i_TE_p0_5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}MCW032i_TE_p0_6_R1_001.fastq.gz ${FQDIR}MCW032i_TE_p0_6_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}MCW032i_TE_p0_6-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}MCW032i_TE_p0_6-star-Aligned.sortedByCoord.out.bam ${FQDIR}MCW032i_TE_p0_6-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}MCW032i_TE_p0_6-star-counts.txt ${FQDIR}MCW032i_TE_p0_6-star-c.bam

