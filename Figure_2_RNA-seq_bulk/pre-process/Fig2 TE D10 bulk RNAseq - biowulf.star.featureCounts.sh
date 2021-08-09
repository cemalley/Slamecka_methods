source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS10_S11_R1_001.fastq.gz ${FQDIR}JS10_S11_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS10_S11-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS10_S11-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS10_S11-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS10_S11-star-counts.txt ${FQDIR}JS10_S11-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS11_S12_R1_001.fastq.gz ${FQDIR}JS11_S12_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS11_S12-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS11_S12-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS11_S12-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS11_S12-star-counts.txt ${FQDIR}JS11_S12-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS12_S13_R1_001.fastq.gz ${FQDIR}JS12_S13_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS12_S13-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS12_S13-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS12_S13-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS12_S13-star-counts.txt ${FQDIR}JS12_S13-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS13_S14_R1_001.fastq.gz ${FQDIR}JS13_S14_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS13_S14-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS13_S14-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS13_S14-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS13_S14-star-counts.txt ${FQDIR}JS13_S14-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS14_S15_R1_001.fastq.gz ${FQDIR}JS14_S15_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS14_S15-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS14_S15-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS14_S15-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS14_S15-star-counts.txt ${FQDIR}JS14_S15-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS15_S16_R1_001.fastq.gz ${FQDIR}JS15_S16_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS15_S16-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS15_S16-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS15_S16-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS15_S16-star-counts.txt ${FQDIR}JS15_S16-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS16_S17_R1_001.fastq.gz ${FQDIR}JS16_S17_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS16_S17-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS16_S17-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS16_S17-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS16_S17-star-counts.txt ${FQDIR}JS16_S17-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS17_S18_R1_001.fastq.gz ${FQDIR}JS17_S18_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS17_S18-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS17_S18-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS17_S18-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS17_S18-star-counts.txt ${FQDIR}JS17_S18-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS18_S19_R1_001.fastq.gz ${FQDIR}JS18_S19_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS18_S19-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS18_S19-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS18_S19-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS18_S19-star-counts.txt ${FQDIR}JS18_S19-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS19_S20_R1_001.fastq.gz ${FQDIR}JS19_S20_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS19_S20-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS19_S20-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS19_S20-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS19_S20-star-counts.txt ${FQDIR}JS19_S20-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS1_S1_R1_001.fastq.gz ${FQDIR}JS1_S1_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS1_S1-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS1_S1-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS1_S1-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS1_S1-star-counts.txt ${FQDIR}JS1_S1-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS20_S21_R1_001.fastq.gz ${FQDIR}JS20_S21_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS20_S21-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS20_S21-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS20_S21-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS20_S21-star-counts.txt ${FQDIR}JS20_S21-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS21_S22_R1_001.fastq.gz ${FQDIR}JS21_S22_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS21_S22-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS21_S22-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS21_S22-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS21_S22-star-counts.txt ${FQDIR}JS21_S22-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS22_S23_R1_001.fastq.gz ${FQDIR}JS22_S23_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS22_S23-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS22_S23-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS22_S23-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS22_S23-star-counts.txt ${FQDIR}JS22_S23-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS23_S24_R1_001.fastq.gz ${FQDIR}JS23_S24_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS23_S24-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS23_S24-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS23_S24-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS23_S24-star-counts.txt ${FQDIR}JS23_S24-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS24_S25_R1_001.fastq.gz ${FQDIR}JS24_S25_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS24_S25-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS24_S25-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS24_S25-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS24_S25-star-counts.txt ${FQDIR}JS24_S25-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS25_S26_R1_001.fastq.gz ${FQDIR}JS25_S26_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS25_S26-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS25_S26-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS25_S26-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS25_S26-star-counts.txt ${FQDIR}JS25_S26-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS26_S27_R1_001.fastq.gz ${FQDIR}JS26_S27_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS26_S27-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS26_S27-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS26_S27-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS26_S27-star-counts.txt ${FQDIR}JS26_S27-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS27_S28_R1_001.fastq.gz ${FQDIR}JS27_S28_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS27_S28-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS27_S28-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS27_S28-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS27_S28-star-counts.txt ${FQDIR}JS27_S28-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS28_S29_R1_001.fastq.gz ${FQDIR}JS28_S29_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS28_S29-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS28_S29-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS28_S29-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS28_S29-star-counts.txt ${FQDIR}JS28_S29-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS29_S30_R1_001.fastq.gz ${FQDIR}JS29_S30_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS29_S30-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS29_S30-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS29_S30-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS29_S30-star-counts.txt ${FQDIR}JS29_S30-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS2_S3_R1_001.fastq.gz ${FQDIR}JS2_S3_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS2_S3-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS2_S3-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS2_S3-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS2_S3-star-counts.txt ${FQDIR}JS2_S3-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS30_S31_R1_001.fastq.gz ${FQDIR}JS30_S31_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS30_S31-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS30_S31-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS30_S31-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS30_S31-star-counts.txt ${FQDIR}JS30_S31-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS31_S32_R1_001.fastq.gz ${FQDIR}JS31_S32_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS31_S32-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS31_S32-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS31_S32-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS31_S32-star-counts.txt ${FQDIR}JS31_S32-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS32_S2_R1_001.fastq.gz ${FQDIR}JS32_S2_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS32_S2-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS32_S2-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS32_S2-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS32_S2-star-counts.txt ${FQDIR}JS32_S2-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS33_S33_R1_001.fastq.gz ${FQDIR}JS33_S33_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS33_S33-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS33_S33-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS33_S33-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS33_S33-star-counts.txt ${FQDIR}JS33_S33-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS34_S34_R1_001.fastq.gz ${FQDIR}JS34_S34_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS34_S34-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS34_S34-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS34_S34-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS34_S34-star-counts.txt ${FQDIR}JS34_S34-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS35_S35_R1_001.fastq.gz ${FQDIR}JS35_S35_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS35_S35-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS35_S35-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS35_S35-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS35_S35-star-counts.txt ${FQDIR}JS35_S35-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS36_S36_R1_001.fastq.gz ${FQDIR}JS36_S36_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS36_S36-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS36_S36-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS36_S36-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS36_S36-star-counts.txt ${FQDIR}JS36_S36-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS37_S37_R1_001.fastq.gz ${FQDIR}JS37_S37_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS37_S37-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS37_S37-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS37_S37-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS37_S37-star-counts.txt ${FQDIR}JS37_S37-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS38_S38_R1_001.fastq.gz ${FQDIR}JS38_S38_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS38_S38-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS38_S38-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS38_S38-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS38_S38-star-counts.txt ${FQDIR}JS38_S38-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS39_S39_R1_001.fastq.gz ${FQDIR}JS39_S39_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS39_S39-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS39_S39-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS39_S39-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS39_S39-star-counts.txt ${FQDIR}JS39_S39-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS3_S4_R1_001.fastq.gz ${FQDIR}JS3_S4_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS3_S4-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS3_S4-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS3_S4-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS3_S4-star-counts.txt ${FQDIR}JS3_S4-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS40_S40_R1_001.fastq.gz ${FQDIR}JS40_S40_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS40_S40-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS40_S40-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS40_S40-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS40_S40-star-counts.txt ${FQDIR}JS40_S40-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS41_S41_R1_001.fastq.gz ${FQDIR}JS41_S41_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS41_S41-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS41_S41-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS41_S41-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS41_S41-star-counts.txt ${FQDIR}JS41_S41-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS42_S42_R1_001.fastq.gz ${FQDIR}JS42_S42_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS42_S42-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS42_S42-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS42_S42-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS42_S42-star-counts.txt ${FQDIR}JS42_S42-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS43_S43_R1_001.fastq.gz ${FQDIR}JS43_S43_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS43_S43-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS43_S43-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS43_S43-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS43_S43-star-counts.txt ${FQDIR}JS43_S43-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS44_S44_R1_001.fastq.gz ${FQDIR}JS44_S44_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS44_S44-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS44_S44-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS44_S44-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS44_S44-star-counts.txt ${FQDIR}JS44_S44-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS45_S45_R1_001.fastq.gz ${FQDIR}JS45_S45_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS45_S45-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS45_S45-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS45_S45-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS45_S45-star-counts.txt ${FQDIR}JS45_S45-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS46_S46_R1_001.fastq.gz ${FQDIR}JS46_S46_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS46_S46-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS46_S46-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS46_S46-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS46_S46-star-counts.txt ${FQDIR}JS46_S46-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS47_S47_R1_001.fastq.gz ${FQDIR}JS47_S47_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS47_S47-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS47_S47-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS47_S47-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS47_S47-star-counts.txt ${FQDIR}JS47_S47-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS48_S48_R1_001.fastq.gz ${FQDIR}JS48_S48_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS48_S48-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS48_S48-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS48_S48-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS48_S48-star-counts.txt ${FQDIR}JS48_S48-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS49_S49_R1_001.fastq.gz ${FQDIR}JS49_S49_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS49_S49-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS49_S49-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS49_S49-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS49_S49-star-counts.txt ${FQDIR}JS49_S49-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS4_S5_R1_001.fastq.gz ${FQDIR}JS4_S5_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS4_S5-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS4_S5-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS4_S5-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS4_S5-star-counts.txt ${FQDIR}JS4_S5-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS50_S50_R1_001.fastq.gz ${FQDIR}JS50_S50_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS50_S50-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS50_S50-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS50_S50-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS50_S50-star-counts.txt ${FQDIR}JS50_S50-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS51_S51_R1_001.fastq.gz ${FQDIR}JS51_S51_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS51_S51-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS51_S51-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS51_S51-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS51_S51-star-counts.txt ${FQDIR}JS51_S51-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS52_S52_R1_001.fastq.gz ${FQDIR}JS52_S52_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS52_S52-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS52_S52-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS52_S52-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS52_S52-star-counts.txt ${FQDIR}JS52_S52-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS53_S53_R1_001.fastq.gz ${FQDIR}JS53_S53_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS53_S53-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS53_S53-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS53_S53-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS53_S53-star-counts.txt ${FQDIR}JS53_S53-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS5_S6_R1_001.fastq.gz ${FQDIR}JS5_S6_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS5_S6-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS5_S6-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS5_S6-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS5_S6-star-counts.txt ${FQDIR}JS5_S6-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS6_S7_R1_001.fastq.gz ${FQDIR}JS6_S7_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS6_S7-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS6_S7-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS6_S7-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS6_S7-star-counts.txt ${FQDIR}JS6_S7-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS7_S8_R1_001.fastq.gz ${FQDIR}JS7_S8_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS7_S8-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS7_S8-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS7_S8-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS7_S8-star-counts.txt ${FQDIR}JS7_S8-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS8_S9_R1_001.fastq.gz ${FQDIR}JS8_S9_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS8_S9-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS8_S9-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS8_S9-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS8_S9-star-counts.txt ${FQDIR}JS8_S9-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}JS9_S10_R1_001.fastq.gz ${FQDIR}JS9_S10_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}JS9_S10-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}JS9_S10-star-Aligned.sortedByCoord.out.bam ${FQDIR}JS9_S10-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}JS9_S10-star-counts.txt ${FQDIR}JS9_S10-star-c.bam

