export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050169.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050169-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050169-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050169-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050169-star-counts.txt ${FQDIR}SRR11050169-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050170.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050170-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050170-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050170-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050170-star-counts.txt ${FQDIR}SRR11050170-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050171.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050171-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050171-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050171-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050171-star-counts.txt ${FQDIR}SRR11050171-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050172.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050172-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050172-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050172-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050172-star-counts.txt ${FQDIR}SRR11050172-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050173.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050173-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050173-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050173-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050173-star-counts.txt ${FQDIR}SRR11050173-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050174.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050174-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050174-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050174-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050174-star-counts.txt ${FQDIR}SRR11050174-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050175.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050175-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050175-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050175-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050175-star-counts.txt ${FQDIR}SRR11050175-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050176.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050176-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050176-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050176-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050176-star-counts.txt ${FQDIR}SRR11050176-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050177.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050177-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050177-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050177-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050177-star-counts.txt ${FQDIR}SRR11050177-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050178.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050178-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050178-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050178-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050178-star-counts.txt ${FQDIR}SRR11050178-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050179.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050179-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050179-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050179-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050179-star-counts.txt ${FQDIR}SRR11050179-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050180.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050180-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050180-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050180-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050180-star-counts.txt ${FQDIR}SRR11050180-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050181.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050181-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050181-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050181-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050181-star-counts.txt ${FQDIR}SRR11050181-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050182.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050182-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050182-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050182-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050182-star-counts.txt ${FQDIR}SRR11050182-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050183.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050183-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050183-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050183-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050183-star-counts.txt ${FQDIR}SRR11050183-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050184.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050184-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050184-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050184-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050184-star-counts.txt ${FQDIR}SRR11050184-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050185.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050185-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050185-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050185-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050185-star-counts.txt ${FQDIR}SRR11050185-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050186.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050186-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050186-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050186-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050186-star-counts.txt ${FQDIR}SRR11050186-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050187.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050187-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050187-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050187-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050187-star-counts.txt ${FQDIR}SRR11050187-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050188.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050188-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050188-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050188-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050188-star-counts.txt ${FQDIR}SRR11050188-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050189.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050189-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050189-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050189-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050189-star-counts.txt ${FQDIR}SRR11050189-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050190.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050190-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050190-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050190-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050190-star-counts.txt ${FQDIR}SRR11050190-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050191.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050191-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050191-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050191-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050191-star-counts.txt ${FQDIR}SRR11050191-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050192.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050192-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050192-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050192-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050192-star-counts.txt ${FQDIR}SRR11050192-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050193.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050193-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050193-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050193-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050193-star-counts.txt ${FQDIR}SRR11050193-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050194.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050194-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050194-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050194-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050194-star-counts.txt ${FQDIR}SRR11050194-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050195.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050195-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050195-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050195-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050195-star-counts.txt ${FQDIR}SRR11050195-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050196.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050196-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050196-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050196-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050196-star-counts.txt ${FQDIR}SRR11050196-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050197.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050197-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050197-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050197-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050197-star-counts.txt ${FQDIR}SRR11050197-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050198.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050198-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050198-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050198-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050198-star-counts.txt ${FQDIR}SRR11050198-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050199.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050199-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050199-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050199-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050199-star-counts.txt ${FQDIR}SRR11050199-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11050200.fastq.gz --outFileNamePrefix ${FQDIR}SRR11050200-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR11050200-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11050200-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11050200-star-counts.txt ${FQDIR}SRR11050200-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR12789776.fastq.gz --outFileNamePrefix ${FQDIR}SRR12789776-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR12789776-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR12789776-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR12789776-star-counts.txt ${FQDIR}SRR12789776-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR12789777.fastq.gz --outFileNamePrefix ${FQDIR}SRR12789777-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR12789777-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR12789777-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR12789777-star-counts.txt ${FQDIR}SRR12789777-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR12789778.fastq.gz --outFileNamePrefix ${FQDIR}SRR12789778-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR12789778-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR12789778-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR12789778-star-counts.txt ${FQDIR}SRR12789778-star-c.bam

export TMPDIR=/lscratch/$SLURM_JOB_ID \
&& source /data/slameckaj2/Io2021/PRJNA605646/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR12789779.fastq.gz --outFileNamePrefix ${FQDIR}SRR12789779-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outTmpDir /lscratch/$SLURM_JOB_ID/temp \
&& mv ${FQDIR}SRR12789779-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR12789779-star-c.bam \
&& $featureCounts -p -s 0 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR12789779-star-counts.txt ${FQDIR}SRR12789779-star-c.bam

