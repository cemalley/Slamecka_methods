source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791780.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791780-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791780-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791780-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791780-star-counts.txt ${FQDIR}SRR11791780-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791781.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791781-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791781-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791781-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791781-star-counts.txt ${FQDIR}SRR11791781-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791782.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791782-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791782-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791782-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791782-star-counts.txt ${FQDIR}SRR11791782-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791783.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791783-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791783-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791783-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791783-star-counts.txt ${FQDIR}SRR11791783-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791784.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791784-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791784-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791784-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791784-star-counts.txt ${FQDIR}SRR11791784-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791785.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791785-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791785-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791785-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791785-star-counts.txt ${FQDIR}SRR11791785-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791786.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791786-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791786-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791786-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791786-star-counts.txt ${FQDIR}SRR11791786-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791787.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791787-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791787-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791787-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791787-star-counts.txt ${FQDIR}SRR11791787-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791788.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791788-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791788-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791788-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791788-star-counts.txt ${FQDIR}SRR11791788-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791789.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791789-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791789-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791789-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791789-star-counts.txt ${FQDIR}SRR11791789-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791790.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791790-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791790-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791790-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791790-star-counts.txt ${FQDIR}SRR11791790-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791791.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791791-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791791-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791791-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791791-star-counts.txt ${FQDIR}SRR11791791-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791792.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791792-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791792-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791792-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791792-star-counts.txt ${FQDIR}SRR11791792-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791793.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791793-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791793-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791793-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791793-star-counts.txt ${FQDIR}SRR11791793-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791794.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791794-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791794-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791794-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791794-star-counts.txt ${FQDIR}SRR11791794-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791795.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791795-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791795-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791795-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791795-star-counts.txt ${FQDIR}SRR11791795-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791796.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791796-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791796-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791796-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791796-star-counts.txt ${FQDIR}SRR11791796-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791797.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791797-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791797-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791797-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791797-star-counts.txt ${FQDIR}SRR11791797-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791798.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791798-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791798-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791798-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791798-star-counts.txt ${FQDIR}SRR11791798-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791799.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791799-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791799-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791799-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791799-star-counts.txt ${FQDIR}SRR11791799-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791800.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791800-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791800-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791800-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791800-star-counts.txt ${FQDIR}SRR11791800-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791801.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791801-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791801-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791801-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791801-star-counts.txt ${FQDIR}SRR11791801-star-c.bam

source /data/slameckaj2/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}SRR11791802.sra.fastq --outFileNamePrefix ${FQDIR}SRR11791802-star- \
--runThreadN $SLURM_CPUS_PER_TASK \
--outSAMtype BAM SortedByCoordinate \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
&& mv ${FQDIR}SRR11791802-star-Aligned.sortedByCoord.out.bam ${FQDIR}SRR11791802-star-c.bam \
&& $featureCounts -s 1 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}SRR11791802-star-counts.txt ${FQDIR}SRR11791802-star-c.bam

