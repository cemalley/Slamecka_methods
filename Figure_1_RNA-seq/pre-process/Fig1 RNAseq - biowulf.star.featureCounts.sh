source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}16_D3_A1_1_S14_R1_001.fastq.gz ${FQDIR}16_D3_A1_1_S14_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}16_D3_A1_1_S14-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}16_D3_A1_1_S14-star-Aligned.sortedByCoord.out.bam ${FQDIR}16_D3_A1_1_S14-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}16_D3_A1_1_S14-star-counts.txt ${FQDIR}16_D3_A1_1_S14-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}17_D3_A1_2_S15_R1_001.fastq.gz ${FQDIR}17_D3_A1_2_S15_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}17_D3_A1_2_S15-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}17_D3_A1_2_S15-star-Aligned.sortedByCoord.out.bam ${FQDIR}17_D3_A1_2_S15-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}17_D3_A1_2_S15-star-counts.txt ${FQDIR}17_D3_A1_2_S15-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}18_D3_A1_3_S16_R1_001.fastq.gz ${FQDIR}18_D3_A1_3_S16_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}18_D3_A1_3_S16-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}18_D3_A1_3_S16-star-Aligned.sortedByCoord.out.bam ${FQDIR}18_D3_A1_3_S16-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}18_D3_A1_3_S16-star-counts.txt ${FQDIR}18_D3_A1_3_S16-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}1_D0_SB_1_S11_R1_001.fastq.gz ${FQDIR}1_D0_SB_1_S11_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}1_D0_SB_1_S11-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}1_D0_SB_1_S11-star-Aligned.sortedByCoord.out.bam ${FQDIR}1_D0_SB_1_S11-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}1_D0_SB_1_S11-star-counts.txt ${FQDIR}1_D0_SB_1_S11-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}28_D7_A_1_S17_R1_001.fastq.gz ${FQDIR}28_D7_A_1_S17_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}28_D7_A_1_S17-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}28_D7_A_1_S17-star-Aligned.sortedByCoord.out.bam ${FQDIR}28_D7_A_1_S17-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}28_D7_A_1_S17-star-counts.txt ${FQDIR}28_D7_A_1_S17-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}29_D7_A_2_S18_R1_001.fastq.gz ${FQDIR}29_D7_A_2_S18_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}29_D7_A_2_S18-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}29_D7_A_2_S18-star-Aligned.sortedByCoord.out.bam ${FQDIR}29_D7_A_2_S18-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}29_D7_A_2_S18-star-counts.txt ${FQDIR}29_D7_A_2_S18-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}2_D0_SB_2_S12_R1_001.fastq.gz ${FQDIR}2_D0_SB_2_S12_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}2_D0_SB_2_S12-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}2_D0_SB_2_S12-star-Aligned.sortedByCoord.out.bam ${FQDIR}2_D0_SB_2_S12-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}2_D0_SB_2_S12-star-counts.txt ${FQDIR}2_D0_SB_2_S12-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}30_D7_A_3_S19_R1_001.fastq.gz ${FQDIR}30_D7_A_3_S19_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}30_D7_A_3_S19-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}30_D7_A_3_S19-star-Aligned.sortedByCoord.out.bam ${FQDIR}30_D7_A_3_S19-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}30_D7_A_3_S19-star-counts.txt ${FQDIR}30_D7_A_3_S19-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}34_D10_A_1_S20_R1_001.fastq.gz ${FQDIR}34_D10_A_1_S20_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}34_D10_A_1_S20-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}34_D10_A_1_S20-star-Aligned.sortedByCoord.out.bam ${FQDIR}34_D10_A_1_S20-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}34_D10_A_1_S20-star-counts.txt ${FQDIR}34_D10_A_1_S20-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}35_D10_A_2_S21_R1_001.fastq.gz ${FQDIR}35_D10_A_2_S21_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}35_D10_A_2_S21-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}35_D10_A_2_S21-star-Aligned.sortedByCoord.out.bam ${FQDIR}35_D10_A_2_S21-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}35_D10_A_2_S21-star-counts.txt ${FQDIR}35_D10_A_2_S21-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}36_D10_A_3_S22_R1_001.fastq.gz ${FQDIR}36_D10_A_3_S22_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}36_D10_A_3_S22-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}36_D10_A_3_S22-star-Aligned.sortedByCoord.out.bam ${FQDIR}36_D10_A_3_S22-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}36_D10_A_3_S22-star-counts.txt ${FQDIR}36_D10_A_3_S22-star-c.bam

source /data/username/PIPE.RNA-seq/INPUTS/biowulf.rnaseq.paths.sh \
&& $star --genomeDir ${REFPATH}star-index --sjdbGTFfile $GTF \
--readFilesIn ${FQDIR}3_D0_SB_3_S13_R1_001.fastq.gz ${FQDIR}3_D0_SB_3_S13_R2_001.fastq.gz --outFileNamePrefix ${FQDIR}3_D0_SB_3_S13-star- \
--runThreadN $SLURM_CPUS_PER_TASK --sjdbOverhang 149 \
--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
&& mv ${FQDIR}3_D0_SB_3_S13-star-Aligned.sortedByCoord.out.bam ${FQDIR}3_D0_SB_3_S13-star-c.bam \
&& $featureCounts -p -s 2 -t exon -g gene_id -T $SLURM_CPUS_PER_TASK --ignoreDup -a $GTF -G $GEN -o ${FQDIR}3_D0_SB_3_S13-star-counts.txt ${FQDIR}3_D0_SB_3_S13-star-c.bam

