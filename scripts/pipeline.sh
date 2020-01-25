

BASEDIR=$HOME/scratch/variant-calling
cd $BASEDIR
ml python/3.7-anaconda
source activate ./SNIPPY_ENV
export PERL5LIB=`which perl`
snippy --help
snippy --mincov 5 --targets sites.bed 
snippy --cpus 16 --outdir mysnps --ref Listeria.gbk --R1 FDA_R1.fastq.gz --R2 FDA_R2.fastq.gz


REF_FA='Campylobacter_jejuni_subsp__jejuni_CG8421_GCA_000171795_2-contigs.fa'
REF_GB='Campylobacter_jejuni_subsp__jejuni_CG8421_GCA_000171795_2-contigs.fa'

bwa index reference.fasta
samtools faidx reference.fasta
java -jar ~/bin/picard-tools-1.8.5/CreateSequenceDictionary.jar REFERENCE=reference.fasta OUTPUT=reference.dict
bwa mem -R “@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:test\tSM:someID” reference.fasta R1.fastq.gz R2.fastq.gz > aln.sam
java -jar ~/bin/picard-tools-1.8.5/SortSam.jar I=aln.sam O=sorted.bam SORT_ORDER=coordinate
java -jar ~/bin/picard-tools-version/MarkDuplicates.jar I=sorted.bam O=dedup.bam METRICS_FILE=metrics.txt
java -jar ~/bin/picard-tools-version/BuildBamIndex.jar INPUT=dedup.bam
java -jar ~/bin/GATK3.4/GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I dedup.bam -o targetintervals.list
java -jar ~/bin/GATK3.4/GenomeAnalysisTK.jar -T IndelRealigner -R reference.fasta -I dedup.bam -targetIntervals targetintervals.list -o realigned.bam
java -jar ~/bin/GATK3.4/GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fasta -I realigned.bam -ploidy 1 –emitRefConfidence GVCF -o raw_gVCF.vcf
java -jar ~/bin/GATK3.4/GenomeAnalysisTK.jar -T GenotypeGVCFs -R reference.fasta –variant raw_gVCF.vcf -o raw.vcf



bwa mem -M -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' ref input_1 input_2 > aligned_reads.sam
java -jar picard.jar SortSam INPUT=aligned_reads.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate
java -jar picard.jar CollectAlignmentSummaryMetrics R=ref I=sorted_reads.bam O=alignment_metrics.txt
java -jar picard.jar CollectInsertSizeMetrics INPUT=sorted_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf
samtools depth -a sorted_reads.bam > depth_out.txt
java -jar picard.jar MarkDuplicates INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt
java -jar picard.jar BuildBamIndex INPUT=dedup_reads.bam
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref -I dedup_reads.bam -o realignment_targets.list
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ref -I dedup_reads.bam -targetIntervals realignment_targets.list -o realigned_reads.bam
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref -I realigned_reads.bam -o raw_variants.vcf
java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants.vcf -selectType SNP -o raw_snps.vcf
java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants.vcf -selectType INDEL -o raw_indels.vcf
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps.vcf
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels.vcf
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ref -I realigned_reads.bam -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -o recal_data.table
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ref -I realigned_reads.bam -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -BQSR recal_data.table -o post_recal_data.table
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ref -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf
java -jar GenomeAnalysisTK.jar -T PrintReads -R ref -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref -I recal_reads.bam -o raw_variants_recal.vcf
java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants_recal.vcf -selectType SNP -o raw_snps_recal.vcf
java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants_recal.vcf -selectType INDEL -o raw_indels_recal.vcf
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps_final.vcf
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels_recal.vcf
java -jar snpEff.jar -v snpeff_db filtered_snps_final.vcf > filtered_snps_final.ann.vcf
bedtools genomecov -bga -ibam recal_reads.bam > genomecov.bedgraph