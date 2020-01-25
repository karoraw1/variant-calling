#!/bin/bash

#SBATCH
#SBATCH --job-name=fsnp_%A-%a
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --partition=shared,lrgmem,parallel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=k.arorawilliams@gmail.com
#SBATCH --array=1-28
#SBATCH --output=fsnp_%A-%a.out
#SBATCH --error=fsnp_%A-%a.err

BASEDIR=$HOME/scratch/variant-calling
SAMPLES=$BASEDIR/misc_data/ref_map_pairs.txt
QCDIR=/scratch/users/karoraw1@jhu.edu/__CA__/CA_V3/data/01_QC_final

REF_MATCH=$(cut -f1 $SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p)
SNAME=$(cut -f2 $SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "$(date)" " :::: Starting analysis on ${REF_MATCH} ${SNAME}"

SNPENV=$BASEDIR/bin/SNIPPY_ENV
ml python/3.7-anaconda
source activate $SNPENV
_P1_=$SNPENV/perl5
_P2_=$SNPENV/lib/site_perl/5.26.2
_P3_=$SNPENV/bin/perl
export PERL5LIB="$_P1_:$_P2_:$_P3_"

OUTMAIN=$BASEDIR/snipout4; mkdir -p $OUTMAIN
PFX=${SLURM_ARRAY_TASK_ID}-${REF_MATCH}-${SNAME}
OUTDIR=${OUTMAIN}/${PFX};
REF_FA=${BASEDIR}/ref_fastas/${REF_MATCH}_bb.fa
TMPDIR=${OUTDIR}-tmp; rm -rf $TMPDIR; mkdir -p $TMPDIR
FWD_READS=${QCDIR}/${SNAME}-QUALITY_PASSED_R1.fastq
REV_READS=${QCDIR}/${SNAME}-QUALITY_PASSED_R2.fastq

snippy --prefix $PFX --tmpdir $TMPDIR --minqual 10 \
--minfrac 0.001 --cpus 24 --ram 80 --outdir $OUTDIR --ref $REF_FA --R1 $FWD_READS --R2 $REV_READS

OUTSPEC=$OUTDIR; SAMP_I=$PFX; cd $OUTSPEC;

echo "$(date)" " :::: Begin calling SNPS on $REF"
freebayes-parallel reference/ref.txt 24 -p 1 -P 0 -C 5 --min-repeat-entropy 1.5 --strict-vcf -q 13 -m 60 --min-coverage 5 -F 0.001  -f reference/ref.fa ${SAMP_I}.bam > ${SAMP_I}.raw2.vcf;
bcftools view --include 'QUAL>=10 && FMT/DP>=5 && (FMT/AO)/(FMT/DP)>=0.001' ${SAMP_I}.raw2.vcf  | vt normalize -r reference/ref.fa - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > ${SAMP_I}.filt2.vcf
cp ${SAMP_I}.filt2.vcf ${SAMP_I}.2.vcf
snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf ${SAMP_I}.2.vcf > ${SAMP_I}.2.tab

cd $BASEDIR/scripts

date
echo "finished against curated match, now against main ref"

if [ "$REF_MATCH" == "CG8421" ]; then
    REF_MATCH2="CG8421_sim"
    PFX2=${SLURM_ARRAY_TASK_ID}-${REF_MATCH2}-${SNAME}
    OUTDIR2=${OUTMAIN}/${PFX2};
    REF_FA2=${BASEDIR}/ref_fastas/${REF_MATCH2}_bb.fa
    TMPDIR2=${OUTDIR2}-tmp; rm -rf $TMPDIR2; mkdir -p $TMPDIR2
    
    snippy --prefix $PFX2 --tmpdir $TMPDIR2 --minqual 10 --minfrac 0.001 --cpus 24 --ram 80 --outdir $OUTDIR2 --ref $REF_FA2 --R1 $FWD_READS --R2 $REV_READS
    
    date
    echo "finished all"

	OUTSPEC=$OUTDIR2; SAMP_I=$PFX2; cd $OUTSPEC

	echo "$(date)" " :::: Begin calling SNPS on CG8421"
	freebayes-parallel reference/ref.txt 24 -p 1 -P 0 -C 5 --min-repeat-entropy 1.5 --strict-vcf -q 13 -m 60 --min-coverage 5 -F 0.001  -f reference/ref.fa ${SAMP_I}.bam > ${SAMP_I}.raw2.vcf;
	bcftools view --include 'QUAL>=10 && FMT/DP>=5 && (FMT/AO)/(FMT/DP)>=0.001' ${SAMP_I}.raw2.vcf  | vt normalize -r reference/ref.fa - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > ${SAMP_I}.filt2.vcf
	cp ${SAMP_I}.filt2.vcf ${SAMP_I}.2.vcf
	snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf ${SAMP_I}.2.vcf > ${SAMP_I}.2.tab
	
	echo "$(date)" " :::: SNP calls complete"
fi
