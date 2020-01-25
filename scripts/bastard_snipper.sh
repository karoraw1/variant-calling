#!/bin/bash

#SBATCH
#SBATCH --job-name=assem_ref
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --partition=shared,lrgmem,parallel
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1
#SBATCH --output=bast_%A-%a.out
#SBATCH --error=bast_%A-%a.err

BASEDIR=$HOME/scratch/variant-calling
SNPENV=$BASEDIR/bin/SNIPPY_ENV
ml python/3.7-anaconda
source activate $SNPENV
_P1_=$SNPENV/perl5
_P2_=$SNPENV/lib/site_perl/5.26.2
_P3_=$SNPENV/bin/perl
export PERL5LIB="$_P1_:$_P2_:$_P3_"

OUTDIR=${BASEDIR}/snipout3
SAMPLES=${BASEDIR}/misc_data/ref_map_pairs.txt

REF=$(cut -f1 $SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p)
QUERY=$(cut -f2 $SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p)

OUTSPEC=$(dirname $OUTDIR/*${REF}-${QUERY}/*.bam)
SAMP_I=$(basename $OUTSPEC)
cd $OUTSPEC

echo "$(date)" " :::: Begin calling SNPS on $REF"

freebayes-parallel reference/ref.txt 24 -p 1 -P 0 -C 5 --min-repeat-entropy 1.5 --strict-vcf -q 13 -m 60 --min-coverage 5 -F 0.001  -f reference/ref.fa ${SAMP_I}.bam > ${SAMP_I}.raw2.vcf;
bcftools view --include 'QUAL>=10 && FMT/DP>=5 && (FMT/AO)/(FMT/DP)>=0.001' ${SAMP_I}.raw2.vcf  | vt normalize -r reference/ref.fa - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > ${SAMP_I}.filt2.vcf
cp ${SAMP_I}.filt2.vcf ${SAMP_I}.2.vcf
snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf ${SAMP_I}.2.vcf > ${SAMP_I}.2.tab

OUTSPEC=$(dirname $OUTDIR/*CG8421-${QUERY}/*.bam)
SAMP_I=$(basename $OUTSPEC)
cd $OUTSPEC

echo "$(date)" " :::: Begin calling SNPS on CG8421"

freebayes-parallel reference/ref.txt 24 -p 1 -P 0 -C 5 --min-repeat-entropy 1.5 --strict-vcf -q 13 -m 60 --min-coverage 5 -F 0.001  -f reference/ref.fa ${SAMP_I}.bam > ${SAMP_I}.raw2.vcf;
bcftools view --include 'QUAL>=10 && FMT/DP>=5 && (FMT/AO)/(FMT/DP)>=0.001' ${SAMP_I}.raw2.vcf  | vt normalize -r reference/ref.fa - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > ${SAMP_I}.filt2.vcf
cp ${SAMP_I}.filt2.vcf ${SAMP_I}.2.vcf
snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf ${SAMP_I}.2.vcf > ${SAMP_I}.2.tab

echo "$(date)" " :::: SNP calls complete"

