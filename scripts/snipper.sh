#!/bin/bash

#SBATCH
#SBATCH --job-name=snp2%A-%a
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --partition=shared,lrgmem,parallel
#SBATCH --mail-type=END
#SBATCH --mail-user=k.arorawilliams@gmail.edu
#SBATCH --array=1-50
#SBATCH --output=snp_%A-%a.out
#SBATCH --error=snp_%A-%a.err

BASEDIR=$HOME/scratch/variant-calling
SAMPLES=$BASEDIR/misc_data/ref_map_pairs.txt
QCDIR=/scratch/users/karoraw1@jhu.edu/__CA__/CA_V3/data/01_QC_final

REF_MATCH=$(cut -f1 $SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p)
SNAME=$(cut -f2 $SAMPLES | sed -n ${SLURM_ARRAY_TASK_ID}p)

ml python/3.7-anaconda-2019.03
source activate anvio5
rm ${QCDIR}/${SNAME}-QUALITY_PASSED_R*.fast*
rm ${QCDIR}/${SNAME}-STATS.txt
iu-filter-quality-minoche -p 0.3 --ignore-deflines ${QCDIR}/${SNAME}.ini

source deactivate 

FWD_READS=${QCDIR}/${SNAME}-QUALITY_PASSED_R1.fastq
if [ -f ${FWD_READS}.bz2 ]; then bzip2 -d ${FWD_READS}.bz2; fi
REV_READS=${QCDIR}/${SNAME}-QUALITY_PASSED_R2.fastq
if [ -f ${REV_READS}.bz2 ]; then bzip2 -d ${REV_READS}.bz2; fi

SNPENV=$BASEDIR/bin/SNIPPY_ENV
ml python/3.7-anaconda
source activate $SNPENV

_P1_=$SNPENV/perl5
_P2_=$SNPENV/lib/site_perl/5.26.2
_P3_=$SNPENV/bin/perl
export PERL5LIB="$_P1_:$_P2_:$_P3_"

OUTMAIN=$BASEDIR/snipout2; mkdir -p $OUTMAIN
OUTDIR=$OUTMAIN/${SLURM_ARRAY_TASK_ID}-${REF_MATCH}-${SNAME}
REF_FA=$BASEDIR/ref_fastas/${REF_MATCH}_bb.fa
TMPDIR=$[OUTDIR]-tmp; mkdir -p $TMPDIR

PFX=${SLURM_ARRAY_TASK_ID}-${SNAME}

snippy --prefix $PFX --tmpdir $TMPDIR --minqual 30 \
--minfrac 0.01 --cpus 24 --ram 80 --outdir $OUTDIR --ref $REF_FA --R1 $FWD_READS --R2 $REV_READS

rm -rf $TMPDIR

bzip2 ${FWD_READS} &
bzip2 ${REV_READS}

