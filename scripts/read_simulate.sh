#!/bin/bash

#SBATCH
#SBATCH --job-name=assem_ref
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --partition=shared,lrgmem,parallel
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --output=assem_ref.out
#SBATCH --error=assem_ref.err

# git clone https://github.com/lh3/wgsim.git
# cd wgsim
# gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
# chmod +x wgsim
# mv wgsim ../wgsim_exe
# cd ..; rm -rf wgsim
# mv wgsim_exe wgsim

#BASEDIR=/Volumes/KeithSSD/variant-calling
#EXE=$BASEDIR/bin/wgsim
#REF=$BASEDIR/refs/Campylobacter_jejuni_subsp__jejuni_CG8421_GCA_000171795_2-contigs.fa
#OUTD=$BASEDIR/simulated_reference
#mkdir -p $OUTD

#$EXE -S42 -e0.003 -r0 -d500 -s50 -N5000000 -1250 -2250 -R0 -X0 -h $REF $OUTD/CG8421_1.fq $OUTD/CG8421_2.fq

SPADES_=/home-3/karoraw1@jhu.edu/scratch/__CA__/bin/SPAdes-3.13.0-Linux/bin/spades.py
OUTD=/home-3/karoraw1@jhu.edu/scratch/variant-calling/simulated_reference
#TMEP=$OUTD/tmep
#OUT1=$OUTD/CG8421_assem

#mkdir $TMEP
#$SPADES_ -o $OUT1 -t 24 -1 $OUTD/CG8421_1.fq -2 $OUTD/CG8421_2.fq -m 74 --tmp-dir $TMEP --cov-cutoff 'auto'
#mv $OUT1/contigs.fasta $OUTD/CG8421.fa
#rm -rf $TMEP

ml python/3.7-anaconda-2019.03
source activate anvio5
anvi-script-reformat-fasta $OUTD/CG8421.fa -o $OUTD/CG8421_contigs_fixed.fa --simplify-names

# make contigs database
anvi-gen-contigs-database -f $OUTD/CG8421_contigs_fixed.fa -o $OUTD/CG8421_contigs.db -n "CG8421"
anvi-export-table -o $OUTD/CG8421_Proteins.txt --table gene_amino_acid_sequences $OUTD/CG8421_contigs.db
anvi-export-table -o CG8421_GenePositions.txt --table genes_in_contigs CG8421_contigs.db



