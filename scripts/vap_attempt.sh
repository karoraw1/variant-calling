REF_=CAMPY28
SAMPLE_=CAMPY10
PFX=${REF_}-${SAMPLE_}
BASEDIR=/Volumes/KeithSSD/variant-calling
IN_FA=$BASEDIR/ref_fastas/${REF_}_bb.fa
VCF=$(ls $BASEDIR/snipout3_vcfs/*${PFX}*)
VEF_OUT=$(dirname $VCF)/$(basename $VCF .vcf).vef
IN_GFF=$BASEDIR/ref_fastas/${REF_}_bb.gff
IN_GFF2=$BASEDIR/ref_fastas/${REF_}_bb.srt.gff.gz

prodigal -f gff -i $IN_FA -o $IN_GFF
grep -v "#" $IN_GFF | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > $IN_GFF2

vep --offline --i $VCF -o $VEF_OUT --fasta $IN_FA --vcf  --sift b --polyphen b --custom ${IN_GFF2},CAMPY28,gff,exact,0

SNPENV=
_P1_=$SNPENV/perl5
_P2_=$SNPENV/lib/site_perl/5.26.2
_P3_=$SNPENV/bin/perl
export PERL5LIB="$_P1_:$_P2_:$_P3_"
