@# Read CG8421 Genbank file to get gene feature positions and annotations and DNA sequence
@# Read in Crofts variation locations 
@# Validate genbank parse by confirming expected reference alleles & nearby gene features

@ # aggregate positions by feature so as not to make multiple anchors per effect

@# Create 3 nucleotide anchor sequences per feature:
@## 250 bp on either side of mutation 
@## Full gene sequences of 2 closest genes on either side (excluding affected gene if in sequence)

@# Create 2-3 protein sequence anchors 
@# left side, right side, & in sequence (if applicable)

@# create a new folder for all data written to disk
@# Write out both anchor files into files called reference_anchors

@# Read in functional enrichment dataframe and extract relapse pairs

@# Create columns of assembly & protein sequence locations

@# Read in all of both types of data

@# Create unique parseable header names 

@# Add CG8421 reference genome & protein sequence to the dataframe  

@# write out all into a file (test) and data for one sample into a file (train)

@# print 'makeblastdb' commands & 'blastp/n' commands to copy into terminal

@# read in the blast hits produced 

@# determine how many unique hits were recorded for each anchor position 

@# cross check protein hits and nucleotide hits 
@# determine if the effected region was hit, and if synteny was conserved in other anchors 
@# group the contigs which contain homolog that exists in > 2 relapse pairs

# write them into a fasta file, taking care to ensure the copy number of each homolog matches CG8421
# ( select the more complete contig if a choice needs to be made )

# convert fasta -> gff w/ prodigal/prokka (if snippy can't handle interpreting mutations w/ just fasta)

# run snippy to compare (relapse v primary) & (primary v primary) 

# read in all variants observed 

# check to see if both reported variants were observed & if coding changes were observed by prodigal:
# 1. innoc v 4R 82R | 1262397^1262398 / 1262399 / 1263109^1263110 | CJ8421_RS06630 | Asn297 or Glu534
# 2. innoc v 31R 4R | 335945^335946 / 335972 | CJ8421_RS01805 | Asn157 or Cys166
# 3. innoc 82R v (81R, 1R, 4R, 9R , 31R) | 629686^629687/ 629687 | CJ8421_RS03350 | Frameshift 61% into gene
# 4. innoc v (31R, 8.2R, 8.1R)| 1361226 or 1361225^1361226 | CJ8421_RS07085 & CJ8421_RS07090 

# check to see if other reported variants were observed programatically 
