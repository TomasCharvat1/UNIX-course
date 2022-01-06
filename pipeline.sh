#!/bin/bash

# First look at the data

zcat luscinia_vars.vcf.gz | less -S
zcat luscinia_vars.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq -c | sort -k1,1
zcat luscinia_vars.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq -c | sort -k1,1 | wc -l

# Clean and prepare the data for sorting
IN= luscinia_vars.vcf.gz

<$IN zcat | grep -v '^#' | grep -v '_random' | sed 's/chr//' | grep -v '^[A-Z]' | cut -f2-8 > 2-8.tsv
<$IN zcat | grep -v '^#' | grep -v '_random' | sed 's/chr//' | grep -v '^[A-Z]' | sed 's/A/.1/' | sed 's/B/.2/' | cut -f1 > 1.tsv
<$IN zcat | grep -v '^#' | grep -e 'chrZ\s' | sed 's/chr//' | sort -k2,2n > Zchrom.tsv

paste 1.tsv 2-8.tsv >fix.tsv  
cat fix.tsv | sort -k1,1n -k2,2n  >sorted.tsv
cat sorted.tsv Zchrom.tsv | sed 's/^/chr/' >clear.tsv

# Extract relevant data and creating a new data frame
IN=clear.tsv

<$IN cut -f1-6 > clearmain.tsv
<$IN egrep -o 'DP=[^;]*' | sed 's/DP=//' > cleardp.tsv
<$IN awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' >cleartype.tsv
paste clearmain.tsv cleardp.tsv cleartype.tsv >all9.tsv

# Extra - Extracting relevant data and creating a new data frame just for chr 1
zcat luscinia_vars.vcf.gz | grep -v '^#' | grep -v '_random' >nohead1.vcf
IN=nohead1.vcf
<$IN cut -f1-6 > main1.tsv
<$IN egrep -o 'DP=[^;]*' | sed 's/DP=//' > dp1.tsv
<$IN awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' >type1.tsv
paste main1.tsv dp1.tsv type1.tsv > all1.tsv
