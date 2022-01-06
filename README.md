# The final excercise - Distribution of read depth (DP) over the whole genome and by chromosome
## Overview
This repository represents my solution of the given task as a part of the final exam from the course Unix and work with genomic data.
## Handling the data in UNIX
### First look at the data
Using the following functions we can see that we have 71 variables in the CHROM column, some of which are marked random and unmapped.
``` 
zcat luscinia_vars.vcf.gz | less -S

zcat luscinia_vars.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq -c | sort -k1,1

zcat luscinia_vars.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq -c | sort -k1,1 | wc -l
```

### Clean and prepare the data for sorting
Because the dataset includes some random chromosomes, unmapped positions, and weird linkage groups, all of small sizes, I decided to delete them for safety reasons and for the sake of practicing with UNIX features. If all given chromosomes should be included, this part can be left out. Also, for the purposes of sorting, Z chromosome is cut separately and will be pasted back after the other chromosomes are sorted.

```
IN= luscinia_vars.vcf.gz

<$IN zcat | grep -v '^#' | grep -v '_random' | sed 's/chr//' | grep -v '^[A-Z]' | cut -f2-8 > 2-8.tsv
<$IN zcat | grep -v '^#' | grep -v '_random' | sed 's/chr//' | grep -v '^[A-Z]' | sed 's/A/.1/' | sed 's/B/.2/' | cut -f1 > 1.tsv
<$IN zcat | grep -v '^#' | grep -e 'chrZ\s' | sed 's/chr//' | sort -k2,2n > Zchrom.tsv
```

Now we have to sort the obtained files and paste them together.

```
paste 1.tsv 2-8.tsv >fix.tsv  
cat fix.tsv | sort -k1,1n -k2,2n  >sorted.tsv
cat sorted.tsv Zchrom.tsv | sed 's/^/chr/' >clear.tsv
 ```
 
 ### Extracting relevant data and creating a new data frame
 
 
 ```
IN=clear.tsv

<$IN cut -f1-6 > clearmain.tsv
<$IN egrep -o 'DP=[^;]*' | sed 's/DP=//' > cleardp.tsv
<$IN awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' >cleartype.tsv
paste clearmain.tsv cleardp.tsv cleartype.tsv >all9.tsv
```
Such prepared table can be used for majority (if not all) of the final exercises. 

## Creating the figures in R

### Upload the data and load packages

```
library(tidyverse)
data<-read_tsv("all9.tsv",
               col_names=c("CHROM", "POS", "DOT", "REF", "ALT", "QUAL", "DP", "TYPE"))
```

### The distribution of DP over the genome and per chromosome

Following script will rond the position of each measured DP by 10000 and subsequently calculate the mean for each set of 1000 bases. We then plot it 

```
data %>%
  group_by(CHROM) %>%
  mutate(POS_block = plyr::round_any(POS, 1e3)) ->
  dc

dc %>%
  group_by(CHROM, POS_block) %>%
  summarise(DP = mean(DP)) -> dcc

dcc %>%
  ggplot(aes(POS_block, DP)) +
  geom_line() +
  facet_wrap(~as.numeric(CHROM), ncol = 8) +
  labs(x="Position on the chromosome", y= "Mean read depth per 10k bases", title = "Read depth of all chromosomes") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

```

![image](https://user-images.githubusercontent.com/95172475/148413142-e5062240-c787-4418-8867-f21c8392db8d.png)




