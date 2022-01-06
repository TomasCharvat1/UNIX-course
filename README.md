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
Such prepared table can be used to solve more of the final exercises. For our purposes are some features of the table redundant.

### *Extra - Extracting relevant data and creating a new data frame just for chr 1*

 ```
zcat luscinia_vars.vcf.gz | grep -v '^#' | grep -v '_random' >nohead1.vcf
IN=nohead1.vcf
<$IN cut -f1-6 > main1.tsv
<$IN egrep -o 'DP=[^;]*' | sed 's/DP=//' > dp1.tsv
<$IN awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' >type1.tsv
paste main1.tsv dp1.tsv type1.tsv > all1.tsv
```


## Creating the figures in R

### Upload the data and load packages

```
library(tidyverse)
data<-read_tsv("all9.tsv",
               col_names=c("CHROM", "POS", "DOT", "REF", "ALT", "QUAL", "DP", "TYPE"))
```

### The distribution of DP over the genome and per chromosome

Following script will round the position of each measured DP by 10000 and subsequently calculate the mean for each set of 1000 bases. We then plot it 

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
  labs(x="Position on the chromosome", y= "Mean read depth per 10k bases", title = "Distribution of read depths of all chromosomes") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ```
 

* Z chromosome is marked as "NA" due to the numeric sorting.

![image](https://user-images.githubusercontent.com/95172475/148428131-71e2b1b7-c67c-498e-879f-438c2fff9615.png)


This figure however does not show the data in detail. Better idea gives a figure of just one chromosome. 

### The distribution of DP on chromosome 1

By simple adjustments of the previous code we get:
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
  labs(x="Position on the chromosome", y= "Mean read depth per 1k bases", title = "Distribution of read depths on chromosome 1") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![image](https://user-images.githubusercontent.com/95172475/148426711-941b7ed7-114f-45c3-861d-860aa1fae03b.png)



