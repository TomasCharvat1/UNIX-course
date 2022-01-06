#R
#Upload the data and load packages
library(tidyverse)
data<-read_tsv("all9.tsv",
               col_names=c("CHROM", "POS", "DOT", "REF", "ALT", "QUAL", "DP", "TYPE"))

#Read depth of the genome and per chromosome
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
  theme(axis.text.x = element_text(angle = 90))


# Read depth of the whole genome

data$CHROMZ<- gsub("Z","29", data$CHROM) 

data %>%
 ggplot(aes(as.numeric(CHROM),DP)) +
 geom_smooth() +
 labs(x="Position in the genome - From chr1 to 29 (Z)", y= "Read depth", title = "Distribution of read depths over whole genome") + 
 theme(axis.text.x = element_text(angle = 90))

#Distribution of the DP per genome

data %>%
  ggplot(aes(DP)) +
  geom_histogram() +
  scale_x_log10() + 
  labs(x="Read depth", y= "N of bases", title = "Distribution of read depth over whole genome")

#Distribution of the DP per chromosome
data %>%
  ggplot(aes(DP)) +
  geom_histogram() +
  scale_x_log10() + 
  facet_wrap(~as.numeric(CHROM), ncol = 8) +
  labs(x="Read depth", y= "N of bases", title = "Distribution of read depth per chromosome") +
  theme(axis.text.x = element_text(angle = 90))

#Distribution of the DP by chromosome
data %>%
  ggplot(aes(CHROM)) +
  geom_bar() +
  labs(x="Chromosome", y= "N of bases", title = "Distribution of read depth by chromosome") +
  theme(axis.text.x = element_text(angle = 90))
