library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)

# use this to reset the location of the pwd if the here command fails
# set_here(path = ".", verbose = TRUE)

# # import per-base coverage info over exons and introns (multiBigWigSummary in deeptools)
# # need to manually insert sample names
# # with these set up correctly, the following commands should work until the inner_join command
exon <- read.csv(here("/res/exon_summary.tab"), header = F, sep = "\t")
intron <- read.csv(here("/res/intron_summary.tab"), header = F, sep = "\t")

# this will just print the sample order to help format column names until I can automate it in
read.csv(here("/res/colnames.tab"), header = F, sep = "\t")
column_names <- c("locus","feature", "chr","start","end", "Rep1_1point5min", "Rep2_1point5min", 
                  "Rep1_2point5min", "Rep2_2point5min",
                  "Rep1_5min", "Rep2_5min", "Rep1_Background" , 
                  "Rep2_Total", "WH3_4tU_Rep1", "WH3_4tU_Rep2", "WH3_Total_Rep1","WH3_Total_Rep2")

number_of_samples <- ncol(exon)-5

names(intron) <- column_names
intron$in_length <- intron$end - intron$start
intron_sum <- cbind(intron[,1:5] , intron[,6:(5+number_of_samples)] * intron$in_length, "in_length" = intron$in_length)
intron_summary <- aggregate(. ~ locus, data=intron_sum, sum)
in_norm_avg <- cbind(intron_summary[,1:5], intron_summary[,6:(5+number_of_samples)] / intron_summary$in_length, intron_summary$in_length) 
in_norm_avg <- in_norm_avg[order(in_norm_avg$locus),]

names(exon) <-  column_names
exon$ex_length <- exon$end - exon$start
exon_sum <- cbind(exon[,1:5] , exon[,6:(5+number_of_samples)] * exon$ex_length, "ex_length" = exon$ex_length)
exon_summary <- aggregate(. ~ locus+feature, data=exon_sum, sum)
ex_norm_avg <- cbind(exon_summary[,1:5], exon_summary[,6:(5+number_of_samples)] / exon_summary$ex_length, exon_summary$ex_length) 
ex_norm_avg <- ex_norm_avg[order(ex_norm_avg$locus),]

setdiff(intron_summary[,1],exon_summary[,1])
setdiff(exon_summary[,1],intron_summary[,1])
# This line is to find any stragglers that are in one data set but not the other 
# summary <- anti_join(in_norm_avg, ex_norm_avg, by = "locus")
# then grep away the offending locus and restart analysis (without re-importing the tab file)
# in_norm_avg <- in_norm_avg[-grep("PUTLOCUSNAMEHERE", in_norm_avg$locus),]

# combines normalized per/base average coverage in exons and introns into a single dataframe
inner_join(ex_norm_avg, in_norm_avg, by = "locus")
IES <- cbind(ex_norm_avg[,1:2], 
             in_norm_avg[,6:(5+number_of_samples)] / ex_norm_avg[,6:(5+number_of_samples)])
colnames(IES)[1] <- "locus"

# replaces NaN w/ 0
IES[is.na(IES)] <- 0

# this line displays the average IES per colum
  apply(IES[3:(2 + number_of_samples)], 2, mean, na.rm = T)
# only counting genes in average
  apply(IES[IES$feature == "gene", 3:(2 + number_of_samples)], 2, mean, na.rm = T)

# this line displays the median IES per colum
  apply(IES[3:(2 + number_of_samples)], 2, median, na.rm = T)
  
  
# counting only genes in median
apply(IES[IES$feature == "gene", 3:(2 + number_of_samples)], 2, median, na.rm = T)
  
# generates row means of nascent and total RNA and puts them in a labled matrix 
### MODIFY LATER
colnames(IES)
aggregateIES <- cbind(IES[,1:2], rowMeans(IES[,3:4]), 
                      rowMeans(IES[,5:6]), rowMeans(IES[,7:8]), IES[,10],  rowMeans(IES[11:12]),  rowMeans(IES[13:14]))
colnames(aggregateIES)[3:8] <- c("min_1.5", "min_2.5", "min_5", "Total", "WH3_4tU", "WH3_Total")

# Plots Histograms of Barass IESs
ggplot(aggregateIES[aggregateIES$feature=="gene",], aes(x = Total)) + geom_histogram(binwidth=.03,color="black", fill="blue", alpha = 0.4) +
  geom_histogram(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_1.5), binwidth=.03 ,color="black", fill="orange", alpha = 0.8) +
  labs(y="# of Observations", x = "Intron Enrichment Score") +
  theme_grey(base_size = 15) +
  geom_vline(xintercept=1) + xlim(0, 2)


colnames(aggregateIES)
methodDif <- cbind(aggregateIES[,1:2], log2(aggregateIES[,3]/aggregateIES[,6]), 
                   log2(aggregateIES[,4]/aggregateIES[,6]), 
                   log2(aggregateIES[,5]/aggregateIES[,6]))
colnames(methodDif)[3:5] <- c("min_1.5vsTotal", "min_2.5vsTotal", "min_5vsTotal")


# Comparing Labeling Times
pdf("imgs/Barass_IESdist.pdf")
ggplot(aggregateIES[aggregateIES$feature=="gene",], aes(x = Total)) + geom_density(color="black", fill="red", alpha = 0.3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_5), color="black", fill="orange", alpha = 0.3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_2.5), color="black", fill="yellow", alpha = .3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_1.5), color="black", fill="green", alpha = .3) +
  labs(y="Density of Observations", x = "Intron Enrichment Score") +
  theme_grey(base_size = 15) +
  geom_vline(xintercept=1) + xlim(0, 1.5)
dev.off()

# smoothed plots
pdf("imgs/deltaIESdistribution_Barass.pdf")
ggplot(methodDif[methodDif$feature=="gene",], aes(x = min_1.5vsTotal)) + geom_density(color="black", fill="green", alpha = 0.3) +
    geom_density(data = methodDif[methodDif$feature=="gene",], aes(x = min_2.5vsTotal) ,color="black", fill="yellow", alpha = 0.3) +
    geom_density(data = methodDif[methodDif$feature=="gene",], aes(x = min_5vsTotal) ,color="black", fill="orange", alpha = 0.3) +
    labs(y="Density of Observations", x = "Log2fc Difference in Intron Enrichment Score") +
    theme_grey(base_size = 15) +
    geom_vline(xintercept=0) + xlim(-2, 6)
dev.off()

# Comparing WH3 Total & Nascent Data to Barass
ggplot(aggregateIES[aggregateIES$feature=="gene",], aes(x = WH3_Total)) + geom_histogram(binwidth=.02,color="black", fill="purple", alpha = 0.1) +
  geom_histogram(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = WH3_4tU), binwidth=.02 ,color="black", fill="blue", alpha = 0.1) +
  geom_histogram(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = Total), binwidth=.02 ,color="black", fill="red", alpha = 1) +
  labs(y="# of Observations", x = "Intron Enrichment Score") +
  theme_grey(base_size = 15) +
  geom_vline(xintercept=1) + xlim(0, 1.5)

# smoothed plots
pdf("IESdistribution_ALL.pdf")
ggplot(aggregateIES[aggregateIES$feature=="gene",], aes(x = Total)) + geom_density(color="black", fill="red", alpha = 0.3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_5), color="black", fill="orange", alpha = 0.3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_2.5), color="black", fill="yellow", alpha = .3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = min_1.5), color="black", fill="green", alpha = .3) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = WH3_4tU), color="black", fill="grey50", alpha = 0.4) +
  geom_density(data = aggregateIES[aggregateIES$feature=="gene",], aes(x = WH3_Total) ,color="black", fill="purple", alpha = .4) +
  labs(y="Density of Observations", x = "Intron Enrichment Score") +
  theme_grey(base_size = 30) +
  geom_vline(xintercept=1) + xlim(0, 1.5)
dev.off()

colnames(IES)
# Comparison Between REPS
Compare_x <- c("Rep2_5min")
Compare_y <- c("Rep1_5min")
limits_x <- c(0,2)
ggplot(data = IES[IES$feature=="gene",], aes_(as.name(Compare_x), as.name(Compare_y), color=quote(feature))) + geom_point() +
  scale_x_continuous(name= as.name(Compare_x), limits = limits_x) +
  scale_y_continuous(name= as.name(Compare_y), limits = limits_x) +
  theme_grey(base_size = 18) +
  theme(legend.position = "none") +
  geom_abline(slope=1, intercept = 0)


Compare_x <- c("WH3_4tU")
Compare_y <- c("WH3_Total")
limits_x <- c(0,2)
ggplot(data = aggregateIES[aggregateIES$feature=="gene",], aes_(as.name(Compare_x), as.name(Compare_y), color=quote(feature))) + geom_point() + 
  scale_x_continuous(name= as.name(Compare_x), limits = limits_x) +
  scale_y_continuous(name= as.name(Compare_y), limits = limits_x) +
  theme_grey(base_size = 18) +
  theme(legend.position = "none") +
  geom_abline(slope=1, intercept = 0)

Compare_x <- c("WH3_Total")
Compare_y <- c("Total")
limits_x <- c(0,1)
ggplot(data = aggregateIES[aggregateIES$feature=="gene",], aes_(as.name(Compare_x), as.name(Compare_y), color=quote(feature))) + geom_point() + 
  scale_x_continuous(name= as.name(Compare_x), limits = limits_x) +
  scale_y_continuous(name= as.name(Compare_y), limits = limits_x) +
  theme_grey(base_size = 18) +
  theme(legend.position = "none") +
  geom_abline(slope=1, intercept = 0)

