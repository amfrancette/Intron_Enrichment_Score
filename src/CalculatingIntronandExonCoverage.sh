#!/usr/bin/env bash
bigwigfolder=Barass_bigwig_good

# using multiBigwigSummary in deeptools to calculate read depth over exons and introns
# sources of example bigwig data is Barass et al. 2015 doi: 10.1186/s13059-015-0848-1 
# and unpublished data from Dr. Karen Arndt's Lab from the University of Pittsburgh

cd ..

# get average read depth for all introns and exons of genes with introns
multiBigwigSummary BED-file \
 --bwfiles "$bigwigfolder"/*.bw \
 --BED exons_formatted.bed \
 --smartLabels \
 -out res/bw_exon_scores.npz --outRawCounts res/bw_exon_scores.tab
# removes header
sed 1d res/bw_exon_scores.tab > bw_exon_scores_decap.tab

multiBigwigSummary BED-file \
 --bwfiles "$bigwigfolder"/*.bw \
 --BED introns_formatted.bed \
 --smartLabels \
 -out res/bw_intron_scores.npz --outRawCounts res/bw_intron_scores.tab
 # removes header
 sed 1d res/bw_intron_scores.tab > bw_intron_scores_decap.tab
 
 # preserves the column names
 head -n 1 res/bw_intron_scores.tab > res/colnames.tab
 
# Maintains gene names for each location
awk '{print $9 "\t" $5}' exons_formatted.bed | paste - bw_exon_scores_decap.tab > res/exon_summary.tab
awk '{print $9 "\t" $5}' introns_formatted.bed | paste - bw_intron_scores_decap.tab > res/intron_summary.tab
rm *decap*

# now move to the R script to finish the analysis 

cd src

echo "now move to the R script to finish the analysis"
