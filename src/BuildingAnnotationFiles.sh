#!/usr/bin/env bash

#building the annotation files for intronic enrichment score
# source of coding UCSC tab (https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/)
# sacCer3AresIntron locations were downloaded from ares intron database and formated into bed 

cd ../AnnotationFileIntermediates

# sacCer3AresIntron locations were downloaded from ares intron database and formated into bed 

# removes chr from chromosome name to be consistent with bam file
sed -i '' 's/chr//' *.tab
sed -i '' 's/%28/(/' *.tab
sed -i '' 's/%29/)/' *.tab
sed -i '' 's/%29/)/' *.tab
# pulls full gene length (TSS to TES from feature table)
grep "\ttranscript\t" sacCer3UCSCallFeatures.tab > sacCer3UCSCgenes.tab

# removes mitochondrial loci from intron features 
grep -v "mt\t" full_introns.tab > introns.tab
grep "intron\t" full_introns.tab > full_intron.tab

bedtools sort -i introns.tab > introns_sorted.tab
bedtools sort -i sacCer3UCSCgenes.tab > sacCer3UCSCgenes_sorted.tab

# extracts gene names from annotation files
awk '{print $9}' introns_sorted.tab | awk -F '[_=;""]' '{print $2}' > inttemp.txt
awk '{print $10}' sacCer3UCSCgenes_sorted.tab | awk -F '[_=;"]' '{print $2}' > exotemp.txt  

# replaces gene identifier info with just their names 
awk 'FNR==NR{a[NR]=$1;next}{$9=a[FNR]}1' inttemp.txt introns_sorted.tab > introns_abbr.tab
awk 'FNR==NR{a[NR]=$1;next}{$9=a[FNR]}1' exotemp.txt sacCer3UCSCgenes_sorted.tab | \
awk '{$10=$11=$12=$13=$14=""; print $0}' > gene_abbr.tab

# finds common genes in annotation files 
# genes in exons that are also in intron
grep -xf inttemp.txt exotemp.txt > commongenes_temp.txt
# genes in introns that are also in exons then pull out unique names
grep -xf commongenes_temp.txt inttemp.txt | sort | uniq > commongenes.txt

# extracts lines with common genes from annotation files BUT SOME GENES HAVE -A at the end and get included on accident
grep -wf commongenes.txt gene_abbr.tab > genes_with_introns.tab

# now I have intron and gene with intron annotation files with identical naming schemes 
# but some slip through that are only found in one or the other 
awk '{print $9}' introns_abbr.tab | sort | uniq > in_genenames.txt
awk '{print $9}' genes_with_introns.tab | sort | uniq > gene_genenames.txt 

# this just tells you what is going to be left out in the next command 
# genes in gene file but not intron file
grep -xvf in_genenames.txt gene_genenames.txt > genesInOneButNotTheOther.txt 
# genes in intron file but not in gene file
grep -xvf gene_genenames.txt in_genenames.txt >> genesInOneButNotTheOther.txt 

# removes genes found in one but not the other
grep -vwf genesInOneButNotTheOther.txt genes_with_introns.tab > genes_with_introns_filtered.tab
grep -vwf genesInOneButNotTheOther.txt introns_abbr.tab > introns_filtered.tab

# awk '{print $9}' introns_filtered.tab | sort | uniq 
# awk '{print $9}' genes_with_introns_filtered.tab | sort | uniq 


awk '{print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' genes_with_introns_filtered.tab > genes_with_introns_formatted.bed
awk '{print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' introns_filtered.tab > ../introns_formatted.bed

# makes sure tRNAs and snoRNAs are labeled as such
sed -i '' 's/	transcript	.	+	.	Y/	gene	.	+	.	Y/' *.bed
sed -i '' 's/	transcript	.	+	.	t/	tRNA	.	+	.	t/' *.bed
sed -i '' 's/	transcript	.	+	.	s/	snoRNA	.	+	.	s/' *.bed
sed -i '' 's/	transcript	.	-	.	Y/	gene	.	-	.	Y/' *.bed
sed -i '' 's/	transcript	.	-	.	t/	tRNA	.	-	.	t/' *.bed
sed -i '' 's/	transcript	.	-	.	s/	snoRNA	.	-	.	s/' *.bed

subtractBed -a genes_with_introns_formatted.bed -b ../introns_formatted.bed > ../exons_formatted.bed

cd ../src