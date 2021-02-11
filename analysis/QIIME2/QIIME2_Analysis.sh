### QIIME 2 scripts for Cleland/Mount Lagoon koala faecal 16S paper
### Raphael Eisenhofer, 10/2/2021

conda activate qiime2-2020.6

#Trim primer sequences from Mount Lagoon data
for i in 0_Reads/R1/ML*.fastq.gz;
  do cutadapt -j 5 -g GTGYCAGCMGCCGCGGTAA -o ${i/.fastq.gz/_trimmed.fastq.gz} $i;
done

#Import reads into qiime2
qiime tools import \
--type SampleData[SequencesWithQuality] \
--input-path manifesto.txt \
--input-format SingleEndFastqManifestPhred33V2 \
--output-path demux-r1.qza

qiime demux summarize \
--i-data demux-r1.qza \
--o-visualization demux-r1.qzv
##Qualities look good, though position 3 looks very dodgy.
##Will trim the first 3 bp -- can do this via deblur

#Run deblur with a 150 bp trim length, trimming the first 3 bp from 5' end
qiime deblur denoise-16S \
--i-demultiplexed-seqs demux-r1.qza \
--p-left-trim-len 3 \
--p-trim-length 150 \
--p-sample-stats \
--o-representative-sequences Koala-rep-seqs.qza \
--o-table Koala-table.qza \
--o-stats Koala-deblur-stats.qza \
--verbose \
--p-jobs-to-start 6

qiime feature-table tabulate-seqs \
--i-data Koala-rep-seqs.qza \
--o-visualization Koala-rep-seqs.qzv

qiime feature-table summarize \
--i-table Koala-table.qza \
--o-visualization Koala-table.qzv

#Assign taxonomy using a naive bayesian classifier against the SILVA 138 v4 database
qiime feature-classifier classify-sklearn \
--i-reads Koala-rep-seqs.qza \
--i-classifier silva-138-99-515-806-nb-classifier.qza \
--o-classification Koala-SILVA-138.qza \
--p-n-jobs 16

qiime metadata tabulate \
--m-input-file Koala-SILVA-138.qza \
--o-visualization Koala-SILVA-138.qzv

#Create a phylogenetic tree using SEPP
qiime fragment-insertion sepp \
--i-representative-sequences Koala-rep-seqs.qza \
--i-reference-database sepp-refs-silva-128.qza \
--o-tree Koala-sepp-tree.qza \
--o-placements Koala-sepp-placements.qza \
--verbose \
--p-threads 16

#Filter out features that couldn't be inserted on tree:
qiime fragment-insertion filter-features \
  --i-table Koala-table.qza \
  --i-tree Koala-sepp-tree.qza \
  --o-filtered-table Koala-filtered_table.qza \
  --o-removed-table Koala-removed_table.qza \
  --verbose

  qiime feature-table summarize \
  --i-table Koala-removed_table.qza \
  --o-visualization Koala-removed_table.qzv
##Looks like nothing was removed, so will just use the unfiltered table:
##(table.qza) for subsequent analyses

#Run alpha-rarefaction analysis to see what a good sampling depth will be
qiime diversity alpha-rarefaction \
--i-table Koala-table.qza \
--m-metadata-file Koala_Metadata.tsv \
--o-visualization table-rarefaction.qzv \
--p-min-depth 500 \
--p-max-depth 50000
##Diversity seems to plateau ~30,000 sequences. The same with the lowset # of
##sequences = Cleland-R46, with 36,622 sequences.
## Let's set the rarefaction depth to 36,622 to keep all of our samples :-)

#Core diversity analyses
qiime diversity core-metrics-phylogenetic \
--i-table Koala-table.qza \
--i-phylogeny Koala-sepp-tree.qza \
--p-sampling-depth 36622 \
--m-metadata-file Koala_Metadata.tsv \
--output-dir Core-metrics-both-populations-36622

#Loop through alpha diversity files and compute significance
for i in Core-metrics-both-populations-36622/*vector.qza; \
  do qiime diversity alpha-group-significance \
     --i-alpha-diversity $i \
     --m-metadata-file Koala_Metadata.tsv \
     --o-visualization ${i/.qza/.qzv};
done

#Make some exploratory taxa bar plots
qiime taxa barplot \
--i-table Koala-table.qza \
--i-taxonomy Koala-SILVA-138.qza \
--m-metadata-file Koala_Metadata.tsv \
--o-visualization Koala-table-BarPlots-by-sample.qzv

#Split table by population
qiime feature-table filter-samples \
--i-table Koala-table.qza \
--o-filtered-table Koala-Cleland-table.qza \
--m-metadata-file Koala_Metadata.tsv \
--p-where "[Population]='Mountain_Lagoon'" \
--p-exclude-ids True

qiime feature-table filter-samples \
--i-table Koala-table.qza \
--o-filtered-table Koala-MtLagoon-table.qza \
--m-metadata-file Koala_Metadata.tsv \
--p-where "[Population]='Mountain_Lagoon'"

#Core diversity analyses on split tables
qiime diversity core-metrics-phylogenetic \
--i-table Koala-Cleland-table.qza \
--i-phylogeny Koala-sepp-tree.qza \
--p-sampling-depth 36622 \
--m-metadata-file Koala_Metadata.tsv \
--output-dir Core-metrics-Cleland-36622

qiime diversity core-metrics-phylogenetic \
--i-table Koala-MtLagoon-table.qza \
--i-phylogeny Koala-sepp-tree.qza \
--p-sampling-depth 36622 \
--m-metadata-file Koala_Metadata.tsv \
--output-dir Core-metrics-MtLagoon-36622

for i in Core-metrics-*-36622/*vector.qza; \
  do qiime diversity alpha-group-significance \
     --i-alpha-diversity $i \
     --m-metadata-file Koala_Metadata.tsv \
     --o-visualization ${i/.qza/.qzv};
done
