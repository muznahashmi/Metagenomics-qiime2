mkdir single_end_qza
# Import data for Sample 1
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_sample1.txt \
  --output-path single_end_qza/single-end_sample1.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Import data for Sample 2
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_sample2.txt \
  --output-path single_end_qza/single-end_sample2.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Import data for Sample 3
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_sample3.txt \
  --output-path single_end_qza/single-end_sample3.qza \
  --input-format SingleEndFastqManifestPhred33V2




  # qza to qzv

  INPUT_DIR="single_end_qza"


  for SAMPLE_FILE in "$INPUT_DIR"/*.qza
do
  # Get the sample ID from the file name
  SAMPLE_ID=$(basename "$SAMPLE_FILE" .qza)

  # Define the output visualization file path
  OUTPUT_VISUALIZATION="$INPUT_DIR/${SAMPLE_ID}.qzv"

  # Run the qiime demux summarize command
  qiime demux summarize \
    --i-data "$SAMPLE_FILE" \
    --o-visualization "$OUTPUT_VISUALIZATION"

  echo "Generated visualization for $SAMPLE_ID"
done


## dada-denoising


INPUT_DIR="single_end_qza"

# Loop through the sample files
for SAMPLE_FILE in "$INPUT_DIR"/*.qza
do
  # Get the sample ID from the file name
  SAMPLE_ID=$(basename "$SAMPLE_FILE" .qza)

  # Define output file names
  REP_SEQS_OUT="$INPUT_DIR/rep-seqs-dada2_${SAMPLE_ID}.qza"
  TABLE_OUT="$INPUT_DIR/table-dada2_${SAMPLE_ID}.qza"
  STATS_OUT="$INPUT_DIR/stats-dada2_${SAMPLE_ID}.qza"

  # Run the qiime dada2 denoise-single command
  qiime dada2 denoise-single \
    --i-demultiplexed-seqs "$SAMPLE_FILE" \
    --p-trunc-len 429 \
    --p-n-threads 32 \
    --p-chimera-method consensus \
    --p-min-fold-parent-over-abundance 8 \
    --o-representative-sequences "$REP_SEQS_OUT" \
    --o-table "$TABLE_OUT" \
    --o-denoising-stats "$STATS_OUT"

  echo "Processed denoising for $SAMPLE_ID"
done

# qza to qzv


# Define the path to the input directory
INPUT_DIR="single_end_qza"

# Loop through the sample files
for SAMPLE_FILE in "$INPUT_DIR"/table-dada2_*.qza
do
  # Get the sample ID from the file name
  SAMPLE_ID=$(basename "$SAMPLE_FILE" .qza | cut -d'_' -f3-)
  
  # Define the output visualization file path
  OUTPUT_VISUALIZATION="$INPUT_DIR/table_${SAMPLE_ID}.qzv"

  # Run the qiime feature-table summarize command
  qiime feature-table summarize \
    --i-table "$SAMPLE_FILE" \
    --o-visualization "$OUTPUT_VISUALIZATION" \
    --m-sample-metadata-file metadata_PE.tsv

  echo "Generated visualization for $SAMPLE_ID"
done


# qza to qzv

# Define the path to the input directory
INPUT_DIR="single_end_qza"

# Loop through the sample files
for SAMPLE_FILE in "$INPUT_DIR"/rep-seqs-dada2_*.qza
do
  # Get the sample ID from the file name
  SAMPLE_ID=$(basename "$SAMPLE_FILE" .qza | cut -d'_' -f3-)
  
  # Define the output visualization file path
  OUTPUT_VISUALIZATION="$INPUT_DIR/rep-seqs-dada2_${SAMPLE_ID}.qzv"

  # Run the qiime feature-table tabulate-seqs command
  qiime feature-table tabulate-seqs \
    --i-data "$SAMPLE_FILE" \
    --o-visualization "$OUTPUT_VISUALIZATION"

  echo "Generated visualization for $SAMPLE_ID"
done



qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path silva_132_97_16S.fna\
  --output-path silva_132_97_16S.qza

  qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --input-format HeaderlessTSVTaxonomyFormat \
   --input-path consensus_taxonomy_7_levels.txt \
   --output-path consensus_taxonomy_7_levels.qza


  #!/bin/bash

# Define the path to the input directory
INPUT_DIR="single_end_qza"

# Define the reference sequences file
REFERENCE_SEQS="silva_132_97_16S.qza"

# Loop through the sample files
for SAMPLE_ID in 1 2 3
do
  # Define the paths for the current sample
  SAMPLE_TABLE="$INPUT_DIR/table-dada2_single-end_sample${SAMPLE_ID}.qza"
  SAMPLE_SEQS="$INPUT_DIR/rep-seqs-dada2_single-end_sample${SAMPLE_ID}.qza"

  # Define the output files
  OUTPUT_TABLE="$INPUT_DIR/table-cr-97_sample${SAMPLE_ID}.qza"
  OUTPUT_SEQS="$INPUT_DIR/rep-seqs-cr-97_sample${SAMPLE_ID}.qza"
  OUTPUT_UNMATCHED="$INPUT_DIR/unmatched-cr-97_sample${SAMPLE_ID}.qza"

  # Run the qiime vsearch cluster-features-closed-reference command
  qiime vsearch cluster-features-closed-reference \
    --i-table "$SAMPLE_TABLE" \
    --i-sequences "$SAMPLE_SEQS" \
    --i-reference-sequences "$REFERENCE_SEQS" \
    --p-perc-identity 0.97 \
    --o-clustered-table "$OUTPUT_TABLE" \
    --o-clustered-sequences "$OUTPUT_SEQS" \
    --o-unmatched-sequences "$OUTPUT_UNMATCHED"

  echo "Processed closed-reference clustering for sample ${SAMPLE_ID}"
done


qiime feature-classifier extract-reads \
  --i-sequences silva_132_97_16S.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTA \
  --p-r-primer N \
  --p-trunc-len 0 \
  --o-reads gene-wiz-v4-ref-seqs.qza

#TRAINING ----omitted for now
 # qiime feature-classifier fit-classifier-naive-bayes \
 #   --i-reference-reads  gene-wiz-v4-ref-seqs.qza \
 #   --i-reference-taxonomy  consensus_taxonomy_7_levels.qza \
 #   --o-classifier  gene-wiz-v4-classifier.qza 

#trained data obtained from: https://docs.qiime2.org/2019.7/data-resources/

# Silva 132 99% OTUs from 515F/806R region of sequences (MD5: a0925e86cda18829f84f03dab01ff589)


# taxonomy assignment

pip install scikit-learn==0.24.1

INPUT_DIR="single_end_qza"

# Define the classifiers for each sample
CLASSIFIER_SAMPLE="silva-138-99-515-806-nb-classifier.qza"

# Define the output directory for classifications
OUTPUT_DIR="taxonomy_outputs"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through the samples
for SAMPLE_ID in 1 2 3
do
  # Define the paths for the current sample
  SAMPLE_SEQS="$INPUT_DIR/rep-seqs-dada2_single-end_sample${SAMPLE_ID}.qza"


  # Define the output file
  OUTPUT_CLASSIFICATION="$OUTPUT_DIR/taxonomy_sample${SAMPLE_ID}.qza"

  # Run the classification command
  qiime feature-classifier classify-sklearn \
    --i-reads "$SAMPLE_SEQS" \
    --i-classifier "$CLASSIFIER_SAMPLE" \
    --o-classification "$OUTPUT_CLASSIFICATION"

  echo "Processed classification for sample ${SAMPLE_ID}"
done


# table-filtering

# Define the path to the input directory
INPUT_DIR_1="single_end_qza"
INPUT_DIR_2="taxonomy_outputs"

# Define the output directory for filtered tables
OUTPUT_DIR="filtered_tables"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through the samples
for SAMPLE_ID in 1 2 3
do
  # Define the paths for the current sample
  SAMPLE_TABLE="$INPUT_DIR_1/table-dada2_single-end_sample${SAMPLE_ID}.qza"
  SAMPLE_TAXONOMY="$INPUT_DIR_2/taxonomy_sample${SAMPLE_ID}.qza"
  
  # Define the output file
  OUTPUT_TABLE="$OUTPUT_DIR/filter-table_sample${SAMPLE_ID}.qza"

  # Run the filter-table command
  qiime taxa filter-table \
    --i-table "$SAMPLE_TABLE" \
    --i-taxonomy "$SAMPLE_TAXONOMY" \
    --p-exclude Eukaryota \
    --o-filtered-table "$OUTPUT_TABLE"

  echo "Processed filter-table for sample ${SAMPLE_ID}"
done


# SEQS-filtering

# Define the path to the input directory
INPUT_DIR_1="single_end_qza"
INPUT_DIR_2="taxonomy_outputs"

# Define the output directory for filtered tables
OUTPUT_DIR="filtered_seqs"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through the samples
for SAMPLE_ID in 1 2 3
do
  # Define the paths for the current sample
  SAMPLE_SEQS="$INPUT_DIR_1/rep-seqs-dada2_single-end_sample${SAMPLE_ID}.qza"
  SAMPLE_TAXONOMY="$INPUT_DIR_2/taxonomy_sample${SAMPLE_ID}.qza"
  
  # Define the output file
  OUTPUT_SEQS="$OUTPUT_DIR/filter-seqs_sample${SAMPLE_ID}.qza"

  # Run the filter-table command
  qiime taxa filter-seqs \
    --i-sequences "$SAMPLE_SEQS" \
    --i-taxonomy "$SAMPLE_TAXONOMY" \
    --p-exclude Eukaryota \
    --o-filtered-sequences "$OUTPUT_SEQS"

  echo "Processed filter-table for sample ${SAMPLE_ID}"
done


#merger filter-dada table

# Create a new directory for the merged table
mkdir merge_folder

# Merge the tables
qiime feature-table merge \
  --i-tables filtered_tables/filter-table_sample1.qza \
  --i-tables filtered_tables/filter-table_sample2.qza \
  --i-tables filtered_tables/filter-table_sample3.qza \
  --o-merged-table merge_folder/merged-table.qza


  #merge filter seqs

qiime feature-table merge-seqs \
  --i-data filtered_seqs/filter-seqs_sample1.qza \
  --i-data filtered_seqs/filter-seqs_sample2.qza \
  --i-data filtered_seqs/filter-seqs_sample3.qza \
  --o-merged-data merge_folder/merged-seqs.qza


#merge taxa

qiime feature-table merge-taxa \
    --i-data taxonomy_outputs/taxonomy_sample1.qza \
    --i-data taxonomy_outputs/taxonomy_sample2.qza \
    --i-data taxonomy_outputs/taxonomy_sample3.qza \
    --o-merged-data merge_folder/merged-taxonomy.qza



#bar blot

qiime taxa barplot \
    --i-table merge_folder/merged-table.qza \
    --i-taxonomy merge_folder/merged-taxonomy.qza \
    --m-metadata-file metadata_PE.tsv \
    --o-visualization merge_folder/merged-taxa-bar-plots.qzv
done

##################__________________________________----___________________----____________

#lets try tree from here:

#  https://john-quensen.com/tutorials/merging-dada2-results-in-qiime2/


#To align and tree the representative sequences, enter the QIIME2 command:


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merge_folder/merged-seqs.qza \
  --o-alignment merge_folder/merged-aligned-seqs.qza \
  --o-masked-alignment merge_folder/merged-masked-aligned-seqs.qza \
  --o-tree merge_folder/merged-unrooted-tree.qza \
  --o-rooted-tree merge_folder/merged-rooted-tree.qza

#Export the ASV table in biom format with the command:

mkdir phyloseq
qiime tools export \
  --input-path merge_folder/merged-table.qza \
  --output-path phyloseq

#Convert the biom format to a tab-separated file easily imported into R with the commands:

biom convert \
  -i phyloseq/feature-table.biom \
  -o phyloseq/otu_table.tsv \
  --to-tsv
cd phyloseq
sed -i '1d' otu_table.tsv 
sed -i 's/#OTU ID//' otu_table.tsv 
cd ..


#Export the representative sequences and tree files with the QIIME2 commands:

qiime tools export \
  --input-path merge_folder/merged-seqs.qza \
  --output-path phyloseq
qiime tools export \
  --input-path merge_folder/merged-unrooted-tree.qza \
  --output-path phyloseq
cd phyloseq
mv tree.nwk merged_unrooted_tree.nwk
cd ../

qiime tools export \
  --input-path merge_folder/merged-rooted-tree.qza \
  --output-path phyloseq
cd phyloseq
mv tree.nwk rooted_tree.nwk

###################################



#Alpha diversity


library(qiime2R)
library(tidyverse)
library(phyloseq)

physeq<-qza_to_phyloseq(
  features="merged-table.qza",
  tree="merged-rooted-tree.qza",
  taxonomy="merged-taxonomy.qza",
  metadata = "metadata_PE.tsv"
)

richness <- estimate_richness(physeq)
head(richness)
# Create a boxplot of Shannon diversity measures for each sample, with the sample IDs on the x-axis
p <- plot_richness(physeq, x="Group", measures="Shannon", color = "Group") +
  geom_boxplot(alpha=0.6) + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# Get the sample IDs from the metadata file
sample_ids <- as.character(physeq@sam_data$`body-site`)

# Modify the x-axis labels to show the sample IDs
p + scale_x_discrete(labels = sample_ids)


##################################################################################

# alpha-beta diversity in R
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(vegan)

physeq <- qza_to_phyloseq(
  features = "merged-table.qza",
  tree = "merged-rooted-tree.qza",
  taxonomy = "merged-taxonomy.qza",
  metadata = "metadata_PE.tsv"
)

# Estimate richness
richness <- estimate_richness(physeq)

# Create a boxplot of Shannon diversity measures for each sample, with the sample IDs on the x-axis
p <- plot_richness(physeq, x = "Group", measures = "Shannon", color = "Group") +
  geom_boxplot(alpha = 0.6) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12))

# Get the sample IDs from the metadata file
sample_ids <- as.character(physeq@sam_data$`body-site`)

# Modify the x-axis labels to show the sample IDs
p <- p + scale_x_discrete(labels = sample_ids)

# Save the boxplot as TIFF
tiff("boxplot_shannon_diversity.tiff", width = 2000, height = 1500, res = 400)
print(p)
dev.off()

# Perform PCoA on weighted UniFrac distance
wunifrac_dist <- phyloseq::distance(physeq, method = "unifrac", weighted = TRUE)
ordu <- ordinate(physeq, "PCoA", "unifrac", distance = wunifrac_dist)

# Plot and save PCoA plot as TIFF
p <- plot_ordination(physeq, ordu, color = "Group", shape = "Group") +
  geom_point(size = 1, alpha = 0.75) + 
  scale_colour_brewer(type = "qual", palette = "Set1") + 
  ggtitle("MDS/PCoA on weighted-UniFrac distance")

tiff("MDS-PcoA-weighted-UniFrac-distance.tiff", width = 2000, height = 1500, res = 400)
print(p)
dev.off()

# Perform PCoA on unweighted UniFrac distance
Unwunifrac_dist <- phyloseq::distance(physeq, method = "unifrac", weighted = FALSE)
ordu <- ordinate(physeq, "PCoA", "unifrac", distance = Unwunifrac_dist)

# Plot and save PCoA plot as TIFF
p <- plot_ordination(physeq, ordu, color = "Group", shape = "Group") +
  geom_point(size = 1, alpha = 0.75) + 
  scale_colour_brewer(type = "qual", palette = "Set1") + 
  ggtitle("MDS/PCoA on Unweighted-UniFrac distance")

tiff("MDS-PcoA-Unweighted-UniFrac-distance.tiff", width = 2000, height = 1500, res = 400)
print(p)
dev.off()

# Perform adonis
adonis_result <- adonis(wunifrac_dist ~ sample_data(physeq)$Group)
print(adonis_result)








#########

#alpha diversity from qiime2

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny merge_folder/merged-rooted-tree.qza \
  --i-table merge_folder/merged-table.qza \
  --p-sampling-depth 401\
  --m-metadata-file merge_folder/metadata_PE.tsv \
  --output-dir core-metrics-results

  ##########################


  #krona plot

  qiime krona collapse-and-plot \
--i-table merge_folder/merged-table.qza \
--i-taxonomy merge_folder/merged-table.qza \
--o-krona-plot krona.qzv