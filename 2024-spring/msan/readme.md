## Identification of Differentially Expressed Genes from Marine Heterotrophic Bacterial RNAseq data
### Introduction
RNAseq data is publicly available for a marine heterotrophic bacteria grown on different carbon sources. For example, the raw reads of *Zobellia galactanivorans* are found in **SRA** (sequence archive) of **NCBI**.
  - where can I get the public data (raw reads)? [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
  - search keyword: `(marine bacteria carbon source) AND "Zobellia galactanivorans"[orgn:__txid63186] `

#### Background
* Purpose: To identify transcriptional changes in gene levels or to understand mechanisms that occur under certain conditions and either suppress them or induce more expression.
-------
![image](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165700031/1ac46226-16a1-4c94-a514-24e8c32f4599)

  - Check Quality : The process of filtering out low-quality parts of something fetched via Rfastp
  - Read Sample : Sorting and calculations for the sequences 
  - Counting in exons models : Differentiate exons by collapsing exons to the nonoverlapping, disjoint regions
  - Check Mean and variance relationship : This will assist us to detect differential expression with a small number of replicates.
  - Multiple testing : If we are 95% confident a gene is differentially expressed when we look at all genes we can expect 5% to be false. When looking across all genes we apply a multiple testing correction to account for this.
  - Genes overlaps : Compare multiple results to clean up overlaps.

------
### Data list
  
  - GSM5699202: Free-living cells with maltose rep 1; Zobellia galactanivorans; RNA-Seq   
 https://www.ncbi.nlm.nih.gov/sra/SRX13191408[accn]
  - GSM5699203: Free-living cells with maltose rep 2; Zobellia galactanivorans; RNA-Seq    
https://www.ncbi.nlm.nih.gov/sra/SRX13191409[accn]
  - GSM5699204: Free-living cells with maltose rep 3; Zobellia galactanivorans; RNA-Seq    
https://www.ncbi.nlm.nih.gov/sra/SRX13191410[accn]
  - GSM5699196: Free-living cells with alginate rep 1; Zobellia galactanivorans; RNA-Seq
  - https://www.ncbi.nlm.nih.gov/sra/SRX13191402[accn]
  - GSM5699197: Free-living cells with alginate rep 2; Zobellia galactanivorans; RNA-Seq   
https://www.ncbi.nlm.nih.gov/sra/SRX13191403[accn]
  - GSM5699198: Free-living cells with alginate rep 3; Zobellia galactanivorans; RNA-Seq   
https://www.ncbi.nlm.nih.gov/sra/SRX13191404[accn]

------
### Setup
 ```R
# Manual file setup
gff_file <- "C:/R/Zob/Zob.gff"
# You can download reference genome.
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=63186
bamFls <- c("C:/R/bam/1176.bam",
            "C:/R/bam/1179.bam",
            "C:/R/bam/1180.bam",
            "C:/R/bam/1182.bam",
            "C:/R/bam/1183.bam",
            "C:/R/bam/1186.bam")
# FILE LINK : https://zenodo.org/records/11187387
root_dir <- 'C:/R/main'
project_name <- '2024_spr'

# Ensure all necessary packages are installed and loaded
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("S4Vectors", "Rsamtools", "edgeR", "GenomicRanges",
                       "GenomicAlignments", "rtracklayer", "DelayedArray", "SummarizedExperiment"))

install.packages(c("lattice", "pheatmap", "ggplot2", "ggExtra", "Matrix"))

# Load libraries
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(edgeR)
library(GenomicRanges)
library(GenomicAlignments)
library(pheatmap)
library(ggplot2)
library(ggExtra)

# Set working directory
setwd(root_dir)  # Set the working directory to root_dir
```
### Read genes & RPKM & MDS
 ```R
# Read genome feature file
gff <- import.gff(gff_file)  # Import the GFF file into gff object

# Extract only gene features
gffgene <- gff[elementMetadata(gff)[,"type"] == "gene"]  # Extract features of type 'gene' into gffgene

# Extract gene IDs
gid <- elementMetadata(gffgene)[,"ID"]  # Store gene IDs in gid

# Determine the locus tag field
if ("locus_tag" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"locus_tag"]  # If locus_tag exists, store in glocus_tag
} else if ("proteinId" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"proteinId"]  # If proteinId exists, store in glocus_tag
}

# Extract IDs and products for genes and CDS
id_parent <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "gene", c("ID", "Parent")]  # Extract ID and Parent for features of type 'gene'
id_product <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "CDS", c("Parent", "product")]  # Extract Parent and product for features of type 'CDS'
id_product <- unique(id_product)  # Remove duplicates
id_product$product <- sapply(id_product$product, paste, collapse=',')  # Collapse product field into a comma-separated string

# Considering RNA genes (rRNA/tRNA)
gffrgene <- gff[elementMetadata(gff)[,"type"] %in% c("rRNA", "tRNA", "tmRNA", "transcript")]  # Extract features of types rRNA, tRNA, tmRNA, and transcript

# RNA gene ID & index in total gid
rgid <- unlist(elementMetadata(gffrgene)[,"Parent"])  # Store Parent IDs of RNA genes in rgid
rna.ind <- match(rgid, gid)  # Match rgid in gid and store indices

# Initialize matrices for counts
all.counts <- matrix(0, nrow=length(gid), ncol=length(bamFls))  # Initialize matrix to store counts for all genes
rownames(all.counts) <- glocus_tag  # Set row names to glocus_tag
all.rcounts <- matrix(0, nrow=length(rgid), ncol=length(bamFls))  # Initialize matrix to store counts for RNA genes
colnames(all.counts) <- sub(".bam", "", basename(bamFls))  # Set column names to BAM file names without extension
colnames(all.rcounts) <- sub(".bam", "", basename(bamFls))  # Set column names to BAM file names without extension

# Read alignments and count overlaps
for (i in seq_along(bamFls)) {
  aln <- readGAlignments(bamFls[i])  # Read alignments from BAM file
  seqlevels(aln) <- seqlevels(gff)  # Set sequence levels of alignments to match gff
  counts <- countOverlaps(gffgene, aln, ignore.strand=TRUE)  # Count overlaps for gene features
  counts.r <- countOverlaps(gffrgene, aln, ignore.strand=TRUE)  # Count overlaps for RNA gene features
  all.counts[, i] <- counts  # Store counts in all.counts matrix
  all.rcounts[, i] <- counts.r  # Store RNA gene counts in all.rcounts matrix
}

# Coverage table preparation
coverage_table <- data.frame(
  name=sub(".bam", "", basename(bamFls)),  # Set file names
  all_counts=colSums(all.counts),  # Total counts for all genes
  ncrna_counts=colSums(all.rcounts),  # Total counts for ncRNA genes
  ncrna_ratio=round(100 * colSums(all.rcounts) / colSums(all.counts), 2)  # Calculate ncRNA ratio
)

# Write the coverage table
write.table(coverage_table, file.path(root_dir, paste(project_name, '_cov_table.txt', sep='')), sep='\t', row.names=FALSE, quote=FALSE)  # Save coverage table to file

# Remove rRNA genes counts for DEG analysis
all.gcounts <- if (length(rna.ind) > 0) all.counts[-rna.ind, ] else all.counts  # Remove counts for rRNA genes for DEG analysis

# Create DGEList object
dge <- DGEList(counts = all.gcounts)  # Create DGEList object
dge <- calcNormFactors(dge)  # Calculate normalization factors

# Save gene counts matrix
save(all.gcounts, file=file.path(root_dir, "all.gcounts.RData"))  # Save gene counts matrix to file

# Sample grouping
sample_names <- sub(".bam", "", basename(bamFls))  # Set sample names
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")  # Set conditions
design <- model.matrix(~0 + factor(conditions))  # Create design matrix
colnames(design) <- levels(factor(conditions))  # Set column names to conditions

# Estimate dispersion
dge <- estimateDisp(dge, design)  # Estimate dispersion

# Create model and perform DE analysis
fit <- glmFit(dge, design)  # Create and fit generalized linear model

# RPKM calculation - gene set subtraction
gsize <- width(gffgene)  # Calculate gene lengths
names(gsize) <- glocus_tag  # Assign names to gene lengths
gsize.c <- gsize[rownames(dge$counts)]  # Match gene lengths to DGEList
geneLengthsInKB <- (gsize.c / 1000)  # Convert gene lengths to kilobases

# Normalization using pseudo coefficients
dge <- estimateCommonDisp(dge, design) # Estimate common dispersion
millionsMapped <- colSums(dge$pseudo.counts) / 1e+06 # Convert mapped reads to millions
rpm <- dge$pseudo.counts / millionsMapped # Count reads per million
tt.rpkm <- round(rpm / geneLengthsInKB, 1) # Calculate RPKM values
colnames(tt.rpkm) <- sub('\\$', '_rpkm', colnames(tt.rpkm)) # Set column names
keep <- rowSums(tt.rpkm >= 2) >= 1 # Keep genes with RPKM >= 2 in at least one sample

# Initialize matrix for result table
tt.raw.counts <- dge$counts # Raw counts
colnames(tt.raw.counts) <- sub('\\$', '_raw_counts', colnames(tt.raw.counts)) # Set column names

# Pseudo count
dge.commondisp <- estimateCommonDisp(dge, design) # Estimate common dispersion
tt.pseudocount <- dge.commondisp$pseudo.counts # Calculate pseudo counts
tt.pseudocount.round <- round(tt.pseudocount) # Round pseudo counts
colnames(tt.pseudocount.round) <- sub('\\$', '_pseudocount', colnames(tt.pseudocount.round)) # Set column names

# Coefficients from fitting object (log 2 values)
tt.coeff <- fit$coefficients # Extract coefficients
tt.coeff.log2 <- tt.coeff / log(2) # Convert to log2 values
tt.coeff.log2.round <- round(tt.coeff.log2, 2) # Round log2 values
colnames(tt.coeff.log2.round) <- sub('\\$', '_coeff', colnames(tt.coeff.log2.round)) # Set column names

# Initialize result table
result.table <- cbind(tt.raw.counts, tt.pseudocount.round, tt.coeff.log2.round, tt.rpkm)

# Set up differential expression analysis
all_comb <- list()
comb_name <- list()
k <- 1
for (i in seq_along(levels(factor(conditions)))) {
  for (j in (i + 1):length(levels(factor(conditions)))) {
    if (i != j) {
      contrast_name <- paste(levels(factor(conditions))[i], "-", levels(factor(conditions))[j])
      all_comb[[k]] <- makeContrasts(contrasts = contrast_name, levels = design) # Create contrast matrix
      comb_name[[k]] <- paste(levels(factor(conditions))[i], "vs", levels(factor(conditions))[j]) # Set contrast name
      k <- k + 1
    }
  }
}

# Check correct contrast settings
print(all_comb)

# Initialize result table list
result.tables <- list()
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]] # Get comparison name
  contrast <- all_comb[[k]] # Get contrast matrix
  
  if (any(is.na(contrast))) next # Skip invalid contrast
  
  lrt <- glmLRT(fit, contrast = contrast) # Perform LRT
  res <- as.data.frame(topTags(lrt, n = nrow(dge))) # Extract results
  res$contrast <- comparison # Add contrast name to result
  
  # Check for consistent column names
  colnames(res) <- make.names(colnames(res), unique = TRUE)
  
  result.tables[[k]] <- res # Save to result list
}

# Check column names before combining
for (i in seq_along(result.tables)) {
  print(colnames(result.tables[[i]]))
}

# Combine all results into one data frame
combined_results <- do.call(rbind, result.tables) # Combine all results into one data frame

# Save combined results
write.table(combined_results, file = file.path(root_dir, paste(project_name, "_combined_DE_results.txt", sep = "")), sep = "\t", row.names = FALSE, quote = FALSE) # Save the combined results

# MDS plot for sample quality
png(file = file.path(root_dir, paste(project_name, '_MDS_all_samples.png', sep = '')), width = 2000, height = 2000, res = 300) # Save MDS plot as PNG file
plotMDS(dge, col = as.numeric(factor(conditions)), pch = 16, cex = 1.5, main = "MDS Plot of all samples") # Create MDS plot
dev.off() # Close PNG device
 ```

### Heatmap & Top 20 DE genes
 ```R
# Heatmap
if (exists("tt.rpkm")) {
  rpkm.avg <- round(rowMeans(tt.rpkm), 2) # Calculate average RPKM values
  
  if (length(unique(conditions)) > 1) {
    pheatmap(mat = tt.rpkm, show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, scale = 'row') # Create heatmap
  } else {
    cat("At least two groups are required for heatmap visualization.\n") # Print message if at least two groups are required
  }
} else {
  cat("Object 'tt.rpkm' not found.\n") # Print message if tt.rpkm object does not exist
}

# Extract top 20 genes from differential expression analysis
sigOE_ordered <- result.tables[[1]][order(result.tables[[1]]$PValue), ] # Sort by p-value
top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ]) # Extract top 20 genes

# Extract normalized count data using cpm (counts per million)
normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE) # Calculate normalized counts

# Extract normalized counts of top 20 significantly differentially expressed genes
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]

# Convert to data frame
melted_top20_sigOE <- melt(top20_sigOE_norm) # Convert to data frame

# Change column names
colnames(melted_top20_sigOE) <- c("gene", "samplename", "normalized_counts") # Set column names

# Add metadata
meta <- data.frame(sampletype = conditions, samplename = colnames(normalized_counts)) # Create metadata
melted_top20_sigOE <- merge(melted_top20_sigOE, meta, by = "samplename") # Merge metadata

# Draw plot
ggplot(melted_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) + # Create scatter plot
  scale_y_log10() + # Set y-axis to log scale
  xlab("Genes") + # Set x-axis label
  ylab("Normalized Counts") + # Set y-axis label
  ggtitle("Top 20 Significant DE Genes") + # Set plot title
  theme_bw() + # Apply theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis text
  theme(plot.title = element_text(hjust = 0.5)) # Center align title
```

### Rplot
 ```R
# Specify conditions for group comparison
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")
group1 <- which(conditions == "Condition1") # Designate group 1 as Condition1
group2 <- which(conditions == "Condition2") # Designate group 2 as Condition2

# Set layout to display four plots
par(mfrow = c(2, 2)) # Set layout to display four plots

# Plot average difference for all samples
plotMD(dge, column = 1:ncol(dge$counts), main = "Total") # Generate mean difference plot for all samples
abline(h = 0, col = "grey") # Add a grey horizontal line representing the median

# Mean difference plot for group 1
plotMD(dge, column = group1, main = "Group 1: Condition1") # Create mean difference plot for group 1
abline(h = 0, col = "blue") # Add a blue horizontal line representing the median

# Mean difference plot for group 2
plotMD(dge, column = group2, main = "Group 2: Condition2") # Create mean difference plot for group 2
abline(h = 0, col = "red") # Add a red horizontal line representing the median

# Reset to default layout
par(mfrow = c(1, 1)) # Reset to default layout
```

------
### Result
 - mds
![2024_spr_MDS_all_samples](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165700031/b22ddc60-f9b5-4628-9584-c2ec3317aef2)
 - heatmap
![Heatmap](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165700031/5b818b78-e235-47a9-a431-51c84c9f7912)
 - Top 20 DE Genes
![Top 20 DE Genes](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165700031/dff00aa9-69d2-4628-9781-2c864691083b)
 - Rplot
  ![Rplot](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165700031/cfed80fd-02be-4b94-942b-7ba47e26ec7d)


------
### Links
1. [RU_RNASeq Workflow](https://rockefelleruniversity.github.io/RU_RNAseq/)
2. [Bioconductor Courses & Workshops](https://www.bioconductor.org/help/course-materials/)
3. https://github.com/alfonsosaera/RNAseq
4. https://www.kola.kr/course/course_view.jsp?id=12434&page=7#cview2
5. https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Quality_control

