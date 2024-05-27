# Manual file setup
gff_file <- "C:/R/Zob/Zob.gff"
bamFls <- c("C:/R/bam/1176.bam",
            "C:/R/bam/1179.bam",
            "C:/R/bam/1180.bam",
            "C:/R/bam/1182.bam",
            "C:/R/bam/1183.bam",
            "C:/R/bam/1186.bam")

root_dir <- 'C:/R/main'
project_name <- '2024_spr'

# Ensure all necessary packages are installed and loaded
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("S4Vectors", "Rsamtools", "edgeR", "GenomicRanges", "GenomicAlignments", "rtracklayer", "DelayedArray", "SummarizedExperiment"))

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
setwd(root_dir)

# Read genome feature file
gff <- import.gff(gff_file)

# Extract only gene features
# The Type part can be changed to mRNA, etc.
gffgene <- gff[elementMetadata(gff)[,"type"] == "gene"]

# Extract gene IDs
gid <- elementMetadata(gffgene)[,"ID"]

# Determine the locus tag field
#as.data.frame(gffgene) converts this object to a dataframe. 
#names(as.data.frame(gffgene)) gets the column names (i.e., metadata field names) of the dataframe.
#The first if conditional checks to see if the locus_tag field exists in the metadata. 
#If it does, elementMetadata(gffgene)[,"locus_tag"] extracts the locus_tag field from gffgene's metadata and stores it in the glocus_tag variable.
if ("locus_tag" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"locus_tag"]
} else if ("proteinId" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"proteinId"]
}
#You can add a variety of field values. For example
#else if ("geneID" %in% names(as.data.frame(gffgene))) 
#{glocus_tag <- elementMetadata(gffgene)[,"geneID"]}

# Extract IDs and products for genes and CDS
id_parent <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "gene", c("ID", "Parent")]
id_product <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "CDS", c("Parent", "product")]
id_product <- unique(id_product)
id_product$product <- sapply(id_product$product, paste, collapse=',')
#The sapply function applies the paste function to each element of id_product$product. 
#The paste function combines all the elements of the vector into one string by default. 
#collapse=',' specifies that each element should be separated by a comma to combine them.


# Considering RNA genes (rRNA/tRNA)
gffrgene <- gff[elementMetadata(gff)[,"type"] %in% c("rRNA", "tRNA", "tmRNA", "transcript")]
#Type, filter out those belonging to %in%, and store the filtered items in gffregene.

# RNA gid & index in total gid
rgid <- unlist(elementMetadata(gffrgene)[,"Parent"])
rna.ind <- match(rgid, gid)

# Initialize matrices for counts
all.counts <- matrix(0, nrow = length(gid), ncol = length(bamFls))
#matrix(0, nrow = length(gid), ncol = length(bamFls)): Create a matrix with all values initialized to zero.
#nrow = length(gid): The number of rows is equal to the length of the gene ID (gid).
#ncol = length(bamFls): The number of columns is equal to the number in the BAM file (bamFls).
rownames(all.counts) <- glocus_tag
all.rcounts <- matrix(0, nrow = length(rgid), ncol = length(bamFls))
colnames(all.counts) <- sub(".bam", "", basename(bamFls))
colnames(all.rcounts) <- sub(".bam", "", basename(bamFls))

# Read alignments and count overlaps
for (i in seq_along(bamFls)) {
  aln <- readGAlignments(bamFls[i])
  seqlevels(aln) <- seqlevels(gff)
#seqlevels(gff): Extract the sequence levels of the GFF object. This is necessary information when calculating overlap with genes and RNA genes.  
  counts <- countOverlaps(gffgene, aln, ignore.strand = TRUE)
  counts.r <- countOverlaps(gffrgene, aln, ignore.strand = TRUE)
#countOverlaps(gffgene, aln, ignore.strand = TRUE): Counts the overlaps between genes and reads. gffgene is a GFF object representing a gene feature.
#countOverlaps(gffrgene, aln, ignore.strand = TRUE): Counts overlaps between RNA genes and reads. gffrgene is a GFF object representing RNA gene features.  
  all.counts[, i] <- counts
  all.rcounts[, i] <- counts.r
}

# Coverage table preparation
coverage_table <- data.frame(
  name = sub(".bam", "", basename(bamFls)),
  all_counts = colSums(all.counts),
  ncrna_counts = colSums(all.rcounts),
  ncrna_ratio = round(100 * colSums(all.rcounts) / colSums(all.counts), 2)
)
#all_counts = colSums(all.counts): sum the number of reads for genes in each column of the all.counts matrix, to calculate the total number of reads for all genes in each BAM file. 
#ncrna_counts = colSums(all.rcounts): sum the number of reads for RNA genes in each column of the all.rcounts matrix, to calculate the total number of reads for all RNA genes in each BAM file.
#ncrna_ratio = round(100 * colSums(all.rcounts) / colSums(all.counts), 2): Calculate the ratio of the number of RNA gene reads. 
#To do this, we divide the number of RNA gene reads in each BAM file by the total number of reads in that BAM file and express it as a percentage. 
#We use the round function to round to two decimal places.

# Write the coverage table
write.table(coverage_table, file.path(root_dir, paste(project_name, '_cov_table.txt', sep='')), sep='\t', row.names = FALSE, quote = FALSE)

# Remove rRNA genes counts for DEG analysis
if (length(rna.ind) > 0) {
  all.gcounts <- all.counts[-rna.ind, ]
} else {
  all.gcounts <- all.counts
}
#If `length(rna.ind) > 0`: If `rna.ind` has the index of an rRNA gene, that is, it contains an rRNA gene 
#`all.gcounts <- all.counts[-rna.ind, ]`: Create a new matrix `all.gcounts` by removing the rows corresponding to rRNA genes. Here `-rna.ind` selects the remaining rows except the index corresponding to the rRNA gene - `else`: If there are no rRNA genes, 
#`all.gcounts` has the same value as `all.counts`.

# Save gene counts matrix
save(all.gcounts, file = file.path(root_dir, "all.gcounts.RData"))

# Sample grouping
sample_names <- sub(".bam", "", basename(bamFls))
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")
design <- model.matrix(~0 + factor(conditions))
colnames(design) <- levels(factor(conditions))

# EdgeR analysis
dge <- DGEList(counts = all.gcounts, group = conditions)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# Fit model and compute contrasts
fit <- glmFit(dge, design)

# Define contrasts for differential expression analysis
all_comb <- list()  # Initialising the list to store comparisons
comb_name <- list()  # Initialising a list to store comparison names
k <- 1  # Initialising the comparison index
for (i in seq_along(levels(factor(conditions)))) {
  for (j in (i+1):length(levels(factor(conditions)))) {
    # Generate comparisons for all combinations
    all_comb[[k]] <- makeContrasts(contrasts = levels(factor(conditions))[c(i, j)], levels = design)
    # Create a comparison name and add it to the list
    comb_name[[k]] <- paste(levels(factor(conditions))[c(i, j)], collapse = "vs")
    k <- k + 1  # Increase comparison index
  }
}
print(all_comb[[2]])

# Example of correcting a problematic contrast vector manually
all_comb[[2]] <- c(1, -1, 0)  # Adjust according to your specific needs

all_comb <- all_comb[-2]  # Removes the second element from all_comb
comb_name <- comb_name[-2]  # Do the same for comb_name if it's used in correspondence with all_comb


result.tables <- list()  # Initialising a list to store results
for (k in seq_along(all_comb)) {
  # Extracting comparison names
  comparison <- comb_name[[k]]
  # Running DEG analysis using LRT
  lrt <- glmLRT(fit, contrast = all_comb[[k]])
  # Extracting top gene tags
  topTags <- topTags(lrt, n = Inf)
  # Add the results table to the results list
  result.tables[[comparison]] <- topTags$table
}

for (k in seq_along(all_comb)) {
  contrast <- all_comb[[k]]
  # Check if the comparison has NA, NaN, or Inf
  if (any(is.na(contrast)) || any(is.nan(contrast)) || any(is.infinite(contrast))) {
    cat("Skipping contrast for all_comb[[", k, "]] because it contains NA, NaN, or Inf\n")
  } else {
    tryCatch({
      # Run DEG analysis
      lrt <- glmLRT(fit, contrast = contrast)
      topTags <- topTags(lrt, n = Inf)
      # Add the results table to the results list
      result.tables[[comb_name[[k]]]] <- topTags$table
    }, error = function(e) {
      # Outputting messages when errors occur
      cat("Error in processing contrast for all_comb[[", k, "]]: ", e$message, "\n")
    })
  }
}

# Check for NA, NaN, or Inf for all elements in the model
for (i in seq_along(fit)) {
  element <- fit[[i]]
  if (is.list(element)) {
    for (j in seq_along(element)) {
      sub_element <- element[[j]]
      # Check for NA, NaN, or Inf in child elements
      if (any(is.na(sub_element)) || any(is.nan(sub_element)) || any(is.infinite(sub_element))) {
        cat("NA, NaN, or Inf found in sub-element", j, "of element", i, "of the list.\n")
      }
    }
  } else {
    # Check for NA, NaN, or Inf in elements
    if (any(is.na(element)) || any(is.nan(element)) || any(is.infinite(element))) {
      cat("NA, NaN, or Inf found in element", i, "of the list.\n")
    }
  }
}



# Analyze differential expression for each comparison
result.tables <- list()  # Initialize list to store results
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]  # Extract comparison name
  lrt <- glmLRT(fit, contrast = all_comb[[k]])  # Perform DEG analysis using LRT
  topTags <- topTags(lrt, n = Inf)  # Extract top gene tags
  result.tables[[comparison]] <- topTags$table  # Add result table to result list
}


# Save DEG results
for (comparison in names(result.tables)) {
  # Save the DEG results for each comparison to a file.
  write.table(result.tables[[comparison]], 
              file.path(root_dir, paste(project_name, '_', comparison, '_DE_results.txt', sep = '')), 
              sep = '\t', quote = FALSE, row.names = TRUE)
}

# MDS plot for sample quality
png(file = file.path(root_dir, paste(project_name, '_MDS_all_samples.png', sep = '')), 
    width = 2000, height = 2000, res = 300)  
plotMDS(dge, col = as.numeric(factor(conditions)), pch = 16, cex = 1.5, main = "MDS Plot of all samples")  # 샘플 품질을 시각화합니다.
dev.off() 


# Raw count
tt.raw.counts <- dge$counts
colnames(tt.raw.counts) <- sub('\\$', '_raw_counts', colnames(tt.raw.counts))

# Pseudocount
dge.commondisp <- estimateCommonDisp(dge, design)  # Estimate the common variance of each gene.
tt.pseudocount <- dge.commondisp$pseudo.counts  # Extract the pseudocount for each gene.
tt.pseudocount.round <- round(tt.pseudocount)  # Rounds the pseudo-count and converts it to an integer.
colnames(tt.pseudocount.round) <- sub('\\$', '_pseudocount', colnames(tt.pseudocount.round))

# Coefficient from fit object (note that it is log2 value)
tt.coeff <- fit$coefficients  # Extract the coefficients from the fit object.
tt.coeff.log2 <- tt.coeff / log(2)  # Transform the coefficient based on log2 to make it a log2 value.
tt.coeff.log2.round <- round(tt.coeff.log2, 2)  # The log2 value is rounded to two decimal places.
colnames(tt.coeff.log2.round) <- sub('\\$', '_coeff', colnames(tt.coeff.log2.round))  

# Get each gene's product
row_tt <- rownames(tt.raw.counts)  # Gets the row names of tt.raw.counts.
tmp_index <- match(row_tt, id_locus_tag$locus_tag)  # Match the locus_tag of each row in tt.raw.counts with the locus_tag in the id_locus_tag dataframe.
tt_gene_id <- id_locus_tag$ID[tmp_index]  # Get the gene ID corresponding to the matched locus_tag.
tmp_index2 <- match(tt_gene_id, unlist(id_parent$ID))  # Match the IDs in the tt_gene_id and id_parent dataframes.
tt_rna_id <- id_parent$ID[tmp_index2]  # Get the ID of the RNA corresponding to the matched gene ID.
tmp_index3 <- match(tt_rna_id, unlist(id_product$Parent))  # Match the tt_rna_id with the Parent in the id_product dataframe.
if (any(!is.na(tmp_index3))) {
  tt_product <- unlist(id_product$product[tmp_index3])  # Get the product name corresponding to the ID of the matched RNA.
} else {
  cat('Something wrong with product name\n')  # Output an error message if there is a problem with the product name.
}


# RPKM calculation - gene set subtraction
gsize <- width(gffgene)  # Get the length of each gene from gffgene.
names(gsize) <- glocus_tag  # Assign a gene locus tag to the length of each gene.
gsize.c <- gsize[rownames(dge$counts)]  # Select only the lengths of the genes in dge$counts.
geneLengthsInKB <- (gsize.c / 1000)  # Converts gene length to kilobases.

dge <- estimateCommonDisp(dge, design)  # Estimate the common variance using the DGEList object.
millionsMapped <- colSums(dge$pseudo.counts) / 1e+06  # Sum the pseudo-count values for each sample and convert them to millions.
rpm <- dge$pseudo.counts / millionsMapped  # Calculate the RPM value by dividing the pseudo-count value by a million.
tt.rpkm <- round(rpm / geneLengthsInKB, 1)  # Calculate and round the RPKM value by dividing the RPM value by the gene length.
colnames(tt.rpkm) <- sub('\\$', '_rpkm', colnames(tt.rpkm))  # Append '_rpkm' to the column name of the RPKM value to display it.

# Combine the raw count data, product information, pseudo counts, count log2 values, and RPKM values in a results table.
result.table <- cbind(tt.raw.counts, tt_product, tt.pseudocount.round, tt.coeff.log2.round, tt.rpkm)

# Perform the glmLRT test for each comparison, analyse the results, and add the results to the results table.
for (i in 1:length(all_comb)) {  
  lrTest <- glmLRT(fit, contrast = all_comb[[i]])  # Use the glmLRT function to perform the test for each comparison.
  tt <- topTags(lrTest, Inf)  # Extract the top statistics.
  tt.tab <- tt$table  # Extract the resulting table.
  if (!is.null(dim(tt.tab))) {  # Performs only if a dimension in the resulting table is non-NULL.
    tt.tab[,'logFC'] <- round(tt.tab[,'logFC'], 2)  # Rounds the ogFC column.
    tt.tab[,'logCPM'] <- round(tt.tab[,'logCPM'], 2)  # Rounds the ogCPM column.
    tt.tab[,'LR'] <- round(tt.tab[,'LR'], 2)  # Rounds the LR columns.
    qcomb_name <- gsub(" ","vs",gsub("1[*]","",gsub("-1[*]","",lrTest$comparison)))  # Refine the name of the comparison.
    idx.coeff <- match(rownames(fit$coefficients), rownames(tt.tab))  # Match the index of a coefficient.
    tt.tab.sorted <- tt.tab[idx.coeff,]  # Sort the resulting table.
    colnames(tt.tab.sorted) <- sub('\\$', paste('_', qcomb_name, sep=''), colnames(tt.tab.sorted))  # Add '_qcomb_name' to the column name.
    result.table <- cbind(result.table, tt.tab.sorted)  # Add a sorted results table to the Results table.
  }
}


# Write result table to a file
# Define the file path for the DEG table
dge_file <- file.path(root_dir, paste(project_name, '_deg_table_20210405.txt', sep=''))

# Write the result table to a tab-separated text file
# Set row names to be included, without quoting them
write.table(result.table, dge_file, sep='\t', row.names=T, quote=F)

# Heatmap visualization
# Calculate the average RPKM values across all samples
rpkm.avg <- round(rowMeans(tt.rpkm), 2)

# Check if there are at least two groups for heatmap visualization
if (length(unique(grp)) > 1) {
  # Plot the heatmap using pheatmap function
  pheatmap(mat = tt.rpkm, show_rownames = F, show_colnames = T, cluster_rows = T, cluster_cols = T, scale = 'row')
} else {
  # Print a message indicating that at least two groups are required for heatmap visualization
  cat("At least two groups are required for heatmap visualization.\n")
}
