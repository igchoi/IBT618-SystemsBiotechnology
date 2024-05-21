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
setwd(root_dir)

# Read genome feature file
gff <- import.gff(gff_file)

# Extract only gene features
gffgene <- gff[elementMetadata(gff)[,"type"] == "gene"]

# Extract gene IDs
gid <- elementMetadata(gffgene)[,"ID"]

# Determine the locus tag field
if ("locus_tag" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"locus_tag"]
} else if ("proteinId" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"proteinId"]
}

# Extract IDs and products for genes and CDS
id_parent <- elementMetadata(gff)[elementMetadata(gff)[,"type"] 
                                  == "gene", c("ID", "Parent")]
id_product <- elementMetadata(gff)[elementMetadata(gff)[,"type"] 
                                   == "CDS", c("Parent", "product")]
id_product <- unique(id_product)
id_product$product <- sapply(id_product$product, paste, collapse=',')

# Considering RNA genes (rRNA/tRNA)
gffrgene <- gff[elementMetadata(gff)[,"type"] 
                %in% c("rRNA", "tRNA", "tmRNA", "transcript")]

# RNA gid & index in total gid
rgid <- unlist(elementMetadata(gffrgene)[,"Parent"])
rna.ind <- match(rgid, gid)

# Initialize matrices for counts
all.counts <- matrix(0, nrow = length(gid), ncol = length(bamFls))
rownames(all.counts) <- glocus_tag
all.rcounts <- matrix(0, nrow = length(rgid), ncol = length(bamFls))
colnames(all.counts) <- sub(".bam", "", basename(bamFls))
colnames(all.rcounts) <- sub(".bam", "", basename(bamFls))

# Read alignments and count overlaps
for (i in seq_along(bamFls)) {
  aln <- readGAlignments(bamFls[i])
  seqlevels(aln) <- seqlevels(gff)
  counts <- countOverlaps(gffgene, aln, ignore.strand = TRUE)
  counts.r <- countOverlaps(gffrgene, aln, ignore.strand = TRUE)
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

# Write the coverage table
write.table(coverage_table, file.path(root_dir, 
paste(project_name, '_cov_table.txt', sep='')), sep='\t', 
row.names = FALSE, quote = FALSE)

# Remove rRNA genes counts for DEG analysis
if (length(rna.ind) > 0) {
  all.gcounts <- all.counts[-rna.ind, ]
} else {
  all.gcounts <- all.counts
}

# Save gene counts matrix
save(all.gcounts, file = file.path(root_dir, "all.gcounts.RData"))

# Sample grouping
sample_names <- sub(".bam", "", basename(bamFls))
conditions <- c("Condition1", "Condition1", "Condition1", 
                "Condition2", "Condition2", "Condition2")
design <- model.matrix(~0 + factor(conditions))
colnames(design) <- levels(factor(conditions))

# Fit model and compute contrasts
fit <- glmFit(dge, design)

# Define contrasts for differential expression analysis
all_comb <- list()
comb_name <- list()
k <- 1
for (i in seq_along(levels(factor(conditions)))) {
  for (j in (i+1):length(levels(factor(conditions)))) {
    all_comb[[k]] <- makeContrasts(contrasts 
    = levels(factor(conditions))[c(i, j)], levels = design)
    comb_name[[k]] <- paste(levels(factor(conditions))[c(i, j)], 
                            collapse = "vs")
    k <- k + 1
  }
}

print(all_comb[[2]])
# Example of correcting a problematic contrast vector manually
all_comb[[2]] <- c(1, -1, 0)  # Adjust according to your specific needs

all_comb <- all_comb[-2]  # Removes the second element from all_comb
comb_name <- comb_name[-2]  
# Do the same for comb_name if it's used in correspondence with all_comb


result.tables <- list()
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]
  lrt <- glmLRT(fit, contrast = all_comb[[k]])
  topTags <- topTags(lrt, n = Inf)
  result.tables[[comparison]] <- topTags$table
}

for (k in seq_along(all_comb)) {
  contrast <- all_comb[[k]]
  if (any(is.na(contrast)) || any(is.nan(contrast)) 
      || any(is.infinite(contrast))) {
    cat("Skipping contrast for all_comb[[", k, "]]
        because it contains NA, NaN, or Inf\n")
  } else {
    tryCatch({
      lrt <- glmLRT(fit, contrast = contrast)
      topTags <- topTags(lrt, n = Inf)
      result.tables[[comb_name[[k]]]] <- topTags$table
    }, error = function(e) {
      cat("Error in processing contrast for all_comb[[", k, "]]: 
          ", e$message, "\n")
    })
  }
}

for (i in seq_along(fit)) {
  element <- fit[[i]]
  if (is.list(element)) {
    for (j in seq_along(element)) {
      sub_element <- element[[j]]
      if (any(is.na(sub_element)) || any(is.nan(sub_element)) 
          || any(is.infinite(sub_element))) {
        cat("NA, NaN, or Inf found in sub-element", 
            j, "of element", i, "of the list.\n")
      }
    }
  } else {
    if (any(is.na(element)) || any(is.nan(element)) 
        || any(is.infinite(element))) {
      cat("NA, NaN, or Inf found in element", i, "of the list.\n")
    }
  }
}


# Analyze differential expression for each comparison
result.tables <- list()
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]
  lrt <- glmLRT(fit, contrast = all_comb[[k]])
  topTags <- topTags(lrt, n = Inf)
  result.tables[[comparison]] <- topTags$table
}

# Save DEG results
for (comparison in names(result.tables)) {
  write.table(result.tables[[comparison]], 
              file.path(root_dir, paste(project_name, '_', comparison, 
  '_DE_results.txt', sep = '')), sep = '\t', quote = FALSE, row.names = TRUE)
}

# MDS plot for sample quality
grp = all.gcounts
png(file = file.path(root_dir, paste(project_name, 
  '_MDS_all_samples.png', sep = '')), width = 2000, height = 2000, res = 300)
plotMDS(dge, col = as.numeric(factor(conditions)), 
        pch = 16, cex = 1.5, main = "MDS Plot of all samples")
dev.off()

# Raw count
tt.raw.counts <- dge$counts
colnames(tt.raw.counts) <- sub('\\$', '_raw_counts', colnames(tt.raw.counts))

# Pseudocount
dge.commondisp <- estimateCommonDisp(dge, design)
tt.pseudocount <- dge.commondisp$pseudo.counts
tt.pseudocount.round <- round(tt.pseudocount)
colnames(tt.pseudocount.round) <- sub('\\$', '_pseudocount', 
                                      colnames(tt.pseudocount.round))

# Coefficient from fit object (note that it is log2 value)
tt.coeff <- fit$coefficients
tt.coeff.log2 <- tt.coeff / log(2)
tt.coeff.log2.round <- round(tt.coeff.log2, 2)
colnames(tt.coeff.log2.round) <- sub('\\$', '_coeff', 
                                     colnames(tt.coeff.log2.round))
# EdgeR analysis
dge <- DGEList(counts = all.gcounts, group = conditions)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# RPKM calculation - gene set subtraction
gsize <- width(gffgene)
names(gsize) <- glocus_tag
gsize.c <- gsize[rownames(dge$counts)]
geneLengthsInKB <- (gsize.c / 1000)
dge <- estimateCommonDisp(dge, design)
millionsMapped <- colSums(dge$pseudo.counts) / 1e+06
rpm <- dge$pseudo.counts / millionsMapped
tt.rpkm <- round(rpm / geneLengthsInKB, 1)
colnames(tt.rpkm) <- sub('\\$', '_rpkm', colnames(tt.rpkm))
result.table <- cbind(tt.raw.counts,tt.pseudocount.round, 
                      tt.coeff.log2.round, tt.rpkm)
keep <- rowSums(tt.rpkm >= 10) >= 1
result.table <- result.table[keep, ]
dim(result.table)

for (i in 1:length(all_comb)) {  
  lrTest <- glmLRT(fit, contrast = all_comb[[i]])
  tt <- topTags(lrTest, Inf)
  tt.tab <- tt$table
  if (!is.null(dim(tt.tab))) {
    tt.tab[,'logCPM'] <- round(tt.tab[,'logCPM'], 2)
    tt.tab[,'LR'] <- round(tt.tab[,'LR'], 2)
    qcomb_name <- gsub(" ","vs",gsub("1[*]",""
                                     ,gsub("-1[*]","",lrTest$comparison)))
    idx.coeff <- match(rownames(fit$coefficients), rownames(tt.tab))
    tt.tab.sorted <- tt.tab[idx.coeff,]
    colnames(tt.tab.sorted) <- sub('\\$', paste('_', qcomb_name, sep=''), 
                                   colnames(tt.tab.sorted))
    result.table <- cbind(result.table, tt.tab.sorted)
  }
}


# Write result table to a file
dge_file <- file.path(root_dir, paste(project_name, 
                                      '_deg_table_20210405.txt', sep=''))
write.table(result.table, dge_file, sep='\t', row.names=T, quote=F)

# Heatmap
library(pheatmap)
rpkm.avg <- round(rowMeans(tt.rpkm), 2)

if (length(unique(grp)) > 1) {
  pheatmap(mat = tt.rpkm, show_rownames = F, show_colnames = T, 
           cluster_rows = T, cluster_cols = T, scale = 'row')
} else {
  cat("At least two groups are required for heatmap visualization.\n")
}

