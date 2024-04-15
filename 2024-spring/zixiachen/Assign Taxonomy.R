# 1. Install and load dada2
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repo=site)
library(BiocManager)
BiocManager::install("dada2",version = 3.18)
library(dada2)
install.packages("Rcpp")
library(Rcpp)
packageVersion("dada2")

# 2. Creation of a path to link the folder containing the fastq files
path <- "/Users/livia/Desktop/Daqu sample/R analysis/Rawdata" 
list.files(path)

# 3. Rename fastq files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq",full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq",full.names = TRUE))

# 4. Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 5. Reads quality check
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# 6. Filter and trimm
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names,"_F_filt.fastq.gz"))
# Remove primers sequences (trimLeft), remove reads with low quality score
out <- filterAndTrim(fnFs,filtFs,truncLen = 150, trimLeft = 17,
                     maxEE = 2,maxN=2,truncQ = 2,rm.phix = TRUE,compress=TRUE, multithread = T)
head(out)

# 7. Check error rates
errF <- learnErrors(filtFs,multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# 8. Sample inference
# Applying the core sample inference algorithm to filtered-trimmed data
dadaFs <- dada(filtFs,err=errF,multithread = TRUE)
dadaFs[[1]]

# 9. Construct sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# 10. Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE,verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# 11. Export the generated ASV table
write.csv(seqtab.nochim,file="/Users/livia/Desktop/Daqu sample/R analysis/Rawdata/filtered/ASV.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

# 12. Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, paste0(path, "/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, paste0(path, "/silva_species_assignment_v132.fa.gz"))
write.csv(taxa,file="/Users/livia/Desktop/Daqu sample/R analysis/Rawdata/taxa.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print,file="/Users/livia/Desktop/Daqu sample/R analysis/Rawdata/taxa.print.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
