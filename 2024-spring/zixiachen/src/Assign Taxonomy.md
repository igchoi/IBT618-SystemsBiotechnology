## Assign Taxonomy
------
### 1. Install and load dada2
```r
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if(!requireNamespace("BiocManager",quietly = TRUE))
install.packages("BiocManager",repo=site)  
library(BiocManager)
BiocManager::install("dada2",version = "3.18")
library(dada2)
packageVersion("dada2")
```
------
### 2. Creation of a path to link the folder containing the fastq files
```r
path <- "/Users/livia/Desktop/kombucha/16s/Rawdata"
list.files(path)
```
------
### 3. Rename fastq files
```r
fnFs <- sort(list.files(path,pattern="_R1.fq",full.names = TRUE))
fnRs <- sort(list.files(path,pattern="_R2.fq",full.names = TRUE))
```
------
### 4. Extract sample names
```r
sample.names <- sapply(strsplit(basename(fnFs),"_"),'[',1)
sample.names
```
------
### 5. Reads quality check
```r
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```
------
### 6. Filter and trimm
```r
filtFs <- file.path(path,"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs <- file.path(path,"filtered",paste0(sample.names,"_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs,truncLen = c(240,160),
                     maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,
                     compress=TRUE, multithread=T)
head(out)
```
------
### 7. Check error rates
```r
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```
------
### 8. Sample inference
```r
dadaFs <- dada(filtFs,err=errF,multithread = TRUE)
dadaRs <- dada(filtRs,err=errR,multithread = TRUE)
dadaFs[[1]]
dadaRs[[1]]
```
------
### 9. Merge paired reads
```r
mergers <- mergePairs(dadaFs,filtFs,dadaRs,filtRs,verbose=TRUE)
head(mergers[[1]])
```
------
### 10. Construct sequence table
```r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
```
------
### 11. Remove chimeras
```r
seqtab.nochim <- removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE,verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```
------
### 12. Export the generated ASV table
```r
write.csv(seqtab.nochim,file="/Users/livia/Desktop/kombucha/16s/Rawdata/filtered/ASV.CSV",
          append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", 
          dec = ".", row.names = TRUE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
```
------
# 13. Assign taxonomy
```r
taxa <- assignTaxonomy(seqtab.nochim, paste0(path, "/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, paste0(path, "/silva_species_assignment_v132.fa.gz"))
write.csv(taxa,file="/Users/livia/Desktop/kombucha/16s/Rawdata/taxa.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", 
          na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print,file="/Users/livia/Desktop/kombucha/16s/Rawdata/taxa.print.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", 
          na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
```
------
