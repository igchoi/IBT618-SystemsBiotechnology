## DADA2+Pacbio
[Using Daqu data from NCBI](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP508091)

[Analysed using DADA2+PacBio](https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html)

------
```r
if(!requireNamespace("BiocManager",quietly = TRUE))
  install.packages("BiocManager")  
library(BiocManager)
BiocManager::install("dada2",version = "3.19")
library(dada2)
library(Biostrings)
library(ShortRead)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(phyloseq)
path1 <- "/Users/livia/Desktop/Study on the contribution of Daqu and environment to the microecology of fermented grains of Baijiu/Test"
path.out <- "Figures/"
path.rds <- "/Users/livia/Desktop/Study on the contribution of Daqu and environment to the microecology of fermented grains of Baijiu/Test/RDS"
fns1 <- list.files(path1, pattern="fastq.gz.1", full.names=TRUE)
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())
nops1 <- file.path(path1, "noprimers", basename(fns1))
prim1 <- removePrimers(fns1, nops1, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)
filts1 <- file.path(path1, "noprimers", "filtered", basename(fns1))
track1 <- filterAndTrim(nops1, filts1, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)
track1
drp1 <- derepFastq(filts1, verbose=TRUE)
err1 <- learnErrors(drp1, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
saveRDS(err1, file.path(path.rds, "Daqu_err1.rds"))
plotErrors(err1)
dd1 <- dada(drp1, err=err1, BAND_SIZE=32, multithread=TRUE)
saveRDS(dd1, file.path(path.rds, "Daqu_dd1.rds"))
cbind(ccs=prim1[,1], primers=prim1[,2], filtered=track1[,2], denoised=sapply(dd1, function(x) sum(x$denoised)))
st1 <- makeSequenceTable(dd1); dim(st1)
bim1 <- isBimeraDenovo(st1, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim1)
sum(st1[,bim1])/sum(st1)
write.csv(st1,file="/Users/livia/Desktop/Study on the contribution of Daqu and environment to the microecology of fermented grains of Baijiu/Test/ASV.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
taxa <- assignTaxonomy(st1, paste0(path1, "/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, paste0(path1, "/silva_species_assignment_v132.fa.gz"))
write.csv(taxa,file="/Users/livia/Desktop/Study on the contribution of Daqu and environment to the microecology of fermented grains of Baijiu/Test/taxa.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", 
          na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print,file="/Users/livia/Desktop/Study on the contribution of Daqu and environment to the microecology of fermented grains of Baijiu/Test/taxa.print.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", 
          na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
```
