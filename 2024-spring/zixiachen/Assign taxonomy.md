# Assign taxonomy
``` r
path <- "/Users/livia/Desktop/Daqu sample/R analysis/Rawdata"
taxa <- assignTaxonomy(seqtab.nochim, paste0(path, "/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, paste0(path, "/silva_species_assignment_v132.fa.gz"))
write.csv(taxa,file="/Users/livia/Desktop/Daqu sample/R analysis/Rawdata/taxa.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print,file="/Users/livia/Desktop/Daqu sample/R analysis/Rawdata/taxa.print.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
 ```
