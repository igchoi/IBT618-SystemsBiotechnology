
### Manual file 
gff_file <- "C:/R/Zob/Zob.gff"
bamFls = c("C:/R/bam/1176.bam",
           "C:/R/bam/1179.bam",
           "C:/R/bam/1180.bam",
           "C:/R/bam/1182.bam",
           "C:/R/bam/1183.bam",
           "C:/R/bam/1186.bam")

root_dir = 'C:/R/main'

project_name = '2024_spr'

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("S4Vectors")
BiocManager::install("Rsamtools")

install.packages("lattice")
BiocManager::install("edgeR")

GenomicAlignments
install.packages("Matrix")
BiocManager::install("DelayedArray")
BiocManager::install("SummarizedExperiment")
BiocManager::install("rtracklayer")

# load library
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(edgeR)
library(GenomicRanges)
library(GenomicAlignments)
library(R.utils)


# Set working directory
setwd(root_dir)

#
# bam to counts (from alignment to gene counts)
#
# 1. genome feature (gff) extraction
#
# read genome feature file
##gff <- import.gff(gff_file,asRangedData=FALSE)
gff <- import.gff(gff_file)

gffgene <- gff[elementMetadata(gff)[,"type"] == "gene"] # extract gene only

gid = elementMetadata(gffgene)[,"ID"]

# Locus tag. If it's JGI gff3, use proteinId instead of locus_tag
if("locus_tag" %in% names(as.data.frame(gffgene))) {
  glocus_tag = elementMetadata(gffgene)[,"locus_tag"]
  id_locus_tag = elementMetadata(gffgene[, c("ID", "locus_tag")])
} else if("proteinId" %in% names(as.data.frame(gffgene))) {
  glocus_tag = elementMetadata(gffgene)[,"proteinId"]
  id_locus_tag = elementMetadata(gffgene[, c("ID", "proteinId")])
  colnames(id_locus_tag)[2] = 'locus_tag'
}

#gDbxref = unlist(elementMetadata(gffgene)[,"Dbxref"])
id_parent = elementMetadata(gff)[elementMetadata(gff)[,"type"] == "gene", c("ID", "Parent")]
id_product = elementMetadata(gff)[elementMetadata(gff)[,"type"] == "CDS", c("Parent", "product")]
id_product = unique(id_product)
id_product$product = sapply(id_product$product, paste, collapse=',')

# considering RNA genes (rRNA/tRNA)
gffrgene = gff[elementMetadata(gff)[,"type"] =="rRNA"|elementMetadata(gff)[,"type"]=="tRNA"|elementMetadata(gff)[,"type"]=="tmRNA"|elementMetadata(gff)[,"type"]=="transcript"]

# RNA gid & index in total gid
rgid = unlist(elementMetadata(gffrgene)[,"Parent"])
rna.ind = match(rgid, gid)

# RNA dbxref & index in total dbxref
#gdbx = unlist(elementMetadata(gffgene)[,"Dbxref"])
#rgdbx = unlist(elementMetadata(gffrgene)[,"Dbxref"])
#rna.ind2 = match(rgdbx, gdbx)
#
# read alignment & count overlaps
# for rRNA removal statistics (all .bam files)
# !! it is QC for NGS library construction !!
#

all.counts = c()
all.rcounts = c()
bam_names = c()
coverage_table = c()

for (i in 1:length(bamFls)) {
  aln = readGAlignments(bamFls[i])
  seqlevels(aln) = seqlevels(gff)    # <--- this part is important!! seqlevels should be same
  hits = countOverlaps(aln,gffgene,ignore.strand=TRUE)    # gene hits
  hits.r = countOverlaps(aln,gffrgene,ignore.strand=TRUE) # RNA gene hits
  counts <- countOverlaps(gffgene, aln[hits==1],ignore.strand=TRUE)          # count gene coverage
  counts.r <- countOverlaps(gffrgene, aln[hits==1],ignore.strand=TRUE)       # count RNA gene coverage
  names(counts) = glocus_tag
  all.counts = cbind(all.counts,counts)
  all.rcounts = cbind(all.rcounts,counts.r)
  # rna gene
  gen = sum(counts)
  rrn = sum(counts.r)
  rratio = ((rrn)/gen)*100
  
  bam_name = unlist(strsplit(bamFls[i], '/'))
  bam_name = tail(bam_name, 1)
  bam_name = sub('.sort.bam', '', bam_name)
  bam_names = c(bam_names, bam_name)
  
  tmp_vector = c(bam_name, gen, rrn, round(rratio,2))# 4825266 4303249 10.8184087675167
  coverage_table = rbind(coverage_table, tmp_vector)
}

## Write table
colnames(coverage_table) = c('name', 'all_counts', 'ncrna_counts', 'ncrna_ratio')
row.names(coverage_table) = 1:nrow(coverage_table)
write.table(coverage_table, file.path(root_dir, paste(project_name, '_cov_table.txt', sep='')), sep='\t', row.names=F, quote=F)

#
# all.counts, all.gcounts, all.rcount
#
colnames(all.counts) = bam_names

# remove rna genes
if(length(rna.ind)>0) {
  all.gcounts = all.counts[-rna.ind,]
} else {
  all.gcounts = all.counts
}

save(all.gcounts, file = "all.gcounts")
#
# sample grouping
#
n = colnames(all.gcounts)
n.grp = sub('[0-9]$', '', n)
n.grp = sub('-$', '', n.grp) # <- you should complete this line by your sample groups

####
#### EdgeR analysis
####
#
# all group analysis
#
grp = all.gcounts
dge = DGEList(grp,group=n.grp)
dge = calcNormFactors(dge)

#head(dge$sample); colnames(dge$counts)

#
# MDS analysis for library quality
# 
grp.col = table(dge$sample$group)
grp.col = rainbow(length(grp.col)) 
grp.cols = grp.col[dge$sample$group]

png(file = file.path(root_dir, paste(project_name, '_mds.png', sep='')),
    width = 3*500, height = 3*500, res = 300)
plotMDS(dge,cex=1,col=grp.cols, method='bcv')
dev.off()

#
# group by group comparison
#
grp = factor(n.grp)
design = model.matrix(~0+grp, data=dge$samples)
colnames(design) = levels(dge$samples$group)

all_comb = list()
comb_name = list()

k = 1
for(i in (1:length(levels(factor(n.grp))))) {
  for(j in 1:length(levels(factor(n.grp)))) {
    if(i < j) {
      all_comb[[k]] = numeric(length(unique(n.grp)))
      all_comb[[k]][i] = -1
      all_comb[[k]][j] = 1
      comb_name[[k]] = paste(unique(n.grp)[i], unique(n.grp)[j], sep='vs')
      k = k + 1
    }
  }
}

# Massage table
result.table = c() # Empty vector