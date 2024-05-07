
### Manual file 
gff_file = '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_ref/Colwellia_echini_GCF_002843355_2_ASM284335v2_genomic.gff'
bamFls = c('/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/aga1.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/aga2.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/aga3.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/glu1.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/glu2.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/glu3.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/car1.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/car2.sort.bam',
           '/eevee/backup1/bhpark/rnaseq/2018_col/20181213_col/new_bowtie/car3.sort.bam')

root_dir = '/home/bgkim/Assembly'

project_name = '2018_col_new'

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

# DEG
dge.glm = estimateGLMCommonDisp(dge, design)
fit = glmFit(dge.glm, design,dispersion=dge.glm$common.dispersion)

# Raw count
tt.raw.counts = dge$counts
colnames(tt.raw.counts) = sub('$', '_raw_counts', colnames(tt.raw.counts))

# Pseudocount
#dge.commondisp = estimateCommonDisp(dge, design, tol=0.0000001)
dge.commondisp = estimateCommonDisp(dge, design)
tt.pseudocount = dge.commondisp$pseudo.counts
tt.pseudocount.round = round(tt.pseudocount)
colnames(tt.pseudocount.round) = sub('$', '_pseudocount', colnames(tt.pseudocount.round))

# Coefficient from fit object (note that it is log2 value)
tt.coeff = fit$coefficients
tt.coeff.log2 = tt.coeff / log(2)
tt.coeff.log2.round = round(tt.coeff.log2, 2)
colnames(tt.coeff.log2.round) = sub('$', '_coeff', colnames(tt.coeff.log2.round))

# Get each gene's product
row_tt = rownames(tt.raw.counts)
tmp_index = match(row_tt, id_locus_tag$locus_tag)
tt_gene_id = id_locus_tag$ID[tmp_index]
tmp_index2 = match(tt_gene_id, unlist(id_parent$ID))
tt_rna_id = id_parent$ID[tmp_index2]
tmp_index3 = match(tt_rna_id, unlist(id_product$Parent))
if(length(which(is.na(tmp_index3) == FALSE)) > 0) {
  tt_product = unlist(id_product$product[tmp_index3])  
} else {
  cat('Something wrong with product name\n')
}
#
# RPKM calculation - gene set subtraction
#
gsize = width(gffgene)
names(gsize) = glocus_tag
gsize.c = gsize[rownames(dge$counts)]

# rpkm
geneLengthsInKB <- (gsize.c/1000)
dge = estimateCommonDisp(dge,design)
millionsMapped <- colSums(dge$pseudo.counts)/1e+06 # Factor for converting to million of mapped reads. 
rpm <- dge$pseudo.counts/millionsMapped
tt.rpkm <- round(rpm/geneLengthsInKB,1)
colnames(tt.rpkm) = sub('$', '_rpkm', colnames(tt.rpkm))

result.table = cbind(tt.raw.counts, tt_product, tt.pseudocount.round, tt.coeff.log2.round, tt.rpkm)



for(i in 1:length(all_comb)) {  
  lrTest = glmLRT(fit, contrast=all_comb[[i]])
  
  tt = topTags(lrTest, Inf)
  tt.tab = tt$table#[which(tt$table$PValue<0.001),]
  tt.tab[,'logFC'] = round(tt.tab[,'logFC'], 2)
  tt.tab[,'logCPM'] = round(tt.tab[,'logCPM'], 2)
  tt.tab[,'LR'] = round(tt.tab[,'LR'], 2)
  
  qcomb_name = gsub(" ","vs",gsub("1[*]","",gsub("-1[*]","",lrTest$comparison)))
  
  idx.coeff = match(rownames(fit$coefficients), rownames(tt.tab))
  tt.tab.sorted = tt.tab[idx.coeff,]
  colnames(tt.tab.sorted) = sub('$', paste('_', qcomb_name, sep=''), colnames(tt.tab.sorted))
  result.table = cbind(result.table, tt.tab.sorted)
}

dge_file = file.path(root_dir, paste(project_name, '_deg_table_20210405.txt', sep=''))
write.table(result.table, dge_file, sep='\t', row.names=T, quote=F)

#
# heatmap
#
library(pheatmap)
rpkm.avg = round(t(sapply(rownames(tt.rpkm),function(x) tapply(tt.rpkm[x,],grp,mean))),2)

tt.pval = result.table[,grepl("PValue", colnames(result.table))]
tt.logfc = result.table[,grepl("logFC", colnames(result.table))]


#tt.pval.min = apply(tt.pval, 1, min)
#tt.logfc.max = apply(abs(tt.logfc),1,max)

#deg.genes = rownames(result.table)[(abs(tt.logfc.max)>1) * (tt.pval.min<0.01)==1]
deg.genes = rownames(result.table)[(abs(tt.logfc)>2) * (tt.pval<0.01)==1]
length(deg.genes)
nrow(result.table)

tt.rpkm.forheatmap = tt.rpkm[rowSums(tt.rpkm) != 0,]
pheatmap(mat = tt.rpkm.forheatmap, show_rownames = F, show_colnames = T, cluster_rows = T, cluster_cols = T, scale = 'row')

tt.rpkm.forheatmap.onlydeg = tt.rpkm.forheatmap[deg.genes,]
pheatmap(mat = tt.rpkm.forheatmap.onlydeg, show_rownames = F, show_colnames = T, cluster_rows = T, cluster_cols = T, scale = 'row')


#
# doubleplot
#
install.packages("ggExtra")
library(ggplot2)
indf = as.data.frame(rpkm.avg)
indf$Expression = "No DEGs"
indf[(tt.logfc>1) * (tt.pval<0.01) == 1,"Expression"] <- "Down DEGs"
indf[(tt.logfc<(-1)) * (tt.pval<0.01) == 1,"Expression"] <- "Up DEGs"
colnames(indf) = gsub("NA", "Control", gsub("LA", "Microgravity", colnames(indf)))
indf$LogRPKM_Microgravity = log(indf$Microgravity+1, base = 10)
indf$LogRPKM_Control = log(indf$Control+1, base = 10)
indf$Expression = factor(indf$Expression, levels = c("Up DEGs", "No DEGs", "Down DEGs"))

p1 <- ggplot(data = indf, aes(x=LogRPKM_Control,y=LogRPKM_Microgravity, color = Expression)) + geom_point() + scale_color_manual(values = c("tomato", "#999999", "steelblue")) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

library(ggExtra)
ggMarginal(p1, type = "histogram")

####GOenrichment####

indf



####DEG getting#### - 20170115
# 
# rpkm.avg = round(t(sapply(rownames(tt.rpkm),function(x) tapply(tt.rpkm[x,],grp,mean))),2)
# deg.rpkm.avg = rpkm.avg[deg.genes,]
# 
# deg.rpkm.avg.scaled = c()
# for(i in 1:nrow(deg.rpkm.avg)){
#   query = deg.rpkm.avg[i,]
#   scaled = (query-mean(query))/sd(query)
#   deg.rpkm.avg.scaled = rbind(deg.rpkm.avg.scaled, scaled)
# }
# rownames(deg.rpkm.avg.scaled) = rownames(deg.rpkm.avg)
# 
# ph = pheatmap(deg.rpkm.avg.scaled, scale = 'row', color = colorRampPalette(c("skyblue2", "gray95", "coral2"))(100),
#               cluster_rows = T, cluster_cols = T, show_rownames =F)
# ph$tree_row$labels
# 
# hclr = cutree(ph$tree_row, k=6)
# hclrano = data.frame(paste("cluster", hclr))
# rownames(hclrano) = names(hclr)
# colnames(hclrano) = 'cluster'
# 
# cnames = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5", "cluster 6")
# 
# colorset = c('#F1948A', '#F7DC6F', '#7DCEA0', '#7FB3D5', '#BB8FCE', '#85929E')
# 
# cluster        <- colorset
# names(cluster) <- cnames
# anno_colors    <- list(cluster = cluster)
# 
# ph = pheatmap(deg.rpkm.avg.scaled, scale = 'row', color = colorRampPalette(c("skyblue2", "gray95", "coral2"))(100), fontsize = 12,
#               cluster_rows = T, cluster_cols = T, annotation_row = hclrano, show_rownames =F, annotation_colors = anno_colors)
# 
# 
# sum((abs(tt.logfc$logFC_TF_ESvsWT_ES) > 1 ) * (tt.pval$PValue_TF_ESvsWT_ES < 0.01))
# sum((abs(tt.logfc$logFC_TF_MLvsWT_ML) > 1 ) * (tt.pval$PValue_TF_MLvsWT_ML < 0.01))
# sum((abs(tt.logfc$logFC_WT_ESvsWT_ML) > 1 ) * (tt.pval$PValue_WT_ESvsWT_ML < 0.01))
# sum((abs(tt.logfc$logFC_TF_ESvsTF_ML) > 1 ) * (tt.pval$PValue_TF_ESvsTF_ML < 0.01))
# sum((abs(tt.logfc$logFC_TF_ESvsWT_ML) > 1 ) * (tt.pval$PValue_TF_ESvsWT_ML < 0.01))
# sum((abs(tt.logfc$logFC_TF_MLvsWT_ES) > 1 ) * (tt.pval$PValue_TF_MLvsWT_ES < 0.01))
# 
# 
# 
# 
# topdegs = names(sort(tt.pval.min)[1:30])
# topdegs = c("ELI_4394", "ELI_3815", "ELI_0016", "ELI_1842")
# tpd = cbind(as.vector(result.table[topdegs,'tt_product']), rpkm.avg[topdegs,], tt.pval.min[topdegs])
# tpd = data.frame(tpd, stringsAsFactors = F)
# for(i in 2:6){
#   tpd[,i] = as.numeric(tpd[,i])
# }
# as.numeric(tpd$TF_ES)
# tpd
