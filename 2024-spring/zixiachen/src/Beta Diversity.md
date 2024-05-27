## Beta Diversity
Beta analysis is a method used to compare the diversity between different habitats or samples.  It focuses on assessing the differences in species composition between environments, rather than within a single environment.  Beta diversity provides insights into the turnover of species from one habitat to another and helps in understanding the patterns and drivers of diversity across spatial scales.

------
### 1. Installation and loading of packages
```r
install.packages("ape")
install.packages("vegan")
library(ape)
library(vegan)
```
------

### 2. Import data
```r
data <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV",row.names = 1)
tdata <- t(data)
head(data)
```
------

### 3. Compute distance matrix and perform hierarchical clustering
```r
data.dist <- vegdist(t(tdata))
hc <- hclust(data.dist,method = "average")
plot(hc)
plot(hc,hang = -1)
```
<img width="358" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/24a76800-0128-4c27-9f03-4fefdc2c3660">
<img width="342" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/8423b8c7-4f4c-4799-a6af-22e01207c17d">
------

### 4. Formatting and plotting of clustering results
```r
tree <- as.dendrogram(hc)
pdf("hcluste.tree.pdf")
plot(tree,type = "rectangle",horiz=T)
dev.off()
```
<img width="611" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/895c3259-efe1-49bd-a4fe-e7dd7afa0b28">
------

### 5. Add groups and save pictures
```r
pdf("hcluster.tree.group.pdf")
hp <- as.phylo(hc)
co <- def(hp$tip.label,"BDC_F_filt.fastq.gz"="darkblue","SQ-2_F_filt.fastq.gz"="darkblue",
          "FJ-4_F_filt.fastq.gz"="red","FJ-7_F_filt.fastq.gz"="red",regexp = T)
plot.phylo(hp,tip.color=co,direction="rightwards")
legend("topright",legend = c("Light","Sauce"),box.lty = 0,fill = c("red","darkblue"),col = c("red","darkblue"))
dev.off()
```
<img width="559" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/12046c48-93ac-4354-96b8-45eb3fc21d2b">
------

### 6. PCA analysis
```r
ASV <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV",row.names = 1)
pca <- prcomp(ASV, scal = TRUE)
## Get an overview of PCA results
pca_sum <- summary(pca)
pc12 <- pca$x[,1:2]
## Calculate the dilution of each principal component
pc <- pca_sum$importance[2,]*100
install.packages("ggplot2")
library(ggplot2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=PC1, y=PC2))+
  geom_point(size=3)+
  theme_bw()+
  geom_text(aes(label=samples, y=PC2),size = 4,vjust=1.5)
p 
group <- read.table("group.txt",sep = "\t",header = T)
colnames(group) <- c("samples","group")
head(group)
df <- merge(pc12,group,by="samples")
p <- ggplot(df,aes(x=PC1,y=PC2,color=group))+
  geom_point(size=3)+
  geom_text(aes(label=samples,y=PC2+0.1),size=4,vjust=0)
p
p <- p+theme_bw()+ #Using the black and white theme
  theme(axis.title.x = element_text(size = 15,family = "sans"),
        axis.title.y = element_text(size = 15,family = "sans",angle=90),
        axis.text.y = element_text(size = 12,family = "sans"),
        axis.text.x = element_text(size = 12,family = "sans"),
        panel.grid = element_blank())
p
p <- p+xlab(paste0("PC1 (",pc[1],"%)"))+ylab(paste0("PC2 (",pc[2],"%)"))
p
ggsave("pca.pdf",p,width = 5,height = 5)
```
<img width="361" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/84ea109e-3340-4658-8d9e-d2f92f254fac">
<img width="354" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/0db997db-289c-46dd-aaf8-c1d161749e5f">
<img width="358" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/eaf825c5-be57-49eb-908a-ced694393c56">
<img width="360" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/70d41cff-5212-45d6-9b7a-b61fe8faf31f">
------

### 7. PCoA analysis
```r
ASV <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV",row.names = 1)
ASV.dist <- vegdist(ASV)
pcoa <- cmdscale(ASV.dist,eig = TRUE)
pc12 <- pcoa$points[,1:2]
pc_importance <- round(pcoa$eig/sum(pcoa$eig)*100,digits = 2)
library(ggplot2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=V1,y=V2))+
  geom_point(size=3)
p
group <- read.table("group.txt",sep = "\t",header = T)
colnames(group) <- c("samples","group")
df <- merge(pc12,group,by="samples")
head(df)
p <- ggplot(df,aes(x=V1,y=V2,colour = group))+
  geom_point(size=3)+
  geom_text(aes(label = samples, y=V2+0.01),size = 4,vjust = 0)+
  guides(color=guide_legend(title = NULL))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 15,family = "sans"),
        axis.title.y = element_text(size = 15,family = "sans",angle=90),
        axis.text.y = element_text(size = 12,family = "sans"),
        axis.text.x = element_text(size = 12,family = "sans"),
        panel.grid = element_blank())
p <- p+
  xlab(paste0("PCo1 (",pc_importance[1],"%)"))+
  ylab(paste0("PCo2 (",pc_importance[2],"%)"))
p
p <- p+stat_ellipse(level = 0.95) ## Add confidence ellipse, but here the ellipse is not calculated because there are too few points.
p 
ggsave("pcoa.pdf",p,width = 6,height = 5)
```
<img width="358" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/ef1a094c-ecb5-46a0-b6ef-fa83106a91e9">
<img width="359" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/35bf66cb-8d62-4f04-abc7-93f416e3c78a">
------
