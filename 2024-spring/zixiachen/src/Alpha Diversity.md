## Alpha Diversity
------
### 1. Installation and loading of packages
```r
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("ggsignif")
install.packages("vegan")
install.packages("ggprism")
install.packages("picante")
install.packages("dplyr")
install.packages("RColorBrewer")
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)
```
------

### 2. Import data
```r
setwd("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata")
ASV<-read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV",row.names = 1)
head(ASV)
tASV <- t(ASV)
tASV <- ceiling(as.data.frame(t(ASV)))
```
------

### 3. Calculation of the Alpha Diversity Index
#### 1) Calculate Shannon, Simpson, Richness indices
```r
Shannon <- diversity(tASV, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(tASV, index = "simpson", MARGIN = 2, base = exp(1))
Richness <- specnumber(tASV, MARGIN = 2)
```
#### 2) Tabulate the above diversity indices
```r
index <- as.data.frame(cbind(Shannon, Simpson, Richness))
```
<img width="325" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/01d21feb-a691-409b-b15d-f7bbb9a15aa4">

#### 3) Calculate obs, chao, ace indices
```r
obs_chao_ace <- t(estimateR(ASV))
obs_chao_ace <- obs_chao_ace[rownames(index),]
```
<img width="419" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/5dfe07ca-f206-4344-a634-322cb3e9b8c2">

#### 4) Combine the obs, chao, and ace indices with the results of the previous index calculations
```r
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
```
<img width="474" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/dea25937-dded-4652-9a28-840fff9ff132">

#### 5) Calculating Pielou and Coverage
```r
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(tASV ==1) / colSums(tASV)
```
<img width="644" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/d308b8b1-6813-4aaf-ab6d-60c2568e3cfa">

#### 6) Export Forms
```r
write.table(cbind(sample=c(rownames(index)),index),'diversity.index.txt',row.names=F,sep = '\t', quote =F)
```
------

### 4. Discrepancy calculation and plotting
#### 1) Read in data and group files
###### Read in the file
```r
index <- read.delim('diversity.index.txt',header = T, row.names = 1)
index$samples <- rownames(index)
```
###### Read in grouped files
```r
groups <- read.delim('group.txt',header=T,stringsAsFactors = F)
```
<img width="213" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/e4eef239-102f-4e65-a0fe-145e4bfce393">

###### Merging grouping information with diversity indices
```r
df1 <- merge(index,groups,by = 'samples')
```
<img width="725" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/c43586f8-dd67-45f0-8349-bbc4961a745a">

#### 2) Plotting
##### Setting color
```r
mycol <- c('#E41A1C','#377EB8')
```
##### Shannon
```r
a <- ggbarplot(df1,x = "group",y="Shannon",fill = "group",
          palette=mycol,legend="none",
          add = "mean_se", label=F,)+
  stat_compare_means(method = "anova",label.y = 4.8,label.x = 0.8,
                     label = 'p.format')+
  stat_compare_means(comparisons = list(c('A','B')),
                     method = 'wilcox',label = "p.label")+
  scale_y_continuous(limits = c(0,5),expand = c(0,0))+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text.x = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2,size = 12, face = "bold"))
```
![Shannon](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/d82beac2-dde8-4db0-8018-3684c63b6626)

##### Simpson
```r
b <- ggbarplot(df1,x = "group",y="Simpson",color = "group",
          palette=mycol,legend="none",
          add = "jitter", add.params=list(size=3,jitter=0.1,alpha=0.5),
          outlier.shape =NA)+
  stat_compare_means(method = "kruskal.test",
                     label.y = 2.5, label.x = 0.8,
                     label = 'p.format')+
  stat_compare_means(comparisons = list(c('A','B')),
                     method = 'wilcox',label = "p.signif")+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text.x = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2,size = 12, face = "bold"))
```
![Simpson](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/a02e9b8e-fd26-431f-9ff4-c8bb2445da7f)

##### Richness
```r
c <- ggbarplot(df1,x = "group",y="Richness",color = "group",
          palette=mycol,legend="none",
          add = "mean_sd", error.plot = 'errorbar', 
          size=1,outlier.shape =NA,binwidth=0.5,
       add.params=list(color='black',size=0.8))+
  stat_summary(fun = mean,aes(ymin=..y..,ymax=..y..),colour='black',geom = 'errorbar',size=1.25,width=0.2)+
  stat_compare_means(method = "anova",
                     label.y = 390, label.x = 0.8,
                     label = 'p.format')+
  stat_compare_means(comparisons = list(c('A','B')),
                     method = 'wilcox',label = "p.signif")+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text.x = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2,size = 12, face = "bold"))
```
![Richness](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/24e13c55-7685-4e64-94f9-930441781b5f)

##### Chao
```r
d <- ggviolin(df1,x = "group",y = "Chao",fill = "group",palette = mycol,
         legend="none", alpha = 0.5,
         add = "boxplot", outlier.shape = NA)+
  stat_compare_means(method = "anova",
                     label.y = 450, label.x = 0.8,
                     label = 'p.format')+
  stat_compare_means(comparisons = list(c('A','B')),
                     method = 'wilcox',label = "p.label")+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text.x = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black",face = "bold",hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2,size = 12, face = "bold"))
```
![Chao](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/33c1b376-7e6c-46a5-9b18-0b730341349a)

------

#### 4) Merging images
```r
install.packages("gridExtra")
install.packages("cowplot")
library("gridExtra")
library("cowplot")
plot_grid(a,b,c,d,labels = c("A","B","C","D"),ncol = 4, nrow = 1)
ggsave('diversity.pdf',width = 12, height = 4)
```
![diversity](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/8bbdaa7d-eb1f-42b1-893b-19a1929bcafe)

