## Relative abundance of bacterial communities
Based on the results of taxonomic analysis, the community structure composition at various taxonomic levels (such as domain, kingdom, phylum, class, order, family, genus, species, and ASV) for different groups (or samples) can be determined. 

The community bar chart intuitively presents two key pieces of information: (1) the types of microorganisms present in each sample at a specific taxonomic level; and (2) the relative abundance (proportion) of each microorganism within the samples.

------
### 1. Installation and loading of packages
```r
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install("dplyr")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("ggplot2")
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
```
------

### 2. Phylum level
#### 1) Obtaining Phylum table
```r
t_species <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV",row.names = 1)
species <- t(t_species)
taxa <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/taxa.CSV",row.names = 1)
dspecies <- data.frame(species)
dtaxa <- data.frame(taxa)
Phylum <- dspecies %>%
  group_by(taxa$Phylum) %>%
  summarise_all(sum)
head(Phylum)
```
------
#### 2) Delete the last line (NA) and save Phylum table
```r
Phylum <- Phylum[-nrow(Phylum), ]
head(Phylum)
write.csv(Phylum,file="/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/Phylum.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", 
          na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
```
------
#### 3) Read Phylum table
```r
data <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/Phylum.csv", header=T, check.names=FALSE, row.names = 1)
head(data)
```
------
#### 4) Calculate the total abundance of species and arrange them in descending order
```r
rowsum <- apply(data,1,sum)
data <- data[order(rowsum, decreasing = TRUE),]
```
------
#### 5) Calculate the relative abundance of species
```r
da <- apply(data, 2, function(x) x/sum(x))
head(da)
```
------
#### 6) Variable format conversion
```r
df <- melt(da)
head(df)
```
------
#### 7) Plotting
```r
mycol <- c(119,132,147,454,89,404,123,529,463,104,
           552,28,54,84,256,100,558,43,652,31,610,
           477,588,99,81,503,562,76,96,495)
mycol <- colors()[rep(mycol,20)]
P <- ggplot(df, aes(x=Var2,y=value))+
  geom_bar(aes(fill=Var1),stat = "identity")+
  scale_fill_manual(values = mycol)+
  guides(fill= guide_legend(reverse = T,title = NULL, ncol = 4))+
  labs(x="Sample",y="Relative abundance")+
  theme_classic()+
  theme(legend.position = "bottom")
P
ggsave("Phylum.pdf",P,width=8,height=7)
```
-----
<img width="753" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/9542cc84-17e3-4092-9277-61e20d9c3420">


### 3. Genus level
#### 1) Obtaining Genus table
```r
t_species <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV",row.names = 1)
species <- t(t_species)
taxa <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/taxa.CSV",row.names = 1)
dspecies <- data.frame(species)
dtaxa <- data.frame(taxa)
Genus <- dspecies %>%
  group_by(taxa$Genus) %>%
  summarise_all(sum)
head(Genus)
```
------
#### 2) Delete the last line (NA) and save Genus table
```r
Genus <- Genus[-nrow(Genus), ]
head(Genus)
write.csv(Genus,file="/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/Genus.CSV",,append = FALSE, quote = FALSE , sep = " ",eol = "\n", 
          na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
```
------
#### 3) Read Genus table
```r
data <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/Genus.csv", header=T, check.names=FALSE, row.names = 1)
head(data)
```
------
#### 4) Calculate the total abundance of species and arrange them in descending order
```r
rowsum <- apply(data,1,sum)
data <- data[order(rowsum, decreasing = TRUE),]
```
------
#### 5) Calculate the relative abundance of species
```r
da <- apply(data, 2, function(x) x/sum(x))
head(da)
```
------
#### 6) Select the top 20 Genus
```r
da <-  da[1:20,]
```
------
#### 7) Calculate the total abundance of the remaining species
```r
da1 <- 1-apply(da, 2, sum)
```
------
#### 8) Combining data
```r
da2 <- rbind(da,da1)
row.names(da2)[21]="Others"
```
------
#### 9) Variable format conversion
```r
df <- melt(da2)
head(df)
```
------
#### 10) Plotting
```r
mycol <- c(119,132,147,454,89,404,123,529,463,104,
           552,28,54,84,256,100,558,43,652,31,610,
           477,588,99,81,503,562,76,96,495)
mycol <- colors()[rep(mycol,20)]
P <- ggplot(df, aes(x=Var2,y=value))+
  geom_bar(aes(fill=Var1),stat = "identity")+
  scale_fill_manual(values = mycol)+
  guides(fill= guide_legend(reverse = T,title = NULL, ncol = 4))+
  labs(x="Sample",y="Relative abundance")+
  theme_classic()+
  theme(legend.position = "bottom")
P
ggsave("Gunes.pdf",P,width=8,height=7)
```
------
<img width="754" alt="image" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/b2525903-0e76-4e78-9626-a6ebd34ffd78">

