## Rarefaction Curve
------
```r
install.packages("vegan")
library(vegan)	
ASV <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV", row.names = 1)
head(ASV)
rarecurve(ASV, step = 2000, col = c('red', 'green', 'blue', 'orange'))
```
------
![Rarefaction curve](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/6b724475-dc68-45a4-b744-31d271e7ad5b)
