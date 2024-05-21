## Rarefaction Curve
생태학에서 희소성 분석은 샘플링 결과로부터 종의 풍부도를 평가하는 기법입니다. 희귀도는 소위 희귀도 곡선의 구성을 기반으로 주어진 수의 개별 샘플에 대한 종의 풍부도를 계산할 수 있습니다. 이 곡선은 샘플 수에 따른 종의 수를 그래프로 나타낸 것입니다. 희귀도 곡선은 일반적으로 가장 흔한 종들이 발견되면서 처음에는 빠르게 증가하지만, 가장 희귀한 종들만 표본으로 남게 되면 곡선이 정체됩니다. -[위키피디아](https://en.wikipedia.org/wiki/Rarefaction_(ecology))

```r
install.packages("vegan")
library(vegan)	
ASV <- read.csv("/Users/livia/Desktop/Daqu sample/R analysis-4.30/Rawdata/filtered/ASV.CSV", row.names = 1)
head(ASV)
rarecurve(ASV, step = 2000, col = c('red', 'green', 'blue', 'orange'))
```
------
![Rarefaction curve](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165759135/6b724475-dc68-45a4-b744-31d271e7ad5b)
