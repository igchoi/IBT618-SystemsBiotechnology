# Class project
## Title : Identification of Differentially Expressed Genes from RNAseq data
-------- 
### Introduction
* Data: publicly available marine heterotrophic bacteria grown on different carbon sources
  - where can I get the public data (raw reads)? [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
  - search keyword: `(marine bacteria carbon source) AND "Zobellia galactanivorans"[orgn:__txid63186] `
#### Background
Why do it?
* To identify changes in genes or mechanisms that occur under certain conditions and either suppress them or induce more expression.
-------
![image](https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/165700031/91bd24e3-56e5-443e-b9a7-ed2d8caa5382)

  - Check Quality : The process of filtering out low-quality parts of something fetched via Rfastp
  - Read Sample : Sorting and calculations for the sequences 
  - Counting in exons models : Differentiate exons by collapsing exons to the nonoverlapping, disjoint regions
  - Check Mean and variance relationship : This will assist us to detect differential expression with a small number of replicates.
  - Multiple testing : If we are 95% confident a gene is differentially expressed when we look at all genes we can expect 5% to be false. When looking across all genes we apply a multiple testing correction to account for this.
  - Genes overlaps : Compare multiple results to clean up overlaps.

------
### Data list
  
  - GSM5699202: Free-living cells with maltose rep 1; Zobellia galactanivorans; RNA-Seq / https://www.ncbi.nlm.nih.gov/sra/SRX13191408[accn]
  - GSM5699203: Free-living cells with maltose rep 2; Zobellia galactanivorans; RNA-Seq / https://www.ncbi.nlm.nih.gov/sra/SRX13191409[accn]
  - GSM5699204: Free-living cells with maltose rep 3; Zobellia galactanivorans; RNA-Seq / https://www.ncbi.nlm.nih.gov/sra/SRX13191410[accn]

#### Links
1. [RU_RNASeq Workflow](https://rockefelleruniversity.github.io/RU_RNAseq/)
2. [Bioconductor Courses & Workshops](https://www.bioconductor.org/help/course-materials/)

