# Manual file setup
gff_file <- "C:/R/Zob/Zob.gff"
#https://github.com/igchoi/IBT618-SystemsBiotechnology/blob/main/2024-spring/msan/file/Zob.gff

bamFls <- c("C:/R/bam/1176.bam",
            "C:/R/bam/1179.bam",
            "C:/R/bam/1180.bam",
            "C:/R/bam/1182.bam",
            "C:/R/bam/1183.bam",
            "C:/R/bam/1186.bam")
#https://zenodo.org/records/11187387

root_dir <- 'C:/R/main'
project_name <- '2024_spr'

# Ensure all necessary packages are installed and loaded
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("S4Vectors", "Rsamtools", "edgeR", "GenomicRanges", "GenomicAlignments", "rtracklayer", "DelayedArray", "SummarizedExperiment"))

install.packages(c("lattice", "pheatmap", "ggplot2", "ggExtra", "Matrix"))

# Load libraries
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(edgeR)
library(GenomicRanges)
library(GenomicAlignments)
library(pheatmap)
library(ggplot2)
library(ggExtra)

# Set working directory
setwd(root_dir)  # 작업 디렉토리를 root_dir로 설정

# Read genome feature file
gff <- import.gff(gff_file)  # GFF 파일을 가져와서 gff 객체에 저장

# Extract only gene features
gffgene <- gff[elementMetadata(gff)[,"type"] == "gene"]  # gene 타입의 특징만 추출하여 gffgene에 저장

# Extract gene IDs
gid <- elementMetadata(gffgene)[,"ID"]  # gene ID를 gid 변수에 저장

# Determine the locus tag field
if ("locus_tag" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"locus_tag"]  # locus_tag가 존재하면 glocus_tag에 저장
} else if ("proteinId" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"proteinId"]  # proteinId가 존재하면 glocus_tag에 저장
}

# Extract IDs and products for genes and CDS
id_parent <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "gene", c("ID", "Parent")]  # gene 타입의 ID와 Parent를 추출
id_product <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "CDS", c("Parent", "product")]  # CDS 타입의 Parent와 product를 추출
id_product <- unique(id_product)  # 중복 제거
id_product$product <- sapply(id_product$product, paste, collapse=',')  # product를 콤마로 연결

# Considering RNA genes (rRNA/tRNA)
gffrgene <- gff[elementMetadata(gff)[,"type"] %in% c("rRNA", "tRNA", "tmRNA", "transcript")]  # rRNA, tRNA, tmRNA, transcript 타입의 특징 추출

# RNA gid & index in total gid
rgid <- unlist(elementMetadata(gffrgene)[,"Parent"])  # RNA gene의 Parent ID 추출
rna.ind <- match(rgid, gid)  # RNA gene ID의 인덱스 추출

# Initialize matrices for counts
all.counts <- matrix(0, nrow = length(gid), ncol = length(bamFls))  # gene ID 수와 bam 파일 수에 따라 0으로 초기화된 매트릭스 생성
rownames(all.counts) <- glocus_tag  # 매트릭스 행 이름을 glocus_tag로 설정
all.rcounts <- matrix(0, nrow = length(rgid), ncol = length(bamFls))  # RNA gene ID 수와 bam 파일 수에 따라 0으로 초기화된 매트릭스 생성
colnames(all.counts) <- sub(".bam", "", basename(bamFls))  # 매트릭스 열 이름을 bam 파일 이름으로 설정
colnames(all.rcounts) <- sub(".bam", "", basename(bamFls))

# Read alignments and count overlaps
for (i in seq_along(bamFls)) {
  aln <- readGAlignments(bamFls[i])  # bam 파일에서 정렬된 읽기를 읽어옴
  seqlevels(aln) <- seqlevels(gff)  # 정렬된 읽기의 시퀀스 레벨을 gff와 일치시킴
  counts <- countOverlaps(gffgene, aln, ignore.strand = TRUE)  # gffgene과 정렬된 읽기 간의 겹치는 부분 계산
  counts.r <- countOverlaps(gffrgene, aln, ignore.strand = TRUE)  # gffrgene과 정렬된 읽기 간의 겹치는 부분 계산
  all.counts[, i] <- counts  # 계산된 gene 카운트를 매트릭스에 저장
  all.rcounts[, i] <- counts.r  # 계산된 RNA gene 카운트를 매트릭스에 저장
}

# Coverage table preparation
coverage_table <- data.frame(
  name = sub(".bam", "", basename(bamFls)),  # bam 파일 이름
  all_counts = colSums(all.counts),  # 모든 gene 카운트 합계
  ncrna_counts = colSums(all.rcounts),  # ncRNA 카운트 합계
  ncrna_ratio = round(100 * colSums(all.rcounts) / colSums(all.counts), 2)  # ncRNA 비율 계산
)

# Write the coverage table
write.table(coverage_table, file.path(root_dir, paste(project_name, '_cov_table.txt', sep='')), sep='\t', row.names = FALSE, quote = FALSE)  # 커버리지 테이블 저장

# Remove rRNA genes counts for DEG analysis
if (length(rna.ind) > 0) {
  all.gcounts <- all.counts[-rna.ind, ]  # RNA gene을 제외한 모든 gene 카운트 저장
} else {
  all.gcounts <- all.counts
}

# Save gene counts matrix
save(all.gcounts, file = file.path(root_dir, "all.gcounts.RData"))  # gene 카운트 매트릭스 저장

# Sample grouping
sample_names <- sub(".bam", "", basename(bamFls))  # 샘플 이름 설정
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")  # 샘플 조건 설정
design <- model.matrix(~0 + factor(conditions))  # 모델 매트릭스 생성
colnames(design) <- levels(factor(conditions))  # 디자인 매트릭스 열 이름 설정

# RPKM calculation - gene set subtraction
gsize <- width(gffgene)  # gene의 길이 계산
names(gsize) <- glocus_tag  # gene 길이의 이름을 glocus_tag로 설정
gsize.c <- gsize[rownames(dge$counts)]  # DGE 카운트의 행 이름에 따른 gene 길이 설정
geneLengthsInKB <- (gsize.c / 1000)  # gene 길이를 KB 단위로 변환
dge <- estimateCommonDisp(dge, design)  # 공통 분산 추정
millionsMapped <- colSums(dge$pseudo.counts) / 1e+06  # 백만 단위로 매핑된 읽기 수 계산
rpm <- dge$pseudo.counts / millionsMapped  # RPM 계산
tt.rpkm <- round(rpm / geneLengthsInKB, 1)  # RPKM 계산 및 반올림
colnames(tt.rpkm) <- sub('\\$', '_rpkm', colnames(tt.rpkm))  # RPKM 열 이름 설정
result.table <- cbind(tt.raw.counts, tt_product, tt.pseudocount.round, tt.coeff.log2.round, tt.rpkm)  # 결과 테이블 생성
keep <- rowSums(tt.rpkm >= 2) >= 1  # RPKM이 2 이상인 행 유지
result.table <- result.table[keep, ]  # 필터링된 결과 테이블 저장
dim(result.table)  # 결과 테이블의 차원 출력

for (i in 1:length(all_comb)) {  
  lrTest <- glmLRT(fit, contrast = all_comb[[i]])  # GLM LRT 수행
  tt <- topTags(lrTest, Inf)  # 상위 태그 추출
  tt.tab <- tt$table  # 태그 테이블 저장
  if (!is.null(dim(tt.tab))) {
    tt.tab[,'logFC'] <- round(tt.tab[,'logFC'], 2)  # logFC 반올림
    tt.tab[,'logCPM'] <- round(tt.tab[,'logCPM'], 2)  # logCPM 반올림
    tt.tab[,'LR'] <- round(tt.tab[,'LR'], 2)  # LR 반올림
    qcomb_name <- gsub(" ","vs",gsub("1[*]","",gsub("-1[*]","",lrTest$comparison)))  # 비교 이름 설정
    idx.coeff <- match(rownames(fit$coefficients), rownames(tt.tab))  # 계수 매칭
    tt.tab.sorted <- tt.tab[idx.coeff,]  # 정렬된 태그 테이블 저장
    colnames(tt.tab.sorted) <- sub('\\$', paste('_', qcomb_name, sep=''), colnames(tt.tab.sorted))  # 정렬된 태그 테이블 열 이름 설정
    result.table <- cbind(result.table, tt.tab.sorted)  # 결과 테이블에 추가
  }
}

# EdgeR analysis
dge <- DGEList(counts = all.gcounts, group = conditions)  # DGEList 생성
dge <- calcNormFactors(dge)  # 정규화 계수 계산
dge <- estimateCommonDisp(dge)  # 공통 분산 추정
dge <- estimateTagwiseDisp(dge)  # 태그별 분산 추정

# Fit model and compute contrasts
fit <- glmFit(dge, design)  # 모델 적합 및 대비 계산

# Define contrasts for differential expression analysis
all_comb <- list()  # 대비 리스트 초기화
comb_name <- list()  # 대비 이름 리스트 초기화
k <- 1
for (i in seq_along(levels(factor(conditions)))) {
  for (j in (i+1):length(levels(factor(conditions)))) {
    all_comb[[k]] <- makeContrasts(contrasts = levels(factor(conditions))[c(i, j)], levels = design)  # 대비 생성
    comb_name[[k]] <- paste(levels(factor(conditions))[c(i, j)], collapse = "vs")  # 대비 이름 생성
    k <- k + 1
  }
}

print(all_comb[[2]])  # 두 번째 대비 출력
# Example of correcting a problematic contrast vector manually
all_comb[[2]] <- c(1, -1, 0)  # 특정 요구에 맞게 조정된 대비 벡터

all_comb <- all_comb[-2]  # 두 번째 요소 제거
comb_name <- comb_name[-2]  # comb_name에서도 동일하게 제거

result.tables <- list()  # 결과 테이블 리스트 초기화
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]  # 비교 이름 설정
  lrt <- glmLRT(fit, contrast = all_comb[[k]])  # GLM LRT 수행
  topTags <- topTags(lrt, n = Inf)  # 상위 태그 추출
  result.tables[[comparison]] <- topTags$table  # 결과 테이블에 저장
}

for (k in seq_along(all_comb)) {
  contrast <- all_comb[[k]]  # 대비 설정
  if (any(is.na(contrast)) || any(is.nan(contrast)) || any(is.infinite(contrast))) {
    cat("Skipping contrast for all_comb[[", k, "]] because it contains NA, NaN, or Inf\n")  # NA, NaN 또는 Inf가 있는 대비 건너뛰기
  } else {
    tryCatch({
      lrt <- glmLRT(fit, contrast = contrast)  # GLM LRT 수행
      topTags <- topTags(lrt, n = Inf)  # 상위 태그 추출
      result.tables[[comb_name[[k]]]] <- topTags$table  # 결과 테이블에 저장
    }, error = function(e) {
      cat("Error in processing contrast for all_comb[[", k, "]]: ", e$message, "\n")  # 대비 처리 중 오류 발생 시 메시지 출력
    })
  }
}

for (i in seq_along(fit)) {
  element <- fit[[i]]  # fit의 각 요소에 대해 반복
  if (is.list(element)) {
    for (j in seq_along(element)) {
      sub_element <- element[[j]]  # 하위 요소에 대해 반복
      if (any(is.na(sub_element)) || any(is.nan(sub_element)) || any(is.infinite(sub_element))) {
        cat("NA, NaN, or Inf found in sub-element", j, "of element", i, "of the list.\n")  # 하위 요소에 NA, NaN 또는 Inf가 있는 경우 메시지 출력
      }
    }
  } else {
    if (any(is.na(element)) || any(is.nan(element)) || any(is.infinite(element))) {
      cat("NA, NaN, or Inf found in element", i, "of the list.\n")  # 요소에 NA, NaN 또는 Inf가 있는 경우 메시지 출력
    }
  }
}

# Analyze differential expression for each comparison
result.tables <- list()  # 결과 테이블 리스트 초기화
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]  # 비교 이름 설정
  lrt <- glmLRT(fit, contrast = all_comb[[k]])  # GLM LRT 수행
  topTags <- topTags(lrt, n = Inf)  # 상위 태그 추출
  result.tables[[comparison]] <- topTags$table  # 결과 테이블에 저장
}

# Save DEG results
for (comparison in names(result.tables)) {
  write.table(result.tables[[comparison]], file.path(root_dir, paste(project_name, '_', comparison, '_DE_results.txt', sep = '')), sep = '\t', quote = FALSE, row.names = TRUE)  # DEG 결과 저장
}

# MDS plot for sample quality
grp = all.gcounts  # 모든 gene 카운트
png(file = file.path(root_dir, paste(project_name, '_MDS_all_samples.png', sep = '')), width = 2000, height = 2000, res = 300)  # MDS 플롯 저장 경로 설정
plotMDS(dge, col = as.numeric(factor(conditions)), pch = 16, cex = 1.5, main = "MDS Plot of all samples")  # MDS 플롯 생성
dev.off()  # 플롯 저장

# Raw count
tt.raw.counts <- dge$counts  # DGE 카운트 저장
colnames(tt.raw.counts) <- sub('\\$', '_raw_counts', colnames(tt.raw.counts))  # 열 이름 설정

# Pseudocount
dge.commondisp <- estimateCommonDisp(dge, design)  # 공통 분산 추정
tt.pseudocount <- dge.commondisp$pseudo.counts  # 의사 카운트 저장
tt.pseudocount.round <- round(tt.pseudocount)  # 의사 카운트 반올림
colnames(tt.pseudocount.round) <- sub('\\$', '_pseudocount', colnames(tt.pseudocount.round))  # 열 이름 설정

# Coefficient from fit object (note that it is log2 value)
tt.coeff <- fit$coefficients  # 계수 저장
tt.coeff.log2 <- tt.coeff / log(2)  # log2 값으로 변환
tt.coeff.log2.round <- round(tt.coeff.log2, 2)  # 반올림
colnames(tt.coeff.log2.round) <- sub('\\$', '_coeff', colnames(tt.coeff.log2.round))  # 열 이름 설정

# Write result table to a file
dge_file <- file.path(root_dir, paste(project_name, '_deg_table_20210405.txt', sep=''))  # 결과 테이블 저장 경로 설정
write.table(result.table, dge_file, sep='\t', row.names=T, quote=F)  # 결과 테이블 저장

# Heatmap
library(pheatmap)
rpkm.avg <- round(rowMeans(tt.rpkm), 2)  # RPKM 평균 계산

if (length(unique(grp)) > 1) {
  pheatmap(mat = tt.rpkm, show_rownames = F, show_colnames = T, cluster_rows = T, cluster_cols = T, scale = 'row')  # 히트맵 생성
} else {
  cat("At least two groups are required for heatmap visualization.\n")  # 히트맵 생성을 위해 최소 두 그룹이 필요하다는 메시지 출력
}
