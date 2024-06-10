# 파일 설정
gff_file <- "C:/R/Zob/Zob.gff"
# 참조 유전체를 다운로드할 수 있습니다.
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=63186
bamFls <- c("C:/R/bam/1176.bam",
            "C:/R/bam/1179.bam",
            "C:/R/bam/1180.bam",
            "C:/R/bam/1182.bam",
            "C:/R/bam/1183.bam",
            "C:/R/bam/1186.bam")
# 파일 링크 : https://zenodo.org/records/11187387
root_dir <- 'C:/R/main'
project_name <- '2024_spr'

# 필요한 패키지가 설치되어 있는지 확인하고 설치
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("S4Vectors", "Rsamtools", "edgeR", "GenomicRanges",
                       "GenomicAlignments", "rtracklayer", "DelayedArray", "SummarizedExperiment"))

install.packages(c("lattice", "pheatmap", "ggplot2", "ggExtra", "Matrix"))

# 라이브러리 로드
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(edgeR)
library(GenomicRanges)
library(GenomicAlignments)
library(pheatmap)
library(ggplot2)
library(ggExtra)

# 작업 디렉토리 설정
setwd(root_dir)  # 작업 디렉토리를 root_dir로 설정

# 유전체 특징 파일 읽기
gff <- import.gff(gff_file)  # GFF 파일을 가져와서 gff 객체에 저장

# gene 특징만 추출
gffgene <- gff[elementMetadata(gff)[,"type"] == "gene"]  # gene 타입의 특징만 추출하여 gffgene에 저장

# gene ID 추출
gid <- elementMetadata(gffgene)[,"ID"]  # gene ID를 gid 변수에 저장

# locus tag 필드 결정
if ("locus_tag" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"locus_tag"]  # locus_tag가 존재하면 glocus_tag에 저장
} else if ("proteinId" %in% names(as.data.frame(gffgene))) {
  glocus_tag <- elementMetadata(gffgene)[,"proteinId"]  # proteinId가 존재하면 glocus_tag에 저장
}

# gene과 CDS에 대한 ID와 product 추출
id_parent <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "gene", c("ID", "Parent")]  # gene 타입의 ID와 Parent를 추출
id_product <- elementMetadata(gff)[elementMetadata(gff)[,"type"] == "CDS", c("Parent", "product")]  # CDS 타입의 Parent와 product를 추출
id_product <- unique(id_product)  # 중복 제거
id_product$product <- sapply(id_product$product, paste, collapse=',')  # product를 콤마로 연결

# RNA gene (rRNA/tRNA) 고려
gffrgene <- gff[elementMetadata(gff)[,"type"] %in% c("rRNA", "tRNA", "tmRNA", "transcript")]  # rRNA, tRNA, tmRNA, transcript 타입의 특징 추출

# RNA gene ID와 전체 gene ID에서의 인덱스 추출
rgid <- unlist(elementMetadata(gffrgene)[,"Parent"])  # RNA gene의 Parent ID를 rgid 변수에 저장
rna.ind <- match(rgid, gid)  # rgid를 gid에서 매칭하여 인덱스 저장

# 카운트를 위한 행렬 초기화
all.counts <- matrix(0, nrow=length(gid), ncol=length(bamFls))  # 모든 gene의 카운트를 저장할 행렬 초기화
rownames(all.counts) <- glocus_tag  # 행 이름을 glocus_tag로 설정
all.rcounts <- matrix(0, nrow=length(rgid), ncol=length(bamFls))  # RNA gene의 카운트를 저장할 행렬 초기화
colnames(all.counts) <- sub(".bam", "", basename(bamFls))  # 열 이름을 BAM 파일 이름으로 설정
colnames(all.rcounts) <- sub(".bam", "", basename(bamFls))  # 열 이름을 BAM 파일 이름으로 설정

# 정렬된 읽기를 가져오고 겹치는 부분을 카운트
for (i in seq_along(bamFls)) {
  aln <- readGAlignments(bamFls[i])  # BAM 파일에서 정렬된 읽기를 가져옴
  seqlevels(aln) <- seqlevels(gff)  # 정렬된 읽기의 시퀀스 레벨을 gff의 시퀀스 레벨로 설정
  counts <- countOverlaps(gffgene, aln, ignore.strand=TRUE)  # gene에 대한 카운트 계산
  counts.r <- countOverlaps(gffrgene, aln, ignore.strand=TRUE)  # RNA gene에 대한 카운트 계산
  all.counts[, i] <- counts  # 계산된 카운트를 all.counts 행렬에 저장
  all.rcounts[, i] <- counts.r  # 계산된 RNA 카운트를 all.rcounts 행렬에 저장
}

# 커버리지 테이블 준비
coverage_table <- data.frame(
  name=sub(".bam", "", basename(bamFls)),  # 파일 이름 설정
  all_counts=colSums(all.counts),  # 모든 gene의 총 카운트
  ncrna_counts=colSums(all.rcounts),  # ncRNA의 총 카운트
  ncrna_ratio=round(100 * colSums(all.rcounts) / colSums(all.counts), 2)  # ncRNA 비율 계산
)

# 커버리지 테이블 파일로 저장
write.table(coverage_table, file.path(root_dir, paste(project_name, '_cov_table.txt', sep='')), sep='\t', row.names=FALSE, quote=FALSE)  # 커버리지 테이블을 파일로 저장

# DEG 분석을 위한 rRNA gene 카운트 제거
all.gcounts <- if (length(rna.ind) > 0) all.counts[-rna.ind, ] else all.counts  # rRNA gene 카운트를 제거한 gene 카운트 행렬

# DGEList 객체 생성
dge <- DGEList(counts = all.gcounts)  # DGEList 객체 생성
dge <- calcNormFactors(dge)  # 정규화 인자 계산

# gene 카운트 행렬 저장
save(all.gcounts, file=file.path(root_dir, "all.gcounts.RData"))  # gene 카운트 행렬 저장

# 샘플 그룹화
sample_names <- sub(".bam", "", basename(bamFls))  # 샘플 이름 설정
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")  # 조건 설정
design <- model.matrix(~0 + factor(conditions))  # 실험 설계 행렬 생성
colnames(design) <- levels(factor(conditions))  # 열 이름 설정

# 이분산성 추정
dge <- estimateDisp(dge, design)  # 이분산성 추정

# 모델 생성 및 DE 분석 수행
fit <- glmFit(dge, design)  # 일반화 선형 모델 생성 및 분석

# RPKM 계산 - gene 세트 제거
gsize <- width(gffgene)  # 유전자 길이 계산
names(gsize) <- glocus_tag  # 유전자 길이 이름 할당
gsize.c <- gsize[rownames(dge$counts)]  # DGEList 유전자 길이
geneLengthsInKB <- (gsize.c / 1000)  # 유전자 길이를 킬로베이스로 변환

# 의사 계수를 사용한 정규화
dge <- estimateCommonDisp(dge, design)  # 공통 분산 추정
millionsMapped <- colSums(dge$pseudo.counts) / 1e+06  # 매핑된 읽기를 백만 단위로 변환
rpm <- dge$pseudo.counts / millionsMapped  # 백만당 읽기 수 계산
tt.rpkm <- round(rpm / geneLengthsInKB, 1)  # RPKM 값 계산
colnames(tt.rpkm) <- sub('\\$', '_rpkm', colnames(tt.rpkm))  # 열 이름 설정
keep <- rowSums(tt.rpkm >= 2) >= 1  # RPKM >= 2 인 유전자 유지

# 결과 테이블을 위한 매트릭스 초기화
tt.raw.counts <- dge$counts  # 원시 카운트
colnames(tt.raw.counts) <- sub('\\$', '_raw_counts', colnames(tt.raw.counts))  # 열 이름 설정

# 의사 카운트
dge.commondisp <- estimateCommonDisp(dge, design)  # 공통 분산 추정
tt.pseudocount <- dge.commondisp$pseudo.counts  # 의사 카운트 계산
tt.pseudocount.round <- round(tt.pseudocount)  # 의사 카운트 반올림
colnames(tt.pseudocount.round) <- sub('\\$', '_pseudocount', colnames(tt.pseudocount.round))  # 열 이름 설정

# 피팅 객체로부터의 계수 (로그2 값임)
tt.coeff <- fit$coefficients  # 계수 추출
tt.coeff.log2 <- tt.coeff / log(2)  # 로그2 값으로 변환
tt.coeff.log2.round <- round(tt.coeff.log2, 2)  # 값 반올림
colnames(tt.coeff.log2.round) <- sub('\\$', '_coeff', colnames(tt.coeff.log2.round))  # 열 이름 설정

# 결과 테이블 초기화
result.table <- cbind(tt.raw.counts, tt.pseudocount.round, tt.coeff.log2.round, tt.rpkm)

# 차등 발현 분석 설정
all_comb <- list()
comb_name <- list()
k <- 1
for (i in seq_along(levels(factor(conditions)))) {
  for (j in (i + 1):length(levels(factor(conditions)))) {
    if (i != j) {
      contrast_name <- paste(levels(factor(conditions))[i], "-", levels(factor(conditions))[j])
      all_comb[[k]] <- makeContrasts(contrasts = contrast_name, levels = design)  # 대비 설정
      comb_name[[k]] <- paste(levels(factor(conditions))[i], "vs", levels(factor(conditions))[j])  # 대비 이름 설정
      k <- k + 1
    }
  }
}

# 올바른 대비 설정 확인
print(all_comb)

# 결과 테이블 리스트 초기화
result.tables <- list()
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]  # 비교 이름 가져오기
  contrast <- all_comb[[k]]  # 대비 매트릭스 가져오기
  
  if (any(is.na(contrast))) next  # 유효하지 않은 대비 건너뛰기
  
  lrt <- glmLRT(fit, contrast = contrast)  # LRT 수행
  res <- as.data.frame(topTags(lrt, n = nrow(dge)))  # 결과 추출
  res$contrast <- comparison  # 결과에 대비 이름 추가
  
  # 일관된 열 이름 확인
  colnames(res) <- make.names(colnames(res), unique = TRUE)
  
  result.tables[[k]] <- res  # 결과 리스트에 저장
}

# 결합하기 전에 열 이름 검사
for (i in seq_along(result.tables)) {
  print(colnames(result.tables[[i]]))
}

# 모든 결과를 하나의 데이터 프레임으로 결합
combined_results <- do.call(rbind, result.tables)  # 모든 결과를 하나의 데이터 프레임으로 결합

# 결합된 결과 저장
write.table(combined_results, file = file.path(root_dir, paste(project_name, "_combined_DE_results.txt", sep = "")), sep = "\t", row.names = FALSE, quote = FALSE)  # 결합된 결과 저장

# 샘플 품질을 위한 MDS 플롯
png(file = file.path(root_dir, paste(project_name, '_MDS_all_samples.png', sep = '')), width = 2000, height = 2000, res = 300)  # MDS 플롯을 PNG 파일로 저장
plotMDS(dge, col = as.numeric(factor(conditions)), pch = 16, cex = 1.5, main = "MDS Plot of all samples")  # MDS 플롯 생성
dev.off()  # PNG 파일 닫기


# 히트맵
if (exists("tt.rpkm")) {
  rpkm.avg <- round(rowMeans(tt.rpkm), 2)  # RPKM 평균 계산
  
  if (length(unique(conditions)) > 1) {
    pheatmap(mat = tt.rpkm, show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, scale = 'row')  # 히트맵 생성
  } else {
    cat("At least two groups are required for heatmap visualization.\n")  # 그룹이 2개 이상 필요하다는 메시지 출력
  }
} else {
  cat("Object 'tt.rpkm' not found.\n")  # tt.rpkm 객체가 없다는 메시지 출력
}

# 차등 발현 분석에서 상위 20개 유전자 추출
sigOE_ordered <- result.tables[[1]][order(result.tables[[1]]$PValue), ] # p-value로 정렬
top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ]) # 상위 20개 유전자 추출

# cpm(백만 당 카운트)을 사용하여 정규화된 카운트 데이터 추출
normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE) # 정규화된 카운트 계산

# 상위 20개 유의미한 차등 발현 유전자의 정규화된 카운트 추출
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]

# 데이터 프레임 변환
melted_top20_sigOE <- melt(top20_sigOE_norm) # 데이터 프레임 변환

# 열 이름 변경
colnames(melted_top20_sigOE) <- c("gene", "samplename", "normalized_counts") # 열 이름 설정

# 메타데이터 추가
meta <- data.frame(sampletype = conditions, samplename = colnames(normalized_counts)) # 메타데이터 생성
melted_top20_sigOE <- merge(melted_top20_sigOE, meta, by = "samplename") # 메타데이터 병합

# 플롯 그리기
ggplot(melted_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) + # 포인트 플롯 생성
  scale_y_log10() + # y축을 로그 스케일로 설정
  xlab("Genes") + # x축 레이블 설정
  ylab("Normalized Counts") + # y축 레이블 설정
  ggtitle("Top 20 Significant DE Genes") + # 플롯 제목 설정
  theme_bw() + # 테마 설정
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # x축 텍스트 회전
  theme(plot.title = element_text(hjust = 0.5)) # 제목 가운데 정렬

# 그룹 비교를 위한 조건 지정
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")
group1 <- which(conditions == "Condition1") # 그룹 1을 Condition1로 지정
group2 <- which(conditions == "Condition2") # 그룹 2를 Condition2로 지정

# 네 개의 플롯을 표시하도록 레이아웃 설정
par(mfrow = c(2, 2)) # 네 개의 플롯을 표시하도록 레이아웃 설정

# 모든 샘플에 대한 평균 차이 플롯
plotMD(dge, column = 1:ncol(dge$counts), main = "Total") # 모든 샘플에 대한 평균 차이 플롯 생성
abline(h = 0, col = "grey") # 중간값을 나타내는 회색 수평선 추가

# 그룹 1에 대한 평균 차이 플롯
plotMD(dge, column = group1, main = "Group 1: Condition1") # 그룹 1에 대한 평균 차이 플롯 생성
abline(h = 0, col = "blue") # 중간값을 나타내는 파란색 수평선 추가

# 그룹 2에 대한 평균 차이 플롯
plotMD(dge, column = group2, main = "Group 2: Condition2") # 그룹 2에 대한 평균 차이 플롯 생성
abline(h = 0, col = "red") # 중간값을 나타내는 빨간색 수평선 추가

# 기본 레이아웃으로 재설정
par(mfrow = c(1, 1)) # 기본 레이아웃으로 재설정

#Reference link : https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Quality_control

