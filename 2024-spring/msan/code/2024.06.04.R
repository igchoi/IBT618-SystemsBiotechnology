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
rgid <- unlist(elementMetadata(gffrgene)[,"Parent"])  # RNA gene의 Parent ID를 rgid 변수에 저장
rna.ind <- match(rgid, gid)  # rgid를 gid에서 매칭하여 인덱스 저장

# Initialize matrices for counts
all.counts <- matrix(0, nrow=length(gid), ncol=length(bamFls))  # 모든 gene의 카운트를 저장할 행렬 초기화
rownames(all.counts) <- glocus_tag  # 행 이름을 glocus_tag로 설정
all.rcounts <- matrix(0, nrow=length(rgid), ncol=length(bamFls))  # RNA gene의 카운트를 저장할 행렬 초기화
colnames(all.counts) <- sub(".bam", "", basename(bamFls))  # 열 이름을 BAM 파일 이름으로 설정
colnames(all.rcounts) <- sub(".bam", "", basename(bamFls))  # 열 이름을 BAM 파일 이름으로 설정

# Read alignments and count overlaps
for (i in seq_along(bamFls)) {
  aln <- readGAlignments(bamFls[i])  # BAM 파일에서 정렬된 읽기를 가져옴
  seqlevels(aln) <- seqlevels(gff)  # 정렬된 읽기의 시퀀스 레벨을 gff의 시퀀스 레벨로 설정
  counts <- countOverlaps(gffgene, aln, ignore.strand=TRUE)  # gene에 대한 카운트 계산
  counts.r <- countOverlaps(gffrgene, aln, ignore.strand=TRUE)  # RNA gene에 대한 카운트 계산
  all.counts[, i] <- counts  # 계산된 카운트를 all.counts 행렬에 저장
  all.rcounts[, i] <- counts.r  # 계산된 RNA 카운트를 all.rcounts 행렬에 저장
}

# Coverage table preparation
coverage_table <- data.frame(
  name=sub(".bam", "", basename(bamFls)),  # 파일 이름 설정
  all_counts=colSums(all.counts),  # 모든 gene의 총 카운트
  ncrna_counts=colSums(all.rcounts),  # ncRNA의 총 카운트
  ncrna_ratio=round(100 * colSums(all.rcounts) / colSums(all.counts), 2)  # ncRNA 비율 계산
)

# Write the coverage table
write.table(coverage_table, file.path(root_dir, paste(project_name, '_cov_table.txt', sep='')), sep='\t', row.names=FALSE, quote=FALSE)  # 커버리지 테이블을 파일로 저장

# Remove rRNA genes counts for DEG analysis
all.gcounts <- if (length(rna.ind) > 0) all.counts[-rna.ind, ] else all.counts  # rRNA gene 카운트를 제거한 gene 카운트 행렬

# Create DGEList object
dge <- DGEList(counts = all.gcounts)  # DGEList 객체 생성
dge <- calcNormFactors(dge)  # 정규화 인자 계산

# Save gene counts matrix
save(all.gcounts, file=file.path(root_dir, "all.gcounts.RData"))  # gene 카운트 행렬 저장

# Sample grouping
sample_names <- sub(".bam", "", basename(bamFls))  # 샘플 이름 설정
conditions <- c("Condition1", "Condition1", "Condition1", "Condition2", "Condition2", "Condition2")  # 조건 설정
design <- model.matrix(~0 + factor(conditions))  # 실험 설계 행렬 생성
colnames(design) <- levels(factor(conditions))  # 열 이름 설정

# Estimate dispersion
dge <- estimateDisp(dge, design)  # 이분산성 추정

# Create model and perform DE analysis
fit <- glmFit(dge, design)  # 일반화 선형 모델 생성 및 분석

# RPKM calculation - gene set subtraction
gsize <- width(gffgene)  # gene의 길이 계산
names(gsize) <- glocus_tag  # gene 길이 이름 설정
gsize.c <- gsize[rownames(dge$counts)]  # DGEList 객체의 gene 길이
geneLengthsInKB <- (gsize.c / 1000)  # gene 길이를 킬로베이스 단위로 변환
millionsMapped <- colSums(dge$counts) / 1e+06  # 매핑된 리드 수를 백만 단위로 변환
rpm <- dge$counts / millionsMapped  # 분당 리드 수 계산
tt.rpkm <- round(rpm / geneLengthsInKB, 1)  # RPKM 값 계산
colnames(tt.rpkm) <- sub('\\$', '_rpkm', colnames(tt.rpkm))  # 열 이름 설정
keep <- rowSums(tt.rpkm >= 2) >= 1  # RPKM이 2 이상인 유전자 유지

# Initialize matrices for the result table
tt.raw.counts <- dge$counts  # 원시 카운트 데이터
colnames(tt.raw.counts) <- sub('\\$', '_raw_counts', colnames(tt.raw.counts))  # 열 이름 설정

# Pseudocount
dge.commondisp <- estimateCommonDisp(dge, design)  # 공통 분산 추정
tt.pseudocount <- dge.commondisp$pseudo.counts  # 가상 카운트 계산
tt.pseudocount.round <- round(tt.pseudocount)  # 가상 카운트 반올림
colnames(tt.pseudocount.round) <- sub('\\$', '_pseudocount', colnames(tt.pseudocount.round))  # 열 이름 설정

# Coefficient from fit object (note that it is log2 value)
tt.coeff <- fit$coefficients  # 회귀 계수 추출
tt.coeff.log2 <- tt.coeff / log(2)  # log2 값으로 변환
tt.coeff.log2.round <- round(tt.coeff.log2, 2)  # 반올림
colnames(tt.coeff.log2.round) <- sub('\\$', '_coeff', colnames(tt.coeff.log2.round))  # 열 이름 설정

# Differential expression analysis setup
all_comb <- list()
comb_name <- list()
k <- 1
for (i in seq_along(levels(factor(conditions)))) {
  for (j in (i + 1):length(levels(factor(conditions)))) {
    all_comb[[k]] <- makeContrasts(contrasts = levels(factor(conditions))[c(i, j)], levels = design)  # 대조군 설정
    comb_name[[k]] <- paste(levels(factor(conditions))[c(i, j)], collapse = "vs")  # 대조군 이름 설정
    k <- k + 1
  }
}

# Example of setting contrast manually
print(all_comb[[2]])
all_comb[[2]] <- c(1, -1, 0)
all_comb <- all_comb[-2]
comb_name <- comb_name[-2]

# Initialize result tables list
result.tables <- list()
for (k in seq_along(all_comb)) {
  comparison <- comb_name[[k]]  # 비교군 이름
  contrast <- all_comb[[k]]  # 대조군 설정
  
  if (any(is.na(contrast))) next  # 대조군이 유효하지 않으면 다음으로 건너뛰기
  
  lrt <- glmLRT(fit, contrast = contrast)  # 대조군에 대해 LRT 수행
  res <- as.data.frame(topTags(lrt, n = nrow(dge)))  # 결과 추출
  res$contrast <- comparison  # 결과에 대조군 이름 추가
  
  result.tables[[k]] <- res  # 결과를 리스트에 저장
}

# Save result tables
for (i in seq_along(result.tables)) {
  write.table(result.tables[[i]], file = file.path(root_dir, paste(project_name, comb_name[[i]], "_DE_results.txt", sep = "")), sep = "\t", row.names = FALSE, quote = FALSE)  # 결과 테이블을 파일로 저장
}

# Combine all results into one dataframe
combined_results <- do.call(rbind, result.tables)  # 모든 결과를 하나의 데이터 프레임으로 결합

# Save combined results
write.table(combined_results, file = file.path(root_dir, paste(project_name, "_combined_DE_results.txt", sep = "")), sep = "\t", row.names = FALSE, quote = FALSE)  # 결합된 결과를 파일로 저장

# MDS plot for sample quality
png(file = file.path(root_dir, paste(project_name, '_MDS_all_samples.png', sep = '')), width = 2000, height = 2000, res = 300)  # MDS 플롯을 PNG 파일로 저장
plotMDS(dge, col = as.numeric(factor(conditions)), pch = 16, cex = 1.5, main = "MDS Plot of all samples")  # MDS 플롯 생성
dev.off()  # PNG 파일 닫기

# Heatmap
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

# Differential expression analysis에서 상위 20개 유전자 추출
sigOE_ordered <- result.tables[[1]][order(result.tables[[1]]$PValue), ]  # p-value 기준으로 정렬
top20_sigOE_genes <- rownames(sigOE_ordered[1:20, ])  # 상위 20개 유전자 추출

# 정규화된 카운트 데이터 추출 using cpm (counts per million)
normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)  # 정규화된 카운트 계산

# Extract normalized counts for the top 20 significant DE genes
top20_sigOE_norm <- normalized_counts[top20_sigOE_genes, ]  # 상위 20개 유전자에 대한 정규화된 카운트 추출

# 데이터 프레임을 변형
melted_top20_sigOE <- data.frame(melt(top20_sigOE_norm))  # 데이터 프레임 변형

# 열 이름 변경
colnames(melted_top20_sigOE) <- c("gene", "samplename", "normalized_counts")  # 열 이름 설정

# 메타데이터 추가
meta <- data.frame(sampletype = conditions, samplename = sample_names)  # 메타데이터 생성
melted_top20_sigOE <- merge(melted_top20_sigOE, meta, by = "samplename")  # 메타데이터 병합

# 플롯 그리기
ggplot(melted_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) +  # 포인트 플롯 생성
  scale_y_log10() +  # y 축을 로그 스케일로 설정
  xlab("Genes") +  # x 축 라벨 설정
  ylab("Normalized Counts") +  # y 축 라벨 설정
  ggtitle("Top 20 Significant DE Genes") +  # 플롯 제목 설정
  theme_bw() +  # 테마 설정
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # x 축 텍스트 회전
  theme(plot.title = element_text(hjust = 0.5))  # 제목 중앙 정렬

# Ensure necessary libraries are loaded
library(ggplot2)  # ggplot2 패키지를 로드
library(ggrepel)  # ggrepel 패키지를 로드

# Order the data frame by adjusted p-value
resOE_df_ordered <- resOE_df[order(resOE_df$FDR), ]  # 조정된 p-value 기준으로 데이터 프레임 정렬

# Create a new column 'genelabels' based on whether the row names are in the top 10 rows
resOE_df_ordered$genelabels <- rownames(resOE_df_ordered) %in% rownames(resOE_df_ordered[1:10, ])  # 상위 10개 유전자에 대한 라벨 추가

# Add a threshold column for volcano plot coloring
resOE_df_ordered$threshold <- resOE_df_ordered$FDR < 0.05  # p-value 기준으로 임계값 컬럼 추가

# Set the floor value for the y-axis
y_floor <- 300  # y축의 바닥값 설정 (데이터에 따라 조정 필요)

# Create a volcano plot
ggplot(resOE_df_ordered) +
  geom_point(aes(x = logFC.Condition1, y = pmax(-log10(FDR), -y_floor), colour = threshold)) +  # Volcano plot 생성
  geom_text_repel(aes(x = logFC.Condition1, y = pmax(-log10(FDR), -y_floor), label = ifelse(genelabels, rownames(resOE_df_ordered), "")), vjust = 1, size = 3) +  # 상위 10개 유전자 라벨 추가
  ggtitle("Mov10 overexpression") +  # 플롯 제목 설정
  xlab("log2 fold change") +  # x축 라벨 설정
  ylab("-log10 adjusted p-value") +  # y축 라벨 설정
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),  # 제목 중앙 정렬 및 크기 설정
        axis.title = element_text(size = rel(1.25)))  # 축 제목 크기 설정

# Barplot of library sizes
barplot(resOE_df_ordered$logCPM, names = rownames(resOE_df_ordered), las = 2,
        main = "Barplot of library sizes", xlab = "Samples", ylab = "Library size")  # 라이브러리 크기에 대한 막대 그래프 생성

# Adjusted labeling
barplot(resOE_df_ordered$logCPM / 1e06, names = rownames(resOE_df_ordered), las = 2, ann = FALSE, cex.names = 0.75,
        main = "Barplot of library sizes", xlab = "Samples", ylab = "Library size (millions)")  # 라벨 조정 후 막대 그래프 생성

# Get log2 counts per million
logcounts <- resOE_df_ordered$logCPM  # logCPM 값 추출
logcounts[is.infinite(logcounts)] <- NA  # 무한대 값을 NA로 대체
logcounts[is.nan(logcounts)] <- NA  # NaN 값을 NA로 대체
logcounts <- logcounts[!is.na(logcounts)]  # NA 값을 제거
logcounts <- logcounts[logcounts > 0]  # 음수 값을 제거
logcounts <- log2(logcounts + 1)  # log2 변환

# Boxplot of logCPMs (unnormalized)
boxplot(logcounts, xlab = "Samples", ylab = "Log2 counts per million", las = 2,
        main = "Boxplots of logCPMs (unnormalized)")  # 정규화되지 않은 logCPM에 대한 박스 플롯 생성

# Add a blue horizontal line representing the median logCPM
abline(h = median(logcounts, na.rm = TRUE), col = "blue")  # 중앙값을 나타내는 파란색 수평선 추가

# Set the layout to display four plots
par(mfrow=c(2,2))  # 4개의 플롯을 표시하기 위해 레이아웃 설정

# Plot mean difference plot for all samples
plotMD(dge, column = 1:ncol(dge$counts))  # 모든 샘플에 대한 평균 차이 플롯 생성
abline(h=0, col="grey")  # 중앙값을 나타내는 회색 수평선 추가

# Plot mean difference plot for group 1 vs group 2
plotMD(dge, column = group1)  # 그룹 1에 대한 평균 차이 플롯 생성
abline(h=0, col="blue")  # 중앙값을 나타내는 파란색 수평선 추가
plotMD(dge, column = group2)  # 그룹 2에 대한 평균 차이 플롯 생성
abline(h=0, col="red")  # 중앙값을 나타내는 빨간색 수평선 추가

# Reset the layout to default
par(mfrow=c(1,1))  # 레이아웃을 기본값으로 재설정

Reference link : https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Quality_control