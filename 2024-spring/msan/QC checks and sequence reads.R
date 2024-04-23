setwd("C:/Users/Desktop/data")
library(Rfastp)
json_report <- rfastp(read1 = "C:/Users/Desktop/data/SRR17001186.fastq.gz", outputFastq = "SRR17001186_fastp")
qcSummary(json_report)
curvePlot(json_report)
curvePlot(json_report, curve = "content_curves")



# 필요한 라이브러리를 불러옵니다.
library(Rfastp)

# Working directory 설정
setwd("C:/Users/Desktop/data")

# SRA accession numbers와 각각의 파일명을 정의합니다.
samples <- c("SRR170011796", "SRR17001179", "SRR17001180", "SRR17001182", "SRR17001183", "SRR17001186")
file_paths <- c("SRR17001176.fastq.gz",
                "SRR17001179.fastq.gz",
                "SRR17001180.fastq.gz",
                "SRR17001182.fastq.gz",
                "SRR17001183.fastq.gz",
                "SRR17001186.fastq.gz")

# QC 수행
for (i in seq_along(file_paths)) {
  # Output directory for QC results
  output_dir <- paste0(samples[i], "_fastp")
  
  # Perform QC using Rfastp
  json_report <- rfastp(
    read1 = file_paths[i],
    outputFastq = output_dir
  )
  
  # Generate QC summary and curve plots
  qcSummary(json_report)
  curvePlot(json_report)
  curvePlot(json_report, curve = "content_curves")
  
  cat("QC completed for sample:", samples[i], "\n")
}

















# Load necessary libraries
library(kallisto)
library(tximport)

# Set kallisto output directory
output_dir_kallisto <- "kallisto_results"

# Run kallisto quantification
kallisto::kallisto(
  index = "SRR17001180.fasta.gz",  # Use the provided fasta.gz file
  output_dir = output_dir_kallisto,
  single = input_fastq
)

# Full path to the abundance.tsv file
abundance_file <- file.path(output_dir_kallisto, "abundance.tsv")

# Check if the abundance.tsv file exists
if (file.exists(abundance_file)) {
  # Run tximport
  txi <- tximport(
    abundance_file,
    type = "kallisto",
    tx2gene = NULL,
    countsFromAbundance = "scaledTPM"
  )
  
  # Calculate TPM values
  tpm_values <- txi$abundance / sum(txi$abundance) * 1e6
  
  # Add TPM values to the txi object
  txi$tpm <- tpm_values
  
  # Sort genes by TPM values
  txi_sorted <- txi[order(txi$tpm, decreasing = TRUE), ]
  
  # Extract gene names and TPM values
  gene_names <- rownames(txi_sorted)
  tpm_sorted <- txi_sorted$tpm
  
  # Create a data frame with gene names and TPM values
  tpm_df <- data.frame(Gene = gene_names, TPM = tpm_sorted)
  
  # Display the sorted TPM values
  print(tpm_df)
  
} else {
  stop("abundance.tsv file not found. Please check the kallisto output directory:", output_dir_kallisto)
}




# Set kallisto output directory
output_dir_kallisto <- "kallisto_results"

# Full path to kallisto.exe
kallisto_path <- "C:\\Users\\AN\\Desktop\\RU_RNAseq-master\\kallisto\\kallisto.exe"

# Create kallisto index
index_command <- paste(
  shQuote(kallisto_path),
  "index",
  "-i transcriptome.idx",
  shQuote("C:/Users/Desktop/data/Download/SRR17001180.fasta.gz"),
  sep = " "
)
system(index_command)

# Run kallisto quantification using system command
kallisto_command <- paste(
  shQuote(kallisto_path),
  "quant",
  "-i transcriptome.idx",
  "-o", shQuote(output_dir_kallisto),
  "--single", shQuote(input_fastq),
  "-l 200",  # Mean fragment length
  "-s 20",   # Fragment length standard deviation
  sep = " "
)

# Run the kallisto command and capture the output
kallisto_output <- system(kallisto_command, intern = TRUE)

# Check for errors
if (length(grep("Error", kallisto_output)) > 0) {
  stop("Error running kallisto:", kallisto_output)
}

# Rest of the code for tximport and further analysis remains the same



# Set kallisto output directory
output_dir_kallisto <- "kallisto_results"

# Full path to kallisto.exe
kallisto_path <- "C:\\Users\\Desktop\\data\\kallisto\\kallisto.exe"

# Run kallisto quantification using system command
kallisto_command <- paste(
  shQuote(kallisto_path),
  "quant",
  "-i transcriptome.idx",
  "-o", shQuote(output_dir_kallisto),
  "--single", shQuote(input_fastq),
  "-l 200",  # Mean fragment length
  "-s 20",   # Fragment length standard deviation
  sep = " "
)

# Run the kallisto command and capture the output
kallisto_output <- system(kallisto_command, intern = TRUE)

# Check for errors
if (length(grep("Error", kallisto_output)) > 0) {
  stop("Error running kallisto:", kallisto_output)
}

# Full path to the abundance.tsv file
abundance_file <- file.path(output_dir_kallisto, "abundance.tsv")

# Check if the abundance.tsv file exists
if (file.exists(abundance_file)) {
  # Run tximport
  txi <- tximport(
    abundance_file,
    type = "kallisto",
    tx2gene = NULL,
    countsFromAbundance = "scaledTPM"
  )
  
  # Calculate TPM values
  tpm_values <- txi$abundance / sum(txi$abundance) * 1e6
  
  # Add TPM values to the txi object
  txi$tpm <- tpm_values
  
  # Sort genes by TPM values
  txi_sorted <- txi[order(txi$tpm, decreasing = TRUE), ]
  
  # Extract gene names and TPM values
  gene_names <- rownames(txi_sorted)
  tpm_sorted <- txi_sorted$tpm
  
  # Create a data frame with gene names and TPM values
  tpm_df <- data.frame(Gene = gene_names, TPM = tpm_sorted)
  
  # Display the sorted TPM values
  print(tpm_df)
  
} else {
  stop("abundance.tsv file not found. Please check the kallisto output directory:", output_dir_kallisto)
}











# Full path to the abundance.tsv file
abundance_file <- file.path(output_dir_kallisto, "abundance.tsv")

# Check if the abundance.tsv file exists
if (file.exists(abundance_file)) {
  # Run tximport
  txi <- tximport(
    abundance_file,
    type = "kallisto",
    tx2gene = NULL,
    countsFromAbundance = "scaledTPM"
  )
  
  # Calculate TPM values
  tpm_values <- txi$abundance / sum(txi$abundance) * 1e6
  
  # Add TPM values to the txi object
  txi$tpm <- tpm_values
  
  # Sort genes by TPM values
  txi_sorted <- txi[order(txi$tpm, decreasing = TRUE), ]
  
  # Extract gene names and TPM values
  gene_names <- rownames(txi_sorted)
  tpm_sorted <- txi_sorted$tpm
  
  # Create a data frame with gene names and TPM values
  tpm_df <- data.frame(Gene = gene_names, TPM = tpm_sorted)
  
  # Display the sorted TPM values
  print(tpm_df)
  
} else {
  stop("abundance.tsv file not found. Please check the kallisto output directory:", output_dir_kallisto)
}

head(tpm_df, 10)

