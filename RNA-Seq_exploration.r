if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


BiocManager::install("airway")
BiocManager::install("SummarizedExperiment")

install.packages("gridExtra")

library(airway)
library(SummarizedExperiment)
library(GenomicRanges) # Needed for the sum() operation on ranges
library(ggplot2)
library(gridExtra)

data("airway")

raw_counts <- assay(airway)

# Extract Gene Lengths, The rowRanges are a 'GRangesList' (lists of exons per gene).
# calculate the width of each exon then SUM them to get total gene length.
# sum() gives one number per gene.
gene_ranges <- rowRanges(airway)
gene_lengths <- sum(width(gene_ranges))


# RPK (Reads Per Kilobase), Now gene_lengths is a simple numeric vector.
rpk <- raw_counts / (gene_lengths / 1000)

# Scaling Factor per sample
scaling_factors <- colSums(rpk) / 1e6

tpm <- sweep(rpk, 2, scaling_factors, "/")

sample_data <- data.frame(TPM = tpm[, 1])
sample_data$Log2_TPM <- log2(sample_data$TPM + 1)

p1 <- ggplot(sample_data, aes(x = TPM)) +
  geom_histogram(fill = "#E74C3C", color = "white", bins = 60) +
  scale_x_continuous(limits = c(0, 500)) +
  theme_minimal() +
  labs(title = "Raw TPM (Right-Skewed)", subtitle = "Original Data")

p2 <- ggplot(sample_data, aes(x = Log2_TPM)) +
  geom_histogram(fill = "#3498DB", color = "white", bins = 60) +
  theme_minimal() +
  labs(title = "Log2(TPM + 1)", subtitle = "Transformed")

p1

"""Zero-Inflation"""

p2

data("airway")


raw_counts <- assay(airway)
gene_ranges <- rowRanges(airway)
gene_lengths <- sum(width(gene_ranges))

rpk <- raw_counts / (gene_lengths / 1000)

#Per Million scaling
scaling_factors <- colSums(rpk) / 1e6

tpm <- sweep(rpk, 2, scaling_factors, "/")

sample_data <- data.frame(
  TPM = tpm[, 1]
)

sample_data$Log2_TPM <- log2(sample_data$TPM + 1)


# PLOT A: All Genes (Includes Zeros)
# This shows the "Zero-Inflation" (the spike on the left)
p1 <- ggplot(sample_data, aes(x = Log2_TPM)) +
  geom_histogram(fill = "#3498DB", color = "white", bins = 50) +
  theme_minimal() +
  labs(title = "A. All Genes (Includes Zeros)",
       subtitle = "Zero-inflation hides the normal distribution",
       x = "Log2(TPM + 1)", y = "Count")

# PLOT B: Expressed Genes Only (TPM > 1)
# Filter the data to remove non-expressed genes
# Filter for ANY expression (TPM > 0) to see the full tail
expressed_genes_full <- sample_data[sample_data$TPM > 0, ]

ggplot(expressed_genes_full, aes(x = Log2_TPM)) +
  geom_histogram(aes(y = ..density..), fill = "#BDC3C7", color = "white", bins = 60, alpha = 0.6) +

  # The Density Line (The Smooth "Bell")
  geom_density(color = "#E74C3C", size = 1.2) +

  # Add a theoretical Normal Distribution curve for comparison (Blue Dashed)
  stat_function(fun = dnorm,
                args = list(mean = mean(expressed_genes_full$Log2_TPM),
                            sd = sd(expressed_genes_full$Log2_TPM)),
                color = "#3498DB", linetype = "dashed", size = 1) +

  theme_minimal() +
  labs(title = "Distribution of Non-Zero Genes (TPM > 0)",
       subtitle = "Red = Real Data | Blue Dashed =  Normal Curve",
       x = "Log2(TPM + 1)", y = "Density")

grid.arrange(p1, p2, ncol = 2)

active_genes <- sample_data[sample_data$TPM > 1, ]

ggplot(active_genes, aes(x = Log2_TPM)) +
  geom_histogram(aes(y = ..density..), fill = "#2ECC71", color = "white", bins = 40, alpha = 0.6) +

  geom_density(color = "#27AE60", size = 1.5) +

  theme_minimal() +
  labs(title = "Distribution of Active Genes (TPM > 1)",
       subtitle = "Removing the 'noise' reveals the true Log-Normal distribution",
       x = "Log2(TPM + 1)", y = "Density")
grid.arrange(p1, p2, ncol = 2)

if (!require("DESeq2")) BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)

data("airway")
dds <- DESeqDataSet(airway, design = ~ cell)

#Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

vst_data <- assay(vsd)
sample_1_vst <- data.frame(Expression = vst_data[, 1])


# Filter VST data to remove the bottom 1% (essentially the zeros)
threshold <- quantile(sample_1_vst$Expression, 0.2) # Remove bottom 20% (zeros + noise)
vst_filtered <- sample_1_vst[sample_1_vst$Expression > threshold, , drop=FALSE]

ggplot(vst_filtered, aes(x = Expression)) +
  geom_histogram(aes(y = ..density..), fill = "#8E44AD", color = "white", bins = 40) +
  geom_density(color = "#4A235A", size = 1.2) +
  theme_minimal() +
  labs(title = "VST Data (Active Genes Only)",
       subtitle = "Removing the 'zero spike' finally reveals the smooth distribution",
       x = "VST Expression Value", y = "Density")

# intgroup: The columns in metadata to color/shape by
# 'dex' is the treatment column in the airway dataset (untreated vs treated)
pcaData <- plotPCA(vsd, intgroup=c("dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=dex, shape=cell)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw() +
  labs(title = "PCA Plot of Airway Data",
       subtitle = "Samples cluster by Treatment (Color) and Cell Line (Shape)")
