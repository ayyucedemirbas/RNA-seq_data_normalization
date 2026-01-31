if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("airway", "DESeq2", "vsn", "GenomicRanges"))

library(airway)
library(DESeq2)
library(vsn)
library(GenomicRanges)

data("airway")
se <- airway

raw_counts <- assay(se)

gene_lengths <- sum(width(reduce(rowRanges(se))))

# TPM = (counts / length) * 10^6 / sum(counts / length)
rpk <- raw_counts / gene_lengths
tpm <- t(t(rpk) * 1e6 / colSums(rpk))

# Log2(TPM + 1)
log2_tpm <- log2(tpm + 1)

dds <- DESeqDataSet(se, design = ~ cell + dex)
vst_data <- vst(dds, blind=TRUE)

rlog_data <- rlog(dds, blind=TRUE)

par(mfrow=c(1,3))

required_packages <- c("airway", "DESeq2", "vsn", "hexbin", "patchwork")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)) {
    BiocManager::install(new_packages)
}

library(ggplot2)
library(patchwork)

#Log2(TPM + 1)
p1 <- meanSdPlot(log2_tpm, plot=FALSE)$gg +
  ggtitle("Log2(TPM + 1)") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- meanSdPlot(assay(vst_data), plot=FALSE)$gg +
  ggtitle("VST (Raw Counts)") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- meanSdPlot(assay(rlog_data), plot=FALSE)$gg +
  ggtitle("rlog (Subset)") +
  theme(plot.title = element_text(hjust = 0.5))

p1

p2

p3

pca_vst <- plotPCA(vst_data, intgroup = c("dex", "cell")) +
  ggtitle("PCA with VST") +
  theme_minimal()


rv <- rowVars(log2_tpm)
pca_tpm_mat <- t(log2_tpm[rv, ])
pca_res <- prcomp(pca_tpm_mat)

pcaData <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  dex = colData(dds)$dex,
  cell = colData(dds)$cell
)

pca_log2 <- ggplot(pcaData, aes(PC1, PC2, color=dex, shape=cell)) +
  geom_point(size=3) +
  ggtitle("PCA with Log2(TPM+1) (Less Stable)") +
  theme_minimal()

pca_log2

pca_vst
