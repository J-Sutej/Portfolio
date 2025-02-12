
# setwd("C:/Users/Victus/Desktop/MS-Bioin/Project")

# Install required library
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("airway")

# script to get data from airway package
library(airway)
data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)] # filter column
sample_info$dex <- gsub('trt', 'treated', sample_info$dex) # change trt to treated
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex) 
names(sample_info) <- c('cellLine', 'dexamethasone') # change column head to cellline and dexamethasone
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

##############################################################################

# script to perform differential gene expression analysis using DESeq2 package
# install DESeq2 library
# BiocManager::install("DESeq2")
# load libraries
library(airway)
library(DESeq2)
library(tidyverse)

#######################################################################
## Step 1: Preparing count data
# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)
# counts_data contains gene counts in each sample. Each column represents the value in each sample test.

# read in sample info
colData <- read.csv('sample_info.csv')

# check column name in counts_data matching to row names in colData, respectively
all(colnames(counts_data) %in% rownames(colData))

# check the order of colnames and rownames
all(colnames(counts_data) == rownames(colData))

####################################################################
## Step 2: construct a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)
dds
counts(dds)
# pre-filtering: removing rows with low gene counts < 10 counts
# keeping rows having at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set factor level 
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

######################################################################
# Step 3: Run DESeq
dds <- DESeq(dds) 
# LFC > 0 (up) means how many genes are upregulated
res <- results(dds)

res

# explore results
summary(res) # having 0.23% outliers 
res0.01 <- results(dds, alpha = 0.01) # set p-value = 0.01 and see the result
summary(res0.01)

# contrasts
resultsNames(dds) # make a comparison to dexamethasone_untreated_vs_treated
#############################################################################
# Visualization
# MA plot
plotMA(res)
# Interpret: the blue dot is a significant gene and the dots above 0 are upregulated gene

# Heat map
# Transform data
vsd <- vst(dds, blind = FALSE)

# Select top 50 genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

# Generate heatmap
library(pheatmap)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)))
# 
