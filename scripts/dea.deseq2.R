library(DESeq2)

# load the count file
args <- commandArgs(TRUE)
countTable <- read.table(args[1], header = T)

# countTable <- read.table('./control_PAH_high', header = T)
rownames(countTable) <- countTable[,1]
countData <- countTable[,-1] # remove Gene_id column to get pure count table

control <- args[2]
treat <- args[3]

colData <- data.frame('condition' = factor(c(rep(1, control), rep(2, treat))))

# creat the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# Filtering
keep <- rowSums(counts(dds)) >= 10  # only keep rows that have at least 10 counts in total
dds <- dds[keep,]

# perform DEA
dea <- DESeq(dds)

# result analysis
dea <- results(dea)

# p-adj
padj <- dea$padj

# significantly differently expressed genes
deg <- dea[which(padj < 0.05), ]

# save the DEA result and DEGs to files
write.csv(dea, args[4])
write.csv(deg, args[5])


