library(edgeR)

# Use edgeR to do DEA using paired test

# load the count file
args <- commandArgs(TRUE)
countTable <- read.table(args[1], header = T)

rownames(countTable) <- countTable[,1]
countData <- countTable[,-1] # remove Gene_id column to get pure count table

# save gene list in geneList for extracting gene names later
geneList <- rownames(countData)

# get the sample id
samples <- colnames(countData)

numDuplicate <- length(samples)/2  # pair set-up, so numControl = numTreat = numDuplicate

subject <- factor(c(seq(1, numDuplicate), seq(1, numDuplicate)))
groups <- factor(c(rep("C", numDuplicate), rep("T", numDuplicate)))

# Put the data into a DGEList object
y <- DGEList(counts = countTable[,2:(length(samples)+1)], genes = countTable[,1])

## Filtering
# countsPerMillion <- cpm(y)
# countCheck <- countsPerMillion > 1
# keep <- which(rowSums(countCheck) >= 2)
# y <- y[keep, ]

# Normalization
y <- calcNormFactors(y, method="TMM")

# save the normalized count table
countData.norm <- cpm(y)
write.csv(countData.norm, args[2])

# define the group
y$samples$subject <- subject
y$samples$group <- groups

# The design matrix
design <- model.matrix(~subject+groups)
rownames(design) <- colnames(y)

# Estimating the dispersion

# estimate the NB dispersion for the dataset
y <- estimateDisp(y, design, robust = TRUE)

# Differential expression

# determine differentially expressed genes
# fit genewise glms
fit <- glmFit(y, design)

# conduct likelihood ratio tests for tumour vs normal tissue differences and show the top genes
lrt <- glmLRT(fit)

# the DEA result for all the genes
# dea <- lrt$table
toptag <- topTags(lrt, n = nrow(y$genes), p.value = 1)
dea <- toptag$table  # just to add one more column of FDR 

# differentially expressed genes
toptag <- topTags(lrt, n = nrow(y$genes), p.value = 0.05)
deg <- toptag$table

# save the DEA result and DEGs to files
write.csv(dea, args[3], row.names = F)
write.csv(deg, args[4], row.names = F)  # remove the row names (index number) just to make it consistent with dea
