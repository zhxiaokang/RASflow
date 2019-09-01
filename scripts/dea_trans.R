# Use edgeR to do DEA on the outputs from Salmon

library(biomaRt)
library(yaml)
library(edgeR)

# define the global variables

args = commandArgs(trailingOnly=TRUE)

## load the config file
yaml.file <- yaml.load_file('configs/config_dea_trans.yaml')
dataset <- yaml.file$EnsemblDataSet
meta.file <- yaml.file$METAFILE
meta.data <- read.csv(meta.file, sep = "\t")

samples <- args[1]
file.output <- args[2]
quant.file <- 'quant.sf'
quant.table <- read.table(quant.file, header = TRUE, stringsAsFactors = FALSE)

write.table(output, file = file.output)

for (i in c(1:nrow(quant.table))) {
  quant.table$Name[i] <- unlist(strsplit(quant.table$Name[i], split = '[.]'))[1]
}

ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
datasets <- listDatasets(ensembl)

attributes <- listAttributes(mart = ensembl)

trans.id <- quant.table[, 1]

tx2gene <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id'),
                           filters = 'ensembl_transcript_id', values = trans.id, mart = ensembl)

txi.tx <- tximport(quant.file, type = "salmon", txOut = TRUE)

txi <- tximport(quant.file, type = "salmon", tx2gene = tx2gene)

cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
library(edgeR)
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

