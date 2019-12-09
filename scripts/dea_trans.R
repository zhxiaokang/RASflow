# Use edgeR to do DEA on the outputs from Salmon

library(biomaRt)
library(yaml)
library(edgeR)
library(tximport)

# ====================== define the function of DEA ======================

DEA <- function(control, treat, file.control, file.treat, output.path.dea) {
  count.control <- read.table(file.control, header = TRUE, row.names = 1)
  count.treat <- read.table(file.treat, header = TRUE, row.names = 1)
  count.table <- cbind(count.control, count.treat)  # merge the control and treat tables together
  
  # number of samples in control and treat groups (should be the same if it's a pair test)
  num.sample.control <- ncol(count.control)
  num.sample.treat <- ncol(count.treat)
  
  # save gene list in gene.list for extracting gene names later
  gene.list <- rownames(count.table)
  
  # get the sample id
  samples <- colnames(count.table)
  
  # Put the data into a DGEList object
  y <- DGEList(counts = count.table, genes = gene.list)
  
  # Filtering
  if (filter.need) {
    tpm <- y$counts
    countCheck <- tpm > tpm.threshold
    keep <- which(rowSums(countCheck) > 1)
    y <- y[keep, ]
  }
  
  # Do DEA !!!
  # define the group
  subject <- factor(subject.all[c(which(group.all == control), which(group.all == treat))])
  group <- factor(group.all[c(which(group.all == control), which(group.all == treat))])
  
  y$samples$subject <- subject
  y$samples$group <- group
  
  # The design matrix
  if (pair.test) {
    design <- model.matrix(~subject+group)
  } else {
    design <- model.matrix(~group)
  }
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
  write.table(dea, paste(output.path.dea, '/dea_', control, '_', treat, '.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')
  write.table(deg, paste(output.path.dea, '/deg_', control, '_', treat, '.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')
}

# ====================== load parameters in config file ======================

# load the config file
yaml.file <- yaml.load_file('configs/config_dea_trans.yaml')

# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name
gene.level <- yaml.file$GENE_LEVEL  # whether to do gene-level analysis
controls <- yaml.file$CONTROL  # all groups used as control
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
filter.need <- yaml.file$FILTER$yesOrNo
tpm.threshold <- yaml.file$FILTER$tpm
pair.test <- yaml.file$PAIR
meta.file <- yaml.file$METAFILE
dataset <- yaml.file$EnsemblDataSet
output.path <- file.path(yaml.file$OUTPUTPATH, project, "trans/dea")

num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

# extract the metadata
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
samples <- meta.data$sample
group.all <- meta.data$group
subject.all <- meta.data$subject

# ====================== Do DEA ======================

for (ith.comparison in c(1:num.comparison)) {
  control <- controls[ith.comparison]
  treat <- treats[ith.comparison]
  
  # --------------------- On transctipt level ---------------------
  file.control <- paste(output.path, '/countGroup/', control, '_trans_norm.tsv', sep = '')
  file.treat <- paste(output.path, '/countGroup/', treat, '_trans_norm.tsv', sep = '')
  output.path.dea <- paste(output.path, '/DEA/transcript-level', sep = '')
  
  DEA(control, treat, file.control, file.treat, output.path.dea)
  
  # --------------------- On gene level ---------------------
  if (gene.level) {
    file.control <- paste(output.path, '/countGroup/', control, '_gene_norm.tsv', sep = '')
    file.treat <- paste(output.path, '/countGroup/', treat, '_gene_norm.tsv', sep = '')
    output.path.dea <- paste(output.path, '/DEA/gene-level', sep = '')
  
    DEA(control, treat, file.control, file.treat, output.path.dea)
  }
}

