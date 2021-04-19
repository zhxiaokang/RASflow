library(yaml)
library(edgeR)
library(DESeq2)

# ====================== define the function of DEA ======================

DEA <- function(control, treat) {
  count.control <- read.table(paste(output.path, '/countGroup/', control, '_gene_count.tsv', sep = ''), header = TRUE, row.names = 1)
  count.treat <- read.table(paste(output.path, '/countGroup/', treat, '_gene_count.tsv', sep = ''), header = TRUE, row.names = 1)
  count.table <- cbind(count.control, count.treat)  # merge the control and treat tables together
  
  # number of samples in control and treat groups (should be the same if it's a pair test)
  num.sample <- ncol(count.table)
  num.sample.control <- ncol(count.control)
  num.sample.treat <- ncol(count.treat)

  # samples of two groups
  sample.control <- colnames(count.control)
  sample.treat <- colnames(count.treat)
  
  # save gene list in gene.list for extracting gene names later
  gene.list <- rownames(count.table)
  
  # get the sample id
  samples <- colnames(count.table)

  # define the subject and group; define the reference group (control)
  subject <- factor(subject.all[c(which(group.all == control), which(group.all == treat))])
  group <- relevel(factor(group.all[c(which(group.all == control), which(group.all == treat))]), ref = control)

  # The design matrix
  if (pair.test) {
    design <- model.matrix(~subject+group)
  } else {
    design <- model.matrix(~group)
  }
 
  if (dea.tool == 'edgeR') {  # use edgeR for DEA
    # normalize the two groups and save the normalized count table
    y.control <- DGEList(counts = count.control, genes = gene.list)
    y.treat <- DGEList(counts = count.treat, genes = gene.list)

    y.control <- calcNormFactors(y.control, method="TMM")
    count.table.control.norm <- cpm(y.control)
    write.table(count.table.control.norm, paste(output.path, '/countGroup/', control, '_gene_norm.tsv', sep = ''), quote = FALSE, sep = "\t")

    y.treat <- calcNormFactors(y.treat, method="TMM")
    count.table.treat.norm <- cpm(y.treat)
    write.table(count.table.treat.norm, paste(output.path, '/countGroup/', treat, '_gene_norm.tsv', sep = ''), quote = FALSE, sep = "\t")

    # Put the data into a DGEList object
    y <- DGEList(counts = count.table, genes = gene.list)
    
    # do DEA

    # Filtering
    if (filter.need) {
      countsPerMillion <- cpm(y)
      countCheck <- countsPerMillion > 1
      keep <- which(rowSums(countCheck) > 1)
      y <- y[keep, ]
    }
    
    # Normalization
    y <- calcNormFactors(y, method="TMM")
    
    y$samples$subject <- subject
    y$samples$group <- group
    
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
    dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

    # differentially expressed genes
    toptag <- topTags(lrt, n = nrow(y$genes), p.value = fdr.thr)
    deg <- toptag$table
    if (!is.null(deg)) {
      deg <- deg[order(deg$FDR, -abs(deg$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    }

    # save the DEA result and DEGs to files
    write.table(dea, paste(output.path, '/DEA/dea_', control, '_', treat, '.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')
    write.table(deg, paste(output.path, '/DEA/deg_', control, '_', treat, '.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t') 
  } else if (dea.tool == "DESeq2") {  # use DESeq2 for DEA

    ## create the DESeqDataSet
    colData = data.frame(samples, subject, group)
    dds <- DESeqDataSetFromMatrix(count.table, colData = colData, design = design)

    # generate normalized counts
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized=TRUE)

    normalized_counts.control <- normalized_counts[, which(colnames(normalized_counts) %in% sample.control)]
    write.table(normalized_counts.control, paste(output.path, '/countGroup/', control, '_gene_norm.tsv', sep = ''), quote = FALSE, sep = "\t")

    normalized_counts.treat <- normalized_counts[, which(colnames(normalized_counts) %in% sample.treat)]
    write.table(normalized_counts.treat, paste(output.path, '/countGroup/', treat, '_gene_norm.tsv', sep = ''), quote = FALSE, sep = "\t")

    ## Filtering
    if (filter.need) {
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep,]
    }
    
    ## perform DEA
    dds <- DESeq(dds)
    
    ## export the results
    res.dea <- results(dds)
    res.dea <- res.dea[complete.cases(res.dea), ]  # remove any rows with NA

    dea <- as.data.frame(res.dea)
    dea <- dea[order(dea$padj, -abs(dea$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    
    deg <- dea[dea$padj < fdr.thr, ]
    if (nrow(deg) > 1) {
      deg <- deg[order(deg$padj, -abs(deg$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
    }

    # save the DEA result and DEGs to files
    write.table(dea, paste(output.path, '/DEA/dea_', control, '_', treat, '.tsv', sep = ''), row.names = T, quote = FALSE, sep = '\t')
    write.table(deg, paste(output.path, '/DEA/deg_', control, '_', treat, '.tsv', sep = ''), row.names = T, quote = FALSE, sep = '\t')
  }
}

# load the config file
yaml.file <- yaml.load_file('configs/config_main.yaml')

# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name of this analysis
dea.tool <- yaml.file$DEATOOL  # tool used for DEA
controls <- yaml.file$CONTROL  # all groups used as control
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
filter.need <- yaml.file$FILTER$yesOrNo
pair.test <- yaml.file$PAIR
fdr.thr <- yaml.file$FDR  # threshold of FDR/adjusted P-value for significantlly differentially expressed genes
meta.file <- yaml.file$METAFILE
output.path <- file.path(yaml.file$FINALOUTPUT, project, "genome/dea")

# extract the metadata
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t')
group.all <- meta.data$group
subject.all <- meta.data$subject

num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

# Do DEA

# the main function
for (ith.comparison in c(1:num.comparison)) {
  control <- controls[ith.comparison]
  treat <- treats[ith.comparison]
  DEA(control, treat)
}
