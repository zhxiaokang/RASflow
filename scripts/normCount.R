# Normalise the count files according to the library size
library(openxlsx)
# Get an example count file
expl.count <- read.table('1.count')
names(expl.count) <- c('Gene_ID', 'length', 'mapped', 'unmapped')

# Define the data frame for storing the final normalised data
origin.count <- data.frame(gene.id = expl.count[-nrow(expl.count), 1])  # the last row is useless, and should be excluded

# Give the samples id
samples <- c(1:32, 36:53)

for (i in samples){
  count.file <- read.table(paste0(i, '.count'))
  origin.count[paste('sample', i, sep = '_')] <- count.file[-nrow(count.file), 3]  # the last row is useless, and should be excluded
}

# The library size for each sample
lib.size <- colSums(origin.count[, -1])

norm.count <- origin.count
for (i in c(2:ncol(norm.count))){
  norm.count[, i] <- origin.count[, i] / lib.size[i - 1] * 1e6
}

# sort the file according to the gene id
norm.count.sort <- norm.count[order(norm.count$gene.id),]

# save the file to an Excel file
write.xlsx(norm.count.sort, 'Normalised_Counts_ALL.xlsx', col.names = TRUE)

write.xlsx(norm.count.sort[,c(1, 2:5)], 'Normalised_Counts_A.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 6:9)], 'Normalised_Counts_B.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 10:13)], 'Normalised_Counts_C.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 14:17)], 'Normalised_Counts_D.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 18:21)], 'Normalised_Counts_E.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 22:25)], 'Normalised_Counts_F.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 26:29)], 'Normalised_Counts_G.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 30:33)], 'Normalised_Counts_H.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 34:37)], 'Normalised_Counts_J.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 38:40)], 'Normalised_Counts_K.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 41:43)], 'Normalised_Counts_L.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 44:47)], 'Normalised_Counts_M.xlsx', col.names = TRUE)
write.xlsx(norm.count.sort[,c(1, 48:51)], 'Normalised_Counts_N.xlsx', col.names = TRUE)



