# EX2

```{r}
library("compGenomRData")
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                          package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")
colData <- read.table(coldata_file, header = T, sep = '\t',
stringsAsFactors = TRUE)

#1 
counts_df <- read.table(counts_file,header = T, sep = '\t') #reading the table
counts_mat <- as.matrix(counts_df) #transfering the data frame to matrix
geneLengths_mat <- subset(counts_mat, select = c(width)) #take the lenghth of the gene which is the "width column"
geneLengths_vector <- as.vector(geneLengths_mat) #transfer the length to vector 
rpkm <- apply(X = subset(counts_mat, select = c(-width)), #include all column but not width
              MARGIN = 2,#this indicates columns 
              FUN = function(x) {
              (x * 10^9) / (geneLengths_vector * sum(as.numeric(x)))
              })

colSums(rpkm)
#normalize by the sample size using rpkm values
tpm <- apply(rpkm, 2, function(x) {(x / sum(as.numeric(x))) * 10^6})

colSums(tpm)

#2

library(pheatmap)
top500genes <- names(sort(apply(tpm,1,var), decreasing = T)[1:500])
pheatmap(tpm[top500genes,], scale = 'row',show_rownames = FALSE)
top100genes <- names(sort(apply(tpm,1,var), decreasing = T)[1:100])
pheatmap(tpm[top100genes,], scale = 'row', show_rownames = FALSE)
#the correlation between the two samples can be seen better when we choose the first 500 variable genes.

#3
pheatmap(tpm[top100genes,], scale = 'none', column = 'none', show_rownames = FALSE)
pheatmap(tpm[top100genes,], scale = 'column', column = 'none', show_rownames = FALSE)

pheatmap(tpm[top500genes,], scale = 'none', column = 'none', show_rownames = FALSE)
pheatmap(tpm[top500genes,], scale = 'column', column = 'none', show_rownames = FALSE)

#4
library(stats)
correlationMatrix <- cor(tpm)
library(corrplot)
corrplot(correlationMatrix , order = 'hclust', addrect = 2, addCoef.col = 'white', number.cex = 0.7, method = 'ellipse', hclust.method = "average")

#5
correlation100genes <- cor(tpm[top100genes,])

pheatmap(correlation100genes, order = 'hlcust', cutree_cols = 2, scale = "row" ,annotation_col = colData)
pheatmap(correlation100genes, order = 'hlcust', cutree_cols = 2, scale = "column" ,annotation_col = colData)
pheatmap(correlation100genes, order = 'hlcust', cutree_cols = 2, scale = "none",annotation_col = colData )

library(stats)
library(ggplot2)
library(ggfortify)#transpose the matrix
library(tidyverse)

M <-t(tpm[top100genes,]) 
# transform the counts to log2 scale 
M <-log2(M + 1) #compute PCA
pcaResults <-prcomp(M) 

autoplot(pcaResults, data = colData , colour = 'group')


#6

annotationofsample <- data.frame(Batch=colData[,1], row.names=rownames(colData))
pheatmap(correlation100genes, order = 'hlcust', cutree_cols = 2, scale = "row" ,annotation_col = annotationofsample)
pheatmap(correlation100genes, order = 'hlcust', cutree_cols = 2, scale = "column" ,annotation_col = annotationofsample)
pheatmap(correlation100genes, order = 'hlcust', cutree_cols = 2, scale = "none" ,annotation_col = annotationofsample)

autoplot(pcaResults, data = colData , colour = 'source_name', shape = 'group')
#7

filteredgenes <- tpm[rowSums(tpm) != 0, ] #we go over all the genes and chech which are not 0.
pheatmap(filteredgenes, order = 'hlcust', cutree_cols = 2, scale = "row", show_rownames = FALSE ,annotation_col = annotationofsample)


#my GitHub link
# https://github.com/LailaB22
```
