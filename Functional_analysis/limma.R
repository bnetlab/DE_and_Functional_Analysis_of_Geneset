library(PEPPER)
library(limma)
out.file <- "de.dat"
adjust.method <- 'BH'
fdr.cutoff <- 0.05
sample.mapping= read.csv2("sampleMapping.csv", check.names = FALSE, sep = ',' )
expr=read.csv("sample.csv", sep = ',')
sample.mapping$sample=gsub('-','.',sample.mapping$sample)
de <- find.de.genes(expr, sample.mapping, c("case", "control"), method='limma', out.file, adjust.method=adjust.method, cutoff=fdr.cutoff) 
de <- de[abs(de$logFC)>=1,]
