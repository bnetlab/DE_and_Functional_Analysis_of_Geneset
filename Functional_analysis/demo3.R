
# Load libraries
library(limma)
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(pathview)
library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)

#
# Limma get global gene
res_tableOE <- read.csv(file="data/de1.dat", sep = '\t')
res_tableOE_tb <- res_tableOE %>% rownames_to_column(var="gene") %>% as_tibble()
# Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = as.character(res_tableOE$GeneID),  # data to use for retrieval
                                           columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retreive for given data
                                           keytype = "SYMBOL") # type of data given in 'keys' argument
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)
# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]
## Merge the annotations with the results 
res_ids <- inner_join(res_tableOE_tb, annotations_orgDb, by=c("GeneID"="SYMBOL")) 



# DTA gene
geneDTA <- read.csv(file="data/List4_50.csv", colClasses=c("character"))
geneDTA <- c(names(geneDTA))
annotations_orgDbDTA <- AnnotationDbi::select(org.Hs.eg.db, # database
                                              keys = as.character(geneDTA),  # data to use for retrieval
                                              columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retreive for given data
                                              keytype = "SYMBOL") # type of data given in 'keys' argument

allOE_genes <- as.character(res_ids$ENSEMBL)
## Extract significant results
sigOE_genesDTA <- as.character(annotations_orgDbDTA$ENSEMBL[1:20])

## Run GO enrichment analysis 
egoDTA_go <- enrichGO(gene = sigOE_genesDTA, 
                   universe = allOE_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "none", 
                   qvalueCutoff = 1,
                   pvalueCutoff = 1,
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summaryDTA_go <- data.frame(egoDTA_go)
dotplot(egoDTA_go, showCategory=50)
emapplot(egoDTA_go, showCategory = 50)

# kegg

allOE_genes <- as.character(res_ids$ENTREZID)
## Extract significant results
sigOE_genes <- as.character(annotations_orgDbDTA$ENTREZID[1:50])

egoDTA_Kegg <- enrichKEGG(gene = sigOE_genes, 
                  universe = allOE_genes,
                  pAdjustMethod = "none", 
                  qvalueCutoff = 1,
                  pvalueCutoff = 1)

cluster_summaryDTA_kegg <- data.frame(egoDTA_Kegg)
dotplot(egoDTA_Kegg, showCategory=50)
emapplot(egoDTA_Kegg, showCategory = 50)


#######FUNCTIONAL CLASS SCORING #######################

# remove NA and duplicate
## Remove any NA values
res_entrez <- res_ids %>% dplyr::filter(ENTREZID %in% annotations_orgDbDTA$ENTREZID)
res_entrez <- dplyr::filter(res_entrez, ENTREZID != "NA")
res_entrez <- res_entrez[which(duplicated(res_entrez$ENTREZID) == F), ]

foldchanges <- res_entrez$logFC
names(foldchanges) <- res_entrez$ENTREZID
foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 10, # default number permutations
                    minGSSize = 1, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 1,
                    pAdjustMethod= "none", # padj cutoff value
                    verbose = FALSE)
#pathview

detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts

## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
         pathway.id = "hsa05223",
         species = "hsa",
         limit = list(gene = 1, # value gives the max/min limit for foldchanges
                      cpd = 1))

