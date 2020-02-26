# load library 
library(tidyverse)

# read dE data
res_tableOE <- read.csv("data/Mov10oe_DE_results.csv", row.names = 1)

res_tableOE_tb <- res_tableOE %>% rownames_to_column(var="gene") %>% as_tibble()

## Load libraries
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(pathview)
library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)

# Optional for the lesson
library(gProfileR)
library(treemap)
library(SPIA)
library(annotables)

# Load libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = res_tableOE_tb$gene,  # data to use for retrieval
                                           columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retreive for given data
                                           keytype = "SYMBOL") # type of data given in 'keys' argument

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

# Check number of NAs returned
is.na(annotations_orgDb$ENSEMBL) %>%
  which() %>%
  length()


# Load the library
library(EnsDb.Hsapiens.v75)

# Check object metadata
EnsDb.Hsapiens.v75

# Explore the fields that can be used as keys
keytypes(EnsDb.Hsapiens.v75)
# Return the Ensembl IDs for a set of genes
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                         keys = res_tableOE_tb$gene,
                                         columns = c("GENEID", "ENTREZID","GENEBIOTYPE"),
                                         keytype = "SYMBOL")

# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_edb <- annotations_edb[non_duplicates_idx, ]

# Check number of NAs returned
is.na(annotations_edb$GENEID) %>%
  which() %>%
  length()


# FUNCCTIONAL

# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)

## Merge the annotations with the results 
res_ids <- inner_join(res_tableOE_tb, annotations_edb, by=c("gene"="SYMBOL")) 

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(res_ids$GENEID)
## Extract significant results
sigOE <- dplyr::filter(res_ids, padj < 0.05)
sigOE_genes <- as.character(sigOE$GENEID)

## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
dotplot(ego, showCategory=50)
emapplot(ego, showCategory = 50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
OE_foldchanges <- ifelse(OE_foldchanges > 3, 3, OE_foldchanges)
OE_foldchanges <- ifelse(OE_foldchanges < -3, -3, OE_foldchanges)

cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 3, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

