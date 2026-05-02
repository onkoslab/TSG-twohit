### Importing packages
library(limma)
library(GSVA)
library(preprocessCore)
library(GSEABase)
library(GSVAdata)
data(c2BroadSets) 
c2BroadSets
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(readxl)

### Importing genes for 50 different hallmark pathways
dat = as.matrix(read.csv("hall_mark_genes_df.csv", header = TRUE,row.names="X"))

geneSets <- list()
ttd<-list()

for(i in 1:50) 
{
  geneSets[[colnames(dat)[i]]] <- dat[,i]
  ttd<-append(ttd,dat[,i])
}

files <- c("gene_exp_LGG.csv", "gene_exp_UCEC.csv")

results_list <- lapply(files, function(f) {
  cts <- as.matrix(read.csv(f, row.names = "X"))
  cts <- subset(cts, rownames(cts) %in% ttd)
  gsva_es <- gsva(cts, geneSets, kcdf = "Poisson", mx.diff = 1, method = "gsva")
  t(gsva_es)
})

# combine outputs for cancer types
final_df <- do.call(rbind, results_list)

final_df <- final_df[!duplicated(rownames(final_df)), ]

# write output

write.table(
  final_df,
  file = "GSVA_pathway_scores_newLGGUCEC.tsv",
  sep = "\t",
  quote = FALSE
)