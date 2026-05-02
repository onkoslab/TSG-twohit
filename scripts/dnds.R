library(dndscv)
library(dplyr)

data_dir <- "data"

# Gene-level dnds

# defining tsgs
generoles = read.table(file.path(data_dir, "generoles.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
tsgs = unique(generoles[which(generoles$role=='TSG'),]$Hugo_Symbol)
tsgs[tsgs == "CDKN2A"] <- "CDKN2A.p16INK4a"

# Pancancer TSG 1-hit mutations
muts1hit = read.table(file.path(data_dir, "mutcat1hit_tsg.5col"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
dout1hit = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgs)
ci1hit = geneci(dout1hit, gene_list = tsgs)

# Pancancer TSG 2-hit mutations
muts2hit = read.table(file.path(data_dir, "mutcat2hit_tsg.5col"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
dout2hit = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgs)
ci2hit = geneci(dout2hit, gene_list = tsgs)

# 1hit

plot1hit_df_long <- ci1hit %>%
  rename(gene_name = gene) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = c("class", "stat"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

plot1hit_df_long$class <- recode(
  plot_df_long$class,
  mis = "Missense",
  tru = "Truncating"
)

# 2hit

plot2hit_df_long <- ci2hit %>%
  rename(gene_name = gene) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = c("class", "stat"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

plot2hit_df_long$class <- recode(
  plot_df_long$class,
  mis = "Missense",
  tru = "Truncating"
)

write.table(dout1hit$sel_cv, file.path(data_dir, "dndscv_tsgs_onehit_pvals.tsv"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(plot1hit_df_long, file.path(data_dir, "dndscv_tsgs_onehit_cis.tsv"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(dout2hit$sel_cv, file.path(data_dir, "dndscv_tsgs_twohit_pvals.tsv"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(plot2hit_df_long, file.path(data_dir, "dndscv_tsgs_twohit_cis.tsv"), sep='\t', row.names=FALSE, quote=FALSE)

# Global dnds

# defining other gene sets
ongs = unique(generoles[which(generoles$role=='ONG'),]$Hugo_Symbol)
ess = unique(generoles[which(generoles$role=='essential'),]$Hugo_Symbol)
noness = unique(generoles[which(generoles$role=='nonessential'),]$Hugo_Symbol)
noness <- setdiff(noness, "GPX6") # To avoid the dndscv error: The following input gene names are not in the RefCDS database: GPX6

# defining high vs low twohitfreq tsgs
twohitfreq = read.table(file.path(data_dir, "pancan_twohitfreq.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
tsgshigh <- twohitfreq %>%
  dplyr::filter(X2hitfreq > 0.4, FDR < 0.05) %>%
  dplyr::pull(Hugo_Symbol) %>%
  unique()
tsgs[tsgshigh == "CDKN2A"] <- "CDKN2A.p16INK4a"
tsgslow <- twohitfreq %>%
  dplyr::filter(!Hugo_Symbol %in% tsgshigh) %>%
  dplyr::pull(Hugo_Symbol) %>%
  unique()

res1hittsg = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgs)
res1hittsghigh = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgshigh)
res1hittsglow = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgslow)
res1hitong = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = ongs)
res1hitess = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = ess)
res1hitnoness = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = noness)
res1hitother = dndscv(muts1hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = other)

res2hittsg = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgs)
res2hittsghigh = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgshigh)
res2hittsglow = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = tsgslow)
res2hitong = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = ongs)
res2hitess = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = ess)
res2hitnoness = dndscv(muts2hit, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample =1000, outmats=T, gene_list = noness)

# combine results

extract_global <- function(res, category, hit) {
  res$globaldnds %>%
    filter(name %in% c("wmis", "wtru")) %>%
    mutate(
      category = category,
      hit = hit
    )
}

df_allhighlow <- bind_rows(
  ## 1-hit
  extract_global(res1hittsg,    "TSG",    "1-hit"),
  extract_global(res1hittsghigh,    "TSGhigh",    "1-hit"),
  extract_global(res1hittsglow,    "TSGlow",    "1-hit"),
  extract_global(res1hitong,    "ONG",    "1-hit"),
  extract_global(res1hitess,    "ESS",    "1-hit"),
  extract_global(res1hitnoness, "Non-ESS","1-hit"),

  ## 2-hit
  extract_global(res2hittsg,    "TSG",    "2-hit"),
  extract_global(res2hittsghigh,    "TSGhigh",    "2-hit"),
  extract_global(res2hittsglow,    "TSGlow",    "2-hit"),
  extract_global(res2hitong,    "ONG",    "2-hit"),
  extract_global(res2hitess,    "ESS",    "2-hit"),
  extract_global(res2hitnoness, "Non-ESS","2-hit"),
)

# write outout

write.table(
  as.data.frame(df_allhighlow),
  file.path(data_dir, "globaldnds.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
