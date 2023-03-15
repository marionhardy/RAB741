
## Verifying that T cell activation is truly up in Selected T96 compared to
# control T96

library(tidyverse)

# I'll check if the genes i got for GO ORA T cell activation are in the UP
# category of the volcanoplot since T cell activation GSEA analysis says 
# NES > 0 so it's up and we thought it would be down.

res_tbl <- read_csv("./data_output/sel_vs_ctrl_T96/res_tbl.csv")
SvsC96 <- read_csv("./data_output/sel_vs_ctrl_T96/Sign_different.csv")

hallmark_tcellact <- c("Dnaja3","Tarm1","Cebpb","Bmi1","Mettl3","Cbfb","Lgals3","Slc7a1","Dlg5","Shb","Glmn",
              "Efnb1","Il6","Lag3","H2-Oa","Ccr7","Tespa1","Ccl5","Malt1","Cd74","Ctla2a","Ctla4",
              "H2-Ab1","Il7r","Jak3","Anxa1","Nlrp3","Ccr2","Cd4","Cd5","Cd44","Ripk3","Pik3r6","Cd81",
              "Runx1","Gpnmb","Ripor2","Ptpn6","H2-DMb2","Lrrc32","Il27ra","Prdm1","Cd55",
              "Tmem131l","H2-Aa","Sdc4","Vsir","H2-Ob","Zc3h12d","Btla","Cd244a","H2-Eb1",
              "Tnfsf14","Tnfrsf13c","H2-DMb1","Ada","Irf1","Icosl","Cd160","Tnfsf9","Cd83","Sox4")

highlight <- SvsC96 %>% filter(gene%in%hallmark_tcellact)

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|
                              padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "lightcoral")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  theme_bw()+
  geom_point(data=highlight, 
             aes(x=log2FoldChange,y=-log10(padj)), 
             color='blue',
             size=0.75)+
  geom_text_repel(data = highlight, size = 4, segment.color = "blue",
                  max.overlaps = 20, min.segment.length = 0, color = "blue",
                  box.padding = 0.5)+
  theme(legend.title= element_blank())

table(highlight$log2FoldChange>0)
# 49 of those genes are down in expression for 13 up
# so it's weird the GSEA enrichment score is positive right?
# No since it only appears in abs(stat) ranked genes

# Metrics I used to do the GSEA test was absolute value of stat test

res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<0.05)
table(duplicated(res_tbl1$ENTREZID))
# take care of those duplicated ENTREZID before continuing

# I'm trying to use log2fc for ranking instead of abs(stat) as in Laurent's course

ordered_genes_fc <- res_tbl1$log2FoldChange
names(ordered_genes_fc) <- res_tbl1$ENTREZID
ordered_genes_fc <- sort(ordered_genes_fc, decreasing = T)

# Let's also try abs(log2fc)

ordered_genes_absfc <- res_tbl1$log2FoldChange
names(ordered_genes_absfc) <- res_tbl1$ENTREZID
ordered_genes_absfc <- sort(abs(ordered_genes_absfc), decreasing = T)

# And stat

ordered_genes_stat <- res_tbl1$stat
names(ordered_genes_stat) <- res_tbl1$ENTREZID
ordered_genes_stat <- sort(ordered_genes_stat, decreasing = T)

# And abs(stat), as before

ordered_genes_absstat <- res_tbl1$stat
names(ordered_genes_absstat) <- res_tbl1$ENTREZID
ordered_genes_absstat <- sort(abs(ordered_genes_absstat), decreasing = T)

# And pval

ordered_genes_pval <- res_tbl1$padj
names(ordered_genes_pval) <- res_tbl1$ENTREZID
ordered_genes_pval <- sort(ordered_genes_pval, decreasing = T)

# GO with genes ranked based on stat metrics

go_gsea <- gseGO(gene = ordered_genes_stat,
                 OrgDb = org.Mm.eg.db,
                 ont          = "ALL",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 nPermSimple = 10000,
                 eps = 0)

dotplot(go_gsea, split = "ONTOLOGY",showCategory=10) + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all")

dotplot(go_gsea, x = "NES", showCategory = 30)
# 1. There is NO difference when doing abs(stat) or stat
# all significant processes are considered enriched
# 2. A negative stat test would mean a high pval so genes that are down significantly
# would still have a positive stat -> nothing signif would be considered downreg

# GO with genes ranked based on absl2fc metrics

go_gsea <- gseGO(gene = ordered_genes_absfc,
                 OrgDb = org.Mm.eg.db,
                 ont          = "ALL",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 nPermSimple = 100000,
                 eps = 0)

dotplot(go_gsea, split = "ONTOLOGY",showCategory=15) + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all")

dotplot(go_gsea, showCategory=20) 
ggtitle("Dotplot for GSEA all")
dotplot(go_gsea, x = "NES", showCategory = 30)

# GO with genes ranked based on l2fc metrics

go_gsea <- gseGO(gene = ordered_genes_fc,
                 OrgDb = org.Mm.eg.db,
                 ont          = "ALL",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 nPermSimple = 10000,
                 eps = 0)

go_gsea_tbl <- as_tibble(go_gsea)

View(go_gsea_tbl)

write.csv(go_gsea_tbl,"./data_output/sel_vs_ctrl_T96/go_gsea_tbl_Log2Fc.csv")
# Notably, no T cell activation GO is considered to be significant because
# scores of NES got close to 0 because of enrichment in both log2fc at the start
# of the gene list and annuled by the enrichment at the bottom of the list


dotplot(go_gsea, split = "ONTOLOGY",showCategory=10) + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA all", plot = last_plot(), device = png, dpi = 400)

g_BP <- go_gsea[go_gsea@result$ONTOLOGY=="BP",]
g_BP <- g_BP %>% filter(p.adjust<=0.05)
g_BP$type[g_BP$NES < 0] = "downregulated"
g_BP$type[g_BP$NES > 0] = "upregulated"

dotplot(go_gsea, x = "NES", showCategory = 30)
ggsave("./figures/sel_vs_ctrl_T96/dotplot_GSEA_LFC_NES", plot = last_plot(), device = png, dpi = 400)


ggplot(g_BP, aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size = setSize, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")+
  facet_grid(.~type)

gseaplot(go_gsea, geneSetID = 1, by = "runningScore", title = "mitochondrial gene expression")
ggsave("./figures/sel_vs_ctrl_T96/GSEA_mitoch_gene_expr", plot = last_plot(), device = png, dpi = 400)

ann <- readRDS("./data_output/210906_ensembl_to_geneName.rds")
GO_mitoch <- data_frame(go_gsea@geneSets$'GO:0140053')
colnames(GO_mitoch) <- "ENTREZID"
GO_mitoch$ENTREZID <- as.numeric(GO_mitoch$ENTREZID)
GO_mitoch <- left_join(GO_mitoch, ann)

write_csv(GO_mitoch, "./data_output/sel_vs_ctrl_T96/GSEA_mitoch_gene_expr")

highlight <- res_tbl %>% filter(gene%in%mit_target)

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|
                              padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "lightcoral")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  theme_bw()+
  geom_point(data=highlight, 
             aes(x=log2FoldChange,y=-log10(padj)), 
             color='blue',
             size=0.75)+
  geom_text_repel(data = highlight, size = 4, segment.color = "blue",
                  max.overlaps = 20, min.segment.length = 0, color = "blue",
                  box.padding = 0.5)+
  theme(legend.title= element_blank())

# When using log2fc as the metric, NES positive scores are made of upregulated genes
# in majority (12 out of 13 genes are up)
# The list gives 109 observations though

highlight %>% 
  group_by(padj, log2FoldChange) %>% 
  count(sign_padj = padj >= 0.05,nsign_padj = )


# ego2 <- simplify(go_gsea)
# cnetplot(ego2, foldChange=ordered_genes, circular = TRUE, colorEdge = TRUE)

kk <- gseKEGG(ordered_genes, organism = "mmu", eps = 0)
ridgeplot(kk)

dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30)+ ggtitle("Dotplot for GSEA BP")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA BP", plot = last_plot(), device = png, dpi = 400)
dotplot(go_gsea %>% filter(ONTOLOGY=="CC"), showCategory=30)+ ggtitle("Dotplot for GSEA CC")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA CC", plot = last_plot(), device = png, dpi = 400)
dotplot(go_gsea %>% filter(ONTOLOGY=="MF"), showCategory=30)+ ggtitle("Dotplot for GSEA MF")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA MF", plot = last_plot(), device = png, dpi = 400)

heatplot(go_gsea, showCategory = 5)

# GO with genes ranked based on padj metrics

go_gsea <- gseGO(gene = ordered_genes_pval,
                 OrgDb = org.Mm.eg.db,
                 ont          = "ALL",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 scoreType = "pos",
                 nPermSimple = 100000,
                 eps = 0)

go_gsea_tbl <- as_tibble(go_gsea)

View(go_gsea_tbl)

write.csv(go_gsea_tbl,"./data_output/sel_vs_ctrl_T96/go_gsea_tbl_padj.csv")


# Trying to reduce redundancy of GO terms: be more stringent about cutoffs






# Try MSigDb sets (reactome, kegg, hallmark)

library(msigdbr)
library(clusterProfiler)
print(msigdbr_collections(), n=23)

# Database containing annotated gene sets that can be used for pathway or gene set analyses 
# They have 9 collections : Hallmark, C1-C8
# https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H

mm_hallmark_sets <- msigdbr(
  species = "Mus musculus", 
  category = "H") # for hallmark collection

mm_reactome_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C2",
  subcategory = "CP:REACTOME") # for reactome collection

mm_kegg_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C2",
  subcategory = "CP:KEGG") # for KEGG collection

mm_biocarta_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C2",
  subcategory = "CP:BIOCARTA") # for Biocarta collection

mm_wiki_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C2",
  subcategory = "CP:WIKIPATHWAYS") # for Wikipathways collection

mm_regulatory_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C3",
  subcategory = "TFT:TFT_Legacy") # for regulatory genes collection

mm_immune_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C7",
  subcategory = "IMMUNESIGDB") # for immune signature genes

mm_celltype_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C8") # for celltype signature


ordered_genes_fc <- res_tbl1$log2FoldChange
names(ordered_genes_fc) <- res_tbl1$gene
ordered_genes_fc <- sort(ordered_genes_fc, decreasing = T)

set.seed(69)
gsea_results_h <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

dotplot(gsea_results_h, x = "NES", showCategory = 30)+ ggtitle("GSEA hallmark sets LFC")
ggsave("./figures/sel_vs_ctrl_T96/MSigDB/dotplot GSEA LFC hallmark sets", plot = last_plot(), device = png, dpi = 400)

gsea_results_react <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_reactome_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_react, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC")
ggsave("./figures/sel_vs_ctrl_T96/MSigDB/dotplot GSEA LFC reactome", plot = last_plot(), device = png, dpi = 400)

gsea_results_kegg <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_kegg_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_kegg, x = "NES", showCategory = 30)+ ggtitle("GSEA KEGG LFC")
ggsave("./figures/sel_vs_ctrl_T96/MSigDB/dotplot GSEA LFC KEGG", plot = last_plot(), device = png, dpi = 400)

gsea_results_bioc <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_biocarta_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_bioc, x = "NES", showCategory = 30)+ ggtitle("GSEA Biocarta LFC")
ggsave("./figures/sel_vs_ctrl_T96/MSigDB/dotplot GSEA LFC Biocarta", plot = last_plot(), device = png, dpi = 400)

gsea_results_wiki <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_wiki_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_wiki, x = "NES", showCategory = 30)+ ggtitle("GSEA Wikipathway LFC")
ggsave("./figures/sel_vs_ctrl_T96/MSigDB/dotplot GSEA LFC Wikipathway", plot = last_plot(), device = png, dpi = 400)

gsea_results_reg <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_regulatory_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_reg, x = "NES", showCategory = 30)+ ggtitle("GSEA Regulatory sets LFC")
ggsave("./figures/sel_vs_ctrl_T96/MSigDB/dotplot GSEA LFC Regulatory sets", plot = last_plot(), device = png, dpi = 400)


# Try with an external package (check for validation)

library(ReactomePA)

# With gsea

ordered_genes_fc <- res_tbl1$log2FoldChange
names(ordered_genes_fc) <- res_tbl1$ENTREZID
ordered_genes_fc <- sort(ordered_genes_fc, decreasing = T)

reactome_fc <- gsePathway(ordered_genes_fc,
                          organism = "mouse",
                          pvalueCutoff=0.05, 
                          pAdjustMethod = "BH",
                          eps = 0)

head(summary(reactome_fc))

dotplot(reactome_fc, showCategory = 30, x = "NES")

ggplot(reactome_fc, 
       aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size = setSize, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.1), low="red") +
  ylab(NULL) +
  ggtitle("Reactome pathway enrichment")

reactome_fc <- gsePathway(ordered_genes_fc,
                          organism = "mouse",
                          pvalueCutoff=0.05, 
                          pAdjustMethod = "BH",
                          eps = 0)

head(summary(reactome_fc))

dotplot(reactome_fc, showCategory = 30, x = "NES")

# Good thing it gives the same results as MsigDB when permutations are increased

gsea_results_imm <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_immune_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_imm, x = "NES", showCategory = 20)

gsea_results_cell <- GSEA(
  geneList = ordered_genes_fc, # Ordered ranked gene list
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_celltype_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_cell, x = "NES", showCategory = 20)


## Query for Raphaël------------------------------------------------------------

# she wants to know about mitophagy, ampk, t cell function etc
# i made a database with known markers for different processes
library(readxl)
query <- read_xlsx("./data/candidates.xlsx", sheet="Database")

# SvsC96$gene contain the > 1 LFC and <-1 LFC, padj <0,05 genes

hits <- data.frame(matrix(NA,ncol=3, nrow=ncol(query)))
colnames(hits)[1] <- "Percent_found_genes"
colnames(hits)[2] <- "Ratio"
colnames(hits)[3] <- "Found_queried_genes"

for(i in 1:ncol(query)){
  col <- query[,i][!is.na(query[,i])]
  tmp <- round(table(col%in%SvsC96$gene)[2]/table(col%in%SvsC96$gene)[1],2)*100
  tmp1 <- paste0(round(table(col%in%SvsC96$gene)[2],2),"/",
                 round(table(col%in%SvsC96$gene)[1],2))
  tmp_genes <- paste(col[col%in%SvsC96$gene], collapse = "/")
  hits$Percent_found_genes[i] <- tmp
  hits$Ratio[i] <- tmp1
  hits$Found_queried_genes[i] <- tmp_genes
  rownames(hits)[i] <- colnames(query[,i])
}

head(hits)

write_csv(hits,"./data_output/sel_vs_ctrl_T96/hits.csv")



