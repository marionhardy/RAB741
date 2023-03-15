## OVER REPRESENTATION ANALYSIS

library(tidyverse)

res_tbl <- read_csv("./data_output/sel_vs_ctrl_T96/res_tbl.csv")
diff <- read_csv("./data_output/sel_vs_ctrl_T96/sign_different.csv")

# ORA
# These terms are classed into three categories, called namespaces:

# Molecular Function (MF): molecular activities of gene products
# Cellular Component (CC): where gene products are active
# Biological Process (BP): pathways and larger processes made up of the 
# activities of multiple gene products

library("GO.db")
library("org.Mm.eg.db")

# Kyoto Encyclopedia of Genes and Genomes (KEGG) for cell cycle in humans
# reactome.db for réactomes
# Molecular Signatures Database (MSigDB) for use with the GSEA software : 
# msigdbr package


# We will focus on the genes that have an adjusted p-value (those that have 
# been tested) and that have unique ENTREZ gene identifiers.

# GO

# NB, first time i had 1900 hits so Imma be more stringent on the signif genes

library("clusterProfiler")

diffs <- diff %>% filter(padj<=1e-10, abs(log2FoldChange)>= 2)

de_genes <- unique(diffs$ENTREZID)
all_genes <- unique(res_tbl$ENTREZID)

go_ora <- enrichGO(gene = de_genes,
                   OrgDb = org.Mm.eg.db,
                   universe = all_genes,
                   ont = "ALL",
                   readable = TRUE) 

go_ora_tbl <- as_tibble(go_ora)

# https://advaitabio.com/faq-items/understanding-gene-ontology/


go_ora_cc <- go_ora %>% filter(ONTOLOGY=="CC")
go_ora_cc_tbl <- as_tibble(go_ora_cc)
go_ora_bp <- go_ora %>% filter(ONTOLOGY=="BP")
go_ora_bp_tbl <- as_tibble(go_ora_bp)
go_ora_mf <- go_ora %>% filter(ONTOLOGY=="MF")
go_ora_mf_tbl <- as_tibble(go_ora_mf)

# ont = CC ou BP ou MF ou ALL

library("enrichplot")

barplot(go_ora, showCategory=40) + ggtitle("barplot for ORA")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA all", plot = last_plot(), device = png, dpi = 400)
barplot(go_ora_cc, showCategory=40) + ggtitle("barplot for ORA cellular component")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA cellular component", plot = last_plot(), device = png, dpi = 400)
barplot(go_ora_bp, showCategory=40) + ggtitle("barplot for ORA biological process")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA biological process", plot = last_plot(), device = png, dpi = 400)
barplot(go_ora_mf, showCategory=40) + ggtitle("barplot for ORA molecular function")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA molecular function", plot = last_plot(), device = png, dpi = 400)

dotplot(go_ora, showCategory=25) + ggtitle("dotplot for ORA")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA all", plot = last_plot(), device = png, dpi = 400)
dotplot(go_ora_bp, showCategory=50) + ggtitle("dotplot for ORA bp")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA bp", plot = last_plot(), device = png, dpi = 400)
dotplot(go_ora_cc, showCategory=25) + ggtitle("dotplot for ORA cc")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA cc", plot = last_plot(), device = png, dpi = 400)
dotplot(go_ora_mf, showCategory=25) + ggtitle("dotplot for ORA mf")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA mf", plot = last_plot(), device = png, dpi = 400)


heatplot(go_ora, showCategory = 10)+ ggtitle("heatplot for ORA")
ggsave("./figures/sel_vs_ctrl_T96/heatplot go ORA all", plot = last_plot(), device = png, dpi = 400)
heatplot(go_ora_cc, showCategory = 10)+ ggtitle("heatplot for ORA cc")
ggsave("./figures/sel_vs_ctrl_T96/heatplot go ORA cc", plot = last_plot(), device = png, dpi = 400)
heatplot(go_ora_bp, showCategory = 10)+ ggtitle("heatplot for ORA bp")
ggsave("./figures/sel_vs_ctrl_T96/heatplot go ORA bp", plot = last_plot(), device = png, dpi = 400)
heatplot(go_ora_mf, showCategory = 10)+ ggtitle("heatplot for ORA mf")
ggsave("./figures/sel_vs_ctrl_T96/heatplot go ORA mf", plot = last_plot(), device = png, dpi = 400)

write.csv(go_ora,"./data_output/sel_vs_ctrl_T96/go_ora_all.csv")


# Datasets to query panther or gorilla or other
# Allows for external validation of our go_ora tibble

write.table(res_tbl$gene,"./data_output/sel_vs_ctrl_T96/backgrd_gene_names.txt", sep = "\t",
            row.names = FALSE, quote = F)
write.table(diff$gene, file = "./data_output/sel_vs_ctrl_T96/signif_gene_names.txt", sep = "\t",
            row.names = FALSE, quote = F)


## KEGG 

go_kegg <- enrichKEGG(gene = de_genes,
                      organism = "mmu",
                      keyType = "kegg",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      universe = as.character(all_genes),
                      qvalueCutoff = 0.1)

go_kegg

barplot(go_kegg, showCategory = 30) # weird and there are the disease terms :/


## GENE SET ENRICHMENT ANALYSIS

# A common approach to analyzing gene expression profiles is identifying 
# differentially expressed genes that are deemed interesting. The ORA enrichment 
# analysis is based on these differentially expressed genes. This approach 
# will find genes where the difference is large and will fail where the 
# difference is small, but evidenced in coordinated way in a set of related 
# genes. Gene Set Enrichment Analysis (GSEA)(Subramanian et al. 2005) 
# directly addresses this limitation. All genes can be used in GSEA; 
# GSEA aggregates the per gene statistics across genes within a gene set, 
# therefore making it possible to detect situations where all genes in a 
# predefined set change in a small but coordinated way.

# GSEA: genes must be classified in decreasing order

res_tbl1 <-  res_tbl %>% filter(!is.na(ENTREZID), padj<=0.05)
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

# GO with genes ranked based on stat metrics

go_gsea <- gseGO(gene = ordered_genes_stat,
                 OrgDb = org.Mm.eg.db,
                 ont          = "ALL",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 nPermSimple = 1000,
                 eps = 0,
                 scoreType="pos")

dotplot(go_gsea, x="NES",split = "ONTOLOGY", showCategory=10) + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all")

table(go_gsea@result$ONTOLOGY)
go_gsea_tbl <- as_tibble(go_gsea)

go_gsea %>% filter(ONTOLOGY=="BP") %>% 
  dotplot( x = "NES", showCategory = 30)+ ggtitle("dotplot for GSEA biological processes" )
ggsave("./figures/sel_vs_ctrl_T96/dotplot stat go GSEA bp", plot = last_plot(), device = png, dpi = 400)

go_gsea %>% filter(ONTOLOGY=="MF") %>% 
  dotplot(x = "NES", showCategory = 25)+ ggtitle("dotplot for GSEA molecular function" )
ggsave("./figures/sel_vs_ctrl_T96/dotplot stat go GSEA mf", plot = last_plot(), device = png, dpi = 400)



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
                 nPermSimple = 100000,
                 eps=0)

go_gsea_tbl <- as_tibble(go_gsea)

table(go_gsea_tbl$ONTOLOGY)

write.csv(go_gsea_tbl,"./data_output/sel_vs_ctrl_T96/go_gsea_tbl_Log2Fc.csv")
# Notably, no T cell activation GO is considered to be significant because
# scores of NES got close to 0 because of enrichment in both log2fc at the start
# of the gene list and annuled by the enrichment at the bottom of the list

dotplot(go_gsea, split = "ONTOLOGY",showCategory=15, x = "NES") + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA all", plot = last_plot(), device = png, dpi = 400)


dotplot(go_gsea, x = "NES", showCategory = 30)
ggsave("./figures/sel_vs_ctrl_T96/dotplot_GSEA_LFC_NES", plot = last_plot(), device = png, dpi = 400)

go_gsea %>% filter(ONTOLOGY=="MF") %>% 
  dotplot( x = "NES", showCategory = 30)+
  ggtitle("Dotplot for GSEA Molecular Function")
ggsave("./figures/sel_vs_ctrl_T96/dotplot_GSEA_LFC_NES_MF", plot = last_plot(), device = png, dpi = 400)

go_gsea %>% filter(ONTOLOGY=="BP") %>% 
  dotplot( x = "NES", showCategory = 75)+
  ggtitle("Dotplot for GSEA Biological process")
ggsave("./figures/sel_vs_ctrl_T96/dotplot_GSEA_LFC_NES_BP", plot = last_plot(), device = png, dpi = 400,
       height = 25, width = 10)

gseaplot(go_gsea, geneSetID = 1, by = "runningScore", title = "mitochondrial gene expression")
ggsave("./figures/sel_vs_ctrl_T96/GSEA_mitoch_gene_expr", plot = last_plot(), device = png, dpi = 400)

ann <- readRDS("./data_output/210906_ensembl_to_geneName.rds")
GO_mitoch <- data_frame(go_gsea@geneSets$'GO:0140053')
colnames(GO_mitoch) <- "ENTREZID"
GO_mitoch$ENTREZID <- as.numeric(GO_mitoch$ENTREZID)
GO_mitoch <- left_join(GO_mitoch, ann)

write_csv(GO_mitoch, "./data_output/sel_vs_ctrl_T96/GSEA_mitoch_gene_expr")

mit_target <- GO_mitoch$gene
highlight <- res_tbl1 %>% filter(gene%in%mit_target)

library(ggrepel)

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
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
  theme(legend.title= element_blank())+
  labs(title = "Sel vs Ctrl, GO:0140053 mitochondrial gene expression")

ggsave("./figures/sel_vs_ctrl_T96/volcanoplot_targets_GO_mitoch", last_plot(), device = png, dpi= 500)

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


# KEGG

gse_kegg <- gseKEGG(ordered_genes,
                    organism = "mmu",
                    keyType = "kegg",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    scoreType = "pos",
                    eps = 0)

gse_kegg

dotplot(gse_kegg, showCategory = 30)

# KEGG modules (more straightforward in their interpretation)

gse_mkegg <- gseMKEGG(ordered_genes,
                    organism = "mmu",
                    keyType = "kegg",
                    pvalueCutoff = 0.1,
                    pAdjustMethod = "BH",
                    scoreType = "pos")

gse_mkegg

dotplot(gse_mkegg, showCategory = 30)

# It's still possible to use pathview to extract the significant paths

library("pathview")

keggresids <- gse_kegg$ID

foldchanges <- diff$log2FoldChange
names(foldchanges) <- diff$ENTREZID
head(foldchanges)
table(is.na(foldchanges))

tmp <- sapply(keggresids, function(pid) pathview(gene.data = foldchanges,
                                                 pathway.id = pid,
                                                 species = "mmu",
                                                 kegg.dir="./figures/sel_vs_ctrl_T96/KEGG_pathways",
                                                 out.suffix= "_Colored",
                                                 kegg.native = TRUE,
                                                 map.null = FALSE))








