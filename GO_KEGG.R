## OVER REPRESENTATION ANALYSIS

# ORA using GO
# KEGG using GO

res_tbl <- read.csv("./data/sel_vs_ctrl/res_tbl.csv")

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

library("clusterProfiler")

de_genes <- res_tbl$ENTREZID[res_tbl$padj < 0.05]
all_genes <- res_tbl$ENTREZID

go_ora_tbl <- enrichGO(gene = de_genes,
                       OrgDb = org.Mm.eg.db,
                       universe = all_genes,
                       ont = "ALL",
                       readable = TRUE) %>%
  as_tibble

go_ora <- enrichGO(gene = de_genes,
                   OrgDb = org.Mm.eg.db,
                   universe = all_genes,
                   ont = "ALL",
                   readable = TRUE) 


# Check for mitochondrial gene enrichment?

go_ora_cc_tbl <- go_ora_tbl %>% filter(ONTOLOGY=="CC")
go_ora_mf_tbl <- go_ora_tbl %>% filter(ONTOLOGY=="MF")

# ont = CC ou BP ou MF ou ALL

library("enrichplot")

barplot(go_ora, showCategory=40) + ggtitle("barplot for ORA")

dotplot(go_ora, showCategory=30) + ggtitle("dotplot for ORA")

heatplot(go_ora, showCategory = 15) # need to select the genes

write.csv(go_ora,"./data/sel_vs_ctrl/go_ora_all.csv")


# Datasets to query panther or gorilla or other
# Allows for external validation of our go_ora tibble

write.table(res_tbl$gene,"./data/sel_vs_ctrl/backgrd_gene_names.txt", sep = "\t",
            row.names = FALSE, quote = F)
write.table(signif$gene, file = "./data/sel_vs_ctrl/signif_gene_names.txt", sep = "\t",
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

barplot(go_kegg, showCategory = 30) # weird and theere are the disease terms :/


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

ordered_genes <- abs(res_tbl$stat) # could also use log2fc to order them
names(ordered_genes) <- res_tbl$ENTREZID
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea_tbl <- gseGO(gene = ordered_genes,
                 OrgDb = org.Mm.eg.db,
                 scoreType = "pos",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE) %>%
  as_tibble()

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Mm.eg.db,
                 scoreType = "pos",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

go_gsea

write.csv(go_gsea_tbl,"./data/sel_vs_ctrl/go_gsea_tbl.csv")

dotplot(go_gsea %>% filter(ONTOLOGY=="CC"), showCategory=30)
dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30)
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

tmp <- sapply(keggresids, function(pid) pathview(gene.data = foldchanges,
                                                 pathway.id = pid,
                                                 species = "mmu",
                                                 kegg.dir="./figures/sel_vs_ctrl/KEGG_pathways/UP_DOWN_both",
                                                 out.suffix= "_Colored",
                                                 kegg.native = TRUE,
                                                 map.null = FALSE))


mmu05417 <- pathview(gene.data  = ordered_genes,
                     pathway.id = "mmu05417",
                     species    = "mmu",
                     limit      = list(gene=max(abs(ordered_genes)), cpd=1))






