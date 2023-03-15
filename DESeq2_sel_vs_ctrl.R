
## Assumes you have run the RSubread script or
# you have a feature count file made from a .bam file
# Counts should be RAW! not tpkm rkm

library(DESeq2)
library(purrr)
library("tidyverse")
library(stringr)
library(dplyr)
library("biomaRt")
library(Hmisc)

setwd("/Users/hardy/Documents/R_2022/Raphaele/RAB741")
dir <- getwd()

D0_1 <- read_csv(paste0(dir,"/data/counts/Day0_1_counts.csv"))
D0_2 <- read_csv(paste0(dir,"/data/counts/Day0_2_counts.csv"))
C48_1 <- read_csv(paste0(dir,"/data/counts/Ctrl_1_48h_counts.csv"))
C48_2 <- read_csv(paste0(dir,"/data/counts/Ctrl_2_48h_counts.csv"))
C96_1 <- read_csv(paste0(dir,"/data/counts/Ctrl_1_96h_counts.csv"))
C96_2 <- read_csv(paste0(dir,"/data/counts/Ctrl_2_96h_counts.csv"))
S48_1 <- read_csv(paste0(dir,"/data/counts/Sel_1_48h_counts.csv"))
S48_2 <- read_csv(paste0(dir,"/data/counts/Sel_2_48h_counts.csv"))
S96_1 <- read_csv(paste0(dir,"/data/counts/Sel_1_96h_counts.csv"))
S96_2 <- read_csv(paste0(dir,"/data/counts/Sel_2_96h_counts.csv"))

datanames <- c("D0_1","D0_2","C48_1","C48_2","C96_1","C96_2",
               "S48_1","S48_2","S96_1","S96_2")
list <- llist(D0_1,D0_2,C48_1,C48_2,C96_1,C96_2,S48_1,S48_2,S96_1,S96_2)

# Make count matrix (keep gene_id as rownames and rawcounts)

fun <- function(i){
  i <- i %>% as.data.frame()
  rownames(i) <- i$ensembl_gene_id
  i <- i %>% dplyr::select(RawCounts, ensembl_gene_id)
}

list <- lapply(list,fun)
list2env(list, envir=.GlobalEnv)

head(D0_1)

# Change rawcounts to sample name

colnames(D0_1)[1] <- datanames[1]
colnames(D0_2)[1] <- datanames[2]
colnames(C48_1)[1] <- datanames[3]
colnames(C48_2)[1] <- datanames[4]
colnames(C96_1)[1] <- datanames[5]
colnames(C96_2)[1] <- datanames[6]
colnames(S48_1)[1] <- datanames[7]
colnames(S48_2)[1] <- datanames[8]
colnames(S96_1)[1] <- datanames[9]
colnames(S96_2)[1] <- datanames[10]

counts <- list(D0_1,D0_2,C48_1,C48_2,C96_1,C96_2,S48_1,S48_2,S96_1,S96_2) %>% 
  purrr::reduce(dplyr::full_join,)

# counts <- read.csv("D:/Documents/BVDE/R/RAB741/RAB741counts.csv")
rownames(counts) <- counts$ensembl_gene_id
counts <- counts %>% dplyr::select(-ensembl_gene_id)
head(counts)

write.csv(counts,"./data/RAB741counts.csv")


# Create the coldata for the summarized experiment

coldata <- data.frame(
  celltype=c(rep("day0",2),rep("ctrl",4),rep("selected",4)),
  timeline = as.factor(c(0,0,48,48,96,96,48,48,96,96)),
  replicate=as.factor(rep(c(1,2),5)))

rownames(coldata) <- colnames(counts)

## Analyse de l'expression

# Créer un objet deseq dataset soit dds
# faire le design pour pouvoir faire le test
# comme il y a des timepoints, c'est un test pairwise qu'il faudrait faire
# avec la formule ci-dessous je compare juste les celltype

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~celltype) # MAKE MY OWN MODEL MATRIX I THINK


dds$celltype
class(dds$replicate)
dds$timeline

# Générer un modèle linéaire (regarder console pour voir ce qu'il fait comme manip)

dds <- DESeq(dds)

# This is because i have more than three different type of samples
dds <- DESeq(dds, test="LRT", reduced=~timeline)
res <- results(dds)


as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

# distribution uniforme donc on transforme en log (nécessaire pour la pca, le reste du temps
# il transformera tout seul en log)

# Préparation PCA

rld <- rlogTransformation(dds)

plotPCA(rld,intgroup="celltype") # design choisi le même si celltype ou timeline -> weird?
plotPCA(rld,intgroup="timeline")
plotPCA(rld,intgroup="replicate") # Great replicates damn

# selected and control 48h have almost no difference!

# Vérification que normalisation ok

sizeFactors(dds) # ne prend en compte que la profondeur de séquençage

plotDispEsts(dds) # vérification de la normalisation, ici graphe allure ok
# compliqué de piger le graphe, maybe demander à Camille ou voir cours WSBIM2122


resultsNames(dds) # pas fait dans le bon sens aurait dû être calu_vs_ctrl
# dds$celltype <- relevel(dds$celltype, ref="ctrl") # donc on lui dit quelle est la ref
# dds <- DESeq(dds)
# resultsNames(dds)


# les résultats peuvent mtn être extraits

## SELECTED VS CTRL # ! timelines not taken into account

results <- results(dds, name = "celltype_selected_vs_ctrl")
res_tbl <- as_tibble(results, rownames="ENSMUG")

# Analyse

selected_genes <- res_tbl %>%
  filter(padj < 0.05 & padj > 1e-5 & abs(log2FoldChange) > 5)

# it's lncRNA, not of interest but just testing the plot function
as_tibble(counts(dds["ENSMUSG00000109644"], normalize = TRUE),
          rownames = 'ENSMUG') %>%
  gather(sample, counts, -ENSMUG) %>%
  left_join(as_tibble(coldata, rownames = "sample")) %>%
  ggplot(aes(x = sample, y = counts, fill = celltype)) +
  geom_bar(stat = 'identity', color = "gray30", width = .75) +
  theme(axis.text.x = element_text(size = 11, angle = 1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))+
  ggtitle("mRNA counts")

ggsave(filename = "./figures/211229 mrna count plot.png", plot = last_plot())


hist(res_tbl$padj) # distribution ok, sinon besoin de filtrer les padj, voir cours

plotMA(res) # joli plot waw

## Annotations

library("biomaRt")

# Annotations des codes ENSMUG

listDatasets(useMart("ensembl"))

## Load mus musculus ensembl dataset
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

# listAttributes(mart)

ensembl_to_geneName <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                            "entrezgene_id","description"),
                             mart = mart)
names(ensembl_to_geneName) <- c("ENSMUG", "gene", "ENTREZID", "description")
head(ensembl_to_geneName)

# Add gene names and entrez id to your tibble

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj)

head(res_tbl) # it worked

write.csv(res_tbl,"./data_output/sel_vs_ctrl/res_tbl.csv")

a <- 
  res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1,
             label= gene)) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)

print(a)

# Volcano plot with labels

library(ggrepel)

# highlight_list <- c("PRKAA1","PRKAA2","PRKAB1","PRKAB2","STK11","TP53")
# res_tbl <- res_tbl %>%
#   mutate(plotname = ifelse(plotname %in% highlight_list, plotname, ""))

overexpr <- res_tbl %>% filter(padj<0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

write_csv(overexpr,"./data_output/sel_vs_ctrl/Sign_overexpr.csv")
write_csv(underexpr,"./data_output/sel_vs_ctrl/Sign_underexpr.csv")
write_csv(signif,"./data_output/sel_vs_ctrl/Sign_different.csv")

library("openxlsx")

write.xlsx(overexpr,"./data_output/sel_vs_ctrl/RAB741_sign_genes.xlsx", sheetName = "UP")
write.xlsx(underexpr,"./data_output/sel_vs_ctrl/RAB741_sign_genes.xlsx", sheetName = "DOWN", append = T)
write.xlsx(signif,"./data_output/sel_vs_ctrl/RAB741_sign_genes.xlsx", sheetName = "ALL_sign", append = T)

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/sel_vs_ctrl/RAB741_sign_genes.xlsx")

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 50)+
  theme_bw()

# save the tibbles etc

saveRDS(dds, file = "./data_output/RAB741_dds.rds")
saveRDS(res_tbl, file = "./data_output/RAB741_res_tbl_ctrlvsselected.rds")
write_tsv(res_tbl, file = "./data_output/210906_res_tbl.tsv")
saveRDS(ensembl_to_geneName, file = "./data_output/210906_ensembl_to_geneName.rds")


















