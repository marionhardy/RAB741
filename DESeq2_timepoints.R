
## This assumes you have run the DESeq2_file_preparation script
# and that you have a counts dataframe and a coldata dataframe

library(DESeq2)
library(tidyverse)
library(cowplot)
library("biomaRt")

## Expression analysis of
# Selected T48 vs ctrl T48
# Then
# Selected T96 vs ctrl T96

# Ctrl 48h vs selected 48h is unpaired
# Ctrl 48h vs ctrl 96h is paired
# However, DESeq2 is capable of doing those analyses by itself
# An alternative to pair-wise comparisons is to analyze all levels of a factor at 
# once. By default the Wald test is used to generate the results table, but 
# DESeq2 also offers the LRT which is used to identify any genes that show 
# change in expression across the different levels. This type of test can be 
# especially useful in analyzing time course experiments.

# LRT is equivalent to Wald but turns off LogFoldChange shrinkage which can
# be beneficial when doing subsequent GSEA analysis apparently
# But no real advantage to one or the other except LRT is way more documented

# Day0 and the 0h timepoint are making analysis complicated so i'll remove
# them for this one

counts <- read.csv("./data/RAB741counts.csv") # and create a new coldata

rownames(counts) <- counts$X

counts <- counts %>% 
  dplyr::select(-D0_1, -D0_2, -X)

coldata <- data.frame(
  celltype=c(rep("ctrl",4),rep("selected",4)),
  timepoint = as.factor(c(48,48,96,96,48,48,96,96)),
  replicate=as.factor(rep(c(1,2),4)))

rownames(coldata) <- colnames(counts)

# Create the full model for comparison of samples

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~celltype+timepoint+celltype:timepoint) 

# Generate a linear model

dds <- DESeq(dds, test="LRT", reduced=~celltype+timepoint)

resultsNames(dds)

# These coefficients will allow us to ask multiple questions about the impact
# of timepoints and celltype on differential expression

# Difference between sel and control at T48 : celltype_selected_vs_ctrl
# Difference between sel and ctrl at T96: celltype_selected_vs_ctrl and celltypeselected.timepoint96
# Difference between T48 and T96 in selected: timepoint_96_vs_48 and celltypeselected.timepoint96
# Difference between T48 and T96 in control: timepoint_96_vs_48

## Checking distribution

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

# Checking PCA

rld <- rlogTransformation(dds)

p1 <- plotPCA(rld,intgroup="celltype") 
p2 <- plotPCA(rld,intgroup="timepoint") # ctrl and sel at 48h are very similar
p3 <- plotPCA(rld,intgroup="replicate") # great replicates damn

plot_grid(p1, p2, p3, nrow = 1, align = "v")

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok


## SELECTED VS CTRL at T48------------------------------------------------------

res <- results(dds, name = "celltype_selected_vs_ctrl", test = "LRT")
res_tbl <- as_tibble(res, rownames="ENSMUG")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv")

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj)

write.csv(res_tbl,"./data_output/sel_vs_ctrl_T48/res_tbl.csv")

# Plot for Slc7a5 

as_tibble(counts(dds["ENSMUSG00000040010"], normalize = TRUE),
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
  ggtitle("Slc7a5 mRNA counts")

ggsave(filename = "./figures/sel_vs_ctrl_T48/Slc7a5_mRNA_count.png", plot = last_plot())

hist(res_tbl$padj) # distribution looks ok, if not filtering is needed

plotMA(res) # looks a bit weird but expected since the conditions don't differ much


# Volcano plot with labels

library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

write_csv(overexpr,"./data_output/sel_vs_ctrl_T48/Sign_overexpr.csv")
write_csv(underexpr,"./data_output/sel_vs_ctrl_T48/Sign_underexpr.csv")
write_csv(signif,"./data_output/sel_vs_ctrl_T48/Sign_different.csv")

library("readxl")

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/sel_vs_ctrl_T48/RAB741_sign_genes.xlsx")

# Volcano plot

target <- read_xlsx("./data/candidates.xlsx", sheet = 1)
target <- target$Gene
highlight <- signif %>% filter(gene%in%target)

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
  labs(title = "Selected vs control at T48")+
  theme_bw()

ggsave("./figures/sel_vs_ctrl_T48/volcanoplot", last_plot(), device = png, dpi= 500)


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
  labs(title = "Selected vs control at T48")

ggsave("./figures/sel_vs_ctrl_T48/volcanoplot_targets", last_plot(), device = png, dpi= 500)

# save the tibbles etc

saveRDS(dds, file = "./data_output/sel_vs_ctrl_T48/RAB741_dds_ctrl_sel_T48.rds")
saveRDS(res_tbl, file = "./data_output/sel_vs_ctrl_T48/RAB741_res_tbl_ctrl_sel_T48.rds")
write_csv(res_tbl, file = "./data_output/sel_vs_ctrl_T48/res_tbl.csv")

## SELECTED VS CTRL at T96------------------------------------------------------

res <- results(dds, contrast = list(c("celltype_selected_vs_ctrl",
                                      "celltypeselected.timepoint96")),test = "LRT")

res_tbl <- as_tibble(res, rownames="ENSMUG")

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj)

write.csv(res_tbl,"./data_output/sel_vs_ctrl_T96/res_tbl.csv")


# Plot for Slc7a5 

as_tibble(counts(dds["ENSMUSG00000040010"], normalize = TRUE),
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
  ggtitle("Slc7a5 mRNA counts")

ggsave(filename = "./figures/sel_vs_ctrl_T96/Slc7a5_mRNA_count.png", plot = last_plot())

hist(res_tbl$padj) # distribution looks ok, if not filtering is needed

plotMA(res) # looks ok


# Volcano plot with labels

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

write_csv(overexpr,"./data_output/sel_vs_ctrl_T96/Sign_overexpr.csv")
write_csv(underexpr,"./data_output/sel_vs_ctrl_T96/Sign_underexpr.csv")
write_csv(signif,"./data_output/sel_vs_ctrl_T96/Sign_different.csv")

library("openxlsx")

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/sel_vs_ctrl_T96/RAB741_sign_genes.xlsx")

# Volcano plot

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
  geom_text_repel(max.overlaps = 10)+
  labs(title = "Selected vs control at T96")+
  theme_bw()

# save the tibbles etc

saveRDS(dds, file = "./data_output/sel_vs_ctrl_T96/RAB741_dds_ctrl_sel_T96.rds")
saveRDS(res_tbl, file = "./data_output/sel_vs_ctrl_T96/RAB741_res_tbl_ctrl_sel_T96.rds")
write_csv(res_tbl, file = "./data_output/sel_vs_ctrl_T96/res_tbl.csv")


## Targeted search based on gene list from Raphaele-----------------------------

library(readxl)
library(ggrepel)
target <- read_xlsx("./data/candidates.xlsx", sheet = 1)
target <- target$Gene

highlight <- signif %>% filter(gene%in%target)
# 10 genes from her list seem to be sign upregulated at 96 in selected cells
# (which corresponds to those without tryptophan)

# Volcano plot

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
  theme(legend.title= element_blank())

ggsave("./figures/sel_vs_ctrl_T96/volcanoplot_targets", last_plot(), device = png, dpi= 500)


## TIMEPOINT IN-BETWEEN SAMPLES ANALYSIS----------------------------------------

# Difference between T48 and T96 in control: timepoint_96_vs_48-----------------

res <- results(dds, name = "timepoint_96_vs_48", test = "LRT")
res_tbl <- as_tibble(res, rownames="ENSMUG")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv")

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj)

write.csv(res_tbl,"./data_output/ctrl_96_vs_48/res_tbl.csv")

hist(res_tbl$padj) # distribution looks ok, if not filtering is needed

plotMA(res) # looks a bit weird but expected since paired conditions?

# Volcano plot with labels

library(ggrepel)

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

write_csv(overexpr,"./data_output/ctrl_96_vs_48/Sign_overexpr.csv")
write_csv(underexpr,"./data_output/ctrl_96_vs_48/Sign_underexpr.csv")
write_csv(signif,"./data_output/ctrl_96_vs_48/Sign_different.csv")

library("openxlsx")

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/ctrl_96_vs_48/RAB741_sign_genes.xlsx")

# Volcano plot

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
  geom_text_repel(max.overlaps = 20)+
  labs(title = "Control cells T96 vs T48")+
  theme_bw()

ggsave("./figures/ctrl_96_vs_48/Volcanoplot_ctrl_96_vs_48", 
       last_plot(), device = png, dpi = 400)

# Highlight targets

highlight <- signif %>% filter(gene%in%target)

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
                  max.overlaps = 30, min.segment.length = 0, color = "blue",
                  box.padding = 0.5)+
  theme(legend.title= element_blank())

ggsave("./figures/ctrl_96_vs_48/Volcanoplot_highlights", last_plot(), device = png, dpi= 500)

# save the tibbles etc

saveRDS(dds, file = "./data_output/ctrl_96_vs_48/RAB741_dds_ctrl_sel_T48.rds")
saveRDS(res_tbl, file = "./data_output/ctrl_96_vs_48/RAB741_res_tbl_ctrl_sel_T48.rds")
write_csv(highlight, file = "./data_output/ctrl_96_vs_48/highlights_ctrl_96vs48.csv")
write_csv(res_tbl, file = "./data_output/ctrl_96_vs_48/res_tbl.csv")


# Difference between T96 and T48 in selected:-----------------------------------
# timepoint_96_vs_48 and celltypeselected.timepoint96

res <- results(dds, contrast = list(c("timepoint_96_vs_48","celltypeselected.timepoint96")), 
               test = "LRT")
res_tbl <- as_tibble(res, rownames="ENSMUG")

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv")

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj)

write.csv(res_tbl,"./data_output/sel_96_vs_48/res_tbl.csv")

hist(res_tbl$padj) # distribution looks ok, if not filtering is needed

plotMA(res) # looks a bit weird but expected since paired conditions?

# Volcano plot with labels

library(ggrepel)

# highlight_list <- c("PRKAA1","PRKAA2","PRKAB1","PRKAB2","STK11","TP53")
# res_tbl <- res_tbl %>%
#   mutate(plotname = ifelse(plotname %in% highlight_list, plotname, ""))

overexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange>=1) 
underexpr <- res_tbl %>% filter(padj<=0.05 & log2FoldChange<=-1)
signif <- full_join(overexpr, underexpr)

write_csv(overexpr,"./data_output/sel_96_vs_48/Sign_overexpr.csv")
write_csv(underexpr,"./data_output/sel_96_vs_48/Sign_underexpr.csv")
write_csv(signif,"./data_output/sel_96_vs_48/Sign_different.csv")

library("openxlsx")

list_of_datasets <- list("UP" = overexpr, "DOWN" = underexpr,"ALL_sign" = signif)
write.xlsx(list_of_datasets, file = "./data_output/sel_96_vs_48/RAB741_sign_genes.xlsx")

# Volcano plot

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
  geom_text_repel(max.overlaps = 20)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()

ggsave("./figures/sel_96_vs_48/Volcanoplot_sel_96_vs_48", 
       last_plot(), device = png, dpi = 400)

# Highlight targets

highlight <- signif %>% filter(gene%in%target)

# save the tibbles etc

saveRDS(dds, file = "./data_output/sel_96_vs_48/RAB741_dds_ctrl_sel_T48.rds")
saveRDS(res_tbl, file = "./data_output/sel_96_vs_48/RAB741_res_tbl_ctrl_sel_T48.rds")
write_csv(res_tbl, file = "./data_output/sel_96_vs_48/res_tbl.csv")







