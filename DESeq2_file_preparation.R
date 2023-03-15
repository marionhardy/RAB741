
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

data()

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

write.csv(coldata,"./data/coldata.csv")


