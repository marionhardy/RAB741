## Test for Venn diagrams to have an overview of what changes from T48 to T96
# in selected cells and in control and whether they have things that change in 
# common and others not

library(tidyverse)
library(ggvenn)
library(ComplexHeatmap)
library(UpSetR)

dir <- getwd()
list.files(paste0(dir,"/data_output/"))

sc48 <- read_csv(paste0(dir,"/data_output/sel_vs_ctrl_T48/Sign_different.csv"))
sc96<- read_csv(paste0(dir,"/data_output/sel_vs_ctrl_T96/Sign_different.csv"))
c9648<- read_csv(paste0(dir,"/data_output/ctrl_96_vs_48/Sign_different.csv"))
s9648<- read_csv(paste0(dir,"/data_output/sel_96_vs_48/Sign_different.csv"))

overlap <- list(
  "control_96vs48" = c9648$gene,
  "selected_96vs48" = s9648$gene,
  "selected_vs_control_48" = sc48$gene,
  "selected_vs_control_96"= sc96$gene)

ggvenn(overlap, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 3, text_size = 3)

ggsave("./figures/Venndiagram_compare_all", last_plot(), device = png, dpi = 400)






































