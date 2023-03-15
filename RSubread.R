###################################
##    FASTQ to BAM  to counts   ##
##################################

library("Rsubread")
library("dplyr")
library("tidyverse")

# RsubreadUsersGuide() # For the complete user guide

dir <- getwd()                    

# Build indexes for mouse and homo sapiens genomes (only once)

# buildindex(basename="mouse_index",
#            reference="/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/GRCm38.primary_assembly.genome.fa",
#            memory = 8000, indexSplit = T)

# buildindex(basename="hsapiens_index",reference="./Genome_ref/GRC38.primary_assembly.genome.fa")

# Quality control of fastq files using trim_galore (cutadapt + fastqc)
# this must be run on the terminal
# open the fastq folder containing all folder (themselves containing the 2 fastq files)
# Finder>Services>Open terminal> then run the following command

# for d in */; do trim_galore --illumina --paired --fastqc -o $d/trim_galore/ $d/*.fq; done

# Import your fastq or fasta files

annot <- read_csv("/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/annotations.csv")

# Add gene annotations

# annot <- read_tsv("/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/Mus_musculus.GRCm38.102.gtf.gz", 
#                          col_names = F, skip = 5, col_types = "ccciicccc") %>% 
#    dplyr::filter(X3 == "gene") %>%
#    separate(X9, into = c("X9", "X10", "X11", "X12", "X13", "X14"), sep = ";") %>%
#    dplyr::select(1,4,5,7,9,11,13)
# 
# annot$X9 <- str_remove(annot$X9, "gene_id ")
# annot$X9 <- gsub ('"', '', annot$X9)
# annot$X11 <- str_remove(annot$X11, " gene_name ")
# annot$X11 <- gsub ('"', '', annot$X11)
# annot$X13 <- str_remove(annot$X13, " gene_biotype ")
# annot$X13 <- gsub ('"', '', annot$X13)
# 
# colnames(annot) <- c("chrom", "start", "end", "strand", "ensembl_gene_id", "hgnc_symbol", "gene_biotype")
# 
# head(annot)

# write_csv(annot, "/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/annotations.csv")
# This file is now available in the mus musculus folder

x <- list.dirs(paste0(dir,"/data/fastq"), recursive = F, full.names = TRUE)

for (i in 1:length(x)){

   reads1 <- list.files(path = paste0(x[i],"/trim_galore"), pattern = "*_1.fq$")
   reads2 <- list.files(path = paste0(x[i],"/trim_galore"), pattern = "*_2.fq$")
   name <- reads1
   name <- substring(reads1,1, nchar(name)-11)
   
   align(index = "/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/mouse_index",
         readfile1 = paste0(x[i],"/trim_galore/",reads1),
         readfile2 = paste0(x[i],"/trim_galore/",reads2),
         output_file = paste0("./data/bam/",name,".bam"),
         nthreads = 6,
         sortReadsByCoordinates = TRUE,
         nBestLocations = 1,
         unique = FALSE,
         annot.ext = "/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/Mus_musculus.GRCm38.102.gtf.gz",
         isGTF = T)

   name <- paste0(x[i])
   name <- substring(paste0(x[i]),58,nchar(name))

   # quantify number of reads mapping each gene

   fc <- featureCounts(files = paste0(dir,"/data/bam/", name, ".bam"),
                       annot.ext = "/Users/hardy/Documents/R_2022/Genome_ref/GRCh38_mus_musculus/Mus_musculus.GRCm38.102.gtf.gz",
                       annot.inbuilt = "mm10",
                       isPairedEnd = TRUE,
                       isGTFAnnotationFile = TRUE,
                       useMetaFeatures = TRUE,
                       nthreads = 6,
                       fraction = TRUE)
   write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],
                            fc$counts,stringsAsFactors=FALSE),file="counts.txt",
               quote=FALSE,sep="\t",row.names=FALSE)

    # create tibble with per-gene raw counts and normalized RPKM and TPM values (not here though)
   
   df <- tibble(fc$annotation[, c("GeneID", "Length")]) %>%
      add_column(RawCounts = fc$count[,1]) %>%
      dplyr::rename(ensembl_gene_id = GeneID)
   
   df$RawCounts <- as.double(str_replace(df$RawCounts, ",", "."))
   
   df$ensembl_gene_id <- as.character(df$ensembl_gene_id)
   
   # No normalization needed if subsequently using DESeq2 
   
   # scale_rpkm <- sum(df$RawCounts) / 1000000
   # 
   # df <- df %>%
   #    mutate(RPKM = RawCounts / scale_rpkm / Length * 1000) %>%
   #    mutate(TPM = RawCounts / Length * 1000)
   # 
   # scale_tpm <- sum(df$TPM) / 1000000
   # 
   # df$TPM <- df$TPM / scale_tpm
   
   df <- left_join(df, annot, by = "ensembl_gene_id") %>%
      arrange(hgnc_symbol)

   # save data table with quantified gene expression levels
   write_csv(df, paste0(name, "_counts.csv"))

}




