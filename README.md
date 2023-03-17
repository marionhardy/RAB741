# RAB741

Repository for RAB in BVDE lab (RNAseq analyses)

## Experimental design

Mice T cells at T0h, in presence or absence of tryptophan in the medium.
Total RNA was extracted using Trizol at T0h, 48h and 96h after treatment.
Sequencing was done by Illumina Seq, 10^6 depth of sequencing? (need to ask for more info)

## Scripts used:

- 01_RSubread
- 02_DESeq2_file_preparation
- 03_DESeq2_timepoints
- 04_RAB741_report

## Supplementary scripts

- 05_1_Biomart_annotations
- 05_2_Improvements
- 05_3_Transform_tsv_col
- 05_4_Venn

## Breakdown of script content

A detailed list of the use of the different scripts used on this project at Raphaële's request.
NB: Parts of these data and scripts were used in the learningR repository as examples for RNAseq analysis using R.

### 01_RSubread

1. Build human and mouse genome index (GRC38)
2. Contains the .exe script for running trim_galore (bash)
3. Align fasta files with RSubread
4. Quantify reads with featureCounts()
5. Annotate the resulting gene counts

### 02_DESeq2_file_preparation

1. Merge all counts from all conditions
2. Create the coldata for the DESeq2 object 

### 03_DESeq2_timepoints

1. Merge all counts from all conditions
2. Create the coldata for the DESeq2 object 

NB: the above part is redundant since it's also comprised in the previous script.

3. Create the DESeq2 object and the linear model:
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~celltype+timepoint+celltype:timepoint) 
dds <- DESeq(dds, test="LRT", reduced=~celltype+timepoint)

Which allows to have information on how timepoints and cell types and treatment impact affect the T cell transcription.
4. Quality control check using size factors, dispersion and PCA.
5. Extracting results for different sample comparisons
6. Volcano plots
7. Annotating genes using biomart

### 04_RAB741_report

GSEA analysis using gene lists ranked on logFC, test statistic or pvalue and querying:
- Gene Onology (MF,CC,BP)
- Biocarta
- Hallmark
- Reactome
- KEGG
- Wikipathways
- TFT (regulatory genes)

## Supplementary scripts

- 05_1_Biomart_annotations: Script to fetch annotations from biomart
- 05_2_Improvements: Script listing some improvements that were implemented or not following a discussion on the RAB741 report
- 05_3_Transform_tsv_col: Script to fetch the genes implicated in a pathway (from GSEA object) and output a list with a gene per row in a txt file
- 05_4_Venn: Script to make a venn diagram with more than three intersecting conditions

