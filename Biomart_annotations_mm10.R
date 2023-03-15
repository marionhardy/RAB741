## Annotations

listDatasets(useMart("ensembl"))

# Load mus musculus ensembl dataset

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

listAttributes(mart)

ensembl_to_geneName <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                            "entrezgene_id","description"),
                             mart = mart)
names(ensembl_to_geneName) <- c("ENSMUG", "gene", "ENTREZID", "description")
head(ensembl_to_geneName)






