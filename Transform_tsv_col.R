

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/HALLMARK_GLYCOLYSIS.v2022.1.Mm.tsv")

hgluc$HALLMARK_GLYCOLYSIS[[19]]
gluc <- hgluc$HALLMARK_GLYCOLYSIS[[19]]

gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/HALLMARK_GLYCOLYSIS.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/HALLMARK_MTORC1_SIGNALING.v2022.1.Mm.tsv")

hgluc$HALLMARK_MTORC1_SIGNALING[[19]]
gluc <- hgluc$HALLMARK_MTORC1_SIGNALING[[19]]

gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/HALLMARK_MTORC1_SIGNALING.txt",
            sep = "", row.names = F, quote = F)

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2022.1.Mm.tsv")

hgluc$HALLMARK_OXIDATIVE_PHOSPHORYLATION[[19]]
gluc <- hgluc$HALLMARK_OXIDATIVE_PHOSPHORYLATION[[19]]

gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/HALLMARK_OXPHOS.txt",
            sep = "", row.names = F, quote = F)

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v2022.1.Mm.tsv")

hgluc$HALLMARK_PI3K_AKT_MTOR_SIGNALING[[19]]
gluc <- hgluc$HALLMARK_PI3K_AKT_MTOR_SIGNALING[[19]]

gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/HALLMARK_PI3K_AKT_MTOR.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_GLUCONEOGENESIS.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_GLUCONEOGENESIS[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_GLUCONEOGENESIS.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_GLUCOSE_METABOLISM.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_GLUCOSE_METABOLISM[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_GLUCOSE_METABOLISM.txt",
            sep = "", row.names = F, quote = F)

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_GLYCOGEN_METABOLISM.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_GLYCOGEN_METABOLISM[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_GLYCOGEN_METABOLISM.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCHONDRIAL_BIOGENESIS.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_MITOCHONDRIAL_BIOGENESIS[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_BIOGENESIS.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_CALCIUM_ION_TRANSPORT.txt",
            sep = "", row.names = F, quote = F)

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_BETA_OX_SATURATED_FAT.txt",
            sep = "", row.names = F, quote = F)

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_BETA_OX.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_PROT_IMPORT.txt",
            sep = "", row.names = F, quote = F)

hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCHONDRIAL_UNCOUPLING.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_MITOCHONDRIAL_UNCOUPLING[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_UNCOUPLING.txt",
            sep = "", row.names = F, quote = F)


hgluc <- read_tsv("//Users/hardy/Documents/R_2022/databases/REACTOME_RELEASE_OF_APOPTOTIC_FACTORS_FROM_THE_MITOCHONDRIA.v2022.1.Mm.tsv")
gluc <- hgluc$REACTOME_RELEASE_OF_APOPTOTIC_FACTORS_FROM_THE_MITOCHONDRIA[[19]]
gluc <- strsplit(gluc, ",")
gluc <- data.frame(gluc)
colnames(gluc) <- "X"

write.table(gluc,"//Users/hardy/Documents/R_2022/databases/REACTOME_MITOCH_RELEASE_APOPTOTIC_FACTORS.txt",
            sep = "", row.names = F, quote = F)




















