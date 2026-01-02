# ssGSEA Analysis
# Peter Fiorica
#9 April 2024 # Edited for ARC Data on 1.2.2025

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(UpSetR))
library(msigdbr)
library(GSVA)
"%&%"=function(a,b) paste(a,b,sep="")
setwd("/projects/rpci/songyao/pnfioric/arc_project/")

arc_dir <- "/projects/rpci/songyao/pnfioric/arc_project/"
tils_dir <- "/projects/rpci/songyao/pnfioric/tils_r01_RNAseq/"
# Read in data...

#expected_count<-fread("count_data/all_samples_expected_count.txt")
#sample_info <- fread("Phenotype_Sample_Info_wchs_nsbcs_pathways_w_tils_and_origsampetype.txt")  %>% filter(stage!="0")
#TPM <- fread("count_data/all_samples_TPM.txt")

arc_dds <- readRDS(arc_dir %&% "AsA_EA_Samples_GongYao.rds")
# Pull Sample Names
coldata <- as.data.table(colData(arc_dds))

# Pull TPM information
tpm <- assays(arc_dds)[["TPM"]]
colnames(tpm) <- coldata$FASTQ_RS
rownames(tpm) <- sub("\\..*", "", rownames(tpm))

#Biomart Ensembl Into
gene_info <- fread(tils_dir %&% "mart_export_gene_type.txt", header =T)

#gene_list <- expected_count_bc_filt$gene_id
ref_data<-fread(tils_dir %&% "NCBI_Genome_Annotation_Homosapiens_GRCh38_04162024.tsv", header = T)

# Log2 Transformation:
TPM_filt <- tpm
log2_tpm_filt <- log2(TPM_filt+1)

log_gene_exp <- log2_tpm_filt

library(GSVA)
library(msigdbr)

msigdb_immune_sig<- msigdbr(species = "Homo sapiens",
                            category = "C7",
                            subcollection = "IMMUNESIGDB")

immune_genesets <- split(msigdb_immune_sig$gene_symbol, msigdb_immune_sig$gs_name)


log2_tpm_names <- log2_tpm_filt %>% as.data.frame() %>%
  rownames_to_column(var = "Gene ENSG") %>%
  left_join(gene_info, by = c("Gene ENSG"="Gene stable ID")) #`Gene name`

log2_tpm_named <- log2_tpm_names %>%
  select(`Gene name`, all_of(colnames(log2_tpm_filt))) %>%  # Keep only gene name + expression columns
  filter(!is.na(`Gene name`)) %>%                           # Drop rows without gene symbol
  distinct(`Gene name`, .keep_all = TRUE) %>%               # Drop duplicated gene symbols
  column_to_rownames(var = "Gene name")                     # Set rownames



# Convert to matrix of numeric values
log2_tpm_mat <- as.matrix(sapply(log2_tpm_named, as.numeric))
rownames(log2_tpm_mat) <- rownames(log2_tpm_named)

print("Beginning ssGSEA")
ssgsea_results <- gsva(expr = log2_tpm_mat,
                       gset.idx.list = immune_genesets,
                       method = "ssgsea",
                       kcdf = "Gaussian",  # log2 TPMs → Gaussian
                       min.sz = 10,
                       max.sz = 500,
                       mx.diff = TRUE,
                       verbose = TRUE)

saveRDS(ssgsea_results, "ARC_ssGSEA.rds")

fwrite(ssgsea_results, "ARC_ssgsea_results.txt", col.names = T, row.names = T, sep = "\t")
