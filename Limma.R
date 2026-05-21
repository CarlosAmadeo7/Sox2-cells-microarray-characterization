#' ---
#' title: "Microarray analysis on Sox2+ vs Sox2- cells"
#' output: html_notebook
#' ---
#' # Microarray analysis on Sox2+ vs Sox2- cells from: https://www.sciencedirect.com/science/article/pii/S1535610814002207?via%3Dihub#fig3
## ----Loading libraries---------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2);library(ggrepel);library(ggpubr);library(data.table);library(RColorBrewer)
  ;library(tidyverse);library(preprocessCore);library(future.apply);library(DESeq2)
  ;library(pheatmap);library(sva);library(viridis);library(limma);library(emmeans);library(org.Mm.eg.db)
  ;library(broom);library(janitor);library(tidyplots);library(writexl);library(scToppR);library(enrichplot)
  ;library(fgsea);library(clusterProfiler);library(org.Hs.eg.db);library(dplyr); library(biomaRt); library(ReactomePA)
  ;library(tibble);library(oligo);library(mogene20sttranscriptcluster.db);library(AnnotationDbi);library(readxl); 
  library(sessioninfo)})

## ----Pre processing------------------------------------------------------------------------------------------------
## Loading files 
celFiles <- list.celfiles("data/", full.names = TRUE)
rawData <- read.celfiles(celFiles)

str(rawData)
pData(rawData) ## Phenotype dataset 
head(exprs(rawData)) ### gene expression raw counts 
experimentData(rawData) ### Experimental info 
featureData(rawData) ## Probeset dataset 

# obtaining raw counts
raw_counts<- exprs(rawData)
boxplot(exprs(rawData))

## Removing zeros on the expression matrix
filter <-apply(exprs(rawData), 1, function(x) all(x[1:4]>=0.5) | all(x[5:8]>=0.5))
sum(filter)
summary(exprs(rawData)[, 1:4])
counts<- exprs(rawData)[filter,]

dim(exprs(rawData))
dim(counts)

#----------------------
## Normalizing dataset
normData<- rma(rawData)
head(exprs(normData))
pData(normData)

#-----------------------
# Metadata
condition <- c("Sox2_positive", "Sox2_positive", "Sox2_positive","Sox2_positive",
               "Sox2_negative","Sox2_negative","Sox2_negative","Sox2_negative")
pData(normData)$Condition <- condition
pData(normData)
metadata<- data.frame(
  Filename = c("GSM1184345_MoGene2_061313H_PD1_P1.CEL","GSM1184346_MoGene2_061313H_PD2_P2.CEL", 
             "GSM1184347_MoGene2_061313H_PD3_P3.CEL", "GSM1184348_MoGene2_061313H_PD4_P4.CEL",
             "GSM1184349_MoGene2_061313H_PD5_N1.CEL", "GSM1184350_MoGene2_061313H_PD6_N2.CEL",
             "GSM1184351_MoGene2_061313H_PD7_N3.CEL", "GSM1184352_MoGene2_061313H_PD8_N4_2.CEL"),
  Condition = c("Sox2+", "Sox2+", "Sox2+","Sox2+",
                "Sox2-","Sox2-","Sox2-","Sox2-"))

write.table(metadata, file = "metadata.txt", sep = "\t", row.names = F)
#metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## PCA on normalized expression 
p <-exprs(normData)
pca.Samples<-prcomp(t(p))
PCi<-data.frame(pca.Samples$x, Condition= metadata$Condition, ID= metadata$Filename )
eig <- (pca.Samples$sdev)^2 
variance <- eig*100/sum(eig) 
#pdf("dge/PCA_cpms.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(PCi,
          x = "PC1",
          y = "PC2",
          color = "Condition", palette=c("red","black"),
          shape = "Condition", size = 4,label = "Condition") + 
  xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
  ylab(paste("PC2 (",round(variance[2],1),"% )"))+theme_classic()
#dev.off()

#' # Differential expression analysis using limma
## ----Limma---------------------------------------------------------------------------------------------------------
metadata<- metadata |> column_to_rownames("Filename")
metadata$Condition<- as.factor(metadata$Condition)
levels(metadata$Condition)<- c("Sox2+","Sox2-")
all(colnames(exprs(normData)) == rownames(metadata))

#------
# Limma 
# Designing the matrix, two intercepts per condition so we have control of the reference- control 
tmp <- pData(normData)
tmp$Condition<- as.factor(tmp$Condition)
#levels(tmp$Condition)<- c("Sox2+","Sox2-")
levels(tmp$Condition)

# Design formula
design <- model.matrix(~tmp$Condition -1)
colnames(design) <- c("Sox2_negative", "Sox2_positive")

fit <- lmFit(normData,design)

contrast.matrix <- makeContrasts("Sox2_positive-Sox2_negative", levels = design)
contrast.matrix

# Fitting contrasts
fitC <- contrasts.fit(fit, contrast.matrix)
all(rownames(contrast.matrix) == colnames(fit$coefficients))
# Empirical Bayes correction 
fitC <- eBayes(fitC)
topTable(fitC)

#' 
## ----ID Probes-----------------------------------------------------------------------------------------------------
# Obtaining full list of probs 
FullList <- topTable(fitC, number = Inf, sort.by = "P")

# Get gene symbols for your probes
FullList$Symbol <- mapIds(mogene20sttranscriptcluster.db,keys = rownames(FullList),
  column = "SYMBOL",keytype = "PROBEID",multiVals = "first")
str(FullList)
head(FullList)

# Filtering protein coding genes 
#-------------------------------
listAttributes(mart)
# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# Query for coding genes
coding_genes <- getBM(attributes = c("external_gene_name", "gene_biotype"),filters = "external_gene_name",
  values = FullList$Symbol,mart = mart)
# Adding a list of protein coding genes in the FullList of symbols
#FullList$ProteinCoding<- FullList$Symbol %in% coding_genes$external_gene_name
colnames(coding_genes) <- c("Symbol", "GeneBiotype")
FullList1<- FullList |> left_join(coding_genes, by= "Symbol")
# Filtering the full list of protein coding genes 
FullList_pc <- FullList1 |>filter(GeneBiotype == "protein_coding") |> dplyr::mutate(ABS = abs(logFC))
# List of DEGs
DEGs<- FullList_pc |> dplyr::filter(adj.P.Val < 0.05 & ABS > 1)
# Saving the matrices, creating volcano plot, do GSEA and gene ontology. 
write_xlsx(FullList,"deg/CTX_FullList.xlsx" )
write_xlsx(FullList_pc, "deg/FullList_protein_coding_genes.xlsx")
write_xlsx(DEGs, "deg/CTX_DGE.xlsx")
save(FullList, FullList_pc, metadata,DEGs, file = "deg/CTX_Dge_Data.RData")

#' 
## ----Volcano plot--------------------------------------------------------------------------------------------------
df <- FullList_pc %>%
  mutate(LOG = -log10(adj.P.Val)) %>%
  mutate(Threshold = if_else(adj.P.Val < 0.05 & ABS > 1, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "UpReg", logFC < -1 & adj.P.Val < 0.05 ~ "DownReg")) 

top_labelled <- df %>%group_by(Direction) %>%na.omit() %>%arrange(adj.P.Val) %>%top_n(n = 7, wt = LOG)

pdf("deg/Volcano_Plot_CTX.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df,x="logFC", y="LOG", color = "Direction",  palette=c("dodgerblue1","firebrick1"),
          size = 2, alpha=0.3, shape=19) +
  xlab("log2(Fold Change)")+ ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 1, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -1, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_point(data = top_labelled, aes(x = logFC, y = LOG), shape = 21, fill = "transparent", color = "black", size = 2, stroke = 0.5) + geom_text_repel(data = top_labelled,
                  mapping = aes(label = Symbol),size = 5, box.padding = unit(0.4, "lines"), point.padding = unit(0.4, "lines")) + theme(legend.position="none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 14))+ ylim(0,10) + xlim(-2.5,+5)
dev.off()

#' 
#' Downstream analysis
#' #---------------------
## ----Gene ontology-------------------------------------------------------------------------------------------------
tmp1<- DEGs |> dplyr::mutate(Direction = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "UpReg", logFC < -1 & adj.P.Val < 0.05 ~ "DownReg"))

# Temporary lists for gene ontology 
Up.Down.GO.list<-tmp1 |> pull(Symbol)
Up.GO.list<-tmp1 |>dplyr::filter(Direction == "UpReg") |> pull(Symbol)
Down.GO.list<-tmp1 |>dplyr::filter(Direction == "DownReg") |> pull(Symbol)

## Gene ontology function
GO_analysis <- function(gene_list, condition, status) {
  ## Creating directories
  base_path<-paste0("dge/GO/", condition, "/", status)
  if (!dir.exists(base_path)) {dir.create(base_path, recursive = TRUE)}
  set.seed(0708)
  # Gene ontology 
  GO <- enrichGO(gene = gene_list,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",
    ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,minGSSize = 10, maxGSSize = 500 )
  # Save results as CSV
  write.csv(as.data.frame(GO), file = paste0(base_path, "/", status, ".BP.GO.csv"))
  # Creating Dotplots for each analysis
  pdf(file = paste0(base_path, "/", status, ".BP.GO.pdf"), width = 9, height = 5, useDingbats = FALSE)
  print(dotplot(GO, label_format = 100, showCategory = 20) + ggtitle(paste(status, condition, sep = " - ")) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  dev.off()

  # creating enrichment plot
  enrich_sim <- pairwise_termsim(GO)
  pdf(file = paste0(base_path, "/", status, ".BP.Enrichment.pdf"), width = 9, height = 5, useDingbats = FALSE)
  print(emapplot(enrich_sim, showCategory = 20) +ggtitle(paste(status, "Enrichment map", condition, sep = " - ")) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  dev.off()
  # Clean memory 
  rm(GO, enrich_sim)
  gc()}

# Running Gene ontology analysis 
# Run enrichment analysis for All, Up, and Down DEGs
GO_analysis(Up.Down.GO.list, "Sox2+", "All")
GO_analysis(Up.GO.list, "Sox2+", "Up")
GO_analysis(Down.GO.list, "Sox2+", "Down")

#' 
## ----Gene Set Enrichment analysis----------------------------------------------------------------------------------
run_bulk_gsea <- function(df,pathway_list,  gene_col = "Symbol",  stat_col = "logFC", output_dir = "deg/GSEA/", 
                          nperm = 10000, label = "bulkRNAseq") {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("Running GSEA for: ", label)
  # Filter and rank genes
  ranks <- df %>%filter(!is.na(!!sym(stat_col))) %>% 
    distinct(!!sym(gene_col), .keep_all = TRUE) %>% arrange(desc(!!sym(stat_col)))
  stat_vector <- ranks[[stat_col]]
  names(stat_vector) <- ranks[[gene_col]]
  # Run fgsea
  set.seed(0708)
  fgsea_res <- fgsea(pathways = pathway_list, stats = stat_vector,nperm = nperm,minSize = 10,maxSize = 500)
  # Save raw results
  fgsea_res <- as_tibble(fgsea_res)
  fgsea_res$Label <- label
  write_csv(fgsea_res, file.path(output_dir, paste0(label, "_GSEA.csv")))
  # Save significant pathways
  up_filt <- fgsea_res %>% filter(padj < 0.05, ES > 0)
  down_filt <- fgsea_res %>% filter(padj < 0.05, ES < 0)
  write_csv(up_filt, file.path(output_dir, paste0(label, "_GSEA_up.csv")))
  write_csv(down_filt, file.path(output_dir, paste0(label, "_GSEA_down.csv")))
  return(fgsea_res)
}

#### Running function 
hallmark_pathways<-gmtPathways("deg/GSEA/mh.all.v2024.1.Mm.symbols.gmt")
res<-run_bulk_gsea(df = FullList_pc,pathway_list = hallmark_pathways,gene_col = "Symbol",
              stat_col = "logFC",  output_dir = "deg/GSEA/",label = "Sox2+")

# Saving res ouput 
res %>%  arrange(padj) %>%  head(50)
save(res, file = "deg/GSEA/res_output.RData")
# Saving GSEA output 
gsea_output1<-res |> as.data.frame()
gsea_output1$leadingEdge <- sapply(gsea_output1$leadingEdge, function(x) paste(x, collapse = ", "))
write.csv(gsea_output1, "deg/GSEA/Full_genes_list_GSEA.csv")

#' 
## ----Plotting GSEA output------------------------------------------------------------------------------------------
plot_bulk_gsea_barplot <- function(gsea_df, output_dir = "deg/GSEA/barplots", label = "bulkRNAseq") {
  # Create output dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Preprocess for plotting
  tpm <- gsea_df %>% mutate(adjPvalue = ifelse(padj <= 0.05, "significant", "non-significant"), pathway = gsub("HALLMARK_", "", pathway))  
  # Define colors
  cols <- c("non-significant" = "black", "significant" = "#CD4071FF") 
  # Create barplot for downregulated pathways (ES < 0)
  p <- ggbarplot(tpm %>% filter(ES > 0),x = "pathway", y = "NES",fill = "adjPvalue",color = "white",
                 palette = cols,sort.val = "asc",sort.by.groups = FALSE,x.text.angle = 90,
                 ylab = "NES",legend.title = "FDR < 0.05",rotate = TRUE,ggtheme = theme_minimal()) +
    theme(
      axis.text.x = element_text(size = 18, color = "black"),
      axis.title.x = element_text(size = 16, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 16, color = "black", face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14))
  # Save plot
  ggsave(filename = file.path(output_dir, paste0(label, "_Upregulated_GSEA.jpeg")),plot = p,
    width = 9.5,height = 9,dpi = 600,units = "in")}

# running function 
plot_bulk_gsea_barplot(res, output_dir = "deg/GSEA/barplots", label = "Sox2+")

#' 
## ----Running GSEA STAT3 significant pathway------------------------------------------------------------------------
load("CTX_Dge_Data.RData")
# Using FULLLIST_Pc
colnames(FullList_pc)
stat_col <- "logFC"
df_no_na <- FullList_pc %>%filter(!is.na(!!sym(stat_col)))
# Keep only one row per gene
df_unique <- df_no_na %>%distinct(Symbol, .keep_all = TRUE)
#Sort genes from highest to lowest statistic
ranks <- df_unique %>%arrange(desc(logFC))
# Extract the statistic column as a numeric vector
stat_vector <- ranks$logFC
names(stat_vector) <- ranks$Symbol

set.seed(7081998)
hallmark_pathways<-gmtPathways("mh.all.v2024.1.Mm.symbols.gmt")
fgsea_res <- fgsea(pathways = hallmark_pathways,stats = stat_vector,nperm = 10000,minSize = 10,maxSize = 500)
fgsea_res %>% filter(grepl("IL6|JAK|STAT3", pathway, ignore.case = TRUE)) %>% dplyr::select(pathway, NES, pval, padj)
il6_pathway <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"
# Extract statistics
nes_value <- fgsea_res$NES[fgsea_res$pathway == il6_pathway]
padj_value <- fgsea_res$padj[fgsea_res$pathway == il6_pathway]
#plotting
gsea_plot <- plotEnrichment(hallmark_pathways[[il6_pathway]],stat_vector) +
  labs(title = "IL6/JAK_STAT3 Signaling",subtitle = paste0("NES = ", round(nes_value, 2)," | FDR = ", signif(padj_value, 3)
    ),x = "Rank in Ordered Gene List",y = "Enrichment Score") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(size = 20,face = "bold",hjust = 0.5),
    plot.subtitle = element_text(size = 13,hjust = 0.5,face = "italic"),axis.title = element_text(size = 16,face = "bold"),
    axis.text = element_text(size = 13,color = "black"),panel.border = element_rect(color = "black",fill = NA,linewidth = 1),axis.line = element_line(linewidth = 0.8),plot.margin = margin(15, 15, 15, 15))

pdf("IL6_STAT3_plotEnrichment.pdf", height = 5)
gsea_plot
dev.off()

#' 
## ----REACTOME analysis---------------------------------------------------------------------------------------------
# REACTOME function
reactome_analysis <- function(gene_list, condition, status) {
  # Create Reactome directory if it doesn't exist
  base_path <- paste0("deg/Reactome/", condition, "/", status)
  if (!dir.exists(base_path)) dir.create(base_path, recursive = TRUE)
  # Convert symbols to gene IDs
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  # Check if conversion returned results
  if (nrow(entrez_ids) == 0) {stop("❌ No valid Entrez IDs found for the input gene symbols.")}
  # Extract only Entrez IDs for Reactome analysis
  entrez_gene_list <- entrez_ids$ENTREZID
  print(head(entrez_ids))  # Print converted IDs to debug
  print(paste("Total Entrez IDs mapped:", nrow(entrez_ids)))
  # Run Reactome pathway enrichment
  reactome <- enrichPathway( gene = entrez_gene_list, organism = "mouse", pAdjustMethod = "BH",pvalueCutoff = 0.05,
    minGSSize = 10,maxGSSize = 500,readable = TRUE)
  # Save Reactome results as CSV
  write.csv(as.data.frame(reactome), file = paste0(base_path, "/", status, ".Reactome.Enrichment.csv"))
  # Plot Reactome Pathway Enrichment
  pdf(file = paste0(base_path, "/", status, ".Reactome.Pathway.pdf"), width = 9, height = 10, useDingbats = FALSE)
  print(dotplot(reactome, showCategory = 20, label_format = 45) +ggtitle(paste("Reactome Pathway Enrichment", status, condition, sep = " - ")) +theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  dev.off()
  return(reactome)}

# Running REACTOME analysis
reactome_analysis(Up.Down.GO.list, "Sox2+", "All")
reactome_analysis(Up.GO.list, "Sox2+", "Up")
reactome_analysis(Down.GO.list, "Sox2+", "Down")

#' 
## ----r session-info------------------------------------------------------------------------------------------------
sessioninfo::session_info()

# Session Information

# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
#   LAPACK version 3.12.1
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] sessioninfo_1.2.3                    readxl_1.5.0                         mogene20sttranscriptcluster.db_8.8.0
#  [4] oligo_1.72.0                         Biostrings_2.78.0                    Seqinfo_1.0.0                       
#  [7] XVector_0.50.0                       oligoClasses_1.70.0                  ReactomePA_1.52.0                   
# [10] biomaRt_2.64.0                       org.Hs.eg.db_3.22.0                  clusterProfiler_4.18.4              
# [13] fgsea_1.36.2                         enrichplot_1.30.4                    scToppR_0.99.4                      
# [16] writexl_1.5.4                        tidyplots_0.4.0                      janitor_2.2.1                       
# [19] broom_1.0.12                         org.Mm.eg.db_3.22.0                  AnnotationDbi_1.72.0                
# [22] emmeans_2.0.3                        limma_3.66.0                         viridis_0.6.5                       
# [25] viridisLite_0.4.3                    sva_3.56.0                           BiocParallel_1.44.0                 
# [28] genefilter_1.90.0                    mgcv_1.9-3                           nlme_3.1-168                        
# [31] pheatmap_1.0.13                      DESeq2_1.48.2                        SummarizedExperiment_1.38.1         
# [34] Biobase_2.68.0                       MatrixGenerics_1.20.0                matrixStats_1.5.0                   
# [37] GenomicRanges_1.60.0                 GenomeInfoDb_1.44.3                  IRanges_2.44.0                      
# [40] S4Vectors_0.48.0                     BiocGenerics_0.54.1                  generics_0.1.4                      
# [43] future.apply_1.20.2                  future_1.70.0                        preprocessCore_1.71.2               
# [46] lubridate_1.9.5                      forcats_1.0.1                        stringr_1.6.0                       
# [49] dplyr_1.2.1                          purrr_1.2.2                          readr_2.2.0                         
# [52] tidyr_1.3.2                          tibble_3.3.1                         tidyverse_2.0.0                     
# [55] RColorBrewer_1.1-3                   data.table_1.18.2.1                  ggpubr_0.6.3                        
# [58] ggrepel_0.9.8                        ggplot2_4.0.2                       
# 
# loaded via a namespace (and not attached):
#   [1] fs_2.0.1                httr_1.4.8              tools_4.5.1             backports_1.5.1        
#   [5] R6_2.6.1                lazyeval_0.2.3          withr_3.0.2             graphite_1.54.0        
#   [9] prettyunits_1.2.0       gridExtra_2.3           cli_3.6.6               scatterpie_0.2.6       
#  [13] sass_0.4.10             mvtnorm_1.3-6           S7_0.2.1                systemfonts_1.3.2      
#  [17] yulab.utils_0.2.4       gson_0.1.0              DOSE_4.4.0              R.utils_2.13.0         
#  [21] dichromat_2.0-0.1       parallelly_1.46.1       rstudioapi_0.18.0       RSQLite_2.4.6          
#  [25] gridGraphics_0.5-1      car_3.1-5               zip_2.3.3               GO.db_3.22.0           
#  [29] Matrix_1.7-3            abind_1.4-8             R.methodsS3_1.8.2       lifecycle_1.0.5        
#  [33] yaml_2.3.12             edgeR_4.8.2             snakecase_0.11.1        carData_3.0-6          
#  [37] qvalue_2.42.0           SparseArray_1.10.8      BiocFileCache_3.0.0     affxparser_1.80.0      
#  [41] grid_4.5.1              blob_1.3.0              crayon_1.5.3            ggtangle_0.1.2         
#  [45] lattice_0.22-7          cowplot_1.2.0           annotate_1.88.0         KEGGREST_1.50.0        
#  [49] pillar_1.11.1           knitr_1.51              estimability_1.5.1      codetools_0.2-20       
#  [53] fastmatch_1.1-8         glue_1.8.0              ggiraph_0.9.6           ggfun_0.2.0            
#  [57] fontLiberation_0.1.0    vctrs_0.7.2             png_0.1-9               treeio_1.34.0          
#  [61] cellranger_1.1.0        gtable_0.3.6            cachem_1.1.0            xfun_0.55              
#  [65] openxlsx_4.2.8.1        S4Arrays_1.10.1         tidygraph_1.3.1         coda_0.19-4.1          
#  [69] survival_3.8-3          iterators_1.0.14        statmod_1.5.1           ggtree_4.0.4           
#  [73] bit64_4.6.0-1           fontquiver_0.2.1        progress_1.2.3          filelock_1.0.3         
#  [77] bslib_0.11.0            affyio_1.78.0           otel_0.2.0              DBI_1.3.0              
#  [81] tidyselect_1.2.1        bit_4.6.0               compiler_4.5.1          curl_7.0.0             
#  [85] httr2_1.2.2             graph_1.88.1            xml2_1.5.2              fontBitstreamVera_0.1.1
#  [89] DelayedArray_0.36.0     scales_1.4.0            rappdirs_0.3.4          digest_0.6.39          
#  [93] rmarkdown_2.31          htmltools_0.5.9         pkgconfig_2.0.3         dbplyr_2.5.2           
#  [97] fastmap_1.2.0           rlang_1.1.7             htmlwidgets_1.6.4       UCSC.utils_1.4.0       
# [101] jquerylib_0.1.4         farver_2.1.2            jsonlite_2.0.0          GOSemSim_2.36.0        
# [105] R.oo_1.27.1             magrittr_2.0.4          Formula_1.2-5           GenomeInfoDbData_1.2.14
# [109] ggplotify_0.1.3         patchwork_1.3.2         Rcpp_1.1.1              ape_5.8-1              
# [113] ggnewscale_0.5.2        gdtools_0.5.0           stringi_1.8.7           ggraph_2.2.2           
# [117] MASS_7.3-65             plyr_1.8.9              parallel_4.5.1          listenv_0.10.1         
# [121] graphlayouts_1.2.3      splines_4.5.1           hms_1.1.4               locfit_1.5-9.12        
# [125] igraph_2.2.3            ggsignif_0.6.4          reshape2_1.4.5          XML_3.99-0.23          
# [129] evaluate_1.0.5          BiocManager_1.30.27     foreach_1.5.2           tzdb_0.5.0             
# [133] tweenr_2.0.3            polyclip_1.10-7         ggforce_0.5.0           xtable_1.8-8           
# [137] ff_4.5.2                reactome.db_1.92.0      tidytree_0.4.7          tidydr_0.0.6           
# [141] rstatix_0.7.3           aplot_0.2.9             memoise_2.0.1           cluster_2.1.8.2        
# [145] timechange_0.4.0        globals_0.19.1         
