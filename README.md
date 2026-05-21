# LIMMA-Microarray Pipeline for Sox2 Cell Comparison
Transcriptomic analysis of Sox2+ vs Sox2- cells obtained from primary Sox2-eGFP mouse tumors.

1. Data was obtained from : "Quiescent Sox2+ Cells Drive Hierarchical Growth and Relapse in Sonic Hedgehog Subgroup Medulloblastoma": https://www.sciencedirect.com/science/article/pii/S1535610814002207?via%3Dihub#fig3.
2. Repository contains ouput from LIMMA workflow .R script:

Files:
1. CTX_FullList: Full list of genes from microarray analysis ( either significant or not ) including protein-coding genes, microRNAs, etc.
2. FullList_protein_coding_genes: Genes( either significant or not) and exclusively "protein coding genes".
3. CTX_DGE: protein-coding Differential Expressed genes. Cutoff is FDR < 0.05 and ABS(logFC > 1).








