###### Documentation so everyone understand what is going on in here:

1. Gene Set Enrichment analysis 
--------------------------------------------------------------------------
The FullList of protein-coding genes was used for this analysis and they were ranked based n their LogFC expression.

1.Sox2+_GSEA: CSV file containing all pathways ( either significnat or not).
2.Sox2+_GSEA_down: CSV file containing all the downregulated pathways that came out as significant : FDR < 0.05.
3.Sox2+_GSEA_up: CSV file containing all the upregulated pathways that came out as significant : FDR < 0.05.
4. Full_genes_list_GSEA: CSV file containing all pathways ( either significnat or not)+ the genes ( from our dataset) that belong to each pathway. 


Other files with extension .gmt ( just ignored those).

--------------------------------------------------------------------------
Barplots:
It contains the signifincatly upregulated and downregulated pathways in Sox2+ compared to the Sox2-
