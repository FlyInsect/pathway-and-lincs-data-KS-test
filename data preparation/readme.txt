1. Using the obtained DEGs enrich KEGGpathways.

2. Get the genes on each pathway from the KEGG file (including all pathways and their information)

3. Extract the logFC value of genes on pathways from nrDEG (gene differential expression result matrix)-----"pathway weight match.ipynb"

4. Use the "geneID.csv" file (including ENSEMBL ID, SYMBOL ID, ENTREZ ID of all genes) to match genes.

Converting gene SYMBOL ID to ENTREZ ID ------"pathway IDtrans SYMBLE_to_ENTREZID.ipynb"
(The gene id in the lincs data is ENTREZ ID)
