1、用得到的差异基因富集pathways
2、从KEGG文件（包含所有pathway及其信息）获取各条pathway上的基因

3、从nrDEG（基因差异表达结果矩阵）中提取pathways上基因的logFC值-----“pathway weight match.ipynb”

4、使用“geneID.cs”文件（包含所有基因的ENSEMBL ID、SYMBOL ID、ENTREZ ID）匹配
将gene SYMBOL ID转为ENTREZ ID ------“pathway IDtrans SYMBLE_to_ENTREZID.ipynb”
（lincs数据中基因id为ENTREZ ID）