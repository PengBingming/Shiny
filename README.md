# Shiny  
Shiny-app,  
RNAseq中为 RNAseq 测序的R语言下游分析：  
1、gene_id_change 为转换基因ID；  
2、DEseq2 为使用 DEseq2 包处理，其表达矩阵为整数型数据；  
3、Limma 为使用 limma 包处理，其表达矩阵为 log 标准化的数据；  
4、gene_lis 为使用 DEseq2 或 limma 的处理结果，提取基因ID为 "ENTREZID" 的 gene_list 用于后续富集分析。
