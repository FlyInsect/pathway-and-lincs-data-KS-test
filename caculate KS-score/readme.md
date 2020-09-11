对LINCS数据中每一种处理(剔除DMSO)的基因差异表达值和pathway上基因表达值进行KS检验，计算score值。
1、计算LINCS数据每个处理（每列）的score
2、合并所有列
3、所有药物处理score归一化后 
※按 |score|  从大到小排序。包含了可能治疗和加剧疾病的药物。
※只取 score<0 数据，按 |score|  从大到小排序。只包含可能治疗的药物。

对一种药物（drug_name）多个处理数据，取score值最小列，按score值从小到大对药物排序。

