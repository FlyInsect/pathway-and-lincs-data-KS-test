#####################################
## 对结果进行两种排�?
## ①对结果的|ks|从大到小排序【这部分包含了治疗和加剧疾病的药物�?
## ②只用ks<0的那部分数据（对其按照|ks|从大到小排序）【这部分是治疗的药物�?
# rm(list = ls())
# size就是之前的top_size
 size = c(25)

for(s in size){
  # s = 200
  f = paste("../pathway6_nor/mid/all_score_", s, ".csv", sep = "")
  data = read.csv(f, header = F)
  data = data[-c(1),-c(1)]
#  data = data[,-c(1)]
  data = as.matrix(data)
  len = length(data[8,])
  score = c(0)
  for(i in 1:len){
    score[i] = as.numeric(as.character(data[8,i]))
  }
  # 所有|ks|排序
  data1 = data[,order(-abs(score))]
  f = paste("../pathway6_nor/mid/all_score_sort/all_normal_abs_sort_", s,".csv", sep = "")
  write.csv(data1, f,row.names = F, col.names = F)
  # 只用ks<0的排�?
  data1 = data[,order(score)]
  data1 = data1[,which(data1[8,] < 0)]
  print(length(data1[3,]))
  f = paste("../pathway6_nor/mid/all_score_sort/all_normal_less0_sort_", s,".csv", sep = "")
  write.csv(data1, f,row.names = F, col.names = F)
  
}
########### 对排序结果提取药物最优排�? #############
## 得到的药物排序就是最终的排名，接下来可以验证看top 10,20,30里有多少与疾病真正有�?
# rm(list = ls())
# size就是top_size
 size = c(25)
 ## 需要分成两次（要修改文件名�?
for(s in size){
  # s = 800
  f = paste("../pathway6_nor/mid/all_score_sort/all_normal_less0_sort_", s,".csv", sep = "")
#   f = paste("E:\\work\\chol\\result\\nodelete\\mid\\all_score_sort\\all_normal_abs_sort_", s,".csv", sep = "")
  data = read.csv(f, header = F)
  data = data[-c(1),]
  data = as.matrix(data)
  ## 一共有多少�?
  drugs = as.character(data[4,])
  drugs = unique(drugs)
  print(length(drugs))
  ## 按顺序提出药物最优排�?
  result = data
  drugs = c("")
  for(i in 1:length(result[4,])){
    drugs[i] = as.character(result[4,i])
  }
  drugs = unique(drugs)
  drugs = tolower(drugs)
  print(length(drugs))
  of = paste("../pathway6_nor/mid/all_score_", s, "_less0_sort_drugs.csv", sep = "")
#  of = paste("E:\\work\\chol\\result\\nodelete\\mid\\sort_drugs\\all_score_", s, "_abs_sort_drugs.csv", sep = "")
  write.csv(drugs, of, col.names = F, row.names = F)
  
}

