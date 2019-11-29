######################################################################################### 
######## collect all sizes` score for all drugs ############
#########################################################################################

################ 先合并（所有结果合并到一个文件中�? ##############
#### 无显著�?
# rm(list = ls())
# ra就是之前的top_size
 ra =c(25)
for(topsize in ra){
  # topsize = 990
  pathname = paste("../pathway6_nor/mid/D_",topsize, sep = "")
  setwd(pathname)
  filename <- list.files(pattern = '*_head_score.csv')
  merge.data = read.csv(file = filename[1], header = FALSE, row.names = 1)
  merge.data = merge.data[-c(1),]
  for (i in 2:length(filename)){
    new.data = read.csv(file = filename[i], header = FALSE, row.names = 1)
    new.data = new.data[-c(1),]
    merge.data = cbind(merge.data, new.data)
  }
  score_fdr = merge.data
  score_fdr = as.matrix(score_fdr)
  ### 归一�?
  head = score_fdr[c(1:7),]
  score = as.numeric(as.character(score_fdr[8,]))
  pos_p = which(score > 0)
  pos_n = which(score < 0)
  max = max(score)
  min = min(score)
  score[pos_p] = score[pos_p] / max
  score[pos_n] = -(score[pos_n] / min)
  head = rbind(head, score)
  # rownames
  name = as.character(row.names(score_fdr))
  row.names(head) = name
  score_fdr = head
  # output
  f = paste("all_score_", topsize, ".csv", sep = "")
  write.csv(score_fdr, f)
  # 下面是另一个目录（都在mid下）
  f = paste("../../../pathway6_nor/mid/all_score_", topsize, ".csv", sep = "")
  write.csv(score_fdr,f)
}

 