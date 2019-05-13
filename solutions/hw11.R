
setwd("/Users/andreidm/ETH/courses/smcb/solutions/")

normal.data = read.csv("../data/NGSPileup/Normal_pileup.txt", sep="\t")
tumor.data = read.csv("../data/NGSPileup/Tumor_pileup.txt", sep="\t")

# some checks
sum(normal.data$Position != tumor.data$Position)
sum(length(normal.data$Position) != length(unique(normal.data$Position)))

# plot the histograms of coverage for each
hist(normal.data$Coverage)
hist(tumor.data$Coverage)

# plot together
library(ggplot2)
normal.data$type = 'normal'
tumor.data$type = 'tumor'
ggplot(rbind(normal.data, tumor.data), aes(Coverage, fill = type)) + geom_density(alpha = 0.4)


library(stringr)

get.contingency.table = function(position){
  
  new.table = matrix(rep(0,4), 2)
  
  # total number of alleles
  new.table[1,2] = tumor.data$Coverage[tumor.data$Position == position]
  new.table[2,2] = normal.data$Coverage[normal.data$Position == position]
  
  # now calculate number of variances
  tumor.variants = 0
  normal.variants = 0
  for (variance in c("A", "C", "G", "T", "a", "c", "g", "t")){
    
    # for tumor
    tumor.variants = tumor.variants + str_count(as.character(
      tumor.data$Read_bases[tumor.data$Position == position]), variance)
    # for control
    normal.variants = normal.variants + str_count(as.character(
      normal.data$Read_bases[normal.data$Position == position]), variance)
  }
  
  new.table[1,1] = tumor.variants
  new.table[2,1] = normal.variants
  
  # subtract variances from total number of alleles
  new.table[1,2] = new.table[1,2] - tumor.variants
  new.table[2,2] = new.table[2,2] - normal.variants
  
  return(new.table)
}

library(parallel)

# get all the contingency matrices
tables = mclapply(normal.data$Position, get.contingency.table, mc.cores = 8)

get.p.value = function(c.matrix){
  return(fisher.test(c.matrix, alternative = 't')$p.value)
}

# get all the p values
p.values = mclapply(tables, get.p.value, mc.cores = 8)

# report positions of significant SNVs
tumor.data$Position[p.values < 0.1]

# adjust p-values
adjusted.p.values = p.adjust(p.values, method = "fdr")

# report updated positions of significant SNVs
tumor.data$Position[adjusted.p.values < 0.1]



