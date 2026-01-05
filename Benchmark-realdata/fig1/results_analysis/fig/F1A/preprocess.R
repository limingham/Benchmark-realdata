
source("~/Benchmark-realdata-main/fig1/results_analysis/readER_psedobulk.R")
source("~/Benchmark-realdata-main/fig1/results_analysis/readER_realbulk.R")
source("~/Benchmark-realdata-main/fig1/results_analysis/readTNBC_psedobulk.R")
source("~/Benchmark-realdata-main/fig1/results_analysis/readTNBC_realbulk.R")


pcc_pseudo_TNBC10 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC10.rds")
pcc_pseudo_TNBC7 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC7.rds")
pcc_pseudo_TNBC8 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC8.rds")
pcc_pseudo_TNBC9 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC9.rds")
pcc_pseudo_TNBC1 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC1.rds")
pcc_pseudo_TNBC2 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC2.rds")
pcc_pseudo_TNBC3 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC3.rds")
pcc_pseudo_TNBC4 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC4.rds")
pcc_pseudo_TNBC5 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC5.rds")
pcc_pseudo_TNBC6 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC6.rds")



pcc_pseudo_TNBC10_df = as.data.frame(pcc_pseudo_TNBC10) 
pcc_pseudo_TNBC9_df = as.data.frame(pcc_pseudo_TNBC9) 
pcc_pseudo_TNBC8_df = as.data.frame(pcc_pseudo_TNBC8) 
pcc_pseudo_TNBC7_df = as.data.frame(pcc_pseudo_TNBC7) 
pcc_pseudo_TNBC6_df = as.data.frame(pcc_pseudo_TNBC6) 
pcc_pseudo_TNBC5_df = as.data.frame(pcc_pseudo_TNBC5) 
pcc_pseudo_TNBC4_df = as.data.frame(pcc_pseudo_TNBC4) 
pcc_pseudo_TNBC3_df = as.data.frame(pcc_pseudo_TNBC3) 
pcc_pseudo_TNBC2_df = as.data.frame(pcc_pseudo_TNBC2) 
pcc_pseudo_TNBC1_df = as.data.frame(pcc_pseudo_TNBC1) 
# 假设所有向量已存在于环境中
# 生成数据框名称向量（从TNBC1到TNBC10）
df_names <- paste0("pcc_pseudo_TNBC", 1:10, "_df")

# 检查名称是否存在于环境中（避免报错）
existing_dfs <- df_names[sapply(df_names, exists)]

# # 批量获取数据框并按元素相加
# pcc_pseudo_TNBC_df <- Reduce("+", mget(existing_dfs))



library(abind)

# 1️⃣ 获取数据框列表
df_list <- mget(existing_dfs)

# 2️⃣ 转为纯矩阵（去掉data.frame属性）
df_list <- lapply(df_list, as.matrix)

# 3️⃣ 堆叠为三维数组
arr <- abind(df_list, along = 3)

# 4️⃣ 按元素取中位数
pcc_pseudo_TNBC_df <- apply(arr, c(1, 2), median, na.rm = TRUE)

# 5️⃣ 转为data.frame（恢复原列名）
pcc_pseudo_TNBC_df <- as.data.frame(pcc_pseudo_TNBC_df)


pcc_real_TNBC <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_real_TNBC.rds")

pcc_real_TNBC_df = as.data.frame(pcc_real_TNBC) 


cor(apply(pcc_pseudo_TNBC_df, 2, mean),apply(pcc_real_TNBC_df, 2, mean),method = 'spearman'
)



pcc_pseudo_ER10 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER10.rds")
pcc_pseudo_ER7 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER7.rds")
pcc_pseudo_ER8 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER8.rds")
pcc_pseudo_ER9 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER9.rds")
pcc_pseudo_ER1 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER1.rds")
pcc_pseudo_ER2 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER2.rds")
pcc_pseudo_ER3 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER3.rds")
pcc_pseudo_ER4 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER4.rds")
pcc_pseudo_ER5 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER5.rds")
pcc_pseudo_ER6 <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_pseudo_ER6.rds")



pcc_pseudo_ER10_df = as.data.frame(pcc_pseudo_ER10) 
pcc_pseudo_ER9_df = as.data.frame(pcc_pseudo_ER9) 
pcc_pseudo_ER8_df = as.data.frame(pcc_pseudo_ER8) 
pcc_pseudo_ER7_df = as.data.frame(pcc_pseudo_ER7) 
pcc_pseudo_ER6_df = as.data.frame(pcc_pseudo_ER6) 
pcc_pseudo_ER5_df = as.data.frame(pcc_pseudo_ER5) 
pcc_pseudo_ER4_df = as.data.frame(pcc_pseudo_ER4) 
pcc_pseudo_ER3_df = as.data.frame(pcc_pseudo_ER3) 
pcc_pseudo_ER2_df = as.data.frame(pcc_pseudo_ER2) 
pcc_pseudo_ER1_df = as.data.frame(pcc_pseudo_ER1) 
# 假设所有向量已存在于环境中
df_names <- paste0("pcc_pseudo_ER", 1:10, "_df")

# 检查名称是否存在于环境中（避免报错）
existing_dfs <- df_names[sapply(df_names, exists)]

# 批量获取数据框并按元素相加
# pcc_pseudo_ER_df <- Reduce("+", mget(existing_dfs))
# pcc_pseudo_ER_df = pcc_pseudo_ER_df[row.names(pcc_pseudo_ER_df) != 'CID3941',]



# 1️⃣ 获取数据框列表
df_list <- mget(existing_dfs)

# 2️⃣ 转为纯矩阵（去掉data.frame属性）
df_list <- lapply(df_list, as.matrix)

# 3️⃣ 堆叠为三维数组
arr <- abind(df_list, along = 3)

# 4️⃣ 按元素取中位数
pcc_pseudo_ER_df <- apply(arr, c(1, 2), median, na.rm = TRUE)

# 5️⃣ 转为data.frame（恢复原列名）
pcc_pseudo_ER_df <- as.data.frame(pcc_pseudo_ER_df)




pcc_real_ER <- readRDS("~/Benchmark-realdata-main/fig1/results_analysis/ER_pseudobulk/pcc_real_ER.rds")

pcc_real_ER_df = as.data.frame(pcc_real_ER) 
# pcc_real_ER_df = pcc_real_ER_df[row.names(pcc_real_ER_df) != 'CID3941',]




pcc_pseudo_ER_df = pcc_pseudo_ER_df[,c(1,3,4,5,6)]
pcc_pseudo_TNBC_df = pcc_pseudo_TNBC_df[,c(1,3,4,5,6)]
pcc_real_ER_df = pcc_real_ER_df[,c(1,3,4,5,6)]
pcc_real_TNBC_df = pcc_real_TNBC_df[,c(1,3,4,5,6)]

# 
pcc_pseudo_BRCA_df = rbind(pcc_pseudo_ER_df,pcc_pseudo_TNBC_df)
pcc_real_BRCA_df = rbind(pcc_real_ER_df,pcc_real_TNBC_df)

cc_vec = c()
for (cc in 1:nrow(pcc_pseudo_BRCA_df)) {
  cc_vec =c(cc_vec,cor(as.numeric(pcc_pseudo_BRCA_df[cc,c(1,2,3,4,5)]),
                       as.numeric(pcc_real_BRCA_df[cc,c(1,2,3,4,5)]),method = 'spearman'))
}


cor1 = cor(t(pcc_pseudo_BRCA_df[,c(1,2,3,4,5)]),method = 'spearman')
cor2 = cor(t(pcc_real_BRCA_df[,c(1,2,3,4,5)]),method = 'spearman')
cor3 = cor(t(pcc_real_BRCA_df[,c(1,2,3,4,5)]),t(pcc_pseudo_BRCA_df[,c(1,2,3,4,5)]),method = 'spearman')

diag(cor1) =NA
diag(cor2) =NA
diag(cor3) =NA

# psedo模态内不同样本排名相关性
mean(cor1,na.rm = TRUE)
median(cor1,na.rm = TRUE)

# real模态内不同样本排名相关性
mean(cor2,na.rm = TRUE)
median(cor2,na.rm = TRUE)

# 跨模态不同样本排名相关性
mean(cor3,na.rm = TRUE)
median(cor3,na.rm = TRUE)

# 跨模态相同样本排名相关性
mean(cc_vec)

sort(apply(cor1, 2, mean,na.rm = TRUE))
sort(apply(cor2, 2, mean,na.rm = TRUE))
sort(apply(cor3, 2, mean,na.rm = TRUE))

sort(apply(pcc_pseudo_BRCA_df, 2, mean,na.rm = TRUE))


apply(pcc_real_BRCA_df, 2, mean)
apply(pcc_pseudo_ER_df, 2, mean)
apply(pcc_pseudo_TNBC_df, 2, mean)


df1 = as.data.frame(apply(-1*pcc_real_BRCA_df, 1, rank))
t(df1)[,sort(row.names(df1))]
rowSums(df1)
apply(t(df1), 2, mean)


df2 = as.data.frame(apply(-1*pcc_pseudo_BRCA_df, 1, rank))
t(df2)[,sort(row.names(df2))]
apply(t(df2), 2, mean)


apply(pcc_real_BRCA_df, 2, mean)

apply(pcc_pseudo_BRCA_df, 2, mean)

c1 = c()
c2 = c()
c3 = c()

c_prop = c()
for (ii in 1:nrow(cor1)) {
  for (jj in 1:ncol(cor1)) {
    if((ii != jj)){
      c1 = c(c1,cor1[ii,jj])
      c2 = c(c2,cor2[ii,jj])
      c3 = c(c3,cor3[ii,jj])
      
      
      c_prop = c(c_prop,cor1[ii,jj])
      
    }
  }
}


wilcox.test(c1,c3,paired = TRUE)
wilcox.test(c2,c3,paired = TRUE)
wilcox.test(c1,cc_vec)
