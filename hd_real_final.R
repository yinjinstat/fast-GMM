library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(transport)
library(LassoSIR)
library(SummarizedExperiment)
library(edgeR)
library(ggplot2)


###Data_processing

query_rna_seq <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

tryCatch({
  GDCdownload(
    query = query_rna_seq,
    method = "api",
    files.per.chunk = 50,
    directory = "GDCdata/COAD/RNA_Seq"
  )
  
  rna_seq_data <- GDCprepare(
    query = query_rna_seq,
    directory = "GDCdata/COAD/RNA_Seq"
  )
  
  cat("RNA-Seq 数据下载并准备完成。\n")
  print(rna_seq_data)
  
}, error = function(e) {
  stop(paste("下载或准备RNA-Seq数据时出错:", e$message))
})

counts_matrix <- assay(rna_seq_data, "unstranded") 
tryCatch({
  atlas_data_list <- TCGAbiolinks::PanCancerAtlas_subtypes()
  coad_subtypes <- subset(atlas_data_list, cancer.type == "COAD")
  cat("COAD 分子亚型数据获取完成。\n")
  print("COAD 亚型信息预览:")
  head(coad_subtypes)
  print("亚型类别及其计数:")
  print(table(coad_subtypes$Subtype_Selected))
  
}, error = function(e) {
  stop(paste("获取PanCancer Atlas亚型数据时出错:", e$message))
})

counts_patient_ids <- substr(colnames(counts_matrix), 1, 12)
subtype_patient_ids <- coad_subtypes$pan.samplesID

common_samples <- intersect(counts_patient_ids, subtype_patient_ids)

unique_counts_indices <- match(common_samples, counts_patient_ids)
counts_matrix_matched <- counts_matrix[, unique_counts_indices]

coad_subtypes_matched <- coad_subtypes[coad_subtypes$pan.samplesID %in% common_samples, ]
coad_subtypes_matched <- coad_subtypes_matched[match(common_samples, coad_subtypes_matched$pan.samplesID), ]
colnames(counts_matrix_matched) <- common_samples

target_subtypes <- c("GI.CIN", "GI.GS", "GI.HM-SNV", "GI.MSI")
valid_indices <- which(coad_subtypes_matched$Subtype_Selected %in% target_subtypes)

coad_subtypes_filtered <- coad_subtypes_matched[valid_indices, ]
counts_matrix_filtered <- counts_matrix_matched[, valid_indices]

Y <- ifelse(coad_subtypes_filtered$Subtype_Selected == "GI.MSI", 1, 0)

keep_genes <- rowSums(counts_matrix_filtered > 9) >= floor(341*0.8)
counts_matrix_low_filtered <- counts_matrix_filtered[keep_genes, ]

NAME=gene_info$gene_name[keep_genes]
setdiff(key_genes,NAME)

colData <- data.frame(
  subtype = coad_subtypes_filtered$Subtype_Selected,
  label = as.factor(Y) # design需要因子类型
)
rownames(colData) <- colnames(counts_matrix_low_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix_low_filtered,
  colData = colData,
  design = ~ label 
)

vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

cat("\nSize factor 归一化与 VST 变换完成。\n")
print("VST 变换后的数据矩阵维度:")
print(dim(vst_matrix))

X <- t(vst_matrix)
X<-scale(X)

#############################sparse_k_means processing#############################################

library(sparcl)

set.seed(123)  # 保证可复现
perm.out <- KMeansSparseCluster.permute(
  X,
  K = 2,                 # 目标聚类数，请按需要修改
  nperms = 25,           # 置换次数，越大越稳但耗时更长，如 50~100
  wbounds = seq(1.5, 8, by = 0.5)  # 候选稀疏度范围，可按列数调大/调小
)

perm.out$bestw
set.seed(123)
sparse.out <- KMeansSparseCluster(
  X,
  K = 2,                    # 同上，与调参时一致
  wbounds = perm.out$bestw,
  nstart = 25,              # 多次随机初始化，提升稳定性
  silent = TRUE
)

km.perm <- KMeansSparseCluster.permute(X,K=2,wbounds=seq(3,7,len=15),nperms=25)
km.out <- KMeansSparseCluster(X,K=2,wbounds=km.perm$bestw)
1-cal_acc(km.out[[1]]$Cs,Y)

###########################DR-GMM processing##################################################

###We run the code separately 

#The selected column in the reported result



p<-ncol(as.matrix(X))
n<-nrow(as.matrix(X))
kernel_matrix<-cbind(third_kernel(X),sin_kernel(X))
##column estimation #######
sfInit(parallel = TRUE,cpus=15)
sfLibrary(glmnet)
#sparse_estimate<-sfApply(kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds,lambda=seq(0.025,0.2,0.025))
sparse_estimate<-sfApply(kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds,lambda=seq(0.1,0.4,0.1))
sfStop()
index<-c()
num=p
for(i in 1:dim(sparse_estimate)[2])
{
  if(all(sparse_estimate[,i]==0))
  {
    index=append(index,i)
  }
}
index<-c(1:num)[-index]
sparse_index<-forward_column_selection(sparse_estimate,D=10)
print(sparse_index)

M<-cov(X)
U<-U_lambda(X,index)
count=1
for(t in index)
{
  sparse_estimate[,t]=(admm(M,U[,t],n,p,0.4,0,1000,d=1)$B[[1]])[,1]
  count=count+1
}
index<-c()
for(i in 1:dim(sparse_estimate)[2])
{
  if(all(sparse_estimate[,i]==0))
  {
    index=append(index,i)
  }
}
index<-c(1:num)[-index]
if(length(index)==0)
{
  index=sample(1:num,1)
}
print(index)
sparse_column=sparse_estimate[,index]
if(length(index)==1)
{
  sparse_index<-1
}else
{
  sparse_index<-forward_column_selection(sparse_column,D=10)
}
print(sparse_index)

Matrix_estimation_init<-function(data,ytilde,nfolds=3,lambda=NULL)
{
  if(is.null(lambda))
  {
    fit<-cv.glmnet(x=data,y=ytilde,family="mgaussian",nfolds=nfolds,lambda.min.ratio=1e-3,nlambda=50,type.measure = "mse",intercept=FALSE)
  }
  else{
    fit<-cv.glmnet(x=data,y=ytilde,family="mgaussian",nfolds=nfolds,type.measure = "mse",intercept=FALSE,lambda=lambda)
  }
  return(coef(fit,s=fit$lambda.min))
}

Matrix_estimation<-function(data,ytilde,nfolds=5,lambda=NULL,weight)
{
  if(is.null(lambda))
  {
    fit<-cv.glmnet(x=data,y=ytilde,family="mgaussian",nfolds=nfolds,lambda.min.ratio=1e-3,nlambda=100,type.measure = "mse",intercept=FALSE,penalty.factor=weight)
  }
  else{
    fit<-cv.glmnet(x=data,y=ytilde,family="mgaussian",nfolds=nfolds,type.measure = "mse",intercept=FALSE,lambda=lambda,penalty.factor=weight)
  }
  return(coef(fit,s=fit$lambda.1se))
}

##The selected column in the reported result
idx<-c(141,998,1125,10413,5517,6387,10727,1332,9111,8129,1946,2533,539,2940,7348)


###We use glmnet rather than ADMM, since the cross validation for ADMM algorithm is slow when p is super large.
kernel<-kernel_matrix[,idx]
result=Matrix_estimation_init(X,kernel,lambda=seq(0.02,0.05,0.04/50))
temp=do.call(cbind,(lapply(result,deal_list)))
weight=apply(as.matrix(temp),MARGIN=1, choose_weight)
result=Matrix_estimation(X,kernel,weight=weight,lambda=NULL)
temp=do.call(cbind,(lapply(result,deal_list)))
which(temp[,1]!=0)

###Figure for IF-PCA
IFPCA_V<-IFPCA_result_V$V1
df <- data.frame(
  IFPCA_V = IFPCA_V,
  Y = factor(Y) # 确保Y是因子变量(Factor)
)

p<-ggplot(df, aes(x = IFPCA_V)) +
  geom_histogram_pattern(
    aes(
      fill = Y,
      pattern_angle = Y
      # 不写 y：默认 after_stat(count) -> frequency
    ),
    pattern = "stripe",
    pattern_fill = "black",
    fill = "white",
    colour = "black",
    position = "identity",
    alpha = 0.5,
    bins = 20,
    
    pattern_density = 0.05,
    pattern_spacing = 0.015,
    pattern_size = 0.2
  ) +
  scale_pattern_angle_manual(values = c(45, 135)) +
  scale_fill_manual(values = c("white", "white")) +
  base_theme +
  labs(
    title = "IF-PCA Clustering Result",
    x = "The reduced predictor from IF-PCA",
    y = "Frequency"
  ) +
  theme(legend.position = "none")


rng <- range(x)
w <- (rng[2] - rng[1]) / 20
# 让 c 成为某个 break：break = c + m*w
# 向左右扩到覆盖数据范围
lo <- c + floor((rng[1] - c) / w) * w
hi <- c + ceiling((rng[2] - c) / w) * w

breaks <- seq(lo, hi, by = w)

# 验证：c 是否在 breaks 上（考虑浮点误差）
any(abs(breaks - c) < 1e-10)
p<-ggplot(df, aes(x = IFPCA_V)) +
  geom_histogram_pattern(
    aes(fill = Y, pattern_angle = Y),
    breaks = breaks,
    pattern = "stripe",
    pattern_fill = "black",
    fill = "white",
    colour = "black",
    position = "identity",
    alpha = 0.5,
    pattern_density = 0.05,
    pattern_spacing = 0.015,
    pattern_size = 0.2
  ) +
  scale_pattern_angle_manual(values = c(45, 135)) +
  scale_fill_manual(values = c("white", "white")) +
  base_theme +
  labs(
    title = "IF-PCA Clustering Result",
    x = "The reduced predictor from IF-PCA",
    y = "Frequency"
  ) +
  theme(legend.position = "none")




###Read gene info
gene_info <- rowData(rna_seq_data)

# 查看有哪些注释信息
colnames(gene_info)

# 通常包含 gene_id, gene_name, gene_type 等
head(gene_info[, c("gene_id", "gene_name", "gene_type")])

# 创建包含基因符号的完整映射表
# 首先筛选出X矩阵中保留的基因
gene_info_filtered <- gene_info[colnames(X), ]

gene_mapping_full <- data.frame(
  column_index = 1:ncol(X),
  gene_id = colnames(X),
  gene_symbol = gene_info_filtered$gene_name,
  gene_type = gene_info_filtered$gene_type
)

head(gene_mapping_full, 10)
