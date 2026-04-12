setwd("/home/xxm_xxm/CJX_workspace/Data validation/GSE290071/")  # 替换为你的文件夹路径
rm(list = ls())  
options(stringsAsFactors = FALSE)
cell=read_excel("GSE290071_Cell.xlsx")
colnames(cell)
dat1=cell[,2:13]
ids=dat1
ids$median=apply(ids,1,median) 
ids$gene_symbol=cell$GeneName
#ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$gene_symbol,ids$median,decreasing = T),]
#对ids$symbol按照ids$median中位数从大到小排列的顺序排序
##即先按symbol排序，相同的symbol再按照中位数从大到小排列，方便后续保留第一个值。
##将对应的行赋值为一个新的ids，这样order()就相当于sort()
ids=ids[!duplicated(ids$gene_symbol),]
dat1=ids[,1:12]

dat1=as.matrix(dat1)
rownames(dat1)=ids$gene_symbol


ev=read_excel("GSE290071_EV.xlsx")
colnames(ev)
dat2=ev[,2:13]
ids=dat2
ids$median=apply(ids,1,median) 
ids$gene_symbol=ev$GeneName
#ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$gene_symbol,ids$median,decreasing = T),]
#对ids$symbol按照ids$median中位数从大到小排列的顺序排序
##即先按symbol排序，相同的symbol再按照中位数从大到小排列，方便后续保留第一个值。
##将对应的行赋值为一个新的ids，这样order()就相当于sort()
ids=ids[!duplicated(ids$gene_symbol),]
dat2=ids[,1:12]

dat2=as.matrix(dat2)
rownames(dat2)=ids$gene_symbol
write.csv(dat1,file="expressionmetrix_GSE1.csv")
write.csv(dat2,file="expressionmetrix_GSE2.csv")
##########ssGSEA########################
###########################流程化############################################
#!!!!!!根据情况修改geo——id及分组信息即可
## ====== 整合流程：GSVA + ssGSEA 打分 & 绘图 ======
rm(list = ls())  
options(stringsAsFactors = FALSE)

library(GSVA)
library(GSEABase)
library(ggplot2)

## ---------- 输入数据 ----------
geo_id <- "GSE290071"  # <=== 在这里手动指定 GEO 号，文件保存时会加上它

expr <- read.csv(file="expressionmetrix_GSE2.csv",
                 row.names = 1)

range(expr)
expr<- log2(expr+1)
expr.mat <- as.matrix(expr)
range(expr.mat)

group_info <- data.frame(
  Sample = colnames(expr.mat),
  Group = rep(c("control", "oxa"), each = 6)
)
#rep(c("BMH7", "BMH9", "BMH15", "BMH21","BMH22","BMH23","Control") ,each = 4)
#rep(c("Act D_0h", "Act D_1h", "Act D_3h", "Act D_6h", "Act D_8h"), times = 4)
# 提前保存分组顺序，绘图时保证一致
group_levels <- unique(group_info$Group)
# 基因集：提前准备好
load(file="/home/xxm_xxm/CJX_workspace/geneset/genesets.Rdata")
geneSets_all  <- list(NuS_Up = NuS_up_all_unique,
                      NuS_Down = NuS_down_all_unique)
geneSets_core <- list(NuS_Up = NuS_up_core_unique,
                      NuS_Down = NuS_down_core_unique)


ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr.mat,
  geneSets = geneSets_all,
  minSize = 5,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)


ssgsea.res <- gsva(ssgsea_param)
NuS_score_ssGSEA <- ssgsea.res["NuS_Up", ] - ssgsea.res["NuS_Down", ]
print(NuS_score_ssGSEA)
scores <- data.frame(
  Sample = colnames(expr.mat),
  ssGSEA = NuS_score_ssGSEA
)
scores <- merge(scores, group_info, by="Sample")
## 保证绘图分组顺序和 group_info 一致

scores$Group <- factor(scores$Group, levels = unique(group_info$Group))
plot_data_ssgsea <- data.frame(Score = scores$ssGSEA, Group = scores$Group)
library(ggplot2)

ggplot(plot_data_ssgsea, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = paste0("NuS Score (ssGSEA) - ", geo_id, " - ", geo_id),
       x = "Group", y = "NuS Score") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
ggsave("GSE290071-ev-all.pdf",  width = 2, height = 4)
