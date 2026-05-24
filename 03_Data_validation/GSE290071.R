setwd("/home/xxm_xxm/CJX_workspace/Data validation/GSE290071/")  # Replace with your folder path
rm(list = ls())  
options(stringsAsFactors = FALSE)
cell=read_excel("GSE290071_Cell.xlsx")
colnames(cell)
dat1=cell[,2:13]
ids=dat1
ids$median=apply(ids,1,median) 
ids$gene_symbol=cell$GeneName
# Add a median column to ids and calculate the row median for each gene.
ids=ids[order(ids$gene_symbol,ids$median,decreasing = T),]
# Sort ids by gene symbol and then by median value in descending order.
## This keeps the first entry for duplicated gene symbols.
## Reassign the sorted rows to ids.
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
# Add a median column to ids and calculate the row median for each gene.
ids=ids[order(ids$gene_symbol,ids$median,decreasing = T),]
# Sort ids by gene symbol and then by median value in descending order.
## This keeps the first entry for duplicated gene symbols.
## Reassign the sorted rows to ids.
ids=ids[!duplicated(ids$gene_symbol),]
dat2=ids[,1:12]

dat2=as.matrix(dat2)
rownames(dat2)=ids$gene_symbol
write.csv(dat1,file="expressionmetrix_GSE1.csv")
write.csv(dat2,file="expressionmetrix_GSE2.csv")
##########ssGSEA########################
########################### Workflow ########################################
# Modify geo_id and group information as needed.
## ====== Integrated workflow: GSVA + ssGSEA scoring and plotting ======
rm(list = ls())  
options(stringsAsFactors = FALSE)

library(GSVA)
library(GSEABase)
library(ggplot2)

## ---------- Input data ----------
geo_id <- "GSE290071"  # <=== Manually specify the GEO ID here; it will be added to saved filenames.

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
# Save group order in advance to keep plotting consistent.
group_levels <- unique(group_info$Group)
# Gene sets: prepare them in advance.
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
## Keep the plotting group order consistent with group_info.

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
