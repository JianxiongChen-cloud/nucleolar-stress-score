#############################################################################
######################## 创建最终基因集 #####################################
#############################################################################
setwd("/home/xxm_xxm/CJX_workspace/geneset")
# 假设前面已经得到：
# final_up_genes
# final_down_genes

# 1. 标准化处理
final_up_genes <- sort(unique(trimws(toupper(final_up_genes))))
final_down_genes <- sort(unique(trimws(toupper(final_down_genes))))

# 2. 再次检查是否仍有上下调冲突
final_conflict_genes <- intersect(final_up_genes, final_down_genes)

if (length(final_conflict_genes) > 0) {
  cat("发现上下调冲突基因，已自动移除：", length(final_conflict_genes), "\n")
  print(final_conflict_genes)
  
  final_up_genes <- setdiff(final_up_genes, final_conflict_genes)
  final_down_genes <- setdiff(final_down_genes, final_conflict_genes)
}

# 3. 创建基因集
geneSets_final <- list(
  Final_Up = final_up_genes,
  Final_Down = final_down_genes
)

# 4. 保存
save(geneSets_final, file = "Final_GeneSets.Rdata")

# 5. 可选：同时导出为gmt风格文本或csv
write.csv(
  data.frame(gene_symbol = final_up_genes, stringsAsFactors = FALSE),
  "Final_Up_Genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = final_down_genes, stringsAsFactors = FALSE),
  "Final_Down_Genes.csv",
  row.names = FALSE
)

# 6. 输出统计
cat("Final up genes:", length(final_up_genes), "\n")
cat("Final down genes:", length(final_down_genes), "\n")
cat("Gene set saved as: Final_GeneSets.Rdata\n")