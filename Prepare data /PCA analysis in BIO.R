#!/usr/bin/env Rscript
# ==========================================
# 简化版 PCA：筛选高贡献度、低冗余的 bio 变量
# ==========================================


setwd("F:/caas/毕业课题/pig data/CHELSA/")
library(tidyverse)
library(factoextra)

# ---- 配置 ----
in_file <- "extracted_CHELSA_1981-2010_location_pig.csv"
out_dir <- "pca_output"
dir.create(out_dir, showWarnings = FALSE)

TOP_N <- 4           # 保留几个变量
COR_THRESHOLD <- 0.8 # 相关性阈值，超过此值视为高度相关

# ---- 读取数据 ----
df <- read_csv(in_file, show_col_types = FALSE)

# 重命名：bio_1 -> bio1, bio_2 -> bio2, ...
bio_data <- df %>% 
    select(starts_with("bio_")) %>% 
    rename_with(~ str_replace(., "bio_", "bio")) %>%
    na.omit()

# ---- PCA 分析 ----
pca <- prcomp(bio_data, center = TRUE, scale. = TRUE)

# 方差解释
var_exp <- summary(pca)$importance %>% 
    t() %>% 
    as_tibble(rownames = "PC") %>%
    write_csv(file.path(out_dir, "variance_explained.csv"))

# 提取PC1和PC2的方差解释率
pc1_var <- round(var_exp$`Proportion of Variance`[1] * 100, 2)
pc2_var <- round(var_exp$`Proportion of Variance`[2] * 100, 2)

# 变量贡献度（前两个主成分）
contrib <- get_pca_var(pca)$contrib[, 1:2]
contrib_sum <- rowSums(contrib)

# ---- 变量筛选策略 ----
# 1. 按贡献度排序
bio_ranked <- tibble(
    variable = names(contrib_sum),
    contrib_PC1_PC2 = contrib_sum
) %>% arrange(desc(contrib_PC1_PC2))

# 2. 剔除高相关变量
bio_cor <- cor(bio_data)
selected <- c()

for (var in bio_ranked$variable) {
    if (length(selected) == 0) {
        selected <- c(selected, var)
    } else {
        max_cor <- max(abs(bio_cor[var, selected]))
        if (max_cor < COR_THRESHOLD) {
            selected <- c(selected, var)
        }
    }
    if (length(selected) >= TOP_N) break
}

# 输出筛选结果
result <- bio_ranked %>%
    mutate(selected = variable %in% selected) %>%
    write_csv(file.path(out_dir, "variable_selection.csv"))

cat("\n=== 筛选结果 ===\n")
cat("保留的变量：", paste(selected, collapse = ", "), "\n")

# ---- 可视化 ----
# 1. 碎石图
png(file.path(out_dir, "scree_plot.png"), width = 2400, height = 1600, res = 300)
fviz_eig(pca, addlabels = TRUE)
dev.off()

# 2. 变量贡献图
png(file.path(out_dir, "variable_contribution.png"), width = 2400, height = 1600, res = 300)
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10)
dev.off()

# 3. 相关性热图
png(file.path(out_dir, "correlation_heatmap.png"), width = 2400, height = 2400, res = 300)
bio_cor %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor") %>%
    ggplot(aes(var1, var2, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Bio变量相关性矩阵")
dev.off()

# 4. PCA 双标图（仅显示筛选的变量）
png(file.path(out_dir, "biplot_selected.png"), width = 2400, height = 1800, res = 300)
fviz_pca_var(pca, select.var = list(name = selected), 
             repel = TRUE, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()

# 5. 仿照论文图的PCA双标图（所有变量）
png(file.path(out_dir, "biplot_all_variables.png"), width = 3200, height = 3200, res = 300)

# 提取数据
ind_coords <- as.data.frame(pca$x[, 1:2])
var_coords <- as.data.frame(pca$rotation[, 1:2] * 10)  # 放大箭头
var_coords$variable <- rownames(var_coords)

# 标注选中的变量
var_coords$is_selected <- var_coords$variable %in% selected

# 绘图
ggplot() +
    # 样本点
    geom_point(data = ind_coords, 
               aes(x = PC1, y = PC2), 
               color = "#66C2A5", 
               alpha = 0.6, 
               size = 3) +
    # 变量箭头
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 linewidth = ifelse(var_coords$is_selected, 1, 0.5),
                 color = ifelse(var_coords$is_selected, "black", "grey60")) +
    # 变量标签
    geom_text(data = var_coords,
              aes(x = PC1, y = PC2, label = variable),
              vjust = -0.5, hjust = 0.5,
              size = ifelse(var_coords$is_selected, 4.5, 3.5),
              fontface = ifelse(var_coords$is_selected, "bold", "plain")) +
    # 坐标轴
    labs(x = paste0("PC1(", pc1_var, "%)"),
         y = paste0("PC2(", pc2_var, "%)")) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_fixed()

dev.off()

cat("\n完成！结果保存在:", normalizePath(out_dir), "\n")
cat("\n主图文件：biplot_all_variables.png\n")
