#!/usr/bin/env Rscript
# ==========================================
# 改进版 PCA：筛选高贡献度、低冗余的 bio 变量
# 新增功能：
# 1. 动态确定PC数量（基于累计方差解释率）
# 2. VIF分析作为补充验证
# 3. 生态学维度平衡性检查
# 4. 敏感性分析
# ==========================================

setwd("F:/caas/毕业课题/pig data/CHELSA/")
library(tidyverse)
library(factoextra)
library(car)  # 用于VIF计算

# ---- 配置 ----
in_file <- "extracted_CHELSA_1981-2010_location_pig.csv"
out_dir <- "pca_output"
dir.create(out_dir, showWarnings = FALSE)

TOP_N <- 4              # 保留几个变量
COR_THRESHOLD <- 0.8    # 相关性阈值
CUM_VAR_THRESHOLD <- 0.80  # 累计方差解释率阈值
VIF_THRESHOLD <- 10     # VIF阈值（>10表示严重多重共线性）

# ---- 读取数据 ----
df <- read_csv(in_file, show_col_types = FALSE)

# 重命名：bio_1 -> bio1, bio_2 -> bio2, ...
bio_data <- df %>% 
    select(starts_with("bio_")) %>% 
    rename_with(~ str_replace(., "bio_", "bio")) %>%
    na.omit()

cat("\n=== 数据概览 ===\n")
cat("样本数：", nrow(bio_data), "\n")
cat("变量数：", ncol(bio_data), "\n")

# ---- PCA 分析 ----
pca <- prcomp(bio_data, center = TRUE, scale. = TRUE)

# 方差解释 - 直接从summary提取importance矩阵
pca_summary <- summary(pca)
importance_mat <- pca_summary$importance

# 提取各行数据
sdev <- importance_mat[1, ]           # Standard deviation
prop_var <- importance_mat[2, ]       # Proportion of Variance  
cum_var <- importance_mat[3, ]        # Cumulative Proportion

# 保存方差解释率表格
var_exp <- tibble(
    PC = colnames(importance_mat),
    Standard_deviation = sdev,
    Proportion_of_Variance = prop_var,
    Cumulative_Proportion = cum_var
)
var_exp %>% write_csv(file.path(out_dir, "variance_explained.csv"))

# 动态确定保留的PC数
n_pc <- which(cum_var >= CUM_VAR_THRESHOLD)[1]
if (is.na(n_pc)) n_pc <- min(5, ncol(bio_data))  # 至少取5个或所有

cat("\n=== PCA结果 ===\n")
cat("保留主成分数：", n_pc, "（累计方差解释率 ≥", CUM_VAR_THRESHOLD*100, "%）\n")
cat("实际累计方差：", round(cum_var[n_pc]*100, 2), "%\n")

# 提取前n个PC的方差解释率
pc1_var <- round(prop_var[1] * 100, 2)
pc2_var <- round(prop_var[2] * 100, 2)

# ---- 变量贡献度（加权方法）----
# 使用加权贡献度：loading^2 * 方差解释率
loadings <- pca$rotation[, 1:n_pc]
var_weights <- prop_var[1:n_pc]

# 计算加权平方载荷
weighted_contrib <- sweep(loadings^2, 2, var_weights, "*")
contrib_sum <- rowSums(weighted_contrib)

# 也保存相对贡献度供参考
contrib_relative <- get_pca_var(pca)$contrib[, 1:n_pc]

# ---- 变量筛选策略 ----
cat("\n=== 开始变量筛选 ===\n")

# 1. 按加权贡献度排序
bio_ranked <- tibble(
    variable = names(contrib_sum),
    contrib_weighted = contrib_sum,
    contrib_relative = rowSums(contrib_relative)
) %>% arrange(desc(contrib_weighted))

cat("\n前10个高贡献度变量：\n")
print(bio_ranked %>% head(10), n = 10)

# 2. 剔除高相关变量
bio_cor <- cor(bio_data)
selected <- c()

for (var in bio_ranked$variable) {
    if (length(selected) == 0) {
        selected <- c(selected, var)
        cat("\n选择:", var, "（首个变量，加权贡献度:", round(bio_ranked$contrib_weighted[bio_ranked$variable == var], 4), "）\n")
    } else {
        max_cor <- max(abs(bio_cor[var, selected]))
        if (max_cor < COR_THRESHOLD) {
            selected <- c(selected, var)
            cat("选择:", var, "（最大相关系数:", round(max_cor, 3), "）\n")
        } else {
            most_corr_var <- names(which.max(abs(bio_cor[var, selected])))
            cat("剔除:", var, "（与", most_corr_var, 
                "相关系数:", round(max_cor, 3), "）\n")
        }
    }
    if (length(selected) >= TOP_N) break
}

# ---- VIF 分析（补充验证）----
cat("\n=== VIF分析（多重共线性检验）===\n")
if (length(selected) >= 2) {
    # 构建线性模型计算VIF
    vif_data <- bio_data[, selected]
    
    # 对每个变量计算VIF
    vif_results <- sapply(selected, function(var) {
        formula_str <- paste(var, "~", paste(setdiff(selected, var), collapse = " + "))
        model <- lm(as.formula(formula_str), data = vif_data)
        1 / (1 - summary(model)$r.squared)
    })
    
    vif_df <- tibble(
        variable = names(vif_results),
        VIF = vif_results
    ) %>% arrange(desc(VIF))
    
    vif_df %>% write_csv(file.path(out_dir, "vif_analysis.csv"))
    
    cat("\nVIF结果：\n")
    print(vif_df, n = Inf)
    
    if (any(vif_results > VIF_THRESHOLD)) {
        cat("\n⚠️  警告：以下变量VIF >", VIF_THRESHOLD, "，存在严重多重共线性：\n")
        print(vif_df %>% filter(VIF > VIF_THRESHOLD))
    } else {
        cat("\n✓ 所有变量VIF <", VIF_THRESHOLD, "，多重共线性可接受\n")
    }
}

# ---- 生态学维度检查 ----
cat("\n=== 生态学维度平衡性检查 ===\n")

# 定义生态学分类
eco_groups <- list(
    temperature = c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11"),
    precipitation = c("bio12", "bio13", "bio14", "bio16", "bio17", "bio18", "bio19"),
    seasonality = c("bio2", "bio3", "bio4", "bio7", "bio15")
)

eco_coverage <- map_dfr(names(eco_groups), function(group) {
    vars_in_group <- eco_groups[[group]]
    selected_in_group <- intersect(selected, vars_in_group)
    tibble(
        dimension = group,
        total_vars = length(vars_in_group),
        selected_vars = length(selected_in_group),
        coverage = paste(selected_in_group, collapse = ", ")
    )
})

eco_coverage %>% write_csv(file.path(out_dir, "ecological_balance.csv"))
print(eco_coverage)

if (any(eco_coverage$selected_vars == 0)) {
    cat("\n⚠️  警告：以下维度没有变量被选中：\n")
    print(eco_coverage %>% filter(selected_vars == 0))
    cat("建议：考虑增加TOP_N或调整COR_THRESHOLD\n")
}

# ---- 输出筛选结果 ----
result <- bio_ranked %>%
    mutate(
        selected = variable %in% selected,
        rank = row_number()
    ) %>%
    left_join(
        tibble(variable = rownames(bio_cor)) %>%
            mutate(max_cor_with_others = apply(bio_cor, 1, function(x) max(abs(x[x != 1])))),
        by = "variable"
    ) %>%
    write_csv(file.path(out_dir, "variable_selection.csv"))

cat("\n=== 最终筛选结果 ===\n")
cat("保留的变量：", paste(selected, collapse = ", "), "\n")

# 输出选中变量的相关矩阵
selected_cor <- bio_cor[selected, selected]
write.csv(selected_cor, file.path(out_dir, "selected_variables_correlation.csv"))

cat("\n选中变量间的相关系数：\n")
print(round(selected_cor, 3))

# ---- 敏感性分析 ----
cat("\n=== 敏感性分析 ===\n")
cat("测试不同参数组合...\n")

sensitivity <- expand_grid(
    top_n = c(3, 4, 5, 6),
    cor_thresh = c(0.7, 0.75, 0.8, 0.85, 0.9)
) %>%
    rowwise() %>%
    mutate(
        selected_vars = {
            sel <- c()
            for (var in bio_ranked$variable) {
                if (length(sel) == 0) {
                    sel <- c(sel, var)
                } else {
                    max_cor <- max(abs(bio_cor[var, sel]))
                    if (max_cor < cor_thresh) {
                        sel <- c(sel, var)
                    }
                }
                if (length(sel) >= top_n) break
            }
            paste(sel, collapse = ",")
        }
    ) %>%
    ungroup()

sensitivity %>% write_csv(file.path(out_dir, "sensitivity_analysis.csv"))
cat("敏感性分析结果已保存\n")

# ---- 可视化 ----
cat("\n=== 生成可视化图表 ===\n")

# 1. 碎石图
png(file.path(out_dir, "scree_plot.png"), width = 2400, height = 1600, res = 300)
fviz_eig(pca, addlabels = TRUE, ncp = 10) +
    geom_hline(yintercept = CUM_VAR_THRESHOLD * 100, linetype = "dashed", color = "red") +
    labs(subtitle = paste0("红线: ", CUM_VAR_THRESHOLD*100, "% 累计方差阈值"))
dev.off()

# 2. 变量贡献图（基于加权贡献度）
png(file.path(out_dir, "variable_contribution.png"), width = 2400, height = 1600, res = 300)

contrib_plot_data <- tibble(
    variable = names(contrib_sum),
    contribution = contrib_sum * 100  # 转换为百分比
) %>% arrange(desc(contribution))

ggplot(contrib_plot_data, aes(x = reorder(variable, contribution), y = contribution)) +
    geom_col(aes(fill = contribution), show.legend = TRUE) +
    scale_fill_gradient(low = "#00AFBB", high = "#FC4E07") +
    coord_flip() +
    geom_hline(yintercept = mean(contrib_plot_data$contribution), 
               linetype = "dashed", color = "red") +
    labs(title = "变量对前n个主成分的加权贡献度",
         subtitle = paste0("基于PC1-PC", n_pc, " (累计方差: ", round(cum_var[n_pc]*100, 1), "%)"),
         x = "变量",
         y = "加权贡献度 (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
dev.off()

# 3. 相关性热图（全部变量）
png(file.path(out_dir, "correlation_heatmap_all.png"), width = 2400, height = 2400, res = 300)
bio_cor %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor") %>%
    ggplot(aes(var1, var2, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Bio变量相关性矩阵（全部变量）")
dev.off()

# 4. 相关性热图（仅选中变量）
png(file.path(out_dir, "correlation_heatmap_selected.png"), width = 1800, height = 1600, res = 300)
selected_cor %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor") %>%
    ggplot(aes(var1, var2, fill = cor)) +
    geom_tile() +
    geom_text(aes(label = round(cor, 2)), size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "选中变量相关性矩阵")
dev.off()

# 5. PCA 双标图（仅显示筛选的变量）
png(file.path(out_dir, "biplot_selected.png"), width = 2400, height = 1800, res = 300)
fviz_pca_var(pca, select.var = list(name = selected), 
             repel = TRUE, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
    labs(title = "PCA Biplot - Selected Variables Only")
dev.off()

# 6. 论文级PCA双标图（所有变量，突出选中的）
png(file.path(out_dir, "biplot_all_variables.png"), width = 3200, height = 3200, res = 300)

# 提取数据
ind_coords <- as.data.frame(pca$x[, 1:2])
var_coords <- as.data.frame(pca$rotation[, 1:2] * 10)  # 放大箭头
var_coords$variable <- rownames(var_coords)
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
                 linewidth = ifelse(var_coords$is_selected, 1.2, 0.5),
                 color = ifelse(var_coords$is_selected, "#FC4E07", "grey60")) +
    # 变量标签
    geom_text(data = var_coords,
              aes(x = PC1, y = PC2, label = variable),
              vjust = -0.5, hjust = 0.5,
              size = ifelse(var_coords$is_selected, 5, 3.5),
              fontface = ifelse(var_coords$is_selected, "bold", "plain"),
              color = ifelse(var_coords$is_selected, "#FC4E07", "black")) +
    # 坐标轴
    labs(x = paste0("PC1 (", pc1_var, "%)"),
         y = paste0("PC2 (", pc2_var, "%)"),
         title = "PCA Biplot - All Variables",
         subtitle = paste("Selected variables (n =", length(selected), ") highlighted in orange")) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_fixed()

dev.off()

# 7. VIF柱状图（如果计算了VIF）
if (exists("vif_df")) {
    png(file.path(out_dir, "vif_barplot.png"), width = 2000, height = 1600, res = 300)
    ggplot(vif_df, aes(x = reorder(variable, VIF), y = VIF)) +
        geom_col(aes(fill = VIF > VIF_THRESHOLD)) +
        geom_hline(yintercept = VIF_THRESHOLD, linetype = "dashed", color = "red", linewidth = 1) +
        scale_fill_manual(values = c("TRUE" = "#FC4E07", "FALSE" = "#00AFBB"), 
                          labels = c("TRUE" = paste0("VIF > ", VIF_THRESHOLD), 
                                     "FALSE" = paste0("VIF ≤ ", VIF_THRESHOLD))) +
        coord_flip() +
        labs(title = "方差膨胀因子 (VIF) - 多重共线性检验",
             subtitle = paste0("红色虚线：VIF = ", VIF_THRESHOLD, " 阈值"),
             x = "变量",
             y = "VIF",
             fill = "状态") +
        theme_minimal() +
        theme(legend.position = "bottom")
    dev.off()
}

# ---- 生成分析报告 ----
report <- paste0(
    "==============================================\n",
    "PCA气候因子筛选分析报告\n",
    "==============================================\n\n",
    "分析时间：", Sys.time(), "\n",
    "数据文件：", in_file, "\n",
    "样本数量：", nrow(bio_data), "\n",
    "原始变量数：", ncol(bio_data), "\n\n",
    
    "--- 参数设置 ---\n",
    "目标变量数：", TOP_N, "\n",
    "相关性阈值：", COR_THRESHOLD, "\n",
    "累计方差阈值：", CUM_VAR_THRESHOLD, "\n",
    "VIF阈值：", VIF_THRESHOLD, "\n\n",
    
    "--- PCA结果 ---\n",
    "保留主成分数：", n_pc, "\n",
    "累计方差解释率：", round(cum_var[n_pc]*100, 2), "%\n",
    "PC1方差解释率：", pc1_var, "%\n",
    "PC2方差解释率：", pc2_var, "%\n\n",
    
    "--- 筛选结果 ---\n",
    "最终保留变量：", paste(selected, collapse = ", "), "\n\n",
    
    "--- 生态学维度 ---\n",
    paste(capture.output(print(eco_coverage)), collapse = "\n"), "\n\n",
    
    "--- VIF分析 ---\n",
    if (exists("vif_df")) {
        paste(capture.output(print(vif_df)), collapse = "\n")
    } else {
        "未计算（需要至少2个变量）"
    }, "\n\n",
    
    "--- 输出文件 ---\n",
    "1. variance_explained.csv - 方差解释率\n",
    "2. variable_selection.csv - 变量筛选详情\n",
    "3. selected_variables_correlation.csv - 选中变量相关矩阵\n",
    "4. vif_analysis.csv - VIF分析结果\n",
    "5. ecological_balance.csv - 生态学维度平衡\n",
    "6. sensitivity_analysis.csv - 敏感性分析\n",
    "7. *.png - 各类可视化图表\n\n",
    
    "==============================================\n"
)

writeLines(report, file.path(out_dir, "analysis_report.txt"))

cat("\n", report)
cat("\n✓ 完成！所有结果保存在:", normalizePath(out_dir), "\n")
cat("\n主要文件：\n")
cat("  - biplot_all_variables.png（论文级双标图）\n")
cat("  - analysis_report.txt（分析报告）\n")
cat("  - variable_selection.csv（详细筛选结果）\n")
