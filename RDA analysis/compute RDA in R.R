library(LEA)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(qvalue)
library(rnaturalearth)
library(rnaturalearthdata)
library(vegan)
library(sf)
library(dplyr)
library(readxl)
library(data.table)

# 安装缺失的包（如果需要）
required_packages <- c("robust", "qvalue")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "qvalue") {
      # qvalue是Bioconductor包
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("qvalue")
    } else {
      install.packages(pkg)
    }
  }
}

# 加载必要的包
library(robust)
library(qvalue)

# ============================================================
# 参数区
# ============================================================

vcf_path <- "/home/lingjiang/supeng/horse/env/chelas/horse_filter_pruned_956_rename.vcf"
climate_file <- "/home/lingjiang/supeng/horse/env/chelas/horse_ind_bio.xlsx"

output_prefix <- "./horse_rda"
geno_file     <- paste0(output_prefix, ".geno")
vcfsnp_file   <- paste0(output_prefix, ".vcfsnp")
sample_file   <- "./sample.txt"
snp_geno_rdata <- "./horse_rda_snp_geno.RData"

# 是否强制重新转换（设为TRUE将忽略已有文件重新处理）
force_convert <- FALSE

# ============================================================
# 工具函数：稳健提取 SNP 名称
# ============================================================

get_snp_names <- function(vcf_path, vcfsnp_file, n_snp_expected) {
  cat("正在提取 SNP 名称...\n")
  cat("期望 SNP 数量（geno 列数）:", n_snp_expected, "\n")

  if (file.exists(vcfsnp_file)) {
    cat("尝试方法1：read.table(sep=\"\") 从 vcfsnp 拆列读取...\n")
    tmp1 <- tryCatch(
      read.table(vcfsnp_file, header = FALSE, sep = "", fill = TRUE,
                 stringsAsFactors = FALSE, comment.char = "", quote = ""),
      error = function(e) NULL
    )

    if (!is.null(tmp1) && ncol(tmp1) >= 3) {
      snp_names <- as.character(tmp1[[3]])
      bad <- is.na(snp_names) | snp_names == "." | snp_names == ""
      if (any(bad) && ncol(tmp1) >= 2) {
        snp_names[bad] <- paste0(tmp1[[1]][bad], ":", tmp1[[2]][bad])
      }
      snp_names <- make.unique(snp_names)

      if (length(snp_names) == n_snp_expected) {
        cat("方法1成功：从 vcfsnp 第3列提取 SNP ID。\n")
        return(snp_names)
      } else {
        cat("方法1读取到 SNP 数:", length(snp_names),
            "但与 geno 列数不一致，继续尝试其他方法。\n")
      }
    } else {
      cat("方法1失败：vcfsnp 未能拆出 >=3 列。\n")
    }
  } else {
    cat("vcfsnp 文件不存在：", vcfsnp_file, "\n")
  }

  # 兜底：从 VCF 导出生成 ID
  cat("尝试方法2（兜底）：bcftools query 从 VCF 导出 SNP 信息生成 ID...\n")
  snpinfo_tsv <- "./tmp_snpinfo_from_vcf.tsv"
  cmd <- paste0(
    "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' ",
    shQuote(vcf_path),
    " > ", shQuote(snpinfo_tsv)
  )
  system(cmd)

  if (!file.exists(snpinfo_tsv)) {
    stop("兜底失败：bcftools query 未生成 ", snpinfo_tsv, "。请确认 bcftools 可用。")
  }

  snpinfo <- fread(snpinfo_tsv, header = FALSE)
  snp_names <- as.character(snpinfo[[3]])
  bad <- is.na(snp_names) | snp_names == "." | snp_names == ""
  if (any(bad)) {
    snp_names[bad] <- paste0(snpinfo[[1]][bad], ":", snpinfo[[2]][bad],
                            "_", snpinfo[[4]][bad], "_", snpinfo[[5]][bad])
  }
  snp_names <- make.unique(snp_names)

  if (length(snp_names) != n_snp_expected) {
    stop("兜底提取的 SNP 数(", length(snp_names),
         ") 与 geno 列数(", n_snp_expected, ") 仍不一致。")
  }

  cat("方法2成功：从 VCF 导出生成 SNP ID。\n")
  return(snp_names)
}

# ============================================================
# 第一步：检查并转换VCF为LFMM geno格式
# ============================================================

cat("=====================================\n")
cat("步骤1: 检查并转换VCF为LFMM格式\n")
cat("=====================================\n\n")

# 检查是否需要转换
if (file.exists(snp_geno_rdata) && !force_convert) {
  cat("检测到已有RData文件:", snp_geno_rdata, "\n")
  cat("加载已保存的snp_geno对象...\n")
  load(snp_geno_rdata)
  cat("已加载基因型矩阵\n")
  cat("基因型矩阵维度:", dim(snp_geno), "\n")
  cat("  行数（样本数）:", nrow(snp_geno), "\n")
  cat("  列数（SNP数）:", ncol(snp_geno), "\n\n")

  # 检查是否已有snp_names和samples
  if (!exists("snp_names") || !exists("samples")) {
    cat("警告: RData中缺少snp_names或samples，需要重新提取...\n")
    need_extract <- TRUE
  } else {
    need_extract <- FALSE
    cat("RData中已包含snp_names和samples\n")
    cat("SNP数量:", length(snp_names), "\n")
    cat("样本数量:", length(samples), "\n")
  }
} else if (file.exists(geno_file) && file.exists(vcfsnp_file) && !force_convert) {
  cat("检测到已有geno和vcfsnp文件，跳过转换步骤...\n")
  cat("  - ", geno_file, "\n", sep = "")
  cat("  - ", vcfsnp_file, "\n\n", sep = "")

  cat("正在读取geno文件...\n")
  snp_geno <- read.geno(geno_file)
  cat("基因型矩阵维度:", dim(snp_geno), "\n")
  cat("  行数（样本数）:", nrow(snp_geno), "\n")
  cat("  列数（SNP数）:", ncol(snp_geno), "\n\n")
  need_extract <- TRUE
} else {
  # 需要转换VCF
  if (!file.exists(vcf_path)) stop("VCF文件不存在: ", vcf_path)

  cat("VCF文件路径:", vcf_path, "\n")
  cat("\n正在转换VCF为LFMM geno格式...\n")
  vcf2geno(vcf_path, output.file = geno_file)

  cat("转换完成！生成文件:\n")
  cat("  - ", geno_file, " (基因型矩阵)\n", sep = "")
  cat("  - ", vcfsnp_file, " (SNP位置信息)\n\n", sep = "")

  cat("正在读取geno文件...\n")
  snp_geno <- read.geno(geno_file)

  cat("基因型矩阵维度:", dim(snp_geno), "\n")
  cat("  行数（样本数）:", nrow(snp_geno), "\n")
  cat("  列数（SNP数）:", ncol(snp_geno), "\n\n")
  need_extract <- TRUE
}

# ============================================================
# 第二步：检查并提取SNP ID和样本名
# ============================================================

cat("=====================================\n")
cat("步骤2: 检查并提取SNP ID和样本名\n")
cat("=====================================\n\n")

if (exists("need_extract") && need_extract) {
  # 2.1 SNP 名
  snp_names <- get_snp_names(vcf_path, vcfsnp_file, ncol(snp_geno))
  cat("\nSNP数量:", length(snp_names), "\n")
  cat("前10个SNP ID:\n")
  print(head(snp_names, 10))

  colnames(snp_geno) <- snp_names
  cat("\n已为SNP矩阵添加列名\n")

  # 2.2 样本名
  cat("\n正在提取样本名称...\n")
  if (!file.exists(sample_file) || force_convert) {
    system(paste0("bcftools query -l ", shQuote(vcf_path), " > ", shQuote(sample_file)))
  }

  if (!file.exists(sample_file)) {
    stop("sample.txt 未生成，请检查 bcftools 是否可用。")
  }

  samples <- readLines(sample_file, warn = FALSE)
  samples <- samples[nzchar(samples)]

  cat("样本数量:", length(samples), "\n")
  cat("前10个样本名:\n")
  print(head(samples, 10))

  if (length(samples) != nrow(snp_geno)) {
    stop("样本名数量(", length(samples), ") 与 geno 行数(", nrow(snp_geno), ") 不一致！")
  }

  rownames(snp_geno) <- samples
  cat("\n已为SNP矩阵添加行名\n")

  # 缺失值比例
  missing_rate <- sum(snp_geno == 9) / (nrow(snp_geno) * ncol(snp_geno))
  cat("\n缺失值比例:", round(missing_rate * 100, 4), "%\n")

  # 保存到RData文件
  save(snp_geno, snp_names, samples, file = snp_geno_rdata)
  cat("已保存: horse_rda_snp_geno.RData\n\n")
} else {
  cat("跳过SNP ID和样本名提取步骤，使用已加载的数据\n")

  # 确保snp_names和samples与snp_geno维度匹配
  if (length(snp_names) != ncol(snp_geno)) {
    warning("snp_names长度(", length(snp_names), ")与SNP矩阵列数(", ncol(snp_geno), ")不匹配，重新提取...")
    snp_names <- get_snp_names(vcf_path, vcfsnp_file, ncol(snp_geno))
    colnames(snp_geno) <- snp_names
  }

  if (length(samples) != nrow(snp_geno)) {
    warning("samples长度(", length(samples), ")与SNP矩阵行数(", nrow(snp_geno), ")不匹配，重新提取...")
    if (!file.exists(sample_file)) {
      system(paste0("bcftools query -l ", shQuote(vcf_path), " > ", shQuote(sample_file)))
    }
    samples <- readLines(sample_file, warn = FALSE)
    samples <- samples[nzchar(samples)]
    rownames(snp_geno) <- samples
  }

  # 缺失值比例
  missing_rate <- sum(snp_geno == 9) / (nrow(snp_geno) * ncol(snp_geno))
  cat("\n缺失值比例:", round(missing_rate * 100, 4), "%\n\n")
}

# ============================================================
# 第三步：检查并读取气候数据
# ============================================================

cat("=====================================\n")
cat("步骤3: 检查并读取气候数据\n")
cat("=====================================\n\n")

if (!file.exists(climate_file)) stop("Excel文件不存在: ", climate_file)
climate_data <- read_excel(climate_file, sheet = "Sheet1")

# 使用您的6个环境变量
required_cols <- c("Sample", "经度", "纬度", "bio_1", "bio_3", "bio_4", "bio_10", "bio_17", "bio_18")
missing_cols <- setdiff(required_cols, colnames(climate_data))
if (length(missing_cols) > 0) stop("缺少以下必需列: ", paste(missing_cols, collapse = ", "))

cat("气候数据样本数量:", nrow(climate_data), "\n")
cat("必需列检查通过！\n\n")

# ============================================================
# 第四步：匹配样本顺序
# ============================================================

cat("=====================================\n")
cat("步骤4: 匹配样本顺序\n")
cat("=====================================\n\n")

geno_samples <- rownames(snp_geno)
climate_samples <- climate_data$Sample

common_samples <- intersect(geno_samples, climate_samples)
cat("基因型矩阵样本数:", length(geno_samples), "\n")
cat("气候数据样本数:", length(climate_samples), "\n")
cat("共同样本数:", length(common_samples), "\n")

missing_in_climate <- setdiff(geno_samples, climate_samples)
missing_in_geno <- setdiff(climate_samples, geno_samples)

if (length(missing_in_climate) > 0) {
  cat("\n警告：基因型有但气候缺失的样本(显示前10):\n")
  print(head(missing_in_climate, 10))
}
if (length(missing_in_geno) > 0) {
  cat("\n警告：气候有但基因型缺失的样本(显示前10):\n")
  print(head(missing_in_geno, 10))
}

# 用 geno 的顺序作为最终顺序
keep <- geno_samples[geno_samples %in% climate_samples]
snp_geno_filtered <- snp_geno[keep, , drop = FALSE]
climate_data_filtered <- climate_data[match(keep, climate_data$Sample), ]

stopifnot(all(rownames(snp_geno_filtered) == climate_data_filtered$Sample))

cat("\n样本匹配成功！\n")
cat("最终分析样本数:", nrow(snp_geno_filtered), "\n")
cat("最终分析SNP数:", ncol(snp_geno_filtered), "\n\n")

# ============================================================
# 第五步：准备RDA输入数据
# ============================================================

cat("=====================================\n")
cat("步骤5: 准备RDA输入数据\n")
cat("=====================================\n\n")

# 使用您的6个环境变量
ind_clim_data <- data.frame(
  longitude = climate_data_filtered$经度,
  latitude  = climate_data_filtered$纬度,
  bio_1  = climate_data_filtered$bio_1,
  bio_3  = climate_data_filtered$bio_3,
  bio_4  = climate_data_filtered$bio_4,
  bio_10 = climate_data_filtered$bio_10,
  bio_17 = climate_data_filtered$bio_17,
  bio_18 = climate_data_filtered$bio_18
)

cat("检查缺失值:\n")
missing_counts <- colSums(is.na(ind_clim_data))
print(missing_counts)

if (any(is.na(ind_clim_data))) {
  cat("\n警告：气候数据中存在缺失值，将移除包含缺失值的样本\n")
  complete_cases <- complete.cases(ind_clim_data)
  cat("包含缺失值的样本数:", sum(!complete_cases), "\n")
  ind_clim_data <- ind_clim_data[complete_cases, ]
  snp_geno_filtered <- snp_geno_filtered[complete_cases, ]
  cat("移除后剩余样本数:", nrow(ind_clim_data), "\n")
}

cat("\n气候变量统计摘要:\n")
print(summary(ind_clim_data[, 3:8]))

# ============================================================
# 第六步：检查RDA分析结果文件
# ============================================================

cat("\n=====================================\n")
cat("步骤6: 检查RDA分析结果文件\n")
cat("=====================================\n\n")

rda_results_file <- "./RDA_analysis_results.RData"
if (file.exists(rda_results_file) && !force_convert) {
  cat("检测到已有RDA分析结果文件:", rda_results_file, "\n")
  cat("是否跳过RDA分析，直接加载已有结果？\n")
  cat("1: 跳过分析，直接加载结果\n")
  cat("2: 重新运行RDA分析\n")

  # 在交互模式下询问用户，非交互模式下默认重新运行
  if (interactive()) {
    choice <- readline(prompt = "请输入选择(1或2): ")
  } else {
    choice <- "2"
    cat("非交互模式，默认重新运行RDA分析\n")
  }

  if (choice == "1") {
    cat("正在加载已有RDA分析结果...\n")
    load(rda_results_file)

    # 检查是否包含所有必要的对象
    required_objects <- c("RDA", "rdadapt_res", "outliers", "ind_clim_data", "snp_geno_filtered", "snp_names")
    missing_objects <- setdiff(required_objects, ls())

    if (length(missing_objects) == 0) {
      cat("成功加载RDA分析结果\n")
      cat("候选基因座数量:", nrow(outliers), "\n")
      skip_rda <- TRUE
    } else {
      cat("警告：RData文件缺少以下对象:", paste(missing_objects, collapse = ", "), "\n")
      cat("将重新运行RDA分析\n")
      skip_rda <- FALSE
    }
  } else {
    skip_rda <- FALSE
  }
} else {
  skip_rda <- FALSE
}

# ============================================================
# 第七步：定义RDA辅助函数（与参考脚本完全一致）
# ============================================================

if (!skip_rda) {
  cat("\n=====================================\n")
  cat("步骤7: 定义RDA分析函数\n")
  cat("=====================================\n\n")

  # 修正的rdadapt函数 - 不再在函数内部加载包
  rdadapt <- function(rda, K){
    # 检查必要的包是否已加载
    if (!("robust" %in% .packages())) {
      stop("请先加载 'robust' 包: library(robust)")
    }
    if (!("qvalue" %in% .packages())) {
      stop("请先加载 'qvalue' 包: library(qvalue)")
    }

    # 参考脚本使用前K个RDA轴的loadings
    zscores <- as.matrix(rda$CCA$v[, 1:as.numeric(K)])

    # 参考脚本对loadings进行标准化
    resscale <- apply(zscores, 2, scale)

    # 参考脚本使用稳健协方差计算Mahalanobis距离
    resmaha <- robust::covRob(resscale, distance = TRUE,
                              na.action = na.omit, estim = "pairwiseGK")$dist

    # 计算lambda（基因组膨胀因子）
    lambda <- median(resmaha) / qchisq(0.5, df = K)

    # 卡方检验
    reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)

    # 多重检验校正
    padj_BH <- p.adjust(reschi2test, method = "BH")           # FDR校正（Benjamini-Hochberg）
    qval <- qvalue::qvalue(reschi2test)$qvalues               # q值
    padj <- p.adjust(reschi2test, method = "bonferroni")      # Bonferroni校正

    return(data.frame(p.values = reschi2test,
                      padj_BH = padj_BH,
                      padj_q = qval,
                      padj = padj))
  }

  cat("RDA函数定义完成（与参考脚本一致）\n")

  # ============================================================
  # 第八步：运行RDA分析（使用您的6个环境变量）
  # ============================================================

  cat("\n=====================================\n")
  cat("步骤8: 运行RDA分析\n")
  cat("=====================================\n\n")

  cat("运行RDA模型...\n")
  cat("使用的气候变量: bio_1, bio_3, bio_4, bio_10, bio_17, bio_18\n\n")

  # 运行RDA
  RDA <- rda(snp_geno_filtered ~ bio_1 + bio_3 + bio_4 + bio_10 + bio_17 + bio_18,
             data = ind_clim_data)

  cat("RDA模型摘要:\n")
  print(summary(RDA))

  # 计算R²
  r2 <- RsquareAdj(RDA)
  cat("\n=== 模型拟合度 ===\n")
  cat("R² (unadjusted):", round(r2$r.squared, 4), "\n")
  cat("R² (adjusted):  ", round(r2$adj.r.squared, 4), "\n")

  # RDA轴方差解释
  cat("\n=== RDA轴方差解释 ===\n")
  eig <- RDA$CCA$eig
  variance <- eig / sum(eig)
  cat("RDA1解释方差:", round(variance[1] * 100, 2), "%\n")
  cat("RDA2解释方差:", round(variance[2] * 100, 2), "%\n")
  cat("RDA3解释方差:", round(variance[3] * 100, 2), "%\n")

  # ============================================================
  # 第九步：识别候选基因座（与参考脚本完全一致）
  # ============================================================

  cat("\n=====================================\n")
  cat("步骤9: 识别候选基因座 (outlier SNPs)\n")
  cat("=====================================\n\n")

  # 参考脚本使用K=2
  K_axes <- 2
  rdadapt_res <- rdadapt(RDA, K_axes)

  # 参考脚本使用FDR校正阈值0.01
  thres_env <- 0.01

  outliers <- data.frame(
    Loci = snp_names[which(rdadapt_res$padj_BH < thres_env)],
    p.value = rdadapt_res$p.values[which(rdadapt_res$padj_BH < thres_env)]
  )

  cat("\n=== 候选基因座统计 ===\n")
  cat("FDR阈值:", thres_env, "\n")
  cat("候选基因座数量:", nrow(outliers), "\n")
  cat("占总SNP比例:", round(nrow(outliers)/ncol(snp_geno_filtered)*100, 2), "%\n")

  if (nrow(outliers) > 0) {
    cat("\n前20个候选SNP:\n")
    print(head(outliers, 20))

    write.table(outliers, "./RDA_outlier_SNPs.txt",
                row.names = FALSE, quote = FALSE)
    cat("\n候选SNP列表已保存到: RDA_outlier_SNPs.txt\n")
  } else {
    cat("\n未检测到显著的候选SNP\n")
  }

  # ============================================================
  # 第十步：保存分析结果
  # ============================================================

  cat("\n=====================================\n")
  cat("步骤10: 保存分析结果\n")
  cat("=====================================\n\n")

  save(RDA, rdadapt_res, outliers, ind_clim_data, snp_geno_filtered, snp_names,
       file = rda_results_file)
  cat("已保存分析对象:", rda_results_file, "\n")

  rda_pvalues <- data.frame(
    Loci = snp_names,
    p.value = rdadapt_res$p.values,
    padj_BH = rdadapt_res$padj_BH,
    padj_q = rdadapt_res$padj_q,
    padj_bonferroni = rdadapt_res$padj
  )
  write.table(rda_pvalues, "./RDA_all_SNPs_pvalues.txt",
              row.names = FALSE, quote = FALSE)
  cat("所有SNP的p值已保存到: RDA_all_SNPs_pvalues.txt\n")
} else {
  cat("\n跳过RDA分析步骤，使用已加载的结果\n")

  if (nrow(outliers) > 0) {
    cat("\n候选基因座数量:", nrow(outliers), "\n")
    cat("前20个候选SNP:\n")
    print(head(outliers, 20))
  }
}

# ============================================================
# 第十一步：生成简单可视化
# ============================================================

cat("\n=====================================\n")
cat("步骤11: 生成可视化图表\n")
cat("=====================================\n\n")

rda_plot_file <- "RDA_ordination_plot.pdf"
manhattan_plot_file <- "RDA_manhattan_plot.pdf"

if (file.exists(rda_plot_file) && file.exists(manhattan_plot_file) && !force_convert) {
  cat("检测到已有可视化文件:\n")
  cat("  -", rda_plot_file, "\n")
  cat("  -", manhattan_plot_file, "\n")
  cat("是否重新生成可视化图表？\n")

  if (interactive()) {
    choice <- readline(prompt = "重新生成请输入y，跳过请输入n: ")
  } else {
    choice <- "n"
    cat("非交互模式，默认跳过可视化生成\n")
  }

  if (tolower(choice) != "y") {
    cat("跳过可视化生成步骤\n")
    skip_plots <- TRUE
  } else {
    skip_plots <- FALSE
  }
} else {
  skip_plots <- FALSE
}

if (!skip_plots) {
  # 1. RDA排序图（你原来的：画的是样本点sites）
  pdf("RDA_ordination_plot.pdf", width = 10, height = 8)
  plot(RDA, type = "n", main = "RDA - Samples and Environmental Variables")
  points(RDA, display = "sites", pch = 20, col = "blue", cex = 1.2)
  text(RDA, display = "bp", col = "red", cex = 1.2, lwd = 2)
  legend("topleft", legend = c("Samples", "Env Variables"),
         col = c("blue", "red"), pch = c(20, NA), lty = c(NA, 1), bty = "n")
  dev.off()
  cat("已保存: RDA_ordination_plot.pdf\n")

  # 2. 曼哈顿图
  pdf("RDA_manhattan_plot.pdf", width = 12, height = 6)
  par(mar = c(5, 5, 4, 2))
  plot(-log10(rdadapt_res$p.values),
       pch = 20, cex = 0.6, col = "grey40",
       xlab = "SNP Index", ylab = "-log10(P-value)",
       main = "RDA Outlier SNPs Detection")
  abline(h = -log10(thres_env), col = "red", lty = 2, lwd = 2)

  if (nrow(outliers) > 0) {
    outlier_indices <- which(rdadapt_res$padj_BH < thres_env)
    points(outlier_indices,
           -log10(rdadapt_res$p.values[outlier_indices]),
           pch = 20, col = "red", cex = 1)
  }

  legend("topright",
         legend = c(paste0("Non-outlier SNPs (n=",
                           ncol(snp_geno_filtered) - nrow(outliers), ")"),
                    paste0("Outlier SNPs (n=", nrow(outliers),
                           ", FDR<", thres_env, ")")),
         col = c("grey40", "red"), pch = 20, bty = "n")
  dev.off()
  cat("已保存: RDA_manhattan_plot.pdf\n")
}

# ============================================================
# 第十二步：生成 Figure S5 风格的 RDA space（SNP点云）
# ============================================================

cat("\n=====================================\n")
cat("步骤12: 生成 Figure S5 风格的 RDA space (All SNPs vs Outliers)\n")
cat("=====================================\n\n")

# --- 取 SNP loadings（species scores） ---
sp <- scores(RDA, display = "species", choices = 1:2, scaling = 2)
sp <- as.data.frame(sp)
sp$Loci <- rownames(sp)

# 对齐检查：确保 outliers$Loci 在 sp 里能匹配到
matched_n <- sum(outliers$Loci %in% sp$Loci)
cat("Outliers 与 RDA species 行名匹配数量:", matched_n, "/", nrow(outliers), "\n")
if (matched_n == 0 && nrow(outliers) > 0) {
  cat("警告：outliers$Loci 可能与 RDA 的species行名不一致。\n")
  cat("通常原因：RDA对象的species行名不是SNP ID。请确认你用于rda()的矩阵列名是 snp_names。\n")
}

sp$group <- ifelse(sp$Loci %in% outliers$Loci, "Outliers", "All SNPs")

# --- 取环境变量箭头（bp vectors） ---
bp <- scores(RDA, display = "bp", choices = 1:2, scaling = 2)
bp <- as.data.frame(bp)
bp$var <- rownames(bp)

# --- 核心：缩放箭头到点云尺度（避免箭头撑大坐标轴） ---
scale_arrows_to_cloud <- function(bp, sp, frac = 0.85) {
  sx <- diff(range(sp$RDA1))
  sy <- diff(range(sp$RDA2))
  bx <- max(abs(bp$RDA1))
  by <- max(abs(bp$RDA2))
  if (bx == 0 || by == 0) return(bp)
  k <- frac * min(sx/(2*bx), sy/(2*by))
  bp$RDA1 <- bp$RDA1 * k
  bp$RDA2 <- bp$RDA2 * k
  bp
}
bp2 <- scale_arrows_to_cloud(bp, sp, frac = 0.85)

# --- 坐标范围：只按 SNP 点云定（避免箭头影响范围） ---
pad <- 0.06
xlim <- range(sp$RDA1); xr <- diff(xlim); xlim <- xlim + c(-1, 1) * xr * pad
ylim <- range(sp$RDA2); yr <- diff(ylim); ylim <- ylim + c(-1, 1) * yr * pad

# --- 解释率 ---
eig <- RDA$CCA$eig
p1 <- round(eig[1] / sum(eig) * 100, 2)
p2 <- round(eig[2] / sum(eig) * 100, 2)

# --- 作图（png + pdf） ---
p <- ggplot(sp, aes(x = RDA1, y = RDA2)) +
  geom_point(data = subset(sp, group == "All SNPs"),
             color = "grey70", size = 0.6, alpha = 0.8) +
  geom_point(data = subset(sp, group == "Outliers"),
             color = "purple4", size = 0.9, alpha = 0.95) +

  geom_segment(data = bp2,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.22, "cm")),
               linewidth = 0.7, color = "grey20") +
  geom_text(data = bp2,
            aes(x = RDA1, y = RDA2, label = var),
            color = "grey10", size = 4, vjust = -0.6) +

  geom_hline(yintercept = 0, linetype = "dashed", color = "grey85") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey85") +

  coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "RDA space",
       x = paste0("RDA 1 (", p1, "%)"),
       y = paste0("RDA 2 (", p2, "%)")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(p)

ggsave("RDA_space_SNPs_outliers.png", p, width = 7, height = 6, dpi = 300)
ggsave("RDA_space_SNPs_outliers.pdf", p, width = 7, height = 6)

cat("已保存: RDA_space_SNPs_outliers.png / .pdf\n")

cat("\n=====================================\n")
cat("RDA分析完成！\n")
cat("=====================================\n")
cat("总结:\n")
cat("- 使用您的6个环境变量: bio_1, bio_3, bio_4, bio_10, bio_17, bio_18\n")
cat("- rdadapt方法与参考脚本一致：RDA loadings + Mahalanobis + BH\n")
cat("- Figure S5风格图：灰点=全部SNP, 紫点=BH<0.01 outliers, 箭头=环境变量\n")
cat("=====================================\n\n")

# 清理临时文件
if (file.exists("./tmp_snpinfo_from_vcf.tsv")) {
  file.remove("./tmp_snpinfo_from_vcf.tsv")
  cat("已清理临时文件: tmp_snpinfo_from_vcf.tsv\n")
}

