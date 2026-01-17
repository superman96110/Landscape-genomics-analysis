library(readr)
library(dplyr)
library(readxl)
library(ggplot2)
library(stringr)

## ---- 输入文件 ----
freq_file   <- "snp_8_20644555_byBreed.frq.strat"
excel_file  <- "F:/caas/毕业课题/第五章_景观基因组/horse/chelsa/horse_ind_bio.xlsxx"
excel_sheet <- "group_bio19+ele"   # 也可以写 "group_bio19+ele"


## ---- freq.strat 中的列 ----
freq_breed_col <- "CLST"     # breed 列
freq_value_col <- "MAF"      # 用于分析的频率列
freq_n_col     <- "NCHROBS"  # 用于显示 n= 的列

## ---- Excel 中的列 ----
excel_breed_col <- "Breed"               # breed 列
excel_value_col <- "Original_Latitude"   # 你想分析的列（bio_1 / Original_Latitude 等）

## ---- 是否对 Excel 的变量取绝对值 ----
use_abs_value <- TRUE   # TRUE = abs(纬度)；FALSE = 原值

## ---- 输出图名 ----
out_png <- "Freq_vs_Env_dual_axis.png"

## =====================================================
## 1️⃣ 读取 freq.strat
## =====================================================
freq <- read_table(freq_file, col_types = cols(.default = col_guess())) %>%
    rename(
        Breed = all_of(freq_breed_col),
        Freq  = all_of(freq_value_col),
        N     = all_of(freq_n_col)
    ) %>%
    mutate(
        Breed = str_squish(as.character(Breed)),
        Freq  = as.numeric(Freq),
        N     = as.integer(N)
    ) %>%
    select(Breed, Freq, N)

## =====================================================
## 2️⃣ 读取 Excel
## =====================================================
env <- read_excel(excel_file, sheet = excel_sheet) %>%
    rename_with(~str_replace_all(.x, "\\s+", "")) %>%
    rename(
        Breed = all_of(excel_breed_col),
        Env   = all_of(excel_value_col)
    ) %>%
    mutate(
        Breed = str_squish(as.character(Breed)),
        Env   = as.numeric(Env),
        Env   = if (use_abs_value) abs(Env) else Env
    ) %>%
    select(Breed, Env) %>%
    distinct()

## =====================================================
## 3️⃣ 合并 & 排序
## =====================================================
dat <- freq %>%
    left_join(env, by = "Breed") %>%
    filter(!is.na(Freq), !is.na(Env)) %>%
    arrange(Freq)

## x轴标签：Breed (n=)
dat <- dat %>%
    mutate(
        Breed_label = sprintf("%s (n=%d)", Breed, N),
        Breed_label = factor(Breed_label, levels = Breed_label)
    )

## =====================================================
## 4️⃣ 相关性
## =====================================================
cor_sp <- cor.test(dat$Freq, dat$Env, method = "spearman")

message(sprintf(
    "Spearman correlation: rho = %.3f, p = %.3g",
    cor_sp$estimate, cor_sp$p.value
))

## =====================================================
## 5️⃣ 双 Y 轴缩放
## =====================================================
freq_min <- min(dat$Freq)
freq_max <- max(dat$Freq)
env_min  <- min(dat$Env)
env_max  <- max(dat$Env)

scale_factor <- (freq_max - freq_min) / (env_max - env_min)
if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1

dat <- dat %>%
    mutate(Env_scaled = (Env - env_min) * scale_factor + freq_min)

## =====================================================
## 6️⃣ 作图
## =====================================================
p <- ggplot(dat, aes(x = Breed_label)) +
    geom_col(aes(y = Freq), width = 0.75) +
    geom_line(aes(y = Env_scaled, group = 1), linewidth = 1) +
    geom_point(aes(y = Env_scaled), size = 2) +
    scale_y_continuous(
        name = freq_value_col,
        sec.axis = sec_axis(
            trans = ~ (. - freq_min) / scale_factor + env_min,
            name = excel_value_col
        )
    ) +
    labs(
        x = "Breed (sorted by frequency; n from freq.strat)",
        title = paste(freq_value_col, "vs", excel_value_col),
        subtitle = sprintf("Spearman rho = %.3f, p = %.3g",
                           cor_sp$estimate, cor_sp$p.value)
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

print(p)
ggsave(out_png, p, width = 13, height = 6, dpi = 300)

## =====================================================
## 7️⃣ 可选：散点图（强烈推荐）
## =====================================================
p2 <- ggplot(dat, aes(x = Env, y = Freq)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw(base_size = 12) +
    labs(
        x = excel_value_col,
        y = freq_value_col,
        title = paste(freq_value_col, "vs", excel_value_col),
        subtitle = sprintf("Spearman rho = %.3f, p = %.3g",
                           cor_sp$estimate, cor_sp$p.value)
    )

print(p2)
