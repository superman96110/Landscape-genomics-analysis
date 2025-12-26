#需要准备刚刚提取好的csv文件
#还需要准备每个个体对应的经纬度坐标的excel文件，再用次R脚本匹配bio数据

# ============================
# 1. 加载需要的包
# ============================
library(readxl)    # 读取 Excel 文件
library(dplyr)     # 数据处理包
library(openxlsx)  # 写 Excel 文件

# 设置工作目录（根据您的路径调整）
setwd("F:/caas/毕业课题/pig data/CHELSA/")

# ============================
# 2. 读取第一个 Excel 文件（品种信息，包含经纬度）
# ============================
file1 <- "F:/caas/毕业课题/pig data/worldclim/更新后的数据.xlsx"
df1 <- read_excel(file1, sheet = "545_ind")

# ============================
# 3. 读取第二个 CSV 文件（包含经纬度和气候因子）
# ============================
file2 <- "F:/caas/毕业课题/pig data/CHELSA/extracted_CHELSA_1981-2010_location_pig.csv"
df2 <- read.csv(file2, header = TRUE, stringsAsFactors = FALSE)

# ============================
# 4. 确保经纬度字段是 numeric 类型（避免字符导致匹配失败）
# ============================
df1$经度 <- as.numeric(df1$经度)
df1$纬度 <- as.numeric(df1$纬度)

df2$Original_Longitude <- as.numeric(df2$Original_Longitude)
df2$Original_Latitude  <- as.numeric(df2$Original_Latitude)

# ============================
# 5. 合并数据：df1 的经纬度 对应 df2 的 Original 经纬度
# ============================
df_merged <- df1 %>%
    left_join(df2, by = c("经度" = "Original_Longitude", "纬度" = "Original_Latitude"))

# ============================
# 6. 检查合并后行数是否正确（应该仍然 545 行）
# ============================
cat("原始数据 df1 行数:", nrow(df1), "\n")
cat("合并后 df_merged 行数:", nrow(df_merged), "\n")

# ============================
# 7. 查看有多少行没有匹配到 bio（bio_1 为 NA）
# ============================
missing_num <- sum(is.na(df_merged$bio_1))
cat("未匹配到 bio 的行数:", missing_num, "\n")

# ============================
# 8. 输出到新的 Excel 文件
# ============================
output_file <- "F:/caas/毕业课题/pig data/CHELSA/merged_with_bio.xlsx"
write.xlsx(df_merged, output_file)

cat("✅ 合并完成！结果已输出到:", output_file, "\n")
