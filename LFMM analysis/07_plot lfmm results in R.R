#提取合并文件的关键列绘图
{bash}
awk '{print 1,$1,$3,$6}' all_chromosomes_LFMM.tsv > lfmm_bio3_results.txt

setwd("F:/caas/毕业课题/第五章_景观基因组/horse/chelsa/lfmm/bio3/")
library(CMplot)
d <- read.table("lfmm_bio3_results.txt", header = TRUE, stringsAsFactors = FALSE)
# 确保列名正确
colnames(d)[1:4] <- c("SNP", "CHR", "BP", "P")

# 确保CHR为数字，避免CMplot把它当成字符串导致上色异常
d$CHR <- gsub("^chr", "", d$CHR, ignore.case = TRUE)
d$CHR <- as.integer(d$CHR)
d$BP  <- as.numeric(d$BP)
d$P   <- as.numeric(d$P)

# 去掉异常行
d <- d[is.finite(d$CHR) & is.finite(d$BP) & is.finite(d$P), ]
d <- d[order(d$CHR, d$BP), ]

m <- d[, c("SNP","CHR","BP","P")]

CMplot(
  m,
  plot.type = "m",
  LOG10 = TRUE,
  type = "h",
  col = c("#7B1E1E", "grey70"),
  chr.den.col = NULL,              # ✅去除染色体密度色块
  band = 0,                        # ✅关键参数:去除染色体间空白(默认=1)
  file = "png",
  file.output = TRUE,
  file.name = "LFMM_BIO1_Manhattan_line2color",
  dpi = 300,
  width = 12,
  height = 6
)
