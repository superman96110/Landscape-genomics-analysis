#plink --file /data/supeng/env/env/horse/lfmm_vcf/chr15 --extract range locus.range --within group1.txt --freq --out snp_15_53591860_byBreed


library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggnewscale)
library(viridis)

# ================== 1. 路径设置 ==================
frq_file   <- "/home/lingjiang/supeng/horse/env/chelas/plot/lfmm_bio1/snp_8_20644555_byBreed.frq.strat"
group_file <- "/home/lingjiang/supeng/horse/env/chelas/plot/lfmm_bio1/group.txt"
bio1_file  <- "/home/lingjiang/supeng/horse/env/chelas/CHELSA_bio01_1981-2010_V.2.1.tif"
out_png    <- "/home/lingjiang/supeng/horse/env/chelas/plot/lfmm_bio1/bio1_landmask_refalt_pies_8_20644555.png"

# ================== 2. 参数设置 ==================
REF_ALLELE <- "G"
ALT_ALLELE <- "A"

# 绘图参数
pie_scale <- 1.5         # 饼图基础大小 (根据你的地图范围适当调大)
n_circle  <- 80          # 圆滑度
xlim_world <- c(-180, 180)
ylim_world <- c(-60, 85)

# 温度配色
temp_colors <- c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
                 "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026")

# ================== 3. 数据读取与处理 ==================

# --- 3.1 读取 Group 并计算经纬度中心 ---
grp <- read_tsv(
  group_file,
  col_names = c("Sample", "Breed", "Lon", "Lat"),
  col_types = cols(Sample="c", Breed="c", Lon="d", Lat="d"),
  show_col_types = FALSE
) %>%
  filter(!is.na(Breed), !is.na(Lon), !is.na(Lat))

# 计算每个 Breed 的中心位置
breed_loc <- grp %>%
  group_by(Breed) %>%
  summarise(
    Lon = mean(Lon),
    Lat = mean(Lat),
    Count = n(),
    .groups = "drop"
  )

# --- 3.2 读取频率数据 ---
frq <- read_table(
  frq_file,
  col_types = cols(
    CHR="c", SNP="c", CLST="c", A1="c", A2="c", 
    MAF="d", MAC="i", NCHROBS="i"
  ),
  show_col_types = FALSE
)

# --- 3.3 合并数据 (修复：之前缺失这一步导致无法画图) ---
freq_breed <- frq %>%
  inner_join(breed_loc, by = c("CLST" = "Breed")) %>%
  mutate(
    # 判定 A1 是否为 REF，计算真实 REF 频率
    p_REF = case_when(
      A1 == REF_ALLELE ~ MAF,
      A1 == ALT_ALLELE ~ 1 - MAF,
      TRUE ~ NA_real_ 
    ),
    p_ALT = 1 - p_REF,
    # 动态调整半径：基础大小 * sqrt(样本量)
    r = pie_scale * sqrt(Count) / 3 
  ) %>%
  filter(!is.na(p_REF))

# ================== 4. 准备饼图多边形 ==================
make_wedge <- function(x0, y0, r, start, end, n = 80) {
  th <- seq(start, end, length.out = n)
  data.frame(
    x = c(x0, x0 + r*cos(th), x0),
    y = c(y0, y0 + r*sin(th), y0)
  )
}

pie_poly <- freq_breed %>%
  select(CLST, Lon, Lat, r, p_REF, p_ALT) %>%
  mutate(
    a0_ref = 0,
    a1_ref = 2*pi*p_REF,
    a0_alt = 2*pi*p_REF,
    a1_alt = 2*pi
  ) %>%
  rowwise() %>%
  do({
    df <- .
    res <- data.frame()
    if(df$p_REF > 0.001) {
      refp <- make_wedge(df$Lon, df$Lat, df$r, df$a0_ref, df$a1_ref, n_circle)
      refp$part <- "REF"
      refp$group <- paste(df$CLST, "REF", sep="__")
      res <- bind_rows(res, refp)
    }
    if(df$p_ALT > 0.001) {
      altp <- make_wedge(df$Lon, df$Lat, df$r, df$a0_alt, df$a1_alt, n_circle)
      altp$part <- "ALT"
      altp$group <- paste(df$CLST, "ALT", sep="__")
      res <- bind_rows(res, altp)
    }
    res
  }) %>%
  ungroup()

# ================== 5. 处理BIO1数据 (修复：白色条带问题) ==================
bio1 <- rast(bio1_file)
if (!is.lonlat(bio1)) bio1 <- project(bio1, "EPSG:4326")

# 1. 先根据绘图范围裁剪 (Crop)，这比处理全球数据更安全且快
ext_crop <- ext(xlim_world[1], xlim_world[2], ylim_world[1], ylim_world[2])
bio1_crop <- crop(bio1, ext_crop)

# 2. 获取陆地矢量 (直接转为 vect，不要做 st_union)
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
land_vect <- vect(world_sf) 

# 3. 遮盖 (Mask)
# 这会把海洋部分设为 NA，保留陆地数据
bio1_land <- mask(bio1_crop, land_vect)

# 计算温度范围，防止 scale 报错
val_range <- minmax(bio1_land)
temp_min <- floor(val_range[1])
temp_max <- ceiling(val_range[2])

# ================== 6. 绘图 ==================
p <- ggplot() +
  # --- 层1: 温度 Raster (增加 maxcell 防止渲染丢失) ---
  geom_spatraster(data = bio1_land, alpha = 0.9, maxcell = 8e5) +
  scale_fill_gradientn(
    name = "Annual Mean\nTemperature (°C)",
    colours = temp_colors,
    na.value = "transparent",
    limits = c(temp_min, temp_max),
    guide = guide_colorbar(
      order = 1,
      barwidth = 15, barheight = 0.8,
      title.position = "top", title.hjust = 0.5,
      frame.colour = "black", ticks.colour = "black"
    )
  ) +
  
  # --- 层2: 新的填充标尺给饼图 ---
  ggnewscale::new_scale_fill() +
  
  # --- 层3: 饼图 ---
  geom_polygon(
    data = pie_poly,
    aes(x = x, y = y, group = group, fill = part),
    color = "black", linewidth = 0.2, alpha = 1
  ) +
  scale_fill_manual(
    name = "Allele Frequency",
    values = c(REF = "#FFFFFF", ALT = "#2C3E50"),
    labels = c(REF = paste0("REF (", REF_ALLELE, ")"), ALT = paste0("ALT (", ALT_ALLELE, ")")),
    guide = guide_legend(
      order = 2,
      override.aes = list(color = "black", linewidth = 0.3),
      keywidth = 1.5, keyheight = 1.5
    )
  ) +
  
  # --- 层4: 地图边框 (fill = NA, 避免遮盖) ---
  geom_sf(data = world_sf, fill = NA, color = "gray30", linewidth = 0.3) +
  
  # --- 坐标系与主题 ---
  coord_sf(xlim = xlim_world, ylim = ylim_world, expand = FALSE) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.background = element_rect(fill = "#E8F4F8", color = NA), # 海洋颜色
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    axis.title = element_blank()
  ) +
  labs(
    title = "Allele Frequency Distribution and BIO1 Temperature",
    subtitle = paste0("SNP: 8:20644555 | Ref: ", REF_ALLELE, " / Alt: ", ALT_ALLELE)
  )

# ================== 7. 保存 ==================
ggsave(out_png, p, width = 14, height = 7.5, dpi = 350, bg = "white")
cat("✓ 图片已保存:", out_png, "\n")
