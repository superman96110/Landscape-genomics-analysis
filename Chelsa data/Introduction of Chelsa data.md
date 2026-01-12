# CHELSA Bioclim (bio01–bio19) 含义与单位

> 说明：
> - “Quarter” 指 **连续 3 个月**的窗口。
> - 最湿/最干季度按 **降水**确定；最暖/最冷季度按 **气温**确定。
> - 降水单位中 **1 kg m⁻² ≈ 1 mm**（水深）。

## 温度相关（bio01–bio11）

| 变量 | 英文名 | 含义（CHELSA 定义/计算口径） | 单位 |
|---|---|---|---|
| bio01 | Mean Annual Near-Surface Air Temperature | 全年 12 个月“月均温”的平均（年平均气温） | °C |
| bio02 | Mean Diurnal Near-Surface Air Temperature Range | 12 个月的 (tasmax − tasmin) 月均值的平均（平均日较差） | °C |
| bio03 | Isothermality | 100 × (bio02 ÷ bio07)，衡量昼夜温差相对年温差的比例 | 无量纲（常按 % 理解） |
| bio04 | Temperature Seasonality | 12 个月“月均温”的标准差（温度季节性） | °C/100 |
| bio05 | Mean Daily Maximum Temperature of the Warmest Month | 全年中“最暖月”的日最高气温（tasmax）的月均值 | °C |
| bio06 | Mean Daily Minimum Temperature of the Coldest Month | 全年中“最冷月”的日最低气温（tasmin）的月均值 | °C |
| bio07 | Annual Temperature Range | bio05 − bio06（年温差幅度） | °C |
| bio08 | Mean Temperature of the Wettest Quarter | 最湿的连续 3 个月里，月均温的平均 | °C |
| bio09 | Mean Temperature of the Driest Quarter | 最干的连续 3 个月里，月均温的平均 | °C |
| bio10 | Mean Temperature of the Warmest Quarter | 最暖的连续 3 个月里，月均温的平均 | °C |
| bio11 | Mean Temperature of the Coldest Quarter | 最冷的连续 3 个月里，月均温的平均 | °C |

## 降水相关（bio12–bio19）

| 变量 | 英文名 | 含义（CHELSA 定义/计算口径） | 单位 |
|---|---|---|---|
| bio12 | Annual Precipitation | 12 个月降水量之和（年降水） | kg m⁻² year⁻¹（≈ mm/year） |
| bio13 | Precipitation of the Wettest Month | 全年中最大月降水量 | kg m⁻² month⁻¹（≈ mm/month） |
| bio14 | Precipitation of the Driest Month | 全年中最小月降水量 | kg m⁻² month⁻¹（≈ mm/month） |
| bio15 | Precipitation Seasonality | 月降水变异系数：100 × SD ÷ mean | 无量纲（常按 % 理解） |
| bio16 | Mean Monthly Precipitation of the Wettest Quarter | 最湿连续 3 个月里，“月降水量”的平均（注意是 mean monthly，不是季度总量） | kg m⁻² month⁻¹（≈ mm/month） |
| bio17 | Mean Monthly Precipitation of the Driest Quarter | 最干连续 3 个月里，“月降水量”的平均 | kg m⁻² month⁻¹（≈ mm/month） |
| bio18 | Mean Monthly Precipitation of the Warmest Quarter | 最暖连续 3 个月里，“月降水量”的平均 | kg m⁻² month⁻¹（≈ mm/month） |
| bio19 | Mean Monthly Precipitation of the Coldest Quarter | 最冷连续 3 个月里，“月降水量”的平均 | kg m⁻² month⁻¹（≈ mm/month） |
