#参考使用对应的github软件https://github.com/simlal/bioclim-data-extraction?utm_source=chatgpt.com
#准备经纬度坐标location.txt
#pwd:/data/supeng/env/chelsa/bioclim/bioclim-data-extraction

from pathlib import Path
import pandas as pd

from scripts.data_extraction import CrsDataPoint, extract_multiple_bioclim_elev

# ====== 你的输入/输出（你给的） ======
LOC_PATH = Path("location_horse.txt")   # 两列：lon lat
OUT_CSV = Path("bio_1981_2010.csv")     # 输出
# ====================================

def main():
    # 读取 lon lat（支持空格/Tab/逗号分隔；支持 # 注释行）
    df = pd.read_csv(LOC_PATH, sep=None, engine="python", header=None, comment="#")
    if df.shape[1] < 2:
        raise ValueError(f"{LOC_PATH} 至少需要两列：lon lat")
    df = df.iloc[:, :2].copy()
    df.columns = ["lon", "lat"]

    # 构造 CrsDataPoint 列表（你的坐标就是 EPSG:4326）
    points = [
        CrsDataPoint(f"pt{i+1}", epsg=4326, x=float(row["lon"]), y=float(row["lat"]))
        for i, row in df.iterrows()
    ]

    # 直接调用仓库函数：提取 CHELSA (1981-2010) bio1~bio19（内部会做 offset+scale 修正）
    out_df = extract_multiple_bioclim_elev(points, dataset="chelsa", trimmed=True)

    # 保存
    out_df.to_csv(OUT_CSV, index=False)
    print("Saved:", OUT_CSV.resolve())

if __name__ == "__main__":
    main()
