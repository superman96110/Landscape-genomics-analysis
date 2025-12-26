#先从数据库中现在bio数据，案例中使用的是CHELSA数据库，参考网站为：https://www.chelsa-climate.org/datasets/chelsa_bioclim
#气候数据的文件格式为tif结尾的文件
#wget 下载19个气候因子后，根据经纬度坐标提取bio数据生成csv文件，相关提取py脚本如下，179服务器

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
从CHELSA BioClim tif (1981-2010, V2.1) 中按经纬度点提取 bio_1..bio_19，
同时输出原始坐标和实际采样的栅格中心坐标，输出为 CSV。
运行位置建议：/data/supeng/env/chelsa/1981-2010
"""
import os
import re
import glob
from typing import List, Tuple
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol

# ========= 配置（按你的目录/文件名写死也可以）=========
WORKDIR = "/data/supeng/env/chelsa/1981-2010"
coordinates_file = os.path.join(WORKDIR, "location_pig.txt")
bio_glob_pattern = os.path.join(WORKDIR, "CHELSA_bio*_1981-2010_V.2.1.tif")
output_csv = os.path.join(WORKDIR, "extracted_CHELSA_1981-2010_location_pig.csv")

# =====================================================

def read_coords(path: str) -> pd.DataFrame:
    """
    读取经纬度文件。兼容：
    - 有表头(如 Lon Lat / Longitude Latitude)
    - 无表头两列数值
    仅保留前两列作为 经度、纬度
    """
    df = pd.read_csv(path, delim_whitespace=True, header=None, engine="python")
    cols_lower = [str(c).strip().lower() for c in df.columns]
    lon_idx, lat_idx = None, None
    
    for i, c in enumerate(cols_lower):
        if c in ("lon", "longitude", "经度"):
            lon_idx = i
        if c in ("lat", "latitude", "纬度"):
            lat_idx = i
    
    if lon_idx is None or lat_idx is None:
        if df.shape[1] < 2:
            raise ValueError(f"{path} 至少需要两列（经度、纬度）。")
        lon_idx, lat_idx = 0, 1
    
    lon = pd.to_numeric(df.iloc[:, lon_idx], errors="coerce")
    lat = pd.to_numeric(df.iloc[:, lat_idx], errors="coerce")
    
    bad = lon.isna() | lat.isna()
    if bad.any():
        print(f"警告：坐标文件中有 {int(bad.sum())} 行无法转换为数值，已丢弃。")
    
    return pd.DataFrame({"Longitude": lon[~bad].values, "Latitude": lat[~bad].values})


def find_bio_tifs(pattern: str, expected=range(1, 20)) -> List[Tuple[int, str]]:
    """
    查找并按 bio 编号排序返回 (bio_num, path) 列表。
    """
    found = []
    for p in glob.glob(pattern):
        base = os.path.basename(p)
        m = re.search(r"CHELSA_bio(\d{1,2})_1981-2010_V\.2\.1\.tif$", base)
        if m:
            num = int(m.group(1))
            if num in expected:
                found.append((num, p))
    
    found.sort(key=lambda x: x[0])
    missing = sorted(set(expected) - {n for n, _ in found})
    
    if missing:
        print(f"警告：缺失以下 bio 文件：{missing}")
    else:
        print("已找到 bio_1 至 bio_19 的全部文件。")
    
    return found


def get_pixel_center_coords(ds: rasterio.io.DatasetReader, 
                            coords: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    根据输入坐标，计算实际采样的栅格中心坐标
    """
    transform = ds.transform
    pixel_centers = []
    
    for lon, lat in coords:
        # 将地理坐标转换为像素行列号
        row, col = rowcol(transform, lon, lat)
        
        # 将像素行列号转换回地理坐标（栅格中心）
        # transform * (col+0.5, row+0.5) 得到像素中心坐标
        sampled_lon = transform.c + (col + 0.5) * transform.a
        sampled_lat = transform.f + (row + 0.5) * transform.e
        
        pixel_centers.append((sampled_lon, sampled_lat))
    
    return pixel_centers


def sample_one_raster(ds: rasterio.io.DatasetReader, 
                      coords: List[Tuple[float, float]]) -> np.ndarray:
    """
    对一个已打开的数据集，使用 dataset.sample 对所有坐标抽样，
    返回 shape=(N,) 的数组（None 表示 nodata / nan / 越界等）
    """
    values = []
    nodata = ds.nodata
    
    for v in ds.sample(coords):
        val = v[0] if len(v) else np.nan
        
        # 处理 nodata / nan
        if nodata is not None:
            try:
                if (isinstance(val, (float, np.floating)) and 
                    isinstance(nodata, (float, np.floating)) and 
                    np.isclose(val, nodata, equal_nan=True)) or (val == nodata):
                    values.append(None)
                    continue
            except Exception:
                pass
        
        if isinstance(val, (float, np.floating)) and np.isnan(val):
            values.append(None)
        else:
            try:
                values.append(float(val))
            except Exception:
                values.append(None)
    
    return np.array(values, dtype=object)


def main():
    os.chdir(WORKDIR)
    print(f"当前工作目录：{os.getcwd()}")
    
    # 1) 读坐标
    coords_df = read_coords(coordinates_file)
    if coords_df.empty:
        raise ValueError("坐标列表为空。")
    
    coords = list(zip(coords_df["Longitude"].tolist(), 
                     coords_df["Latitude"].tolist()))
    print(f"坐标点数量：{len(coords)}")
    print(f"坐标文件：{coordinates_file}")
    
    # 2) 找到 bio_1..bio_19
    bio_list = find_bio_tifs(bio_glob_pattern)
    if not bio_list:
        raise FileNotFoundError(
            f"没有找到任何 CHELSA bio tif。检查 pattern：{bio_glob_pattern}"
        )
    
    # 3) 初始化输出表（包含原始坐标和实际采样坐标）
    out_df = pd.DataFrame({
        "Original_Longitude": coords_df["Longitude"],
        "Original_Latitude": coords_df["Latitude"]
    })
    
    # 添加实际采样坐标列（将在第一个tif处理时填充）
    out_df["Sampled_Longitude"] = None
    out_df["Sampled_Latitude"] = None
    
    # 添加bio数据列
    for i in range(1, 20):
        out_df[f"bio_{i}"] = None
    
    # 4) 逐个 tif 抽样
    with rasterio.Env(GDAL_CACHEMAX=1024):
        for idx, (num, path) in enumerate(bio_list):
            print(f"抽样 {os.path.basename(path)} -> bio_{num} ...")
            
            with rasterio.open(path) as ds:
                # 第一个文件时，计算并保存实际采样坐标
                if idx == 0:
                    print(f"  栅格分辨率: {abs(ds.transform.a):.6f}° x {abs(ds.transform.e):.6f}°")
                    print(f"  (约 {abs(ds.transform.a)*111:.2f} km)")
                    
                    pixel_centers = get_pixel_center_coords(ds, coords)
                    out_df["Sampled_Longitude"] = [pc[0] for pc in pixel_centers]
                    out_df["Sampled_Latitude"] = [pc[1] for pc in pixel_centers]
                
                # 抽样数据
                vals = sample_one_raster(ds, coords)
                out_df[f"bio_{num}"] = vals
            
            print(f"完成 bio_{num}")
    
    # 5) 计算坐标偏移距离
    out_df["Offset_Longitude"] = out_df["Sampled_Longitude"] - out_df["Original_Longitude"]
    out_df["Offset_Latitude"] = out_df["Sampled_Latitude"] - out_df["Original_Latitude"]
    out_df["Offset_Distance_deg"] = np.sqrt(
        out_df["Offset_Longitude"]**2 + out_df["Offset_Latitude"]**2
    )
    
    # 6) 写 CSV
    cols = (["Original_Longitude", "Original_Latitude", 
             "Sampled_Longitude", "Sampled_Latitude",
             "Offset_Longitude", "Offset_Latitude", "Offset_Distance_deg"] + 
            [f"bio_{i}" for i in range(1, 20)])
    
    out_df.to_csv(output_csv, index=False, columns=cols, float_format="%.8f")
    
    print(f"\n完成！输出：{output_csv}")
    print(f"\n坐标偏移统计：")
    print(f"  平均偏移: {out_df['Offset_Distance_deg'].mean():.6f}°")
    print(f"  最大偏移: {out_df['Offset_Distance_deg'].max():.6f}°")
    print(f"  (1° ≈ 111 km)")


if __name__ == "__main__":
    main()
