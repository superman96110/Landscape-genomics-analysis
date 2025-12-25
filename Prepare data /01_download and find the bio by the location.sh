#先从数据库中现在bio数据，案例中使用的是CHELSA数据库，参考网站为：https://www.chelsa-climate.org/datasets/chelsa_bioclim
#气候数据的文件格式为tif结尾的文件
#wget 下载19个气候因子后，根据经纬度坐标提取bio数据生成csv文件，相关提取py脚本如下，179服务器

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
从CHELSA BioClim tif (1981-2010, V2.1) 中按经纬度点提取 bio_1..bio_19，
输出为 CSV。

运行位置建议：/data/supeng/env/chelsa/1981-2010
"""

import os
import re
import glob
from typing import List, Tuple

import numpy as np
import pandas as pd
import rasterio


# ========= 配置（按你的目录/文件名写死也可以）=========
WORKDIR = "/data/supeng/env/chelsa/1981-2010"

coordinates_file = os.path.join(WORKDIR, "location_pig.txt")

# 你的 tif 文件名是：CHELSA_bio01_1981-2010_V.2.1.tif 这样
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
    # 读取坐标文件并确保保留所有行
    df = pd.read_csv(path, delim_whitespace=True, header=None, engine="python")

    # 如果文件有表头，默认跳过第一行；若无表头，直接读取所有行
    cols_lower = [str(c).strip().lower() for c in df.columns]
    lon_idx, lat_idx = None, None

    for i, c in enumerate(cols_lower):
        if c in ("lon", "longitude", "经度"):
            lon_idx = i
        if c in ("lat", "latitude", "纬度"):
            lat_idx = i

    # 没识别到列名：当作无表头，取前两列
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
    兼容 bio01 / bio1 两种写法，这里你的文件是 bio01..bio19。
    """
    found = []
    for p in glob.glob(pattern):
        base = os.path.basename(p)
        # 匹配 CHELSA_bio01_1981-2010_V.2.1.tif 中的 01
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


def sample_one_raster(ds: rasterio.io.DatasetReader, coords: List[Tuple[float, float]]) -> np.ndarray:
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
                # 浮点 nodata 用 isclose 判断
                if (isinstance(val, (float, np.floating)) and isinstance(nodata, (float, np.floating))
                        and np.isclose(val, nodata, equal_nan=True)) or (val == nodata):
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
    # 0) 保险：切到工作目录（这样相对路径也能用）
    os.chdir(WORKDIR)
    print(f"当前工作目录：{os.getcwd()}")

    # 1) 读坐标
    coords_df = read_coords(coordinates_file)
    if coords_df.empty:
        raise ValueError("坐标列表为空。")

    coords = list(zip(coords_df["Longitude"].tolist(), coords_df["Latitude"].tolist()))
    print(f"坐标点数量：{len(coords)}")
    print(f"坐标文件：{coordinates_file}")

    # 2) 找到 bio_1..bio_19
    bio_list = find_bio_tifs(bio_glob_pattern)

    if not bio_list:
        raise FileNotFoundError(
            f"没有找到任何 CHELSA bio tif。检查 pattern：{bio_glob_pattern}"
        )

    # 3) 初始化输出表
    out_df = coords_df.copy()
    for i in range(1, 20):
        out_df[f"bio_{i}"] = None

    # 4) 逐个 tif 抽样
    #    GDAL 缓存适当加大一点能提升随机读性能（单位MB）
    with rasterio.Env(GDAL_CACHEMAX=1024):
        for num, path in bio_list:
            print(f"抽样 {os.path.basename(path)} -> bio_{num} ...")
            with rasterio.open(path) as ds:
                vals = sample_one_raster(ds, coords)
                out_df[f"bio_{num}"] = vals
            print(f"完成 bio_{num}")

    # 5) 写 CSV
    cols = ["Longitude", "Latitude"] + [f"bio_{i}" for i in range(1, 20)]
    out_df.to_csv(output_csv, index=False, columns=cols, float_format="%.6f")
    print(f"完成！输出：{output_csv}")


if __name__ == "__main__":
    main()
