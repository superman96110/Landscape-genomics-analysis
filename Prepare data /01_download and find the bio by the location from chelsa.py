#先从数据库中现在bio数据，案例中使用的是CHELSA数据库，参考网站为：https://www.chelsa-climate.org/datasets/chelsa_bioclim
#气候数据的文件格式为tif结尾的文件
#wget 下载19个气候因子后，根据经纬度坐标提取bio数据生成csv文件，相关提取py脚本如下，179服务器
#由于bio因子存在scales和offsets，需要进行矫正后再输出正确数值
#python find_bio_chelsa.py

#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
从CHELSA BioClim tif (1981-2010, V2.1) 中按经纬度点提取 bio_1..bio_19，
并对每个 tif 分别应用其自身的 scale/offset：real = raw * scale + offset。

额外增强：
1) 运行开始会扫描并打印每个 bio 的 dtype/nodata/scale/offset（方便核对）
2) 抽样时显式传入该 bio 对应的 scale/offset（避免混用）
3) 输出原始坐标、采样像元中心坐标、偏移统计、bio_1..bio_19（还原后）

运行目录建议：/data/supeng/env/chelsa/1981-2010
"""

import os
import re
import glob
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol


# ========= 配置 =========
WORKDIR = "/data/supeng/env/chelsa/1981-2010"
coordinates_file = os.path.join(WORKDIR, "location_horse.txt")
bio_glob_pattern = os.path.join(WORKDIR, "CHELSA_bio*_1981-2010_V.2.1.tif")
output_csv = os.path.join(WORKDIR, "extracted_CHELSA_1981-2010_location_horse.csv")
# 可选：把每个 bio 的 scale/offset 也写一份表出来
meta_csv = os.path.join(WORKDIR, "CHELSA_bio_scaling_metadata.csv")
# =======================


def read_coords(path: str) -> pd.DataFrame:
    """读取坐标文件：默认取前两列为 Longitude Latitude。"""
    df = pd.read_csv(path, delim_whitespace=True, header=None, engine="python")
    if df.shape[1] < 2:
        raise ValueError(f"{path} 至少需要两列（经度、纬度）。")
    lon = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    lat = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    bad = lon.isna() | lat.isna()
    if bad.any():
        print(f"警告：坐标文件中有 {int(bad.sum())} 行无法转换为数值，已丢弃。")
    return pd.DataFrame({"Longitude": lon[~bad].values, "Latitude": lat[~bad].values})


def find_bio_tifs(pattern: str, expected=range(1, 20)) -> List[Tuple[int, str]]:
    """查找并按 bio 编号排序返回 (bio_num, path)。支持 bio01/bio1 两种。"""
    found = []
    for p in glob.glob(pattern):
        base = os.path.basename(p)
        m = re.search(r"CHELSA_bio(\d{1,2})_1981-2010_V\.2\.1\.tif$", base)
        if m:
            num = int(m.group(1))  # bio01 -> 1
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
    """根据输入坐标计算实际采样的像元中心坐标。"""
    transform = ds.transform
    out = []
    for lon, lat in coords:
        row, col = rowcol(transform, lon, lat)
        sampled_lon = transform.c + (col + 0.5) * transform.a
        sampled_lat = transform.f + (row + 0.5) * transform.e
        out.append((sampled_lon, sampled_lat))
    return out


def read_scale_offset(ds: rasterio.io.DatasetReader) -> Tuple[float, float]:
    """
    读取 scale/offset。
    优先用 rasterio 的 ds.scales/ds.offsets（你机器上 bio01 已验证有值）。
    若缺失则默认 (1.0, 0.0)。
    """
    scale = 1.0
    offset = 0.0

    try:
        if getattr(ds, "scales", None) and ds.scales[0] is not None:
            scale = float(ds.scales[0])
        if getattr(ds, "offsets", None) and ds.offsets[0] is not None:
            offset = float(ds.offsets[0])
    except Exception:
        pass

    return scale, offset


def scan_all_metadata(bio_list: List[Tuple[int, str]]) -> Dict[int, Dict]:
    """扫描每个 bio 的 dtype/nodata/scale/offset，形成字典并打印汇总。"""
    meta: Dict[int, Dict] = {}

    print("\n===== 扫描每个 bio 的 nodata/scale/offset（每个 tif 各自独立）=====")
    rows = []
    for num, path in bio_list:
        with rasterio.open(path) as ds:
            scale, offset = read_scale_offset(ds)
            info = {
                "bio": num,
                "file": os.path.basename(path),
                "dtype": ds.dtypes[0],
                "nodata": ds.nodata,
                "scale": scale,
                "offset": offset,
                "crs": str(ds.crs),
                "res_deg_x": float(abs(ds.transform.a)),
                "res_deg_y": float(abs(ds.transform.e)),
            }
            meta[num] = info
            rows.append(info)

            print(f"bio_{num:02d}: dtype={info['dtype']}, nodata={info['nodata']}, "
                  f"scale={info['scale']}, offset={info['offset']}  ({info['file']})")

    # 可选写出一份 meta csv，便于留档
    pd.DataFrame(rows).to_csv(meta_csv, index=False)
    print(f"已写出元数据汇总：{meta_csv}")
    print("===========================================================\n")

    return meta


def sample_with_scaling(ds: rasterio.io.DatasetReader,
                        coords: List[Tuple[float, float]],
                        nodata,
                        scale: float,
                        offset: float) -> np.ndarray:
    """抽样并应用 real = raw*scale + offset。"""
    out = []
    for v in ds.sample(coords):
        if not len(v):
            out.append(None)
            continue
        raw = v[0]

        # nodata
        if nodata is not None:
            try:
                if raw == nodata:
                    out.append(None)
                    continue
            except Exception:
                pass

        # nan
        if isinstance(raw, (float, np.floating)) and np.isnan(raw):
            out.append(None)
            continue

        try:
            real = float(raw) * float(scale) + float(offset)
        except Exception:
            out.append(None)
            continue

        if isinstance(real, float) and np.isnan(real):
            out.append(None)
        else:
            out.append(real)

    return np.array(out, dtype=object)


def main():
    os.chdir(WORKDIR)
    print(f"当前工作目录：{os.getcwd()}")

    # 1) 坐标
    coords_df = read_coords(coordinates_file)
    if coords_df.empty:
        raise ValueError("坐标列表为空。")
    coords = list(zip(coords_df["Longitude"].tolist(), coords_df["Latitude"].tolist()))
    print(f"坐标点数量：{len(coords)}")
    print(f"坐标文件：{coordinates_file}")

    # 2) bio 文件
    bio_list = find_bio_tifs(bio_glob_pattern)
    if not bio_list:
        raise FileNotFoundError(f"没有找到任何 CHELSA bio tif。检查 pattern：{bio_glob_pattern}")

    # 3) 扫描每个 bio 的 scale/offset/nodata（关键：每个都独立）
    meta = scan_all_metadata(bio_list)

    # 4) 输出表初始化
    out_df = pd.DataFrame({
        "Original_Longitude": coords_df["Longitude"],
        "Original_Latitude": coords_df["Latitude"],
        "Sampled_Longitude": None,
        "Sampled_Latitude": None,
    })
    for i in range(1, 20):
        out_df[f"bio_{i}"] = None

    # 5) 抽样
    with rasterio.Env(GDAL_CACHEMAX=1024):
        for idx, (num, path) in enumerate(bio_list):
            info = meta[num]
            print(f"抽样 {info['file']} -> bio_{num} (scale={info['scale']}, offset={info['offset']})")

            with rasterio.open(path) as ds:
                # 第一个文件：计算采样像元中心坐标
                if idx == 0:
                    pixel_centers = get_pixel_center_coords(ds, coords)
                    out_df["Sampled_Longitude"] = [x for x, _ in pixel_centers]
                    out_df["Sampled_Latitude"] = [y for _, y in pixel_centers]

                vals = sample_with_scaling(
                    ds=ds,
                    coords=coords,
                    nodata=info["nodata"],
                    scale=info["scale"],
                    offset=info["offset"],
                )
                out_df[f"bio_{num}"] = vals

    # 6) 偏移统计
    out_df["Offset_Longitude"] = out_df["Sampled_Longitude"] - out_df["Original_Longitude"]
    out_df["Offset_Latitude"] = out_df["Sampled_Latitude"] - out_df["Original_Latitude"]
    out_df["Offset_Distance_deg"] = np.sqrt(out_df["Offset_Longitude"]**2 + out_df["Offset_Latitude"]**2)

    # 7) 写 CSV
    cols = ([
        "Original_Longitude", "Original_Latitude",
        "Sampled_Longitude", "Sampled_Latitude",
        "Offset_Longitude", "Offset_Latitude", "Offset_Distance_deg",
    ] + [f"bio_{i}" for i in range(1, 20)])

    out_df.to_csv(output_csv, index=False, columns=cols, float_format="%.8f")

    print(f"\n完成！输出：{output_csv}")
    print("坐标偏移统计：")
    print(f"  平均偏移: {out_df['Offset_Distance_deg'].mean():.6f}°")
    print(f"  最大偏移: {out_df['Offset_Distance_deg'].max():.6f}°")
    print("  (1° ≈ 111 km)")


if __name__ == "__main__":
    main()
