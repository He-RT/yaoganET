# ET0 计算说明

本项目后端使用 **FAO-56 Penman-Monteith 参考蒸散公式**（日尺度），所有气象输入来自 Google Earth Engine 数据集 `ECMWF/ERA5_LAND/DAILY_AGGR`。计算在 GEE 服务端逐日、逐像元完成，最后在用户绘制的田块多边形上做空间平均、在回溯窗口上做时间平均，返回 mm/天。

- 实现文件：[`backend/app.py`](../backend/app.py) 中的 `_era5_et0`
- 数据集文档：[ERA5-Land Daily Aggregated on GEE](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR)
- 理论参考：[FAO Irrigation and Drainage Paper 56, Chapter 4](https://www.fao.org/3/x0490e/x0490e00.htm)

---

## 1. 公式

FAO-56 Penman-Monteith：

$$
ET_0 = \frac{0.408\,\Delta\,(R_n - G) + \gamma\,\dfrac{900}{T+273}\,u_2\,(e_s - e_a)}{\Delta + \gamma\,(1 + 0.34\,u_2)}
$$

| 符号 | 含义 | 单位 |
|---|---|---|
| $ET_0$ | 参考蒸散（短草面） | mm·day⁻¹ |
| $R_n$ | 净辐射 | MJ·m⁻²·day⁻¹ |
| $G$ | 土壤热通量，日尺度取 0 | MJ·m⁻²·day⁻¹ |
| $T$ | 日平均气温 | °C |
| $u_2$ | 2 m 处风速 | m·s⁻¹ |
| $e_s$ | 饱和水汽压 | kPa |
| $e_a$ | 实际水汽压 | kPa |
| $\Delta$ | 饱和水汽压-温度曲线斜率 | kPa·°C⁻¹ |
| $\gamma$ | 干湿计常数 | kPa·°C⁻¹ |

> **关于 0.408**：把净辐射的能量通量换成蒸发的水等效厚度。
> $\dfrac{1}{\lambda \cdot \rho_w} = \dfrac{1}{2.45\text{ MJ/kg} \times 1000\text{ kg/m}^3} = 0.408\text{ mm/(MJ/m}^2\text{)}$。
> 当 $R_n$ 已是 MJ/m²/day 时，直接乘 0.408 就得到 mm/day。
>
> **不要混淆 0.0864**：那是 W/m² → MJ/m²/day 的时间尺度换算（86400 s × 10⁻⁶ MJ/J）。只在原始输入是瞬时功率 W/m² 时才需要先乘 0.0864。ERA5-Land 的 `*_radiation_sum` 波段已经是 J/m²/day，除以 10⁶ 直接得到 MJ/m²/day，不需要 0.0864。

---

## 2. 输入波段与换算

| 公式量 | ERA5-Land 波段 | 原始单位 → 目标单位 |
|---|---|---|
| $T_{max}$ | `temperature_2m_max` | K → °C （− 273.15） |
| $T_{min}$ | `temperature_2m_min` | K → °C |
| $T$ | — | $(T_{max}+T_{min})/2$ |
| $T_{dew}$ | `dewpoint_temperature_2m` | K → °C |
| 10 m 风 $u,v$ | `u_component_of_wind_10m`, `v_component_of_wind_10m` | m/s |
| $u_2$ | — | $\sqrt{u^2+v^2}\times 0.748$（10 m → 2 m 对数廓线缩减因子） |
| $P$ | `surface_pressure` | Pa → kPa （/1000） |
| $R_{ns}$ | `surface_net_solar_radiation_sum` | J/m²/day → MJ/m²/day （/1e6） |
| $R_{nl}$ | `surface_net_thermal_radiation_sum` | J/m²/day → MJ/m²/day，**通常为负**（向外长波辐射） |
| $R_n$ | — | $R_{ns} + R_{nl}$ |

---

## 3. 中间量

**饱和水汽压**（Tetens 近似）：
$$
e^\circ(T) = 0.6108\,\exp\!\left(\frac{17.27\,T}{T + 237.3}\right)
$$

**饱和水汽压与实际水汽压**：
$$
e_s = \frac{e^\circ(T_{max}) + e^\circ(T_{min})}{2},\qquad
e_a = e^\circ(T_{dew})
$$

**曲线斜率**：
$$
\Delta = \frac{4098 \cdot e^\circ(T)}{(T + 237.3)^2}
$$

**干湿计常数**（简化式）：
$$
\gamma = 0.000665\,P
$$

**10 m → 2 m 风速换算**（FAO-56 eq. 47）：
$$
u_2 = u_{10} \cdot \frac{4.87}{\ln(67.8\cdot 10 - 5.42)} \approx u_{10}\cdot 0.748
$$

---

## 4. 代码片段（节选自 `backend/app.py`）

### 4.1 加载数据集

```python
coll = (
    ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
    .filterDate(start, end)
    .filterBounds(aoi)
)
```

### 4.2 单日逐像元计算 ET0

```python
def compute(img: ee.Image) -> ee.Image:
    # Temperatures (Kelvin -> Celsius).
    tmax = img.select("temperature_2m_max").subtract(273.15)
    tmin = img.select("temperature_2m_min").subtract(273.15)
    tmean = tmax.add(tmin).divide(2)
    tdew = img.select("dewpoint_temperature_2m").subtract(273.15)

    # Wind speed at 2m from 10m components (log-wind reduction factor 0.748).
    u10 = img.select("u_component_of_wind_10m")
    v10 = img.select("v_component_of_wind_10m")
    u2 = u10.pow(2).add(v10.pow(2)).sqrt().multiply(0.748)

    # Pressure kPa.
    p = img.select("surface_pressure").divide(1000)

    # Net radiation Rn = Rns (net shortwave) + Rnl (net longwave, typically negative).
    # ERA5-Land provides both as daily sums in J/m^2; convert to MJ/m^2/day.
    rns = img.select("surface_net_solar_radiation_sum").divide(1_000_000)
    rnl = img.select("surface_net_thermal_radiation_sum").divide(1_000_000)
    rn = rns.add(rnl)

    # Saturation & actual vapor pressure (FAO-56).
    def es(t: ee.Image) -> ee.Image:
        return t.multiply(17.27).divide(t.add(237.3)).exp().multiply(0.6108)

    es_tmax = es(tmax)
    es_tmin = es(tmin)
    es_mean = es_tmax.add(es_tmin).divide(2)
    ea = es(tdew)

    # Slope of vapor pressure curve.
    delta = (
        es(tmean)
        .multiply(4098)
        .divide(tmean.add(237.3).pow(2))
    )

    # Psychrometric constant.
    gamma = p.multiply(0.000665)

    # FAO-56 reference eq. (G ~ 0 at daily step).
    # 0.408 converts MJ/m^2/day of net radiation to mm/day of evaporated water.
    num = (
        delta.multiply(rn).multiply(0.408)
        .add(
            gamma
            .multiply(900)
            .divide(tmean.add(273))
            .multiply(u2)
            .multiply(es_mean.subtract(ea))
        )
    )
    den = delta.add(gamma.multiply(u2.multiply(0.34).add(1)))
    et0 = num.divide(den).rename("ET0")
    return et0.copyProperties(img, ["system:time_start"])
```

### 4.3 时间 + 空间聚合

```python
et0_coll = coll.map(compute)
size = et0_coll.size().getInfo()
if size == 0:
    return {"available": False}

mean_img = et0_coll.mean()                          # 时间平均
stats = mean_img.reduceRegion(                      # 空间平均
    reducer=ee.Reducer.mean(),
    geometry=aoi,
    scale=1000,
    maxPixels=1e9,
).getInfo()
return {
    "available": True,
    "days": size,
    "et0_mean_mm_per_day": stats.get("ET0"),
}
```

---

## 5. 聚合流程

1. 在 `[start, end]` 窗口（默认今日往前回溯 30 天）取 ERA5-Land 每日一张影像
2. `coll.map(compute)` 对每天每像元独立计算 ET0，形成日 ET0 影像序列
3. `et0_coll.mean()` 对时间维度取均值，得到一张"窗口期日均 ET0"栅格
4. `reduceRegion(mean, scale=1000)` 在田块多边形上做空间平均，`scale` 对齐 ERA5-Land 原生 ~9 km 分辨率（1 km 已足够，更小不增加精度）
5. 返回标量 `et0_mean_mm_per_day` 和参与计算的有效天数

---

## 6. 已知精度折中

| 简化 | 影响 | 原因 |
|---|---|---|
| 露点用日均 `dewpoint_temperature_2m` | $e_s-e_a$ 略偏大 → ET0 小幅偏高 | ERA5-Land 日产品未提供与 $T_{max}/T_{min}$ 对应的露点 |
| 风速用 $\sqrt{\bar u^2+\bar v^2}$ | 风向乱的日子偏低 | 日聚合只给 u、v 分量的均值，非标量风速均值 |
| $G=0$ | 日尺度几乎无误差 | FAO-56 推荐 |
| `scale=1000` | 田块小于 9 km 时内部无空间差异 | ERA5-Land 原生分辨率限制 |

若要追求更高精度，可改用 `ECMWF/ERA5_LAND/HOURLY` 在小时步长计算（$T$、$e_a$、$u_2$、$R_n$ 完全对时），再按天求和。对灌溉推荐而言当前日尺度精度已足够。

---

## 7. 验证样例

| 位置 | 日期窗口 | ET0（mm/day） | 备注 |
|---|---|---|---|
| 北京城区 (116.40, 39.90) ~4 km² | 2026-03-22 – 2026-04-21 | 3.12 | 4 月北京常年 3–4 mm/day，合理 |
| 河南郑州 (113.65, 34.75) ~1 km² | 2026-03-22 – 2026-04-21 | ~4 | 春季华北平原典型值 |

---

## 8. 后续路线

- [ ] 作物系数 Kc × ET0 → 实际作物需水量 ETc
- [ ] 扣除 ERA5 的 `total_precipitation_sum` → 净灌溉需水量
- [ ] 结合土壤含水量（SMAP / ERA5 soil moisture）得到推荐灌水量（mm 或 m³）
- [ ] 可选：切换到小时级 ERA5 提高精度
