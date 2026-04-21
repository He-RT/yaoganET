# 农田灌溉助手（MVP）

地图选田块 → 从 Google Earth Engine 拉取最新无云 Sentinel-2 影像 + 计算 ET0。

GEE 项目：`chuang-yaogan`

## 运行

```bash
cd backend
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 首次：登录 GEE（或设置服务账号，见下）
earthengine authenticate --project=chuang-yaogan

uvicorn app:app --reload --port 8000
```

打开 http://localhost:8000 ，在地图上画多边形并点"分析"。

### 使用服务账号（可选，部署推荐）

```bash
export GEE_PROJECT=chuang-yaogan
export GEE_SERVICE_ACCOUNT_EMAIL=xxx@chuang-yaogan.iam.gserviceaccount.com
export GEE_SERVICE_ACCOUNT_KEY=/path/to/key.json
```

## 接口

`POST /api/analyze`

```json
{ "geometry": { "type": "Polygon", "coordinates": [...] }, "days": 30 }
```

返回：
- `sentinel2`: 最新一景的 ID、成像日期、云量、区域平均 NDVI、可叠加的 RGB tile URL
- `et0`: 回溯窗口内基于 ERA5-Land 的 FAO-56 Penman-Monteith 日均 ET0 (mm/day)

## 路线图
- [ ] Kc × ET0 → 作物实际需水量（ETc）
- [ ] 根据土壤持水量 & 近期降水推荐灌水量
- [ ] 地块保存、历史曲线
