"""Irrigation helper backend.

Endpoints:
  POST /api/analyze  { geometry: GeoJSON Polygon, days: int? }
    -> latest cloud-free Sentinel-2 scene over the field, NDVI mean,
       and ET0 (mm/day) from ERA5-Land for the same date window.
"""
from __future__ import annotations

import os
from datetime import date, timedelta
from pathlib import Path
from typing import Any

import ee
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

GEE_PROJECT = os.environ.get("GEE_PROJECT", "chuang-yaogan")

_ee_ready = False


def ensure_ee() -> None:
    global _ee_ready
    if _ee_ready:
        return
    sa_key = os.environ.get("GEE_SERVICE_ACCOUNT_KEY")
    sa_email = os.environ.get("GEE_SERVICE_ACCOUNT_EMAIL")
    if sa_key and sa_email and Path(sa_key).exists():
        creds = ee.ServiceAccountCredentials(sa_email, sa_key)
        ee.Initialize(credentials=creds, project=GEE_PROJECT)
    else:
        # falls back to the locally cached user credentials from `earthengine authenticate`
        ee.Initialize(project=GEE_PROJECT)
    _ee_ready = True


class AnalyzeRequest(BaseModel):
    geometry: dict[str, Any] = Field(..., description="GeoJSON geometry (Polygon)")
    days: int = Field(30, ge=1, le=120, description="look-back window in days")


def _latest_s2(aoi: ee.Geometry, start: str, end: str) -> ee.Image | None:
    coll = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
        .sort("system:time_start", False)
    )
    size = coll.size().getInfo()
    if size == 0:
        return None
    return ee.Image(coll.first())


def _era5_et0(aoi: ee.Geometry, start: str, end: str) -> dict[str, Any]:
    """FAO-56 Penman-Monteith ET0 from ERA5-Land daily aggregates.

    Uses the daily-aggregated collection which exposes Tmax, Tmin, dewpoint,
    wind components, surface solar radiation, and surface pressure.
    """
    coll = (
        ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
        .filterDate(start, end)
        .filterBounds(aoi)
    )

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

        # Net shortwave proxy: ERA5-Land provides surface_net_solar_radiation_sum (J/m^2/day).
        rn = img.select("surface_net_solar_radiation_sum").divide(1_000_000)  # MJ/m^2/day

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
        num = (
            delta.multiply(rn)
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

    et0_coll = coll.map(compute)
    size = et0_coll.size().getInfo()
    if size == 0:
        return {"available": False}

    mean_img = et0_coll.mean()
    stats = mean_img.reduceRegion(
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


app = FastAPI(title="Irrigation Advisor")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.post("/api/analyze")
def analyze(req: AnalyzeRequest) -> dict[str, Any]:
    try:
        ensure_ee()
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, f"GEE init failed: {exc}")

    try:
        aoi = ee.Geometry(req.geometry)
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(400, f"invalid geometry: {exc}")

    end = date.today()
    start = end - timedelta(days=req.days)
    start_s, end_s = start.isoformat(), end.isoformat()

    s2 = _latest_s2(aoi, start_s, end_s)
    s2_info: dict[str, Any] = {"available": False}
    if s2 is not None:
        ndvi = s2.normalizedDifference(["B8", "B4"]).rename("NDVI")
        ndvi_stats = ndvi.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=aoi,
            scale=10,
            maxPixels=1e9,
        ).getInfo()
        props = s2.getInfo().get("properties", {})
        tile_url = s2.visualize(
            bands=["B4", "B3", "B2"], min=0, max=3000
        ).getMapId()
        s2_info = {
            "available": True,
            "image_id": s2.get("system:index").getInfo(),
            "acquired": props.get("system:time_start"),
            "cloud_pct": props.get("CLOUDY_PIXEL_PERCENTAGE"),
            "ndvi_mean": ndvi_stats.get("NDVI"),
            "tile_url_template": tile_url["tile_fetcher"].url_format,
        }

    et0 = _era5_et0(aoi, start_s, end_s)

    return {
        "window": {"start": start_s, "end": end_s},
        "sentinel2": s2_info,
        "et0": et0,
    }


# serve the static frontend from ../frontend
_frontend = Path(__file__).resolve().parent.parent / "frontend"
if _frontend.exists():
    app.mount("/", StaticFiles(directory=str(_frontend), html=True), name="frontend")
