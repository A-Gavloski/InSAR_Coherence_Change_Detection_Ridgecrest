"""
Coherence change detection (CCD) utilities for mapping surface rupture.

Implements the core steps of the InSAR-correlation differencing approach of
Fielding et al. (2005), applied to Sentinel-1 coherence maps produced in SNAP:

1. Load the area of interest (AOI).
2. Put the pre-event and co-event coherence rasters on one common, metric grid.
3. Difference them (pre - co): positive values = coherence lost during the event.
4. Compare the result with mapped surface rupture (USGS ground truth),
   using the question numbering of notebook 03:
   Q1 distance decay - does coherence change fall off with distance from the
                       rupture, and over what distance? The peak half-width
                       sets the tolerance w used by Q2 and Q3.
   Q2 concentration  - fraction of mapped rupture length with a detection
                       within w (called "completeness" in the code, after
                       Wiedemann et al., 1998)
   Q3 correctness    - fraction of detections within w of mapped rupture
                       (a lower bound, because the USGS map is incomplete)
   Each is tested against copies of the rupture lines shifted several km away.

Notebook 02 uses steps 1-3; notebook 03 uses all four.

"""

from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import shapely
from rasterio import features
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from scipy import ndimage
from shapely.ops import transform as shapely_transform

# UTM zone 11N covers Ridgecrest and gives distances in metres.
DEFAULT_CRS = "EPSG:32611"


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_aoi(path: str | Path) -> gpd.GeoDataFrame:
    """Read an AOI polygon (KML, GeoJSON, shapefile...) and drop any Z values."""
    gdf = gpd.read_file(path)
    gdf["geometry"] = gdf.geometry.apply(
        lambda g: shapely_transform(lambda x, y, z=None: (x, y), g)
    )
    return gdf


def find_coherence_band(src: rasterio.io.DatasetReader) -> int:
    """Return the 1-based index of the coherence band in a SNAP GeoTIFF.

    SNAP's Terrain-Correction keeps every band of the interferogram product
    (i, q, coherence ...). If band descriptions are present, pick the one that
    contains "coh"; otherwise fall back to band 1 and warn.
    """
    for i, desc in enumerate(src.descriptions, start=1):
        if desc and "coh" in desc.lower():
            return i
    if src.count > 1:
        warnings.warn(
            f"{Path(src.name).name}: {src.count} bands and no 'coh' description; "
            "using band 1. Pass `band=` explicitly if coherence is a different band."
        )
    return 1


def read_coherence(path: str | Path, band: int | None = None):
    """Read a coherence band and its georeferencing.

    Returns (array as float32 with no-data set to NaN, profile).
    Coherence is defined on [0, 1]; exact zeros are SNAP's no-data fill.
    """
    with rasterio.open(path) as src:
        band = band or find_coherence_band(src)
        arr = src.read(band).astype("float32")
        nodata = src.nodata if src.nodata is not None else 0.0
        arr[(arr == nodata) | ~np.isfinite(arr)] = np.nan
        profile = src.profile.copy()
    return arr, profile


# --------------------------------------------------------------------------- #
# Gridding and differencing
# --------------------------------------------------------------------------- #
def common_grid(aoi: gpd.GeoDataFrame, crs: str = DEFAULT_CRS, resolution: float = 10.0):
    """Define a north-up grid covering the AOI: returns (transform, (rows, cols), crs)."""
    minx, miny, maxx, maxy = aoi.to_crs(crs).total_bounds
    cols = int(np.ceil((maxx - minx) / resolution))
    rows = int(np.ceil((maxy - miny) / resolution))
    return from_origin(minx, maxy, resolution, resolution), (rows, cols), crs


def to_grid(path: str | Path, transform, shape, crs, band: int | None = None) -> np.ndarray:
    """Resample one coherence raster onto the common grid (bilinear, NaN-aware)."""
    src_arr, prof = read_coherence(path, band)
    dst = np.full(shape, np.nan, dtype="float32")
    reproject(
        source=src_arr,
        destination=dst,
        src_transform=prof["transform"],
        src_crs=prof["crs"],
        src_nodata=np.nan,
        dst_transform=transform,
        dst_crs=crs,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
    )
    return dst


def mask_outside(arr: np.ndarray, aoi: gpd.GeoDataFrame, transform, crs) -> np.ndarray:
    """Set pixels outside the AOI polygon to NaN."""
    inside = features.geometry_mask(
        aoi.to_crs(crs).geometry, out_shape=arr.shape, transform=transform, invert=True
    )
    out = arr.copy()
    out[~inside] = np.nan
    return out


def coherence_difference(pre: np.ndarray, co: np.ndarray, normalized: bool = False) -> np.ndarray:
    """Coherence change, pre minus co-event.

    Positive values mean coherence was lost across the earthquake.
    normalized=True returns (pre - co) / (pre + co), which ranges -1..1 and
    down-weights areas that were already incoherent before the event.
    """
    diff = pre - co
    if normalized:
        with np.errstate(divide="ignore", invalid="ignore"):
            diff = diff / (pre + co)
    return diff


def write_geotiff(path: str | Path, arr: np.ndarray, transform, crs) -> None:
    """Save a single-band float32 GeoTIFF (NaN = no data)."""
    profile = dict(
        driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
        dtype="float32", crs=crs, transform=transform, nodata=np.nan, compress="lzw",
    )
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype("float32"), 1)


# --------------------------------------------------------------------------- #
# Validation against mapped rupture
#
#   dist    = distance_to_lines(lines, transform, diff.shape, crs)
#   prof    = distance_profile(diff, dist, step=20)            # Q1 profile
#   w       = decay_distance(diff, dist, prof, fraction=0.5)    # tolerance (half-width)
#   tail    = decay_distance(diff, dist, prof, fraction=0.1)    # Q1 finding, not a tolerance
#   points  = sample_lines(lines, crs)
#   sweep   = threshold_sweep(smooth(diff), dist, points, transform, w, thresholds)
#   f1      = 2 * c * r / (c + r)  with c = completeness, r = correctness
#   threshold = sweep.threshold[f1.idxmax()]                    # operating threshold
#   obs, null, p = shifted_null(smooth(diff), lines, transform, crs, threshold, w)
#
# Smooth before thresholding: a single noise pixel within w of a line would
# otherwise count up to 2*w of rupture length as detected.
# --------------------------------------------------------------------------- #
def distance_to_lines(lines: gpd.GeoDataFrame, transform, shape, crs) -> np.ndarray:
    """Distance (m) from every pixel to the nearest mapped rupture line."""
    on_line = features.rasterize(
        ((g, 1) for g in lines.to_crs(crs).geometry if g is not None and not g.is_empty),
        out_shape=shape, transform=transform, fill=0, all_touched=True, dtype="uint8",
    )
    if not on_line.any():
        raise ValueError("No rupture lines fall inside the grid - check CRS and AOI.")
    return ndimage.distance_transform_edt(on_line == 0) * abs(transform.a)


def distance_profile(diff: np.ndarray, dist: np.ndarray, max_dist: float = 3000, step: float = 100) -> pd.DataFrame:
    """Q1: median / mean coherence change in distance bins away from the rupture.

    If the method works, the change should peak at 0 m and decay with distance.
    Notebook 03 uses step=20 m out to 3 km.
    """
    edges = np.arange(0, max_dist + step, step)
    valid = np.isfinite(diff)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = valid & (dist >= lo) & (dist < hi)
        vals = diff[sel]
        rows.append({
            "dist_m": (lo + hi) / 2,
            "median": np.median(vals) if vals.size else np.nan,
            "mean": vals.mean() if vals.size else np.nan,
            "p25": np.percentile(vals, 25) if vals.size else np.nan,
            "p75": np.percentile(vals, 75) if vals.size else np.nan,
            "n": int(vals.size),
        })
    return pd.DataFrame(rows)


def decay_distance(diff, dist, profile: pd.DataFrame, background_min: float = 3000,
                   fraction: float = 0.1) -> float:
    """Q1: distance at which the coherence change has decayed by a given fraction.

    Background = median change beyond `background_min` metres. Returns the first
    bin centre past the peak where the median excess over background falls to
    `fraction` of the peak excess, or NaN if it never does.

    Notebook 03 calls it twice:
      fraction=0.5 - peak half-width, the rupture itself; used as the tolerance w
      fraction=0.1 - tail width, the zone of distributed damage; reported only.
    The tail is too wide to use as w: almost any line would count as detected.
    """
    background = np.nanmedian(diff[dist > background_min])
    excess = profile["median"].to_numpy() - background
    peak = np.nanargmax(excess)
    after = np.flatnonzero((np.arange(len(excess)) > peak) & (excess <= fraction * excess[peak]))
    return float(profile["dist_m"].iloc[after[0]]) if after.size else np.nan


def smooth(diff: np.ndarray, size: int = 5) -> np.ndarray:
    """NaN-tolerant median filter to suppress coherence speckle before thresholding."""
    filled = np.where(np.isfinite(diff), diff, 0)
    out = ndimage.median_filter(filled, size=size)
    out[~np.isfinite(diff)] = np.nan
    return out


def detection_metrics(
    diff: np.ndarray, dist: np.ndarray, threshold: float,
    pixel_size: float, buffer_m: float = 300,
) -> dict:
    """Legacy precision / recall score with a fixed buffer. Not used by the notebooks.

    detected  = diff > threshold
    precision = share of detected pixels within `buffer_m` of a mapped rupture
    recall    = share of rupture-line pixels with a detection within `buffer_m`

    Superseded by validation_scores, which uses the data-derived tolerance w
    and measures recall (concentration) by rupture length rather than by
    rasterized line pixels.
    """
    valid = np.isfinite(diff)
    detected = valid & (diff > threshold)
    near = dist <= buffer_m

    precision = (detected & near).sum() / max(detected.sum(), 1)

    on_line = (dist == 0) & valid
    if detected.any() and on_line.any():
        dist_to_detection = ndimage.distance_transform_edt(~detected) * pixel_size
        recall = (dist_to_detection[on_line] <= buffer_m).mean()
    else:
        recall = 0.0

    return {
        "threshold": threshold,
        "buffer_m": buffer_m,
        "detected_pct_of_area": 100 * detected.sum() / max(valid.sum(), 1),
        "precision": float(precision),
        "recall": float(recall),
    }


def sample_lines(lines: gpd.GeoDataFrame, crs, step: float = 5.0) -> np.ndarray:
    """Points every `step` metres along the lines, so point counts measure length.

    MultiLineStrings are split into parts; non-line geometries are skipped.
    """
    pts = []
    for g in lines.to_crs(crs).geometry.explode(index_parts=False):
        if g is None or g.is_empty or g.geom_type != "LineString":
            continue
        d = np.arange(step / 2, g.length, step)
        pts.append(shapely.get_coordinates(shapely.line_interpolate_point(g, d)))
    return np.vstack(pts) if pts else np.empty((0, 2))


def _at_points(arr: np.ndarray, points: np.ndarray, transform) -> np.ndarray:
    """Raster values at map coordinates (NaN outside the grid)."""
    inv = ~transform
    x, y = points[:, 0], points[:, 1]
    cols = inv.a * x + inv.b * y + inv.c
    rows = inv.d * x + inv.e * y + inv.f
    rows, cols = np.floor(rows).astype(int), np.floor(cols).astype(int)
    inside = (rows >= 0) & (rows < arr.shape[0]) & (cols >= 0) & (cols < arr.shape[1])
    out = np.full(len(points), np.nan)
    out[inside] = arr[rows[inside], cols[inside]]
    return out


def validation_scores(diff, dist, points, transform, threshold: float, w: float,
                      valid=None, background_min: float = 3000) -> dict:
    """Q1, Q2 and Q3 scores for one threshold and tolerance w.

    near_excess  - Q1: median change within w minus median change beyond background_min
    completeness - Q2 concentration: share of rupture length (sample points on
                   valid data) with a detection within w
    correctness  - Q3: share of detected pixels within w of a rupture line (lower bound)
    enrichment   - correctness / share of the scene within w (1 = chance level)
    """
    valid = np.isfinite(diff) if valid is None else valid
    detected = valid & (diff > threshold)
    if detected.any():
        to_detection = ndimage.distance_transform_edt(~detected) * abs(transform.a)
    else:
        to_detection = np.full(diff.shape, np.inf)

    on_valid = _at_points(valid, points, transform) == 1  # skip points on no-data
    near = valid & (dist <= w)
    far = valid & (dist > background_min)
    chance = near.sum() / max(valid.sum(), 1)

    completeness = np.mean(_at_points(to_detection, points[on_valid], transform) <= w) if on_valid.any() else np.nan
    correctness = (detected & near).sum() / detected.sum() if detected.any() else np.nan
    return {
        "threshold": threshold,
        "w": w,
        "completeness": float(completeness),
        "correctness": float(correctness),
        "enrichment": float(correctness / chance) if chance else np.nan,
        "near_excess": float(np.median(diff[near]) - np.median(diff[far])) if near.any() and far.any() else np.nan,
        "detected_pct_of_area": 100 * detected.sum() / max(valid.sum(), 1),
    }


def threshold_sweep(diff, dist, points, transform, w: float, thresholds) -> pd.DataFrame:
    """Q2 (completeness) and Q3 (correctness) across thresholds.

    df.attrs["area"] is the step-wise area under the correctness-vs-completeness
    curve. Notebook 03 picks the operating threshold as the one with the highest
    harmonic mean of the two scores.
    """
    df = pd.DataFrame([validation_scores(diff, dist, points, transform, t, w) for t in thresholds])
    s = df.dropna(subset=["completeness", "correctness"]).sort_values("completeness")
    df.attrs["area"] = float(np.sum(np.diff(np.r_[0, s["completeness"]]) * s["correctness"]))
    return df


def shifted_null(diff, lines: gpd.GeoDataFrame, transform, crs, threshold: float, w: float,
                 shifts_m=(3000, 5000, 7000), n_directions: int = 8,
                 exclusion_m: float | None = None, step: float = 5.0):
    """Test Q1-Q3 against copies of the rupture lines shifted away from the rupture.

    Each copy is the whole line set moved `shift` metres in one of `n_directions`
    azimuths. Pixels within `exclusion_m` (default 2*w) of the real lines are
    masked, so a copy is never scored on the real rupture's signal. Copies that
    keep less than half as much line length on valid data as the real lines
    have are skipped (lines on no-data, e.g. outside the swath, don't count).

    Returns (observed scores, DataFrame of null scores, dict of one-sided p-values).
    With n copies the smallest possible p is 1 / (n + 1).
    """
    dist = distance_to_lines(lines, transform, diff.shape, crs)
    points = sample_lines(lines, crs, step)
    observed = validation_scores(diff, dist, points, transform, threshold, w)
    n_real = int((_at_points(np.isfinite(diff), points, transform) == 1).sum())

    exclusion_m = 2 * w if exclusion_m is None else exclusion_m
    away = np.isfinite(diff) & (dist > exclusion_m)
    base = lines.to_crs(crs).geometry

    rows = []
    for shift in shifts_m:
        for az in np.arange(n_directions) * 360 / n_directions:
            dx, dy = shift * np.sin(np.radians(az)), shift * np.cos(np.radians(az))
            moved = gpd.GeoDataFrame(geometry=base.translate(dx, dy), crs=crs)
            pts = sample_lines(moved, crs, step)
            pts = pts[_at_points(away, pts, transform) == 1]
            if len(pts) < 0.5 * n_real:
                continue
            try:
                moved_dist = distance_to_lines(moved, transform, diff.shape, crs)
            except ValueError:
                continue
            scores = validation_scores(diff, moved_dist, pts, transform, threshold, w, valid=away)
            rows.append({"shift_m": shift, "azimuth": az, **scores})

    null = pd.DataFrame(rows)
    keys = ("completeness", "correctness", "near_excess")
    p = {k: (1 + (null[k] >= observed[k]).sum()) / (1 + len(null)) if len(null) else np.nan
         for k in keys}
    return observed, null, p
