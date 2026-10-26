# Data

Only the AOI polygon and the USGS rupture traces are committed. Everything else is downloaded or generated, and is git-ignored.

```
data/
├── aoi/AOI_Ridgecrest_CA.kml     committed: 0.5° × 0.5° box, -117.85 to -117.35 E, 35.55 to 36.05 N
├── raw/                          Sentinel-1 SLC .zip files (notebook 01), about 4–5 GB each
├── processed/                    SNAP coherence GeoTIFFs (PreEventCoherence.tif, PostEventCoherence.tif)
│                                 and validation_results.csv (notebook 03)
└── usgs/USGS_rupturetrace_v4.kml USGS rupture traces
```

## Sentinel-1

Downloaded by `notebooks/01_search_download_sentinel1.ipynb` from the Alaska Satellite Facility. This needs a free [NASA Earthdata Login](https://urs.earthdata.nasa.gov/) token, set as the `EARTHDATA_TOKEN` environment variable.

## USGS surface rupture (ground truth)

Ponti et al. (2020), U.S. Geological Survey data release, https://doi.org/10.5066/P9BZ5IJ9

1. Open the data release page and download the rupture-map shapefile (zip).
2. Unzip it into `data/usgs/`.
3. Set `RUPTURE_PATH` in notebook 02 to the `.shp` file name.

USGS data releases are public domain.
