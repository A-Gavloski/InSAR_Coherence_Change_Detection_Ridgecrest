# SNAP processing: Sentinel-1 coherence

Each interferometric pair (pre-event and co-event) is processed with the same graph, [`coherence_graph.xml`](coherence_graph.xml). The output is a terrain-corrected coherence GeoTIFF for each pair, which is the input to notebooks 02 and 03.

## Processing chain

```
Read (reference) ─► TOPSAR-Split ─► Apply-Orbit-File ─┐
                                                      ├─► Back-Geocoding ─► Enhanced-Spectral-Diversity ─► Interferogram (+ coherence)
Read (secondary) ─► TOPSAR-Split ─► Apply-Orbit-File ─┘                                                          │
                                                                                                                 ▼
                                              Write (GeoTIFF) ◄─ Terrain-Correction ◄─ TopoPhaseRemoval ◄─ TOPSAR-Deburst
```

## Step-by-step

| # | Operator | Key settings | Why |
|---|---|---|---|
| 1 | **Read** ×2 | Reference = earlier date, secondary = later date | Back-Geocoding coregisters the secondary to the reference geometry |
| 2 | **TOPSAR-Split** | Subswath **IW2**, polarization **VV**, bursts **1–9** | IW2 covers the rupture zone. Processing one subswath cuts the data volume by about two thirds. VV has higher SNR over land than VH. |
| 3 | **Apply-Orbit-File** | Sentinel Precise (auto download), polynomial degree 3 | Precise orbits (about 5 cm accuracy) are needed for accurate coregistration |
| 4 | **Back-Geocoding** | DEM: SRTM 3Sec. DEM resampling: bicubic. Image resampling: bisinc 5-point | Uses the DEM and orbits to coregister the two SLCs to sub-pixel accuracy |
| 5 | **Enhanced-Spectral-Diversity** | Defaults: 512 × 512 windows, coherence threshold 0.3 | TOPS mode needs azimuth coregistration to about 0.001 pixel. Otherwise phase jumps appear at burst boundaries. |
| 6 | **Interferogram** | Subtract flat-earth phase. **Include coherence. Window 10 (range) × 3 (azimuth)**. Square pixel. | Estimates coherence γ = \|⟨s₁s₂*⟩\| / √(⟨\|s₁\|²⟩⟨\|s₂\|²⟩) over the window. 10 × 3 IW pixels is about 40 m × 40 m on the ground. A larger window reduces estimation bias but blurs narrow features. |
| 7 | **TOPSAR-Deburst** | Defaults | Merges the bursts into one continuous image |
| 8 | **TopoPhaseRemoval** | SRTM 3Sec | Removes topographic phase from the interferogram |
| 9 | **Terrain-Correction** | SRTM 3Sec, 10 m pixel spacing | Geocodes the coherence map so it can be gridded to UTM 11N in notebooks 02 and 03 |
| 10 | **Write** | GeoTIFF | Read by `rasterio` in notebooks 02 and 03 |

## Running the graph

**Graph Builder (GUI):** *Tools → GraphBuilder → Load* `coherence_graph.xml`, choose the two input `.zip` files in the Read tabs, set the output name, then *Run*.

**Command line:** the graph uses `${reference}`, `${secondary}` and `${output}` placeholders, so the same file processes both pairs:

```bash
# pre-event pair
gpt snap/coherence_graph.xml \
  -Preference=data/raw/S1A_IW_SLC__1SDV_20190610T015047_20190610T015114_027611_031DB6_D8FA.zip \
  -Psecondary=data/raw/S1A_IW_SLC__1SDV_20190622T015048_20190622T015115_027786_0322F1_7A8E.zip \
  -Poutput=data/processed/coherence_pre_20190610_20190622.tif

# co-event pair
gpt snap/coherence_graph.xml \
  -Preference=data/raw/S1A_IW_SLC__1SDV_20190704T015049_20190704T015116_027961_03283A_191E.zip \
  -Psecondary=data/raw/S1A_IW_SLC__1SDV_20190716T015049_20190716T015116_028136_032D86_6DF2.zip \
  -Poutput=data/processed/coherence_co_20190704_20190716.tif
```

## Notes

- **Band order:** Terrain-Correction keeps every band (i, q, intensity, coherence). Notebook 02 looks for a band named `coh`. If the GeoTIFF has no band names, set `COH_BAND` in the notebook to the coherence band's position, as shown in SNAP's Band Maths or Information window. Another option is to list only the `coh_*` band under *Source Bands* in Terrain-Correction, which also makes the files much smaller.
- **File size:** outputs larger than 4 GB need the `GeoTIFF-BigTIFF` writer.
- **Possible refinement:** coherence is estimated before topographic phase removal, so on steep slopes residual topographic fringes inside the window can slightly lower coherence. SNAP's standalone *Coherence* operator has a "subtract topographic phase" option that avoids this. The effect is small in the flat Indian Wells and Searles valleys.
