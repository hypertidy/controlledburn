# Architecture Analysis: exactextract Ecosystem

## 1. Core C++ Library (`exactextract/`)

The core is a well-layered C++ library with a single hard dependency: **GEOS ≥ 3.5**. GDAL is optional (only needed for the CLI binary). The architecture follows a clean pipeline:

```
FeatureSource → Processor → RasterCellIntersection → RasterStats → Operation → OutputWriter
```

### Key abstractions

| Abstraction | Purpose |
|---|---|
| **`RasterSource`** | Abstract I/O: `read_box(bbox)` → typed raster tile. No GDAL dependency. |
| **`FeatureSource`** | Abstract iterator over vector features with GEOS geometries |
| **`Operation`** | Defines what to compute (stat name + value raster + weight raster + options) |
| **`RasterStats<T>`** | Templated accumulator: incrementally computes ~30 statistics from coverage fractions + values + weights |
| **`Processor`** | Orchestrates the loop. Two strategies: **feature-sequential** (iterate features, read rasters per-bbox) and **raster-sequential** (load all features into R-tree, iterate raster chunks) |
| **`OutputWriter`** | Abstract sink for result features |

The computational core — `raster_cell_intersection()` — computes **exact sub-pixel coverage fractions** by tracing polygon edges across grid cells (not sampling). It supports both polygons (area fractions) and lines (length fractions).

### Extension points

All five major interfaces (`RasterSource`, `FeatureSource`, `Feature`, `OutputWriter`, `Operation`) are abstract/virtual and designed for subclassing. This is how both the Python and R packages plug in without requiring GDAL.

---

## 2. Python Package — Binding Mechanism

### How it wraps the core

| Layer | Technology |
|---|---|
| **C++→Python bridge** | **pybind11** (`_exactextract` compiled module) |
| **Build system** | scikit-build-core + CMake; links against `libexactextract` + GEOS |
| **Trampoline pattern** | Each core abstract class has a pybind11 trampoline (e.g., `PyRasterSourceBase`) that delegates C++ virtual calls back to Python methods |

The Python package **uses the full C++ pipeline** — `Processor`, `Operation`, `StatsRegistry`, `OutputWriter` — all run in C++. Python implements the I/O interfaces:

- **Raster data** flows as: Python raster object → Python `read_window()` → numpy array → C++ `NumPyRaster<T>` (zero-copy where possible)
- **Geometry** flows as: Python feature → WKB bytes or GeoJSON string → C++ GEOS parse
- **Results** flow as: C++ `Feature` → Python `Writer.write()` → dict/DataFrame/OGR layer

### What the Python side adds

The Python layer is primarily a **format-detection and convenience wrapper**:

- Auto-detects input formats: GDAL datasets, rasterio, xarray, geopandas, fiona, GeoJSON dicts, QGIS layers, file paths
- Wraps each in the appropriate `RasterSource`/`FeatureSource` subclass
- Handles CRS mismatch warnings (via osgeo.osr or pyproj)
- Provides multiple output backends: GeoJSON dicts, pandas DataFrame, GeoDataFrame, GDAL file output, QGIS layers
- Supports tqdm progress bars

### User-defined functions

Python supports **UDFs** via `PythonOperation`: any callable taking `(values_masked, coverage_fractions[, weights_masked])` is called back from C++ with numpy masked arrays.

---

## 3. R Package (`exactextractr/`) — Binding Mechanism

### How it wraps the core

| Layer | Technology |
|---|---|
| **C++→R bridge** | **Rcpp** (`.Call()` interface) |
| **Build system** | autoconf `configure` + `Makevars`; cherry-picks ~13 core `.cpp` files to compile directly |
| **Core library inclusion** | **Vendored as a subtree** at `src/exactextract/` |

Critically, the R package **does NOT use** the core's `Processor`, `Operation`, `StatsRegistry`, or `OutputWriter` pipeline. Instead it has its own simpler architecture:

```
R lapply() over features → per-feature .Call() into C++ →
  → raster_cell_intersection() → RasterStats<double> → return Rcpp matrix/list to R
```

- **Raster data** flows as: R calls `.getValuesBlock()` on terra/raster objects → `Rcpp::NumericVector` → `NumericVectorRaster` (implements `AbstractRaster<double>`)
- **Geometry** flows as: `sf::st_as_binary()` → `Rcpp::RawVector` → `GEOSWKBReader_read_r()`
- **Results** flow as: C++ fills `Rcpp::NumericMatrix` or `Rcpp::List` → returned to R → assembled into data frames by R code

The `S4RasterSource` adapter class (R-package-specific) wraps R raster objects and calls back into R for I/O, with per-band caching.

---

## 4. Capability Comparison

### Processing Architecture

| Aspect | Python | R |
|---|---|---|
| Uses core Processor pipeline | ✅ Full C++ pipeline | ❌ Own R-level loop + direct `RasterStats` |
| Feature-sequential strategy | ✅ | ✅ (only strategy available) |
| Raster-sequential strategy | ✅ (R-tree index, chunk iteration) | ❌ |
| Parallel processing (TBB) | ✅ (if built with TBB) | ❌ |

### Geometry Support

| Geometry type | Python | R |
|---|---|---|
| Polygon / MultiPolygon | ✅ | ✅ |
| LineString / MultiLineString | ✅ (length fractions) | ❌ (blocked at R level; C++ core supports it) |

### Input Formats

| Format | Python | R |
|---|---|---|
| GDAL datasets | ✅ | ❌ (uses terra/raster, not GDAL directly) |
| rasterio | ✅ | — |
| xarray/rioxarray | ✅ | — |
| numpy arrays | ✅ | — |
| terra `SpatRaster` | — | ✅ |
| raster `RasterLayer/Stack/Brick` | — | ✅ |
| File paths (auto-detect) | ✅ | ❌ |
| sf / sfc | — | ✅ |
| geopandas GeoDataFrame | ✅ | — |
| fiona | ✅ | — |
| GeoJSON dicts | ✅ | — |
| QGIS layers | ✅ | — |

### Output Formats

| Output | Python | R |
|---|---|---|
| In-memory list/dict (GeoJSON) | ✅ | — |
| pandas DataFrame / GeoDataFrame | ✅ | — |
| R data.frame / list of data.frames | — | ✅ |
| Write to file (GDAL) | ✅ (`output="gdal"`) | ❌ |
| QGIS vector layer | ✅ | — |
| Append input feature columns | ✅ (`include_cols`) | ✅ (`append_cols` / `include_cols`) |

### Statistics & Operations

Both packages expose the same ~30 built-in C++ statistics (mean, sum, min, max, median, mode, minority, variety, quantile, variance, stdev, coefficient_of_variation, frac, weighted_frac, weighted_mean, weighted_sum, weighted_stdev, weighted_variance, values, weights, coverage, center_x/y, min/max_center_x/y, cell_id, unique, count).

| Feature | Python | R |
|---|---|---|
| Built-in C++ stats | ✅ All ~30 | ✅ All ~30 |
| User-defined functions | ✅ Python callables (numpy masked arrays) | ✅ R functions (data.frame or vectors) |
| `summarize_df` mode (pass full df to UDF) | — | ✅ |
| `stack_apply` (per-layer vs all-layers) | — | ✅ |
| Operation arguments (`min_coverage_frac`, `default_value`, `coverage_weight`, etc.) | ✅ | ✅ (via `default_value`, `default_weight`, `coverage_area` params) |

### Additional Functions

| Function | Python | R |
|---|---|---|
| `exact_extract` | ✅ | ✅ |
| `exact_resample` | ❌ | ✅ |
| `rasterize_polygons` | ❌ | ✅ |
| `coverage_fraction` (return raster) | ❌ | ✅ |

### Other Differences

| Aspect | Python | R |
|---|---|---|
| CRS auto-reprojection | ⚠️ Warns only | ✅ Auto-reprojects polygons |
| Raster type flexibility | All values through variant (int8→float64) | Always `double` (via `NumericVector`) |
| GDAL block-size diagnostics | ❌ | ✅ (warns about poor caching) |
| Eager-load optimization | ❌ (handled by processor chunking) | ✅ (pre-reads small raster extents) |

---

## 5. Architectural Summary

The two packages represent **different generations of integration** with the same core algorithm:

**R package (`exactextractr`)**: An earlier, simpler integration. Vendors the core library source, uses only the computational primitives (`raster_cell_intersection` + `RasterStats`), and implements its own feature-iteration loop in R. This gives it tight integration with R's spatial ecosystem (sf, terra, raster) but means it doesn't benefit from the core library's `Processor` pipeline, raster-sequential strategy, or parallel processing.

**Python package**: A newer, deeper integration. Uses pybind11 to expose the full `Processor`/`Operation`/`OutputWriter` pipeline, meaning the C++ core drives the entire computation loop. Python only provides I/O adapters (via trampoline classes) and format-detection convenience. This gives it access to raster-sequential processing, potential TBB parallelism, and line geometry support.

The core C++ library is the clear **centre of gravity** — it contains the algorithmic sophistication (exact cell intersection, online statistical accumulators, grid algebra, chunked processing) while both language packages focus on bridging their respective spatial data ecosystems to the core's abstract interfaces.
