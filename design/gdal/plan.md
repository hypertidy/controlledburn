# Fast and exact rasterization in GDAL: a staged plan

Status: planning draft, 2026-11-14. A private working document for the
controlledburn maintainer, not text to send to GDAL (see "Constraints",
first item). The evidence is in [experiments/](experiments/README.md).

## Summary

GDAL has two rasterizers and no shared design.

- **`gdal_rasterize` / `gdal vector rasterize`**
  (`alg/gdalrasterize.cpp`, `alg/llrasterize.cpp`):
  - a cell-centre rule or ALL_TOUCHED;
  - REPLACE or ADD;
  - output to a dense band.
  - It is O(rows × edges) per polygon: 100× slower than a scanline with
    an active edge list on one large polygon (experiment 05).
  - ALL_TOUCHED has open edge-case bugs that Even Rouault has said need a
    substantial rework (OSGeo/gdal#10606, #14615).
- **`gdal raster zonal-stats`** (`alg/zonal.cpp`, Dan Baston, GDAL 3.12):
  - exact fractions from GEOS ≥ 3.14 `GEOSGridIntersectionFractions_r`,
    which is exactextract inside GEOS, into a dense float window per
    feature;
  - otherwise it calls the rasterizer above.

controlledburn is a model for one engine under both, and its two modes are
the two things asked for:

- **approx** (fasterize-style): cell-centre runs from a scanline. The cells
  are identical to GDAL's on the test corpora (experiment 05).
- **coverage** (exactextract-style): exact fractions for boundary cells,
  plus interior runs.

Both are sparse: runs and edges, O(perimeter) memory, never a dense buffer.
Index space carries both modes to rotated, rectilinear and curvilinear grids
without changing the engine (controlledburn-rs `design/grid-frontends.md`).

The plan has eight stages:

1. **Stage 0:** reports and a conversation with Dan Baston and the PSC.
2. **Stage 1:** small fixes to wrong answers found here.
3. **Stage 2:** a faster scanline with output identical to today's.
4. **Stage 3:** fractional coverage in rasterize, using the vocabulary
   zonal-stats already has.
5. **Stage 4:** a sparse burn API that rasterize, zonal-stats and cutline
   all use.
6. **Stage 5:** grids beyond the lattice (affine, geolocation, rectilinear,
   multidim).
7. **Stage 6:** behaviour changes saved for GDAL 4.0.
8. **Stage 7:** conservative resampling in warp, as a separate track.

Each stage stands on its own and is worth having if the later ones never
happen.

## Where things stand

### GDAL (master at ccd3811, 2026-11-13; local build 3.14.0dev)

| area | today | evidence |
|---|---|---|
| polygon fill | `GDALdllImageFilledPolygon`: for every row, loops over every edge; even-odd fill; intercepts rounded with `floor(x + 0.5)` | `alg/llrasterize.cpp`; exp. 05 |
| options | `ALL_TOUCHED`, `BURN_VALUE_FROM=Z`, `MERGE_ALG=REPLACE/ADD`, `CHUNKYSIZE`, `OPTIM=AUTO/RASTER/VECTOR` | `alg/gdalrasterize.cpp` |
| internal seam | the low-level fill already calls back once per span: `llScanlineFunc(info, y, x0, x1, value)`, and `llPointFunc` per pixel | `alg/gdal_alg_priv.h` |
| transform | vertices only, through `GDALCreateGenImgProjTransformer2(nullptr, dst)` or a caller's transformer; exact for affine | exp. 03 |
| geolocation targets | **not supported**: with no geotransform, GCPs or RPCs, `DST_METHOD` is forced to `NO_GEOTRANSFORM`, even over `-to DST_METHOD=GEOLOC_ARRAY`; nothing is burned | `apps/gdal_rasterize_lib.cpp`; exp. 04 |
| zonal-stats | `--pixels default/all-touched/fractional`; fractional needs GEOS ≥ 3.14; one dense coverage window per feature; for lines, coverage is a length | `alg/zonal.cpp`; #15148 (open) |
| zonal-stats, rotated raster | **silently wrong**: count 0, exit 0; the coverage window is built with the rotation terms zeroed | exp. 03 |
| warp | destination geolocation arrays supported (GDAL ≥ 3.5); cutline masks built with `GDALRasterizeGeometries` | `alg/gdalcutline.cpp` |
| multidim | `GetResampled` (2D coordinates through geolocation), `GetGridded`, `AsClassicDataset`, `GetMeshGrid`; `gdal mdim reproject` (3.14); no rasterize, no zonal-stats | `gcore/gdal.h`; `doc/source/programs` |
| open requests | faster rasterize (#7200); ALL_TOUCHED edges (#10606, #14615); MEAN and other merge algorithms (#9601, #7480); lines and points in zonal-stats (#15148) | triage snapshot 2026-10-29 |
| releases | 3.13.3 is current; 3.14.0 was due 2026-11-08; 3.15.0 is due 2027-05-03; a 4.0 milestone exists with a behaviour-break meta-ticket (#8440) | GitHub milestones |

### GEOS and Dan Baston

- **What exists:**
  - `GEOSGridIntersectionFractions_r(geom, xmin, ymin, xmax, ymax, nx, ny, float* buf)`
    and `GEOSSubdivideByGrid`, from GEOS 3.14 (`operation/grid`, ported
    from exactextract).
  - The API is dense and axis-aligned.
  - Since #1532 (2026-09) cells outside the geometry's box are zeroed.
- **Robustness:**
  - GEOS segfaults on valid polygons: 3.14.2 on an ordinary triangle,
    and **3.15.0, the current release**, on another case (exp. 02).
  - Dan merged robustness fixes on main on 2026-11-05 (#1589), following
    test cases from #1581 and #1593 (outside contributor). Both landed
    after 3.15.0.
  - Two fully covered interior cells still come back 0 on main (exp. 01,
    02).
  - On self-intersecting input GEOS returns fractions outside [0, 1]
    (exp. 01).
  - Speed: GEOS is 2.4× to 32× slower than the sparse burn on one large
    polygon and needs `4 × nx × ny` bytes per feature window (exp. 01).
- **Direction, as far as the record shows:**
  - Dan put the exact algorithm in GEOS and built GDAL's zonal-stats on
    it, with GDAL's rasterizer as the fallback. His issue #13059 names
    both paths.
  - His open PR #15148 adds lines and points to zonal-stats, and he is
    fixing GEOS robustness.
  - Nothing in the record says whether he wants fractional rasterize in
    GDAL, a sparse API in GEOS, or a GDAL-native engine. That is the
    first question to ask (Stage 0).
- **The PSC:**
  - adopted RFC 111 (LLM tool policy) in 2026, with +1 from DanB and
    EvenR among others;
  - Even's comment on #14615 invites a rasterization rework ("a reliable
    fix for that will involve reworking substantially the rasterization
    process"), addressed to Dan.
  - The controlledburn maintainer has pointed GDAL at controlledburn on
    #7200 (2023) and #14615 (2026).

## Constraints that shape the plan

1. **GDAL's LLM tool policy (RFC 111, adopted).**
   - The policy:
     - humans must be the primary authors;
     - "vibe-coded" contributions and agent-made submissions are banned;
     - LLMs may be used only for autocompletion or mechanical
       refactoring;
     - any use must be disclosed;
     - verbose text is read as LLM output and may be closed.
   - What that means here:
     - controlledburn's C++ core and its 0.3.1 port were written largely
       with Claude, and this document was written by Claude. Neither can
       go to GDAL as it stands.
     - The code is a **reference design, a test oracle and a source of
       fixtures**.
     - GDAL patches, RFC text, issues and comments must be written by
       their human author, short and in their own words, with any LLM
       use disclosed.
     - The drafts in [upstream-drafts.md](upstream-drafts.md) are only
       lists of facts to report.
   - **Route C below falls under the policy as well.** An LLM-written
     external library is not a contribution to the GDAL repository, but
     maintainers can still decline to depend on it.
2. **Licences.**
   - GDAL is MIT, GEOS is LGPL-2.1, controlledburn is Apache-2.0.
   - controlledburn's `cpp/src/ee/` and parts of the crate are derived
     from exactextract (Apache-2.0, ISciences, the company behind
     exactextract).
   - New GDAL code must be MIT. A clean implementation of a classic
     scanline (Stage 2) raises no question.
   - A GDAL-native exact engine that follows exactextract's design would
     need ISciences' agreement. Dan is ISciences, and `alg/zonal.cpp` is
     already "Copyright ISciences LLC, MIT", so this is a question, not a
     blocker.
3. **No Rust, optional GEOS.**
   - GDAL cannot require a Rust toolchain.
   - GEOS is an optional dependency, so anything built on GEOS is absent
     from builds without it, and fractional zonal-stats already is.
4. **Output compatibility.**
   - Existing rasterize output must not change outside a major release.
   - Each stage keeps today's cells by default and is checked
     cell for cell.
   - The tie convention (`floor(x + 0.5)`), the even-odd fill and the
     ALL_TOUCHED tolerance (1e-4 in places, per Even on #10606) are all
     observable.
5. **Maintenance.**
   - The scarce resource, as RFC 111 says.
   - Prefer changes that delete code (one engine under rasterize,
     zonal-stats and cutline) over changes that add a parallel one.

## Three routes for exact coverage in GDAL

| route | what | for | against |
|---|---|---|---|
| **A. GEOS** | keep `GEOSGridIntersectionFractions`; improve it in GEOS: robustness, a sparse or row-run output, valid fractions on invalid input; use it in rasterize too | already in GDAL; Dan owns it; LGPL dependency already accepted | dense API; optional dependency; GEOS release cycle; slower; crash history (exp. 01, 02) |
| **B. GDAL-native** | a sparse burn in `alg/`, written by a human, with approx and coverage sharing one boundary walk; checked against controlledburn, GEOS and exactextract | no dependency; fits the existing `llScanlineFunc` seam; threads; one engine for all consumers | review and maintenance load in GDAL; ISciences agreement for the exactextract-derived design; the most human authorship |
| **C. External library** | an optional dependency (a C++ library, or a C ABI over the crate) | least GDAL code | RFC 111 concerns; a new dependency for packagers; a second implementation to keep in step |

**Recommendation:**
- **Do Stage 2 (approx speed) as route B** in any case. It is small,
  self-contained, MIT-clean and gives identical output.
- **Start exact coverage on route A** (Stage 3), because it is what GDAL
  and Dan already chose. Bring Dan the evidence for a sparse output API in
  GEOS.
- **Leave B versus A for the long-term engine to Dan and Even** at Stage 4.
- **Treat C as a fallback only.**

## Stages

Each stage is a separate PR, or an RFC where noted, with its own tests. Dates
follow GDAL's six-month releases. Feature freeze is roughly a month before each
release.

### Stage 0: reports and conversation (now to January 2027)

No GDAL code. Deliverables, all human-written:

1. **Bug reports** (facts in [upstream-drafts.md](upstream-drafts.md)):
   - GDAL: zonal-stats on a rotated raster gives count 0 silently (exp. 03);
   - GDAL: rasterize ignores `DST_METHOD=GEOLOC_ARRAY` and geolocation
     metadata (exp. 04);
   - GEOS: interior cells come back 0 on main (exp. 02);
   - GEOS: ask whether #1589/#1593 can be backported to 3.14/3.15, since
     GDAL 3.12+ links those versions and 3.15.0 still crashes on a valid
     polygon (exp. 02);
   - GEOS: fractions outside [0, 1] on invalid input (exp. 01). Expected
     behaviour to agree first: winding, even-odd, or an error.
2. **A conversation with Dan** (and Even), not an RFC. Questions:
   - Where should exact coverage live long term: GEOS (route A) or
     `alg/` (route B)? Would ISciences relicense the exactextract-derived
     design for B?
   - Would GEOS take a sparse output, e.g. a callback per row of
     `(col0, col1, full)` runs and `(col, fraction)` boundary cells, next
     to the dense buffer?
   - Should rasterize and zonal-stats share one vocabulary
     (`--pixels default|all-touched|fractional`) and one engine?
   - Appetite for 4.0 changes to the tie rule and the meaning of
     ALL_TOUCHED (#8440)?
3. **A shared test corpus**, offered to GDAL autotest and GEOS tests:
   - the burnscars oracle cases (Apache-2.0, already published);
   - a subset of controlledburn's canon cases with GEOS-intersection
     truth per cell;
   - the cases from experiments 02 to 04.
   - The data is generated and its provenance is stated; licence and
     policy are fine to raise in the open.

Exit criterion: an answer on route A or B, and on the shared vocabulary.

### Stage 1: correctness fixes (GDAL 3.15, May 2027)

Small human-written PRs. No new features.

- **zonal-stats on rotated or sheared rasters.**
  - Minimum: refuse with an error when the geotransform has rotation
    terms.
  - Better: move each feature into pixel space with the inverse
    geotransform and compute coverage on the unit lattice. Exact for any
    affine map (exp. 03; crate `book/src/grids.md`). GEOS fractions are
    then fractions of the rotated cell's area.
  - Line lengths need rescaling per segment direction.
- **rasterize onto geolocation, GCP and RPC targets.**
  - Stop overriding an explicit `DST_METHOD`.
  - Add a densify-in-pixel-space option, e.g. a `SEGMENTIZE_PIXELS`
    transformer option or a rasterize option, because vertices alone are
    wrong by whole cells (exp. 04: 129 of 132 sliver cells; about one
    cell's step is exact for the cell-centre rule).
- **Document the conventions as they are:**
  - the cell-centre tie rule;
  - even-odd fill;
  - the ALL_TOUCHED tolerance;
  - coverage of lines as length in zonal-stats.

### Stage 2: a faster polygon scanline, identical output (GDAL 3.15 or 3.16)

The fasterize request (#7200), done the way GDAL can take it.

- **What:** replace the "every edge, every row" loop in
  `GDALdllImageFilledPolygon` with an edge table sorted by start row and an
  active edge list, which is the classic scanline. Keep the even-odd fill,
  the `floor(x + 0.5)` rounding, the horizontal-edge handling for
  ALL_TOUCHED and the `llScanlineFunc` callbacks unchanged.
- **Proof:**
  - autotest's rasterize tests;
  - a differential test against the old function, cell for cell, on the
    shared corpus and on random polygons with vertices steered onto
    centres, walls and corners (as controlledburn's canon generator
    does);
  - a benchmark from experiment 05: one 100k-vertex polygon on a
    10000² grid is about 1.5 s today, and controlledburn's sweep does it
    in 15 ms.
- **Optional:** threads across geometries in `GDALRasterizeGeometries` for
  the in-memory case. controlledburn's `burn_par` chunking (contiguous,
  work-balanced, concatenated in order) gives identical output. Writes to
  the band stay on one thread.
- **Size:** a few hundred lines in one file; no API change.

### Stage 3: fractional coverage in rasterize (GDAL 3.16, November 2027)

An RFC if the PSC wants one, since it adds user-facing options.

- **Options:**
  - `gdal vector rasterize --pixels default|all-touched|fractional`, the
    zonal-stats vocabulary;
  - `GDALRasterizeGeometries` option `PIXELS=FRACTIONAL` (ALL_TOUCHED
    stays as an alias).
- **Output:**
  - with REPLACE, the fraction, or burn value × fraction;
  - with ADD, burn × fraction summed, which gives area-weighted totals
    for overlapping features;
  - Float32 or Float64 target bands.
- **Merge algorithms:**
  - add MIN, MAX, FIRST, LAST and COUNT (fasterize's set; #7480);
  - MEAN needs a second accumulator (#9601) and is only well defined for
    a fresh output, not an existing raster (jratike80's point on #9601).
    Offer it only when creating the output.
- **Engine:** GEOS `GridIntersectionFractions` per feature window
  (route A), as zonal-stats does today. The feature is absent without
  GEOS ≥ 3.14, and needs the Stage 0 fixes.
- **Lines:** coverage is length in the cell, as in zonal-stats with
  #15148. Document how that interacts with `--add`.

### Stage 4: one sparse burn under rasterize, zonal-stats and cutline (GDAL 3.16 to 3.17)

The structural change. An RFC.

- **The API:**
  - a C++ iterator or callback in `alg/` that emits, per geometry, rows
    of interior runs `(row, col0, col1)` and boundary cells
    `(row, col, fraction)`, with the feature index;
  - this is `llScanlineFunc` generalised to carry fractions;
  - a C wrapper for the bindings.
- **Engine:**
  - route B, a GDAL-native sweep, if Stage 0 agreed;
  - otherwise route A, with GEOS's sparse output if GEOS added one, or
    GEOS's dense window turned into runs.
  - **Approx and coverage share one boundary walk**, as in controlledburn:
    ALL_TOUCHED becomes "fraction > threshold", with the threshold a
    documented parameter. Experiment 05: controlledburn's 1e-6 drops
    real slivers that GDAL's ALL_TOUCHED keeps.
  - Or use an exact predicate: controlledburn's `relate_cells` gives the
    DE-9IM per cell, and "interiors meet" is what #14615 asks
    ALL_TOUCHED to mean.
- **Consumers:**
  - rasterize: dense, materialised by span;
  - zonal-stats: no per-feature dense window, no MEM dataset round trip,
    weights in O(perimeter);
  - warp cutline masks: fractions give an anti-aliased cutline edge for
    free;
  - **a new sparse output:** `gdal vector rasterize --output x.parquet`
    writing the runs and edges as a table (row, col_start, col_end, fid,
    fraction), plus a "materialise" step that turns such a table into
    raster;
  - this is controlledburn's model, and GDAL gains a sparse raster
    representation through OGR without a new raster driver.
- **Threads:** across geometries, as in controlledburn's `burn_par`. The
  output order is deterministic. RFC 101 (raster thread safety) is not
  needed, because the burn is pure computation.

### Stage 5: grids beyond the lattice (GDAL 3.17+)

Everything works in index space, so no engine change is needed (crate
`design/grid-frontends.md`, rwarp `design/curvilinear-frontend.md`). Each
grid type needs a front end that moves geometry into index space:

| grid | GDAL source of the grid | front end | exact? |
|---|---|---|---|
| rotated or sheared affine | geotransform with rotation terms | inverse geotransform, vertices only | yes (fractions of the cell's area) |
| rectilinear | multidim 1D coordinate variables with uneven spacing; netCDF | cut edges at cell walls, map each axis | yes, given walls; with centres only, walls at midpoints by convention |
| curvilinear, bilinear cells | `GEOLOCATION` metadata; multidim 2D coordinates (`GetResampled` already reads them); swaths | geolocation inverse plus cutting edges into short pieces in index space | close (about 1% of a cell on a polar sector) |
| curvilinear, equal-area cells | same, given cell corners (CF bounds) | split each cell into four triangles at the mean of its corners | yes (crate note, "Equal-area cells") |
| GCP and RPC | GCP polynomial, TPS, RPC | densify, as for bilinear | to the tolerance |

What GDAL would gain:

- **rasterize onto geolocated datasets**, finishing the Stage 1 fix with
  edge cutting in index space rather than world-unit densification;
- **multidim:**
  - rasterize into a `GDALMDArray` slice, e.g. `gdal mdim rasterize`, or
    `--like` taking an mdim array;
  - zonal-stats over an mdim array: one burn, then every extra
    dimension (time, depth) reuses the same runs and edges, because the
    sparse result is independent of the band count;
  - this is where sparse output pays most: one burn, thousands of
    slices.
- **The cell definition is a parameter.** GDAL's geolocation model is
  bilinear. Exact fractions need straight-edged cells with the
  equal-area split. That is a documented choice, not a hidden one.
- **Wrap and folds:**
  - periodic longitude and tripolar seams are handled by unwrapping
    index space and folding columns back, then merging one record per
    cell (controlledburn's `merge_cells`);
  - GDAL's geoloc code already normalises longitude to ±180 (#12651),
    which helps here.

### Stage 6: GDAL 4.0 behaviour (through #8440)

Proposals for a major release only:

- **one documented tie rule** for the cell-centre test (half-open,
  left-inclusive, as in controlledburn and burnscars; or keep
  `floor(x + 0.5)`, but document and test it);
- **ALL_TOUCHED defined** as "the cell's interior meets the polygon's
  interior" (DE-9IM `T********`), with an explicit tolerance;
- **a documented fill rule**: GDAL is even-odd today, while
  controlledburn and GEOS-valid input assume non-zero winding. Choose
  per call, defaulting to today's.

### Stage 7: conservative resampling in warp (separate track, flag only)

- What it would add: `-r average` and `-r sum` with exact overlap weights
  between source and destination cells, i.e. first-order conservative
  remapping, ESMF's `conserve`.
- How it would work: in source index space, it is a coverage burn of
  each destination cell's outline.
- rwarp's `design/conservative.md` stages it, with an ESMF oracle.
- It depends on Stages 4 and 5 and is a large piece of work. It is listed
  so the burn API is designed with it in mind (many small convex
  polygons, scratch reuse), not proposed now.

## Other implications to flag

| topic | implication |
|---|---|
| invalid polygons | GDAL fills by even-odd; controlledburn by non-zero winding; GEOS fractions are undefined and come out < 0 or > 1 (exp. 01). Zonal-stats weights from invalid input are wrong today. Decide and document; consider `MakeValid` as an opt-in. |
| fraction precision | controlledburn and GEOS use f32 fractions. Summing weights over a large zone in f32 loses precision; accumulate in f64 (zonal.cpp already uses f64). |
| thresholds | every engine has one: controlledburn drops ≤ 1e-6, GDAL ALL_TOUCHED has 1e-4 in places, controlledburn snaps vertices within 1e-9 of a cell to walls. Make them parameters or documented constants, not accidents. |
| reproducibility | controlledburn turns off FMA contraction by pragma to get the same bits across compilers and architectures. GDAL issue #9822 (Intel vs M1 rasterize results) shows the same problem in GDAL. A differential test should run on every CI architecture. |
| lines | coverage of a line in a cell is a length (GEOS, controlledburn), not 0/1. `--add` then sums lengths. Document per mode (#15148). |
| points | in a cell or not; a point on a wall goes to one cell by a tie rule that must be documented. |
| curves | rasterize linearises curves since #14509; coverage should linearise to a tolerance in cells, not world units. |
| 3D / `BURN_VALUE_FROM=Z` | GDAL interpolates Z along spans; with fractions, Z per boundary cell needs a definition (at the cell's clipped centroid?). Out of scope until asked. |
| nodata and init | fractional output into an existing band: REPLACE by fraction is ambiguous where features overlap; ADD is well defined. |
| memory | GEOS's dense window costs 4·nx·ny bytes per feature (1.6 GB for one 20000² feature, exp. 01); sparse output is O(perimeter). This decides large-feature zonal-stats. |
| speed | exact coverage costs 5 to 10× the approx burn (exp. 05); on many small polygons it is slower than GDAL's ALL_TOUCHED on one thread. Keep approx the default; threads close the gap. |
| antimeridian and geodesic edges | out of scope for the engine. A front end splits or unwraps; spherical cell areas are weights (crate `Measure::cos_latitude`); great-circle edges are bigcurve's job. |
| Python and numpy | a sparse result maps naturally to numpy arrays of runs and edges, and to Arrow through OGR. Stage 4's table output reaches rasterio, xarray and rasterix with no new binding. |
| ecosystem | exactextract (Python and R), exactextractr, burnscars (Earthmover, Rust), controlledburn-rs, rasterix. A shared oracle corpus (Stage 0) is more valuable than any single implementation, and is the least contentious contribution. |

## What the maintainer can do next

1. Read the experiment results, re-run the scripts, and write the Stage 0
   reports yourself from [upstream-drafts.md](upstream-drafts.md),
   disclosing that the investigation used an LLM.
2. Open the conversation with Dan: route A or B, the sparse GEOS output,
   one shared vocabulary, ISciences relicensing.
3. If Stage 2 is wanted early: write the active-edge-table fill by hand
   (it is a textbook algorithm). Use controlledburn approx and the old
   function as the two oracles for a differential test.
4. Keep controlledburn as the place to try designs: a front end per grid
   type, threshold behaviour, `relate_cells` as an exact ALL_TOUCHED.
