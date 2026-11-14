# Loose ends

Updated 2026-11-14. Kept by hand; tick items off or delete them as they close.

## PR #39: align with crate 0.3.1 (merged 2026-11-14)

- [x] Review and merge. CI was green on 2026-11-14 on all jobs (R on three OSes
      plus oldrel, ctest on Linux and macOS, pytest on three OSes × 2 Pythons,
      pkgdown) except R-devel, which was still running.
- [ ] `CRAN-SUBMISSION` is left untracked on purpose (it is in `.Rbuildignore`).
      Commit or delete it as you prefer.

## Crate 0.3.2 (PR #40, branch `align-crate-0.3.2`)

- [x] `burn_wkb_par` parses each blob on the thread that burns it; work read
      from the bytes (`detail::wkb_sample_extent`), tested equal to the
      parsed work. Approx, 100k small polygons, 32 threads: 52 -> 12 ms.
- [x] `tools/canon-ops` pinned `=0.3.2`; `fixtures/canon-ops.txt` identical
      but for its header. Python builds from 0.3.2.
- [x] C++ `burn_wkb_par` had no test; `test_canon_ops` now checks it against
      `burn_wkb` (byte orders, ISO and EWKB, Z/M, empty, null, truncated).
- [ ] PR #40 (opened 2026-11-14 after #39 merged): review and merge. ctest
      and R tests passed locally; pytest not run locally, so check CI.
- [ ] Crate: `burn_wkb_par` checks `id(len - 1)` up front, so with an id
      overflow among 3+ blobs its error names a different blob than
      `burn_wkb`'s (and a trailing empty blob can fail it where serial
      would not). The C++ port names the first non-empty overflowing blob,
      as serial does. Worth a crate issue or fix.
- [ ] Python: expose `burn_wkb_par` if Python gets a `threads` argument.
- [ ] Crate "Unreleased": the book now calls Natural Earth's Sudan a
      self-intersecting ring. Worth a line in the R docs on invalid input
      (winding fill), since `relate_cells` and coverage both follow the
      winding rule.

## BURN format (`design/burn-format/`)

- [ ] Notes towards a spec for a `.burn` file (sparse tables, grid, latent
      indexes, features for DE-9IM and zonal stats) and GDAL drivers. Open
      questions listed there; the name is a placeholder.

## Smaller items from the 0.3.1 work

- [ ] `merge_cells` is in the C++ core (`ops.hpp`) but not exported to R.
- [ ] `relate_cells()` on a cropped or sliced burn uses the window's own grid.
      It was tested only on a grid whose cell size is a power of two. On other
      grids the window's walls can be off by a rounding step. Test it, or
      document it.
- [ ] Old compiler warning: unused parameter `row_off` in `walk_ring`
      (`cpp/src/scanline.cpp`).

## GDAL planning (`design/gdal/`)

- [ ] Read `design/gdal/plan.md` and the experiments.
- [ ] Stage 0, all human-written per GDAL RFC 111 (facts in
      `design/gdal/upstream-drafts.md`):
  - [ ] GDAL: zonal-stats silently ignores rotation (count 0);
  - [ ] GDAL: rasterize ignores geolocation targets, even with `-to DST_METHOD`;
  - [ ] GEOS: backport the crash fixes to 3.14 and 3.15 (3.15.0 still crashes);
  - [ ] GEOS: interior cells returned as 0 on main;
  - [ ] GEOS: fractions on invalid input (discussion);
  - [ ] the conversation with Dan Baston: route A (GEOS) or B (GDAL-native),
        sparse output, shared `--pixels` vocabulary, ISciences relicensing.
- [ ] GDAL source: `~/gdal` is the primary tree (at 9bb42f5a7b, 2026-11-11).
      `/perm_storage/home/gdalbuilder/gdal` and the copies under
      `/rdsi/PRIVATE/raad/` and old session scratchpads need cleaning up.
- [ ] The GEOS builds used by experiments 01 and 02 were in the session
      scratchpad and will be gone. Rebuild with
      `cmake -S geos -B build -DBUILD_TESTING=OFF && ninja -C build install`
      (branches `3.14`, tag `3.15.0`, `main`).
