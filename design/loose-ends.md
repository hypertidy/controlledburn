# Loose ends

Updated 2026-11-14. Kept by hand; tick items off or delete them as they close.

## PR #39: align with crate 0.3.1 (branch `align-crate-0.3.1`)

- [ ] Review and merge. CI was green on 2026-11-14 on all jobs (R on three OSes
      plus oldrel, ctest on Linux and macOS, pytest on three OSes × 2 Pythons,
      pkgdown) except R-devel, which was still running.
- [ ] `CRAN-SUBMISSION` is left untracked on purpose (it is in `.Rbuildignore`).
      Commit or delete it as you prefer.

## Crate 0.3.2 (released 2026-11-12): the next code work

- [ ] `burn_wkb_par` in 0.3.2 parses each blob on the thread that burns it, and
      balances work read from the bytes without building geometries. The C++
      `burn_wkb_par` (`cpp/src/scanline.cpp`) parses every blob first and then
      runs `burn_par`, which is the slow path that 0.3.2 replaced (2 to 11× at
      128 threads, per the crate CHANGELOG). R's `burn(threads = )` goes
      through it (`src/scanline_shim.cpp`).
- [ ] Records are unchanged from 0.3.1, so `tools/canon-ops` needs only the pin
      bumped to `=0.3.2` and a regenerate. The thread checks in
      `cpp/tests/test_canon_ops.cpp` cover the rest.
- [ ] Python: `python/rust/Cargo.toml` to `controlledburn = "0.3.2"`. Expose
      `burn_wkb_par` if Python gets a `threads` argument (it has none yet).
- [ ] Crate "Unreleased": the book now calls Natural Earth's Sudan a
      self-intersecting ring. Worth a line in the R docs on invalid input
      (winding fill), since `relate_cells` and coverage both follow the
      winding rule.

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
- [ ] The GEOS builds used by experiments 01 and 02 were in the session
      scratchpad and will be gone. Rebuild with
      `cmake -S geos -B build -DBUILD_TESTING=OFF && ninja -C build install`
      (branches `3.14`, tag `3.15.0`, `main`).
