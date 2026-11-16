//! canon-ops: the controlledburn crate's answers for the operations the
//! C++ core ports from it, on the canonical inputs.
//!
//! Reads the inputs (grid and WKB) of `fixtures/canon-tables.txt` and
//! writes `fixtures/canon-ops.txt`: for each case, the records of burns
//! with the default and other options, burns onto windows, tiles,
//! `merge_cells` of a parity burn, `crop`, and `relate_cells`. The C++
//! core checks itself against the file (`cpp/tests/test_canon_ops.cpp`).
//!
//!     cargo run --release -- ../../fixtures/canon-tables.txt > ../../fixtures/canon-ops.txt
//!
//! Format, one record per line: `case`, `grid` (four f64 bit patterns,
//! ncol, nrow) and `wkb` lines as in canon-tables.txt, then for each
//! operation an `op ...` line followed by its records (`run`, `edge`,
//! `line`, `point` as in canon-tables.txt, fractions and lengths as f32
//! bits; `note <id>`; `tile <tile_row> <tile_col> <row_off> <col_off>
//! <grid bits> <ncol> <nrow>` before each tile's records; `subgrid ...`
//! for a crop; `rel <row> <col> <id> <matrix>`), and `end`.

use std::fmt::Write as _;

use controlledburn::{
    burn, burn_window, crop, merge_cells, parse_wkb, relate_cells, slice_tiles, BurnMode, BurnOptions,
    BurnResult, Geometry, GridSpec,
};

fn bits(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn grid_line(g: &GridSpec) -> String {
    format!("{} {} {} {} {} {}", bits(g.xmin), bits(g.ymin), bits(g.xmax), bits(g.ymax), g.ncol, g.nrow)
}

fn records(out: &mut String, r: &BurnResult, notes: bool) {
    for x in &r.runs {
        writeln!(out, "run {} {} {} {}", x.row, x.col_start, x.col_end, x.id).unwrap();
    }
    for x in &r.edges {
        writeln!(out, "edge {} {} {:08x} {}", x.row, x.col, x.fraction.to_bits(), x.id).unwrap();
    }
    for x in &r.lines {
        writeln!(out, "line {} {} {:08x} {}", x.row, x.col, x.length.to_bits(), x.id).unwrap();
    }
    for x in &r.points {
        writeln!(out, "point {} {} {}", x.row, x.col, x.id).unwrap();
    }
    if notes {
        for n in &r.notes {
            writeln!(out, "note {}", n.geom_index).unwrap();
        }
    }
}

fn mode_name(m: BurnMode) -> &'static str {
    match m {
        BurnMode::Approx => "approx",
        _ => "coverage",
    }
}

fn unhex(h: &str) -> Vec<u8> {
    (0..h.len() / 2).map(|i| u8::from_str_radix(&h[2 * i..2 * i + 2], 16).unwrap()).collect()
}

fn from_bits(h: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(h, 16).unwrap())
}

struct Case {
    name: String,
    grid: GridSpec,
    wkb: Vec<String>,
}

fn read_cases(path: &str) -> Vec<Case> {
    let text = std::fs::read_to_string(path).expect("read canon-tables.txt");
    let mut cases: Vec<Case> = Vec::new();
    for line in text.lines() {
        let mut it = line.split_whitespace();
        match it.next() {
            Some("case") => cases.push(Case {
                name: it.next().unwrap().to_string(),
                grid: GridSpec::new(0.0, 0.0, 1.0, 1.0, 1, 1),
                wkb: Vec::new(),
            }),
            Some("grid") => {
                let v: Vec<&str> = it.collect();
                cases.last_mut().unwrap().grid = GridSpec::new(
                    from_bits(v[0]),
                    from_bits(v[1]),
                    from_bits(v[2]),
                    from_bits(v[3]),
                    v[4].parse().unwrap(),
                    v[5].parse().unwrap(),
                );
            }
            Some("wkb") => cases.last_mut().unwrap().wkb.push(it.next().unwrap().to_string()),
            _ => {}
        }
    }
    cases
}

/// The windows checked for a grid: a middle block, one row, one column,
/// and the bottom-right corner cell.
fn windows(g: &GridSpec) -> Vec<(u32, u32, u32, u32)> {
    let (nr, nc) = (g.nrow, g.ncol);
    vec![
        (nr / 3, (2 * nr / 3 + 1).min(nr), nc / 4, (3 * nc / 4 + 1).min(nc)),
        (nr / 2, nr / 2 + 1, 0, nc),
        (0, nr, nc / 2, nc / 2 + 1),
        (nr - 1, nr, nc - 1, nc),
    ]
}

/// WKB (little endian) of a polygon or multipolygon given as rings of
/// (x, y), each ring closed by the writer.
fn wkb_polys(polys: &[Vec<Vec<(f64, f64)>>]) -> String {
    let mut b: Vec<u8> = Vec::new();
    let body = |b: &mut Vec<u8>, rings: &Vec<Vec<(f64, f64)>>| {
        b.push(1);
        b.extend_from_slice(&3u32.to_le_bytes());
        b.extend_from_slice(&(rings.len() as u32).to_le_bytes());
        for r in rings {
            b.extend_from_slice(&((r.len() + 1) as u32).to_le_bytes());
            for &(x, y) in r.iter().chain(std::iter::once(&r[0])) {
                b.extend_from_slice(&x.to_le_bytes());
                b.extend_from_slice(&y.to_le_bytes());
            }
        }
    };
    if polys.len() == 1 {
        body(&mut b, &polys[0]);
    } else {
        b.push(1);
        b.extend_from_slice(&6u32.to_le_bytes());
        b.extend_from_slice(&(polys.len() as u32).to_le_bytes());
        for p in polys {
            body(&mut b, p);
        }
    }
    b.iter().map(|x| format!("{x:02x}")).collect()
}

/// Fixed cases for relate_cells' corners: groups of candidate cells no
/// neighbour reaches (a one-cell grid, a polygon whose boundary touches
/// every cell), polygons inside one cell, a hole that is one cell, and
/// boundaries along walls and through corners.
fn extra_cases() -> Vec<Case> {
    let sq = |x0: f64, y0: f64, x1: f64, y1: f64| vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
    let mk = |name: &str, grid: GridSpec, polys: Vec<Vec<Vec<Vec<(f64, f64)>>>>| Case {
        name: name.to_string(),
        grid,
        wkb: polys.iter().map(|p| wkb_polys(p)).collect(),
    };
    vec![
        mk("x_one_cell_grid", GridSpec::new(0.0, 0.0, 1.0, 1.0, 1, 1), vec![vec![vec![vec![(0.2, 0.2), (0.8, 0.3), (0.5, 0.9)]]]]),
        mk(
            "x_one_cell_grid_cover",
            GridSpec::new(0.0, 0.0, 1.0, 1.0, 1, 1),
            vec![vec![vec![vec![(-1.0, -1.0), (2.0, -1.0), (0.5, 0.5), (2.0, 2.0), (-1.0, 2.0)]]]],
        ),
        mk(
            "x_two_by_two_hole",
            GridSpec::new(0.0, 0.0, 2.0, 2.0, 2, 2),
            vec![vec![vec![sq(0.1, 0.1, 1.9, 1.9), sq(0.7, 0.7, 1.3, 1.3)]]],
        ),
        mk(
            "x_three_by_one_all_touched",
            GridSpec::new(0.0, 0.0, 3.0, 1.0, 3, 1),
            vec![vec![vec![vec![(-0.5, 0.5), (1.5, -0.5), (3.5, 0.5), (1.5, 1.5)]]]],
        ),
        mk(
            "x_inside_one_cell",
            GridSpec::new(0.0, 0.0, 4.0, 4.0, 4, 4),
            vec![vec![vec![vec![(1.2, 1.2), (1.8, 1.3), (1.5, 1.9)]]], vec![vec![sq(2.0, 2.0, 3.0, 3.0)]]],
        ),
        mk(
            "x_hole_is_a_cell",
            GridSpec::new(0.0, 0.0, 5.0, 5.0, 5, 5),
            vec![vec![vec![sq(1.0, 1.0, 4.0, 4.0), sq(2.0, 2.0, 3.0, 3.0)]]],
        ),
        mk(
            "x_walls_and_corners",
            GridSpec::new(0.0, 0.0, 4.0, 4.0, 4, 4),
            vec![
                vec![vec![sq(1.0, 1.0, 3.0, 3.0)]],
                vec![vec![vec![(2.0, 0.0), (4.0, 2.0), (2.0, 4.0), (0.0, 2.0)]]],
                vec![vec![sq(0.0, 0.0, 1.0, 1.0)], vec![sq(1.0, 1.0, 2.0, 2.0)]],
            ],
        ),
        mk(
            "x_ring_of_cells",
            GridSpec::new(0.0, 0.0, 3.0, 3.0, 3, 3),
            vec![vec![vec![sq(0.5, 0.5, 2.5, 2.5), sq(0.9, 0.9, 2.1, 2.1)]]],
        ),
    ]
}

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| "fixtures/canon-tables.txt".into());
    let mut cases = read_cases(&path);
    cases.extend(extra_cases());
    let mut out = String::new();
    out.push_str("# controlledburn crate 0.3.2 answers on the canon-tables.txt inputs, from tools/canon-ops; do not edit\n");
    for c in &cases {
        let blobs: Vec<Vec<u8>> = c.wkb.iter().map(|h| unhex(h)).collect();
        let geoms: Vec<Geometry> = blobs.iter().map(|b| parse_wkb(b).expect("canon WKB parses")).collect();
        let g = &c.grid;
        writeln!(out, "case {}", c.name).unwrap();
        writeln!(out, "grid {}", grid_line(g)).unwrap();
        for h in &c.wkb {
            writeln!(out, "wkb {h}").unwrap();
        }
        for mode in [BurnMode::Coverage, BurnMode::Approx] {
            let m = mode_name(mode);
            let variants = [
                ("default", BurnOptions::default().with_mode(mode)),
                ("merge_only", BurnOptions::parity(mode).with_merge_parts(true)),
                ("coalesce_only", BurnOptions::parity(mode).with_coalesce_runs(true).with_id_offset(7)),
            ];
            for (name, opts) in variants {
                writeln!(out, "op burn {m} {name}").unwrap();
                records(&mut out, &burn(&geoms, g, opts).unwrap(), true);
            }
            let opts = BurnOptions::default().with_mode(mode);
            let whole = burn(&geoms, g, opts).unwrap();
            for (r0, r1, c0, c1) in windows(g) {
                writeln!(out, "op window {m} {r0} {r1} {c0} {c1}").unwrap();
                records(&mut out, &burn_window(&geoms, g, r0..r1, c0..c1, opts).unwrap(), true);
            }
            writeln!(out, "op tiles {m} 3 2").unwrap();
            for t in slice_tiles(&whole, g, 3, 2) {
                writeln!(out, "tile {} {} {} {} {}", t.tile_row, t.tile_col, t.row_off, t.col_off, grid_line(&t.grid))
                    .unwrap();
                records(&mut out, &t.result, false);
            }
            writeln!(out, "op merge_cells {m}").unwrap();
            let mut parity = burn(&geoms, g, BurnOptions::parity(mode)).unwrap();
            merge_cells(&mut parity);
            records(&mut out, &parity, false);
            let (dx, dy) = (g.dx(), g.dy());
            let target = [g.xmin + 1.3 * dx, g.ymin + 0.5 * dy, g.xmax - 2.7 * dx, g.ymax - 1.0 * dy];
            writeln!(out, "op crop {m} {} {} {} {}", bits(target[0]), bits(target[1]), bits(target[2]), bits(target[3]))
                .unwrap();
            if let Some((sub, sg)) = crop(&whole, g, target) {
                writeln!(out, "subgrid {}", grid_line(&sg)).unwrap();
                records(&mut out, &sub, false);
            }
            writeln!(out, "op relate {m}").unwrap();
            for rel in relate_cells(&whole, &geoms, g, opts).unwrap() {
                writeln!(out, "rel {} {} {} {}", rel.row, rel.col, rel.id, rel.matrix).unwrap();
            }
        }
        writeln!(out, "end").unwrap();
    }
    print!("{out}");
}
