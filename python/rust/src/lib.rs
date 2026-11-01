//! Python bindings for controlledburn, built from the Rust crate.
//!
//! Thin marshalling only: WKB bytes in, flat column arrays out (assembled
//! into structured arrays by the Python wrapper), plus a `materialize`
//! that writes into a numpy buffer. All engine logic lives in the
//! `controlledburn` crate. Array coercion (dtype, contiguity) is done by
//! the `controlledburn._core` shim before calling in here.

use controlledburn as cb;
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::ffi;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySequence, PyTuple};

/// A held buffer view: released on drop.
struct View(ffi::Py_buffer);

impl Drop for View {
    fn drop(&mut self) {
        // The GIL is held whenever a View is dropped (views never cross
        // into `detach`), so releasing here is safe.
        unsafe { ffi::PyBuffer_Release(&mut self.0) }
    }
}

/// Raw (pointer, length) spans into held views, safe to read without the
/// GIL while the views are alive.
struct Spans(Vec<(*const u8, usize)>);
unsafe impl Send for Spans {}
unsafe impl Sync for Spans {}

impl Spans {
    fn slices(&self) -> impl Iterator<Item = &[u8]> {
        self.0.iter().map(|&(p, n)| {
            if p.is_null() || n == 0 {
                &[][..]
            } else {
                unsafe { std::slice::from_raw_parts(p, n) }
            }
        })
    }
}

/// Take a contiguous byte view of any buffer-protocol object. Its length is
/// the full byte length, whatever its item size.
fn view(obj: &Bound<'_, PyAny>, i: usize) -> PyResult<View> {
    let ptr = obj.as_ptr();
    if unsafe { ffi::PyObject_CheckBuffer(ptr) } == 0 {
        let name = obj
            .get_type()
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|_| "object".into());
        return Err(PyTypeError::new_err(format!(
            "WKB entry {i} is {name}, expected bytes-like or None"
        )));
    }
    let mut buf = std::mem::MaybeUninit::<ffi::Py_buffer>::uninit();
    // PyBUF_SIMPLE: the exporter must hand over one contiguous block, or fail.
    if unsafe { ffi::PyObject_GetBuffer(ptr, buf.as_mut_ptr(), ffi::PyBUF_SIMPLE) } != 0 {
        // Replace the exporter's BufferError with the binding's own message.
        let _ = PyErr::take(obj.py());
        return Err(PyValueError::new_err(format!(
            "WKB entry {i} is not a contiguous buffer"
        )));
    }
    Ok(View(unsafe { buf.assume_init() }))
}

fn list_spans(seq: &Bound<'_, PySequence>) -> PyResult<(Vec<View>, Spans)> {
    let n = seq.len()?;
    let mut keep = Vec::with_capacity(n);
    let mut spans = Vec::with_capacity(n);
    for i in 0..n {
        let h = seq.get_item(i)?;
        if h.is_none() {
            spans.push((std::ptr::null(), 0));
            continue;
        }
        let v = view(&h, i)?;
        spans.push((v.0.buf as *const u8, v.0.len as usize));
        keep.push(v);
    }
    Ok((keep, Spans(spans)))
}

/// Spans over Arrow (large_)binary buffers: a values buffer, int64 offsets
/// (length n + 1) and an optional uint8 validity mask (1 = valid). Zero copy.
fn arrow_spans(
    values: &Bound<'_, PyAny>,
    offsets: &PyReadonlyArray1<'_, i64>,
    valid: Option<&PyReadonlyArray1<'_, u8>>,
) -> PyResult<(View, Spans)> {
    let mut buf = std::mem::MaybeUninit::<ffi::Py_buffer>::uninit();
    if unsafe { ffi::PyObject_GetBuffer(values.as_ptr(), buf.as_mut_ptr(), ffi::PyBUF_SIMPLE) } != 0
    {
        let _ = PyErr::take(values.py());
        return Err(PyValueError::new_err("values must be a contiguous buffer"));
    }
    let keep = View(unsafe { buf.assume_init() });
    let base = keep.0.buf as *const u8;
    let data_len = keep.0.len as usize;
    let off = offsets.as_slice()?;
    if off.is_empty() {
        return Err(PyValueError::new_err(
            "offsets must have length n + 1 (>= 1)",
        ));
    }
    let n = off.len() - 1;
    let mask = match valid {
        Some(v) => {
            let m = v.as_slice()?;
            if m.len() != n {
                return Err(PyValueError::new_err(
                    "validity mask length must equal offsets length - 1",
                ));
            }
            Some(m)
        }
        None => None,
    };
    let mut spans = Vec::with_capacity(n);
    for i in 0..n {
        if mask.is_some_and(|m| m[i] == 0) {
            spans.push((std::ptr::null(), 0));
            continue;
        }
        let (s, e) = (off[i], off[i + 1]);
        if s < 0 || e < s || e as usize > data_len {
            return Err(PyValueError::new_err(
                "offset out of range for values buffer",
            ));
        }
        spans.push((unsafe { base.add(s as usize) }, (e - s) as usize));
    }
    Ok((keep, Spans(spans)))
}

fn parse_mode(mode: &str) -> PyResult<cb::BurnMode> {
    match mode {
        "coverage" => Ok(cb::BurnMode::Coverage),
        "approx" => Ok(cb::BurnMode::Approx),
        _ => Err(PyValueError::new_err("mode must be 'coverage' or 'approx'")),
    }
}

fn grid(
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
    ncol: i64,
    nrow: i64,
) -> PyResult<cb::GridSpec> {
    if ncol <= 0 || nrow <= 0 {
        return Err(PyValueError::new_err(
            "invalid grid: ncol and nrow must be positive",
        ));
    }
    let (ncol, nrow) = match (u32::try_from(ncol), u32::try_from(nrow)) {
        (Ok(c), Ok(r)) => (c, r),
        _ => {
            return Err(PyValueError::new_err(
                "invalid grid: ncol and nrow must fit in 32 bits",
            ))
        }
    };
    Ok(cb::GridSpec::new(xmin, ymin, xmax, ymax, ncol, nrow))
}

fn burn_spans(
    py: Python<'_>,
    spans: &Spans,
    grid: &cb::GridSpec,
    mode: cb::BurnMode,
    coalesce: bool,
) -> PyResult<cb::BurnResult> {
    let opts = cb::BurnOptions::parity(mode).with_coalesce_runs(coalesce);
    py.detach(|| cb::burn_wkb(spans.slices(), grid, opts))
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

fn col<'py, S, T: numpy::Element>(
    py: Python<'py>,
    v: &[S],
    f: impl Fn(&S) -> T,
) -> Bound<'py, PyArray1<T>> {
    PyArray1::from_vec(py, v.iter().map(f).collect())
}

fn assemble<'py>(py: Python<'py>, r: &cb::BurnResult) -> PyResult<Bound<'py, PyDict>> {
    let runs = PyDict::new(py);
    runs.set_item("row", col(py, &r.runs, |x| x.row))?;
    runs.set_item("col_start", col(py, &r.runs, |x| x.col_start))?;
    runs.set_item("col_end", col(py, &r.runs, |x| x.col_end))?;
    runs.set_item("id", col(py, &r.runs, |x| x.id))?;
    let edges = PyDict::new(py);
    edges.set_item("row", col(py, &r.edges, |x| x.row))?;
    edges.set_item("col", col(py, &r.edges, |x| x.col))?;
    edges.set_item("fraction", col(py, &r.edges, |x| x.fraction))?;
    edges.set_item("id", col(py, &r.edges, |x| x.id))?;
    let lines = PyDict::new(py);
    lines.set_item("row", col(py, &r.lines, |x| x.row))?;
    lines.set_item("col", col(py, &r.lines, |x| x.col))?;
    lines.set_item("length", col(py, &r.lines, |x| x.length))?;
    lines.set_item("id", col(py, &r.lines, |x| x.id))?;
    let points = PyDict::new(py);
    points.set_item("row", col(py, &r.points, |x| x.row))?;
    points.set_item("col", col(py, &r.points, |x| x.col))?;
    points.set_item("id", col(py, &r.points, |x| x.id))?;
    let notes = PyList::empty(py);
    for n in &r.notes {
        notes.append(PyTuple::new(
            py,
            [
                n.geom_index.into_pyobject(py)?.into_any(),
                n.message.as_str().into_pyobject(py)?.into_any(),
            ],
        )?)?;
    }
    let out = PyDict::new(py);
    out.set_item("runs", runs)?;
    out.set_item("edges", edges)?;
    out.set_item("lines", lines)?;
    out.set_item("points", points)?;
    out.set_item("notes", notes)?;
    Ok(out)
}

/// Envelope of WKB blobs; best effort: unparseable or empty blobs are skipped.
fn bbox(spans: &Spans) -> Option<(f64, f64, f64, f64)> {
    let mut b = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    let mut any = false;
    let mut add = |c: &cb::Coord| {
        if c.x.is_finite() && c.y.is_finite() {
            b.0 = b.0.min(c.x);
            b.1 = b.1.min(c.y);
            b.2 = b.2.max(c.x);
            b.3 = b.3.max(c.y);
            any = true;
        }
    };
    for s in spans.slices() {
        if s.is_empty() {
            continue;
        }
        let Ok(g) = cb::parse_wkb(s) else { continue };
        match &g {
            cb::Geometry::Point(c) => add(c),
            cb::Geometry::MultiPoint(cs) | cb::Geometry::LineString(cs) => {
                cs.iter().for_each(&mut add)
            }
            cb::Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(&mut add),
            cb::Geometry::Polygon(p) => p.rings.iter().flatten().for_each(&mut add),
            cb::Geometry::MultiPolygon(ps) => ps
                .iter()
                .flat_map(|p| p.rings.iter().flatten())
                .for_each(&mut add),
        }
    }
    any.then_some(b)
}

#[pyfunction]
#[pyo3(signature = (wkb_list, xmin, ymin, xmax, ymax, ncol, nrow, mode = "coverage", coalesce = false))]
#[allow(clippy::too_many_arguments)]
fn burn_wkb<'py>(
    py: Python<'py>,
    wkb_list: &Bound<'py, PySequence>,
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
    ncol: i64,
    nrow: i64,
    mode: &str,
    coalesce: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let g = grid(xmin, ymin, xmax, ymax, ncol, nrow)?;
    let m = parse_mode(mode)?;
    let (_keep, spans) = list_spans(wkb_list)?;
    let r = burn_spans(py, &spans, &g, m, coalesce)?;
    assemble(py, &r)
}

#[pyfunction]
#[pyo3(signature = (values, offsets, valid, xmin, ymin, xmax, ymax, ncol, nrow, mode = "coverage", coalesce = false))]
#[allow(clippy::too_many_arguments)]
fn burn_wkb_arrow<'py>(
    py: Python<'py>,
    values: &Bound<'py, PyAny>,
    offsets: PyReadonlyArray1<'py, i64>,
    valid: Option<PyReadonlyArray1<'py, u8>>,
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
    ncol: i64,
    nrow: i64,
    mode: &str,
    coalesce: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let g = grid(xmin, ymin, xmax, ymax, ncol, nrow)?;
    let m = parse_mode(mode)?;
    let (_keep, spans) = arrow_spans(values, &offsets, valid.as_ref())?;
    let r = burn_spans(py, &spans, &g, m, coalesce)?;
    assemble(py, &r)
}

#[pyfunction]
fn bbox_wkb(
    py: Python<'_>,
    wkb_list: &Bound<'_, PySequence>,
) -> PyResult<Option<(f64, f64, f64, f64)>> {
    let (_keep, spans) = list_spans(wkb_list)?;
    Ok(py.detach(|| bbox(&spans)))
}

#[pyfunction]
#[pyo3(signature = (values, offsets, valid))]
fn bbox_wkb_arrow<'py>(
    py: Python<'py>,
    values: &Bound<'py, PyAny>,
    offsets: PyReadonlyArray1<'py, i64>,
    valid: Option<PyReadonlyArray1<'py, u8>>,
) -> PyResult<Option<(f64, f64, f64, f64)>> {
    let (_keep, spans) = arrow_spans(values, &offsets, valid.as_ref())?;
    Ok(py.detach(|| bbox(&spans)))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn materialize<'py>(
    py: Python<'py>,
    runs_row: PyReadonlyArray1<'py, i32>,
    runs_col_start: PyReadonlyArray1<'py, i32>,
    runs_col_end: PyReadonlyArray1<'py, i32>,
    runs_id: PyReadonlyArray1<'py, i32>,
    edges_row: PyReadonlyArray1<'py, i32>,
    edges_col: PyReadonlyArray1<'py, i32>,
    edges_fraction: PyReadonlyArray1<'py, f32>,
    edges_id: PyReadonlyArray1<'py, i32>,
    ncol: i64,
    nrow: i64,
    values: Vec<f64>,
    fn_: &str,
    edge_policy: &str,
    threshold: f64,
    background: f64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let (rr, rs, re, ri) = (
        runs_row.as_slice()?,
        runs_col_start.as_slice()?,
        runs_col_end.as_slice()?,
        runs_id.as_slice()?,
    );
    if rs.len() != rr.len() || re.len() != rr.len() || ri.len() != rr.len() {
        return Err(PyValueError::new_err(
            "materialize: runs columns differ in length",
        ));
    }
    let (er, ec, ef, ei) = (
        edges_row.as_slice()?,
        edges_col.as_slice()?,
        edges_fraction.as_slice()?,
        edges_id.as_slice()?,
    );
    if ec.len() != er.len() || ef.len() != er.len() || ei.len() != er.len() {
        return Err(PyValueError::new_err(
            "materialize: edges columns differ in length",
        ));
    }
    let fn_ = match fn_ {
        "first" => cb::PixelFn::First,
        "last" => cb::PixelFn::Last,
        "sum" => cb::PixelFn::Sum,
        "min" => cb::PixelFn::Min,
        "max" => cb::PixelFn::Max,
        "count" => cb::PixelFn::Count,
        "any" => cb::PixelFn::Any,
        _ => {
            return Err(PyValueError::new_err(
                "fn must be one of first/last/sum/min/max/count/any",
            ))
        }
    };
    let edge_policy = match edge_policy {
        "threshold" => cb::EdgePolicy::Threshold,
        "fraction" => cb::EdgePolicy::Fraction,
        _ => {
            return Err(PyValueError::new_err(
                "edge_policy must be 'threshold' or 'fraction'",
            ))
        }
    };
    let g = grid(0.0, 0.0, 1.0, 1.0, ncol, nrow)?;
    let res = cb::BurnResult {
        runs: (0..rr.len())
            .map(|i| cb::GridRun {
                row: rr[i],
                col_start: rs[i],
                col_end: re[i],
                id: ri[i],
            })
            .collect(),
        edges: (0..er.len())
            .map(|i| cb::GridEdge {
                row: er[i],
                col: ec[i],
                fraction: ef[i],
                id: ei[i],
            })
            .collect(),
        ..Default::default()
    };
    let opts = cb::MaterializeOptions {
        fn_,
        edge_policy,
        threshold,
    };
    let mut buf = vec![background; g.ncol as usize * g.nrow as usize];
    let vals = (!values.is_empty()).then_some(values.as_slice());
    py.detach(|| cb::materialize(&res, &mut buf, g.ncol, g.nrow, vals, &opts))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    PyArray1::from_vec(py, buf).reshape([g.nrow as usize, g.ncol as usize])
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__doc__", "controlledburn bindings over the Rust crate (internal; use the controlledburn package API)")?;
    m.add_function(wrap_pyfunction!(burn_wkb, m)?)?;
    m.add_function(wrap_pyfunction!(burn_wkb_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(bbox_wkb, m)?)?;
    m.add_function(wrap_pyfunction!(bbox_wkb_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(materialize, m)?)?;
    Ok(())
}
