// _core.cpp -- pybind11 bindings for the controlledburn C++ core
//
// Thin marshalling only: WKB bytes in, flat column arrays out (assembled
// into structured arrays by the Python wrapper), plus a materialize()
// that writes into a numpy buffer. All engine logic lives in the core.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstring>
#include <vector>

#include "controlledburn/controlledburn.hpp"

namespace py = pybind11;
namespace cb = controlledburn;

namespace {

template <typename T, typename Struct, typename Member>
py::array_t<T> column(const std::vector<Struct>& v, Member Struct::* m) {
    py::array_t<T> out(static_cast<py::ssize_t>(v.size()));
    auto buf = out.template mutable_unchecked<1>();
    for (size_t i = 0; i < v.size(); i++) buf(static_cast<py::ssize_t>(i)) = static_cast<T>(v[i].*m);
    return out;
}

cb::BurnMode parse_mode(const std::string& mode) {
    if (mode == "coverage") return cb::BurnMode::Coverage;
    if (mode == "approx")   return cb::BurnMode::Approx;
    throw py::value_error("mode must be 'coverage' or 'approx'");
}

py::dict assemble_result(const cb::BurnResult& res) {
    py::dict runs;
    runs["row"]       = column<int32_t>(res.runs, &cb::GridRun::row);
    runs["col_start"] = column<int32_t>(res.runs, &cb::GridRun::col_start);
    runs["col_end"]   = column<int32_t>(res.runs, &cb::GridRun::col_end);
    runs["id"]        = column<int32_t>(res.runs, &cb::GridRun::id);

    py::dict edges;
    edges["row"]      = column<int32_t>(res.edges, &cb::GridEdge::row);
    edges["col"]      = column<int32_t>(res.edges, &cb::GridEdge::col);
    edges["fraction"] = column<float>(res.edges, &cb::GridEdge::fraction);
    edges["id"]       = column<int32_t>(res.edges, &cb::GridEdge::id);

    py::dict lines;
    lines["row"]    = column<int32_t>(res.lines, &cb::GridLine::row);
    lines["col"]    = column<int32_t>(res.lines, &cb::GridLine::col);
    lines["length"] = column<float>(res.lines, &cb::GridLine::length);
    lines["id"]     = column<int32_t>(res.lines, &cb::GridLine::id);

    py::dict points;
    points["row"] = column<int32_t>(res.points, &cb::GridPoint::row);
    points["col"] = column<int32_t>(res.points, &cb::GridPoint::col);
    points["id"]  = column<int32_t>(res.points, &cb::GridPoint::id);

    py::list notes;
    for (const auto& n : res.notes) {
        notes.append(py::make_tuple(n.geom_index, n.message));
    }

    py::dict out;
    out["runs"] = runs;
    out["edges"] = edges;
    out["lines"] = lines;
    out["points"] = points;
    out["notes"] = notes;
    return out;
}

// burn_wkb over a list of bytes-like objects. Returns a dict of dicts of
// column arrays plus the notes list; the Python wrapper assembles
// structured arrays and issues warnings.
py::dict py_burn_wkb(py::sequence wkb_list,
                     double xmin, double ymin, double xmax, double ymax,
                     int ncol, int nrow,
                     const std::string& mode = "coverage") {
    cb::GridSpec grid{xmin, ymin, xmax, ymax, ncol, nrow};

    cb::BurnMode burn_mode = parse_mode(mode);

    // Hold buffer views for the duration of the call; WKBSpan is
    // non-owning.
    std::vector<py::buffer_info> keep;
    std::vector<cb::WKBSpan> spans;
    keep.reserve(py::len(wkb_list));
    spans.reserve(py::len(wkb_list));

    for (py::handle h : wkb_list) {
        if (h.is_none()) {
            spans.push_back({nullptr, 0});
            continue;
        }
        py::buffer buf = py::reinterpret_borrow<py::buffer>(h);
        keep.push_back(buf.request());
        const py::buffer_info& info = keep.back();
        spans.push_back({
            static_cast<const uint8_t*>(info.ptr),
            static_cast<size_t>(info.size)
        });
    }

    cb::BurnResult res;
    {
        py::gil_scoped_release release;
        res = cb::burn_wkb(spans, grid, burn_mode);
    }

    return assemble_result(res);
}

// Build non-owning WKBSpans over an Arrow (large_)binary array's buffers:
// a contiguous values buffer, an int64 offsets array (length n + 1), and
// an optional per-element validity mask (uint8, 1 = valid). Returns the
// values buffer_info, which the caller must keep alive for as long as the
// spans are used (the spans point straight into it -- zero copy).
py::buffer_info build_arrow_spans(
        py::buffer values,
        const py::array_t<int64_t, py::array::c_style | py::array::forcecast>& offsets,
        py::object valid,
        std::vector<cb::WKBSpan>& spans) {
    py::buffer_info vinfo = values.request();
    const uint8_t* base = static_cast<const uint8_t*>(vinfo.ptr);
    const size_t data_len =
        static_cast<size_t>(vinfo.size) * static_cast<size_t>(vinfo.itemsize);

    auto off = offsets.unchecked<1>();
    if (off.shape(0) < 1) {
        throw py::value_error("offsets must have length n + 1 (>= 1)");
    }
    const py::ssize_t n = off.shape(0) - 1;

    const uint8_t* validp = nullptr;
    py::array_t<uint8_t, py::array::c_style | py::array::forcecast> valid_arr;
    if (!valid.is_none()) {
        valid_arr = valid.cast<
            py::array_t<uint8_t, py::array::c_style | py::array::forcecast>>();
        if (valid_arr.shape(0) != n) {
            throw py::value_error(
                "validity mask length must equal offsets length - 1");
        }
        validp = valid_arr.data();
    }

    spans.clear();
    spans.reserve(static_cast<size_t>(n));
    for (py::ssize_t i = 0; i < n; i++) {
        if (validp && validp[i] == 0) {
            spans.push_back({nullptr, 0});
            continue;
        }
        const int64_t start = off(i);
        const int64_t end = off(i + 1);
        if (start < 0 || end < start ||
            static_cast<size_t>(end) > data_len) {
            throw py::value_error("offset out of range for values buffer");
        }
        spans.push_back({base + static_cast<size_t>(start),
                         static_cast<size_t>(end - start)});
    }
    return vinfo;
}

// burn from Arrow WKB buffers -- the zero-copy Arrow path, no per-geometry
// Python objects. The Python wrapper (via nanoarrow) supplies the buffers.
py::dict py_burn_wkb_arrow(
        py::buffer values,
        py::array_t<int64_t, py::array::c_style | py::array::forcecast> offsets,
        py::object valid,
        double xmin, double ymin, double xmax, double ymax,
        int ncol, int nrow,
        const std::string& mode = "coverage") {
    cb::GridSpec grid{xmin, ymin, xmax, ymax, ncol, nrow};
    cb::BurnMode burn_mode = parse_mode(mode);

    std::vector<cb::WKBSpan> spans;
    py::buffer_info keep = build_arrow_spans(values, offsets, valid, spans);

    cb::BurnResult res;
    {
        py::gil_scoped_release release;
        res = cb::burn_wkb(spans, grid, burn_mode);
    }
    return assemble_result(res);
}

// Envelope (bounding box) of Arrow WKB buffers via the core's WKB scan.
// Returns (xmin, ymin, xmax, ymax), or None when no finite coordinate is
// found. Lets the Python side derive a default extent without shapely.
py::object py_bbox_wkb_arrow(
        py::buffer values,
        py::array_t<int64_t, py::array::c_style | py::array::forcecast> offsets,
        py::object valid) {
    std::vector<cb::WKBSpan> spans;
    py::buffer_info keep = build_arrow_spans(values, offsets, valid, spans);

    cb::BBox bb;
    {
        py::gil_scoped_release release;
        bb = cb::bbox_wkb(spans);
    }
    if (!bb.valid) return py::none();
    return py::make_tuple(bb.xmin, bb.ymin, bb.xmax, bb.ymax);
}

// Materialize runs/edges column arrays into a (nrow, ncol) float64 array.
py::array_t<double> py_materialize(
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> runs_row,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> runs_col_start,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> runs_col_end,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> runs_id,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> edges_row,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> edges_col,
    py::array_t<float,   py::array::c_style | py::array::forcecast> edges_fraction,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> edges_id,
    int ncol, int nrow,
    std::vector<double> values,
    const std::string& fn,
    const std::string& edge_policy,
    double threshold,
    double background) {

    cb::BurnResult res;
    auto rr = runs_row.unchecked<1>();
    auto rs = runs_col_start.unchecked<1>();
    auto re = runs_col_end.unchecked<1>();
    auto ri = runs_id.unchecked<1>();
    res.runs.reserve(static_cast<size_t>(rr.shape(0)));
    for (py::ssize_t i = 0; i < rr.shape(0); i++) {
        res.runs.push_back({rr(i), rs(i), re(i), ri(i)});
    }
    auto er = edges_row.unchecked<1>();
    auto ec = edges_col.unchecked<1>();
    auto ef = edges_fraction.unchecked<1>();
    auto ei = edges_id.unchecked<1>();
    res.edges.reserve(static_cast<size_t>(er.shape(0)));
    for (py::ssize_t i = 0; i < er.shape(0); i++) {
        res.edges.push_back({er(i), ec(i), ef(i), ei(i)});
    }

    cb::MaterializeOptions opts;
    if      (fn == "first") opts.fn = cb::PixelFn::First;
    else if (fn == "last")  opts.fn = cb::PixelFn::Last;
    else if (fn == "sum")   opts.fn = cb::PixelFn::Sum;
    else if (fn == "min")   opts.fn = cb::PixelFn::Min;
    else if (fn == "max")   opts.fn = cb::PixelFn::Max;
    else if (fn == "count") opts.fn = cb::PixelFn::Count;
    else if (fn == "any")   opts.fn = cb::PixelFn::Any;
    else throw py::value_error("fn must be one of first/last/sum/min/max/count/any");

    if      (edge_policy == "threshold") opts.edge_policy = cb::EdgePolicy::Threshold;
    else if (edge_policy == "fraction")  opts.edge_policy = cb::EdgePolicy::Fraction;
    else throw py::value_error("edge_policy must be 'threshold' or 'fraction'");
    opts.threshold = threshold;

    py::array_t<double> out({static_cast<py::ssize_t>(nrow), static_cast<py::ssize_t>(ncol)});
    double* buf = out.mutable_data();
    std::fill(buf, buf + static_cast<size_t>(nrow) * ncol, background);

    {
        py::gil_scoped_release release;
        cb::materialize(res, buf, ncol, nrow, values, opts);
    }
    return out;
}

} // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "controlledburn C++ core bindings (internal; use the controlledburn package API)";
    m.def("burn_wkb", &py_burn_wkb,
          py::arg("wkb_list"),
          py::arg("xmin"), py::arg("ymin"), py::arg("xmax"), py::arg("ymax"),
          py::arg("ncol"), py::arg("nrow"),
          py::arg("mode") = "coverage");
    m.def("burn_wkb_arrow", &py_burn_wkb_arrow,
          py::arg("values"), py::arg("offsets"), py::arg("valid"),
          py::arg("xmin"), py::arg("ymin"), py::arg("xmax"), py::arg("ymax"),
          py::arg("ncol"), py::arg("nrow"),
          py::arg("mode") = "coverage");
    m.def("bbox_wkb_arrow", &py_bbox_wkb_arrow,
          py::arg("values"), py::arg("offsets"), py::arg("valid"));
    m.def("materialize", &py_materialize,
          py::arg("runs_row"), py::arg("runs_col_start"),
          py::arg("runs_col_end"), py::arg("runs_id"),
          py::arg("edges_row"), py::arg("edges_col"),
          py::arg("edges_fraction"), py::arg("edges_id"),
          py::arg("ncol"), py::arg("nrow"),
          py::arg("values"), py::arg("fn"), py::arg("edge_policy"),
          py::arg("threshold"), py::arg("background"));
}
